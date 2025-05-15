/**
 * @file BufrObsIO.hpp
 * @brief BUFR format backend for observation I/O operations
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains a backend implementation for handling observations in
 * BUFR (Binary Universal Form for the Representation of meteorological data)
 * format.
 */

#pragma once

#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "BufrFortranWrapper.hpp"
#include "DateTime.hpp"
#include "Duration.hpp"
#include "ObsIOConcepts.hpp"

namespace metada::backends::io {

using ObsRecord = framework::ObsRecord;

/**
 * @brief Helper function for comparing doubles with BUFR missing value
 *
 * @param value The value to check
 * @return true if the value is not a BUFR missing value or NaN
 */
inline bool isValidValue(double value) {
  // Use machine epsilon from standard library
  const double epsilon = std::numeric_limits<double>::epsilon() * 100.0;
  // Scale factor applied to handle typical meteorological data precision needs
  return (std::abs(value - R8BFMS) > std::abs(R8BFMS) * epsilon) &&
         !std::isnan(value);
}

/**
 * @brief Backend implementation for handling BUFR format observation I/O
 *
 * @details This class provides functionality for reading and writing
 * observation data in BUFR format, which is commonly used in meteorology
 * and weather forecasting. It satisfies the ObservationIOBackendImpl concept
 * required by the ObservationIO adapter.
 *
 * BUFR is a binary format standardized by the World Meteorological Organization
 * (WMO) for the representation and exchange of meteorological data.
 *
 * @tparam ConfigBackend Configuration backend type that provides access to
 * configuration parameters
 */
template <typename ConfigBackend>
class BufrObsIO {
 public:
  /**
   * @brief Default constructor deleted
   *
   * @details BUFR I/O backends must be initialized with specific parameters.
   */
  BufrObsIO() = delete;

  /**
   * @brief Constructor with initialization parameters
   *
   * @details Initializes the BUFR I/O backend with the provided parameters,
   * which may include configuration options specific to BUFR handling.
   *
   * @param config Configuration parameters for BUFR processing
   */
  explicit BufrObsIO(ConfigBackend&& config) : config_(std::move(config)) {
    // Extract filename from config
    filename_ = config_.Get("filename").asString();
    bufrWrapper_ = std::make_unique<BufrFortranWrapper>();
  }

  /**
   * @brief Copy constructor deleted
   *
   * @details BUFR I/O backends are not meant to be copied.
   */
  BufrObsIO(const BufrObsIO&) = delete;

  /**
   * @brief Copy assignment deleted
   *
   * @details BUFR I/O backends are not meant to be copied.
   */
  BufrObsIO& operator=(const BufrObsIO&) = delete;

  /**
   * @brief Move constructor
   *
   * @details Allows efficient transfer of resources between BUFR I/O backends.
   *
   * @param other The BUFR I/O backend to move from
   */
  BufrObsIO(BufrObsIO&& other) noexcept
      : config_(std::move(other.config_)),
        filename_(std::move(other.filename_)),
        bufrWrapper_(std::move(other.bufrWrapper_)) {}

  /**
   * @brief Move assignment
   *
   * @details Allows efficient transfer of resources between BUFR I/O backends.
   *
   * @param other The BUFR I/O backend to move from
   * @return Reference to this backend after assignment
   */
  BufrObsIO& operator=(BufrObsIO&& other) noexcept {
    if (this != &other) {
      config_ = std::move(other.config_);
      filename_ = std::move(other.filename_);
      bufrWrapper_ = std::move(other.bufrWrapper_);
    }
    return *this;
  }

  /**
   * @brief Destructor
   */
  ~BufrObsIO() = default;

  /**
   * @brief Read observations from the configured data source
   *
   * @details Parses the BUFR data source and extracts observation records.
   * Uses the Fortran BUFR API to decode the binary format.
   *
   * @return Vector of observation records read from the data source
   * @throws std::runtime_error If the data cannot be read or parsed
   */
  std::vector<ObsRecord> read() {
    std::vector<ObsRecord> records;
    records.reserve(100);  // Pre-allocate space to avoid reallocations

    try {
      openBufrFile('r');
      records = readBufrRecords();
      closeBufrFile();
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to read BUFR data: " +
                               std::string(e.what()));
    }

    return records;
  }

  /**
   * @brief Write observations to the configured data destination
   *
   * @details Encodes the provided observation records into BUFR format
   * and writes them to the specified destination.
   *
   * @param records Vector of observation records to write
   * @throws std::runtime_error If writing fails
   *
   * @note Current implementation is a placeholder. Actual BUFR encoding
   * requires implementation of the appropriate BUFR encoding API calls.
   */
  void write(const std::vector<ObsRecord>& records) {
    try {
      openBufrFile('w');
      writeBufrRecords(records);
      closeBufrFile();
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to write BUFR data: " +
                               std::string(e.what()));
    }
  }

 private:
  /**
   * @brief Open the BUFR file for reading or writing
   *
   * @param mode File access mode: 'r' for reading, 'w' for writing
   */
  void openBufrFile(char mode) { bufrWrapper_->open(filename_, mode); }

  /**
   * @brief Close the BUFR file
   */
  void closeBufrFile() { bufrWrapper_->close(); }

  /**
   * @brief Read and process all records from a BUFR file
   *
   * @return Vector of observation records
   */
  std::vector<ObsRecord> readBufrRecords() {
    std::vector<ObsRecord> records;

    // Read subsets until the end of file is reached
    std::string subset;
    int date;
    int ret;
    int nlev;

    // Reuse the same ObsRecord to avoid constructor calls
    ObsRecord record;
    char stid[9] = {0};  // 8 chars + null terminator

    // Header values array for station information
    double header[NHR8PM] = {0.0};

    // Allocate an array for events data
    const int evnsSize = MXR8PM * MXR8LV * MXR8VN * MXR8VT;
    std::unique_ptr<double[]> events(new double[evnsSize]());

    // Process BUFR file using readpb
    while ((ret = bufrWrapper_->readPrepbufr(subset, date, header, events.get(),
                                             &nlev)) >= 0) {
      // Populate base record information
      populateBaseRecord(record, subset, date, header, stid);

      // Process each level of data
      processLevels(records, record, events.get(), nlev);

      // If this is the last subset, break
      if (ret == 1) break;
    }

    return records;
  }

  /**
   * @brief Populate basic record information from BUFR header
   *
   * @param record The record to populate
   * @param subset The BUFR subset identifier
   * @param date The observation date
   * @param header The BUFR header array
   * @param stid Buffer to store station ID
   */
  void populateBaseRecord(ObsRecord& record, const std::string& subset,
                          int date, const double header[], char* stid) {
    // Basic observation info
    record.type = subset;

    // Station information (HDR indices 1-4)
    std::memcpy(stid, &header[0], 8);
    record.station_id = stid;
    record.longitude = header[1];  // XOB
    record.latitude = header[2];   // YOB
    record.elevation = header[3];  // ELV

    // Time information - use DateTime constructor for date integer
    record.datetime = DateTime(date);

    // If there's a time offset, adjust the datetime using Duration
    if (header[4] != 0.0) {
      record.datetime += Duration::fromHoursF(header[4]);
    }

    // Report metadata (HDR indices 6-8)
    record.report_type = std::to_string(static_cast<int>(header[5]));
    record.input_report_type = std::to_string(static_cast<int>(header[6]));
    record.instrument_type = std::to_string(static_cast<int>(header[7]));
  }

  /**
   * @brief Process all observation levels from BUFR data
   *
   * @param records Vector to store the processed records
   * @param baseRecord Base record with common information
   * @param events Events data array
   * @param nlev Number of levels in the data
   */
  void processLevels(std::vector<ObsRecord>& records,
                     const ObsRecord& baseRecord, double* events, int nlev) {
    // Process each level of data
    for (int lv = 1; lv <= nlev; ++lv) {
      processObservationType(records, baseRecord, events, lv, 1, "PRES");
      processObservationType(records, baseRecord, events, lv, 2, "HUMID");
      processObservationType(records, baseRecord, events, lv, 3, "TEMP");
      processObservationType(records, baseRecord, events, lv, 4, "HEIGHT");
      processObservationType(records, baseRecord, events, lv, 5, "UWIND");
      processObservationType(records, baseRecord, events, lv, 6, "VWIND");
    }
  }

  /**
   * @brief Process a specific observation type from BUFR data
   *
   * @param records Vector to store the processed records
   * @param baseRecord Base record with common information
   * @param events Events data array
   * @param lv Level index
   * @param kk Variable type index
   * @param typeName Name of the observation type
   */
  void processObservationType(std::vector<ObsRecord>& records,
                              const ObsRecord& baseRecord, double* events,
                              int lv, int kk, const std::string& typeName) {
    double value = BufrFortranWrapper::getEvnsValue(events, 1, lv, 1, kk);
    if (isValidValue(value)) {
      ObsRecord levelRecord = baseRecord;
      levelRecord.type = typeName;
      levelRecord.value = value;
      double qm = BufrFortranWrapper::getEvnsValue(events, 2, lv, 1, kk);
      levelRecord.qc_marker = static_cast<std::size_t>(qm);
      records.emplace_back(std::move(levelRecord));
    }
  }

  /**
   * @brief Write records to a BUFR file
   *
   * @param records Records to write
   *
   * @note Current implementation is a placeholder
   */
  void writeBufrRecords(const std::vector<ObsRecord>& records) {
    // In a real implementation, would encode each record into BUFR format
    // using the appropriate BUFR API functions
  }

  ConfigBackend config_;
  std::string filename_;  // Store the filename from config
  std::unique_ptr<BufrFortranWrapper> bufrWrapper_;
};

/**
 * @brief Static assertion to verify the backend meets the concept requirements
 *
 * @tparam T Configuration backend type
 */
template <typename T>
struct BufrObsIOConceptCheck {
  static_assert(framework::ObsIOBackendImpl<BufrObsIO<T>, T>,
                "BufrObsIO must satisfy the ObsIOBackendImpl concept");
};

}  // namespace metada::backends::io