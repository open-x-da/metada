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
#include <optional>
#include <string>
#include <vector>

#include "BufrFortranWrapper.hpp"
#include "DateTime.hpp"
#include "Duration.hpp"
#include "ObsIOConcepts.hpp"

namespace metada::backends::io {

using ObsRecord = framework::ObsRecord;

// Helper function for comparing doubles with BUFR missing value
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
    // Parse and store configuration parameters
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
    // Create a vector to hold the observation records
    std::vector<ObsRecord> records;
    records.reserve(100);  // Pre-allocate space to avoid reallocations

    try {
      // Open the BUFR file using the Fortran wrapper
      bufrWrapper_->open(filename_, 'r');

      // Read subsets until the end of file is reached
      std::string subset;
      int date;
      int ret;
      int nlev;

      // Reuse the same ObsRecord to avoid constructor calls
      ObsRecord record;
      char stid[9] = {0};  // 8 chars + null terminator

      // Header values array for station information
      constexpr int NHR8PM = 8;  // NHR8PM from readpb.prm
      double header[NHR8PM] = {0.0};

      // Allocate an array for events data
      const int evnsSize = MXR8PM * MXR8LV * MXR8VN * MXR8VT;
      std::unique_ptr<double[]> events(new double[evnsSize]());

      // Process BUFR file using readpb
      while ((ret = bufrWrapper_->readPrepbufr(subset, date, header,
                                               events.get(), &nlev)) >= 0) {
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

        // Process each level of data
        for (int lv = 1; lv <= nlev; ++lv) {
          // Get observation values (first event, ii=1) for this level
          // kk=1: Pressure, kk=2: Humidity, kk=3: Temperature, kk=4: Height,
          // kk=5: U-wind, kk=6: V-wind
          double pressure =
              BufrFortranWrapper::getEvnsValue(events.get(), 1, lv, 1, 1);
          if (isValidValue(pressure)) {
            ObsRecord levelRecord = record;
            levelRecord.type = "PRES";
            levelRecord.value = pressure;
            double pqm =
                BufrFortranWrapper::getEvnsValue(events.get(), 2, lv, 1, 1);
            levelRecord.qc_marker = static_cast<std::size_t>(pqm);
            records.emplace_back(std::move(levelRecord));
            continue;
          }
          double humidity =
              BufrFortranWrapper::getEvnsValue(events.get(), 1, lv, 1, 2);
          if (isValidValue(humidity)) {
            ObsRecord levelRecord = record;
            levelRecord.type = "HUMID";
            levelRecord.value = humidity;
            double qqm =
                BufrFortranWrapper::getEvnsValue(events.get(), 2, lv, 1, 2);
            levelRecord.qc_marker = static_cast<std::size_t>(qqm);
            records.emplace_back(std::move(levelRecord));
            continue;
          }
          double temp =
              BufrFortranWrapper::getEvnsValue(events.get(), 1, lv, 1, 3);
          if (isValidValue(temp)) {
            ObsRecord levelRecord = record;
            levelRecord.type = "TEMP";
            levelRecord.value = temp;
            double tqm =
                BufrFortranWrapper::getEvnsValue(events.get(), 2, lv, 1, 3);
            levelRecord.qc_marker = static_cast<std::size_t>(tqm);
            records.emplace_back(std::move(levelRecord));
            continue;
          }
          double height =
              BufrFortranWrapper::getEvnsValue(events.get(), 1, lv, 1, 4);
          if (isValidValue(height)) {
            ObsRecord levelRecord = record;
            levelRecord.type = "HEIGHT";
            levelRecord.value = height;
            double zqm =
                BufrFortranWrapper::getEvnsValue(events.get(), 2, lv, 1, 4);
            levelRecord.qc_marker = static_cast<std::size_t>(zqm);
            records.emplace_back(std::move(levelRecord));
            continue;
          }
          double uwind =
              BufrFortranWrapper::getEvnsValue(events.get(), 1, lv, 1, 5);
          if (isValidValue(uwind)) {
            ObsRecord levelRecord = record;
            levelRecord.type = "UWIND";
            levelRecord.value = uwind;
            double wqm =
                BufrFortranWrapper::getEvnsValue(events.get(), 2, lv, 1, 5);
            levelRecord.qc_marker = static_cast<std::size_t>(wqm);
            records.emplace_back(std::move(levelRecord));
            continue;
          }
          double vwind =
              BufrFortranWrapper::getEvnsValue(events.get(), 1, lv, 1, 6);
          if (isValidValue(vwind)) {
            ObsRecord levelRecord = record;
            levelRecord.type = "VWIND";
            levelRecord.value = vwind;
            double wqm =
                BufrFortranWrapper::getEvnsValue(events.get(), 2, lv, 1, 6);
            levelRecord.qc_marker = static_cast<std::size_t>(wqm);
            records.emplace_back(std::move(levelRecord));
            continue;
          }
        }

        // If this is the last subset, break
        if (ret == 1) break;
      }

      // Close the BUFR file
      bufrWrapper_->close();
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
   */
  void write(const std::vector<ObsRecord>& records) {
    // Placeholder implementation
    try {
      // Open the BUFR file for writing
      bufrWrapper_->open(filename_, 'w');

      // In a real implementation, would encode each record into BUFR format
      // using the appropriate BUFR API functions

      // Close the BUFR file
      bufrWrapper_->close();
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to write BUFR data: " +
                               std::string(e.what()));
    }
  }

 private:
  ConfigBackend config_;
  std::string filename_;  // Store the filename from config
  std::unique_ptr<BufrFortranWrapper> bufrWrapper_;
};

// Static assertion to verify the backend meets the concept requirements
template <typename T>
struct BufrObsIOConceptCheck {
  static_assert(framework::ObsIOBackendImpl<BufrObsIO<T>, T>,
                "BufrObsIO must satisfy the ObsIOBackendImpl concept");
};

}  // namespace metada::backends::io