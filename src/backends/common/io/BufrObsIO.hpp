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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "BufrFortranAPI.h"
#include "DateTime.hpp"
#include "Duration.hpp"
#include "ObsIOConcepts.hpp"
#include "UnitNumberManager.hpp"
#include "isValidValue.hpp"

namespace metada::backends::io {

using ObsRecord = framework::ObsRecord;
using ObsRecordShared = framework::ObsRecordShared;
using ObsLevelRecord = framework::ObsLevelRecord;

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
    filename_ = config_.Get("filename").asString();
    // Load station_ids from config if present
    if (config_.HasKey("station_ids")) {
      auto station_ids_val = config_.Get("station_ids");
      if (station_ids_val.isVectorString()) {
        station_ids_ = station_ids_val.asVectorString();
      } else if (station_ids_val.isVectorInt()) {
        const auto& int_ids = station_ids_val.asVectorInt();
        station_ids_.clear();
        station_ids_.reserve(int_ids.size());
        for (const auto& id : int_ids) {
          station_ids_.emplace_back(std::to_string(id));
        }
      } else {
        throw std::runtime_error(
            "station_ids must be a vector of strings or ints");
      }
    }
    // Load data_types (subsets) from config if present
    if (config_.HasKey("data_types")) {
      auto data_types_val = config_.Get("data_types");
      if (data_types_val.isVectorString()) {
        data_types_ = data_types_val.asVectorString();
      }
    }
    // Initialize Fortran buffer fields
    memset(subsetBuffer_, ' ', 8);
    subsetBuffer_[8] = '\0';
    memset(tempHeader_, 0, sizeof(tempHeader_));
    evnsSize_ = MXR8PM * MXR8LV * MXR8VN * MXR8VT;
    tempEvns_ = new double[evnsSize_]();
    tempNlev_ = 0;
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
        unitNumber_(other.unitNumber_),
        tableUnit_(other.tableUnit_),
        isOpen_(other.isOpen_),
        tempEvns_(other.tempEvns_),
        evnsSize_(other.evnsSize_),
        tempNlev_(other.tempNlev_),
        station_ids_(std::move(other.station_ids_)),
        data_types_(std::move(other.data_types_)) {
    std::memcpy(subsetBuffer_, other.subsetBuffer_, sizeof(subsetBuffer_));
    std::memcpy(tempHeader_, other.tempHeader_, sizeof(tempHeader_));
  }

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
      unitNumber_ = other.unitNumber_;
      tableUnit_ = other.tableUnit_;
      isOpen_ = other.isOpen_;
      std::memcpy(subsetBuffer_, other.subsetBuffer_, sizeof(subsetBuffer_));
      std::memcpy(tempHeader_, other.tempHeader_, sizeof(tempHeader_));
      tempEvns_ = other.tempEvns_;
      evnsSize_ = other.evnsSize_;
      tempNlev_ = other.tempNlev_;
      station_ids_ = std::move(other.station_ids_);
      data_types_ = std::move(other.data_types_);
    }
    return *this;
  }

  /**
   * @brief Destructor
   */
  ~BufrObsIO() {
    if (isOpen_) {
      close_bufr_file_(unitNumber_);
    }
    delete[] tempEvns_;
  }

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
      open(filename_);
      records = readBufrRecords();
      close();
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
      open(filename_);
      writeBufrRecords(records);
      close();
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to write BUFR data: " +
                               std::string(e.what()));
    }
  }

 private:
  // Constants for BUFR array dimensions from readpb.prm
  static constexpr int MXR8PM = 10;   // Number of event data types
  static constexpr int MXR8LV = 400;  // Maximum number of levels
  static constexpr int MXR8VN = 10;   // Maximum number of event stacks
  static constexpr int MXR8VT = 6;    // Number of variable types (P,Q,T,Z,U,V)
  static constexpr int NHR8PM = 8;    // Number of header elements

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

    // Process BUFR file using readpb
    while ((ret = readPrepbufr(subset, date, header, tempEvns_, &nlev)) >= 0) {
      // Filter by data_types (subset) if not empty
      if (!data_types_.empty() &&
          std::find(data_types_.begin(), data_types_.end(), trim(subset)) ==
              data_types_.end()) {
        if (ret == 1) break;
        continue;
      }
      // Populate base record information
      populateBaseRecord(record.shared, date, header, stid);

      // Filter by station_id if station_ids_ is not empty
      if (!station_ids_.empty() &&
          std::find(station_ids_.begin(), station_ids_.end(),
                    trim(record.shared.station_id)) == station_ids_.end()) {
        if (ret == 1) break;
        continue;
      }

      // Process each level of data
      record.levels.clear();
      processLevels(record.levels, nlev);

      // Add the record to our collection
      records.push_back(record);

      // If this is the last subset, break
      if (ret == 1) break;
    }

    return records;
  }

  /**
   * @brief Populate basic record information from BUFR header
   *
   * @param shared The shared record to populate
   * @param date The observation date
   * @param header The BUFR header array
   * @param stid Buffer to store station ID
   */
  void populateBaseRecord(ObsRecordShared& shared, int date,
                          const double header[], char* stid) {
    std::memcpy(stid, &header[0], 8);
    shared.station_id = stid;
    shared.longitude = header[1];  // XOB
    shared.latitude = header[2];   // YOB
    shared.elevation = header[3];  // ELV
    shared.datetime = DateTime(date);
    if (header[4] != 0.0) {
      shared.datetime += Duration::fromHoursF(header[4]);
    }
    shared.report_type = std::to_string(static_cast<int>(header[5]));
    shared.input_report_type = std::to_string(static_cast<int>(header[6]));
    shared.instrument_type = std::to_string(static_cast<int>(header[7]));
  }

  /**
   * @brief Process all observation levels from BUFR data
   *
   * @param records Vector to store the processed records
   * @param nlev Number of levels in the data
   */
  void processLevels(std::vector<std::vector<ObsLevelRecord>>& levels,
                     int nlev) {
    static const char* types[] = {"PRES",   "HUMID", "TEMP",
                                  "HEIGHT", "UWIND", "VWIND"};
    for (int lv = 1; lv <= nlev; ++lv) {
      std::vector<ObsLevelRecord> level_records;
      for (int kk = 1; kk <= 6;
           ++kk) {  // Iterate through variable types (P, Q, T, Z, U, V)
        double value = getEvnsValue(tempEvns_, 1, lv, 1, kk);
        if (isValidValue(value)) {
          // Create a new ObsLevelRecord for each valid variable at this level
          ObsLevelRecord level_record;
          level_record.type = types[kk - 1];
          level_record.value = value;
          double qm = getEvnsValue(tempEvns_, 2, lv, 1, kk);
          level_record.qc_marker = static_cast<std::size_t>(qm);
          // Add this variable's data to the current level's records
          level_records.push_back(std::move(level_record));
        }
      }
      // Only add the level records if we have data
      if (!level_records.empty()) {
        levels.push_back(std::move(level_records));
      }
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
  int unitNumber_ = -1;
  int tableUnit_ = -1;
  bool isOpen_ = false;
  std::string filename_;
  char subsetBuffer_[9];
  double tempHeader_[NHR8PM];
  double* tempEvns_ = nullptr;
  int evnsSize_ = 0;
  int tempNlev_ = 0;
  std::vector<std::string> station_ids_;
  std::vector<std::string> data_types_;  // Store allowed data_types (subsets)

  void open(const std::string& filename) {
    if (isOpen_) {
      close_bufr_file_(unitNumber_);
    }
    unitNumber_ = UnitNumberManager::getInstance().allocate();
    tableUnit_ = unitNumber_;
    int status = 0;
    open_bufr_file_(filename.c_str(), unitNumber_, &status);
    if (status != 0) {
      UnitNumberManager::getInstance().release(unitNumber_);
      throw std::runtime_error("Cannot open BUFR file: " + filename);
    }
    isOpen_ = true;
    filename_ = filename;
  }

  int readPrepbufr(std::string& subset, int& date, double* header = nullptr,
                   double* events = nullptr, int* nlev = nullptr) {
    if (!isOpen_) {
      throw std::runtime_error("BUFR file not open");
    }
    int idate = 0;
    int iret = 0;
    double* headerPtr = header ? header : tempHeader_;
    double* eventsPtr = events ? events : tempEvns_;
    int* nlevPtr = nlev ? nlev : &tempNlev_;
    readpb_(&unitNumber_, subsetBuffer_, &idate, headerPtr, eventsPtr, nlevPtr,
            &iret, 8);
    subset = std::string(subsetBuffer_);
    date = idate;
    return iret;
  }

  void close() {
    if (isOpen_) {
      close_bufr_file_(unitNumber_);
      UnitNumberManager::getInstance().release(unitNumber_);
      unitNumber_ = -1;
      tableUnit_ = -1;
      isOpen_ = false;
      filename_.clear();
    }
  }

  static double getEvnsValue(const double* events, int ii, int lv, int jj,
                             int kk) {
    int idx = (ii - 1) + (lv - 1) * MXR8PM + (jj - 1) * MXR8PM * MXR8LV +
              (kk - 1) * MXR8PM * MXR8LV * MXR8VN;
    return events[idx];
  }

  static std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\n\r");
    auto end = s.find_last_not_of(" \t\n\r");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
  }
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