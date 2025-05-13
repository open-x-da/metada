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

    try {
      // Open the BUFR file using the Fortran wrapper
      bufrWrapper_->open(filename_, 'r');

      // Read subsets until the end of file is reached
      std::string subset;
      int date;
      int ret;

      // Header values array for station information
      constexpr int NHR8PM = 8;  // NHR8PM from readpb.prm
      double header[NHR8PM] = {0.0};

      // Process BUFR file using readpb
      while ((ret = bufrWrapper_->readPrepbufr(subset, date, header)) >= 0) {
        // Convert the subset data to ObsRecord objects
        // Create record with data from BUFR header
        ObsRecord record;

        // Basic observation info
        record.type = subset;
        record.value = 0.0;  // Placeholder

        // Station information (HDR indices 1-4)
        record.station_id =
            std::string(reinterpret_cast<char*>(&header[0])).c_str();
        record.longitude = header[1];  // XOB
        record.latitude = header[2];   // YOB
        record.elevation = header[3];  // ELV

        // Time information - use DateTime constructor for date integer
        record.datetime = DateTime(date);

        // If there's a time offset, adjust the datetime using Duration
        if (header[4] != 0.0) {
          // Convert float DHR (hours) to Duration and adjust datetime
          Duration dhrOffset = Duration::fromHoursF(header[4]);
          record.datetime += dhrOffset;
        }

        // Report metadata (HDR indices 6-8)
        record.report_type =
            std::to_string(static_cast<int>(header[5]));  // TYP
        record.input_report_type =
            std::to_string(static_cast<int>(header[6]));  // T29
        record.instrument_type =
            std::to_string(static_cast<int>(header[7]));  // ITP

        // Quality control
        record.qc_marker = 0;  // Default QC marker

        records.push_back(record);

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