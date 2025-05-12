#pragma once

#include <fstream>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>

#include "BufrFortranAPI.h"

namespace metada::backends::io {

/**
 * @brief Manages Fortran logical unit numbers to avoid conflicts
 *
 * This class handles allocation and deallocation of Fortran logical unit
 * numbers in a thread-safe manner, preventing conflicts when multiple files are
 * opened.
 */
class UnitNumberManager {
 public:
  /**
   * @brief Get the singleton instance of the unit number manager
   */
  static UnitNumberManager& getInstance() {
    static UnitNumberManager instance;
    return instance;
  }

  /**
   * @brief Allocate a unit number
   *
   * @return int An available unit number
   * @throws std::runtime_error If no unit numbers are available
   */
  int allocateUnit() {
    std::lock_guard<std::mutex> lock(mutex_);

    // Look for an unused unit number
    for (int unit = MIN_UNIT; unit <= MAX_UNIT; ++unit) {
      if (usedUnits_.find(unit) == usedUnits_.end()) {
        usedUnits_.insert(unit);
        return unit;
      }
    }

    throw std::runtime_error("No available Fortran logical unit numbers");
  }

  /**
   * @brief Release a unit number
   *
   * @param unit The unit number to release
   */
  void releaseUnit(int unit) {
    std::lock_guard<std::mutex> lock(mutex_);
    usedUnits_.erase(unit);
  }

  /**
   * @brief Check if a unit number is in use
   *
   * @param unit The unit number to check
   * @return true if the unit is in use, false otherwise
   */
  bool isUnitInUse(int unit) {
    std::lock_guard<std::mutex> lock(mutex_);
    return usedUnits_.find(unit) != usedUnits_.end();
  }

  /**
   * @brief Manually register a unit number
   *
   * @param unit The unit number to register
   */
  void manuallyRegisterUnit(int unit) {
    std::lock_guard<std::mutex> lock(mutex_);
    usedUnits_.insert(unit);
  }

 private:
  // Private constructor to ensure singleton pattern
  UnitNumberManager() = default;

  // Prevent copying
  UnitNumberManager(const UnitNumberManager&) = delete;
  UnitNumberManager& operator=(const UnitNumberManager&) = delete;

  // Typical range of unit numbers (avoid 0, 5, 6 which are often reserved)
  static constexpr int MIN_UNIT = 10;
  static constexpr int MAX_UNIT = 99;

  std::set<int> usedUnits_;  // Set of currently used unit numbers
  std::mutex mutex_;         // For thread safety
};

/**
 * @brief C++ wrapper for the Fortran BUFR API
 *
 * This class provides a C++ interface to the underlying Fortran BUFR library
 * functions, handling string conversions and memory management.
 */
class BufrFortranWrapper {
 public:
  /**
   * @brief Constructor
   */
  BufrFortranWrapper() = default;

  /**
   * @brief Destructor - ensures all files are closed
   */
  ~BufrFortranWrapper() {
    if (isOpen_) {
      close_bufr_file_(unitNumber_);
    }
  }

  /**
   * @brief Open a BUFR file for processing
   *
   * @param filename Path to the BUFR file
   * @param mode I/O mode ('r' for read, 'w' for write)
   * @param requestedUnit Specific unit number to use (if < 0, auto-allocate)
   * @throws std::runtime_error If file opening fails or requested unit is in
   * use
   */
  void open(const std::string& filename, char mode = 'r',
            int requestedUnit = -1) {
    // Close any previously open file
    if (isOpen_) {
      close_bufr_file_(unitNumber_);
    }

    // Check if we need to allocate or use a specific unit number
    if (requestedUnit < 0) {
      // Auto-allocate a unit number
      unitNumber_ = UnitNumberManager::getInstance().allocateUnit();
    } else {
      // Use the requested unit number if available
      if (UnitNumberManager::getInstance().isUnitInUse(requestedUnit)) {
        throw std::runtime_error("Unit number " +
                                 std::to_string(requestedUnit) +
                                 " is already in use");
      }
      // Register the requested unit number
      UnitNumberManager::getInstance().manuallyRegisterUnit(requestedUnit);
      unitNumber_ = requestedUnit;
    }

    tableUnit_ = unitNumber_;  // Use same unit for table

    // Use our Fortran wrapper to open the file
    int status = 0;
    open_bufr_file_(filename.c_str(), unitNumber_, &status);

    if (status != 0) {
      // Release the allocated unit number if opening fails
      UnitNumberManager::getInstance().releaseUnit(unitNumber_);
      throw std::runtime_error("Cannot open BUFR file: " + filename);
    }

    isOpen_ = true;
    filename_ = filename;
  }

  /**
   * @brief Read the next subset from a BUFR file
   *
   * @param subset Output parameter that will contain the subset name
   * @param date Output parameter that will contain the date
   * @return true if a subset was read, false if end of file
   */
  bool readNextSubset(std::string& subset, int& date) {
    if (!isOpen_) {
      throw std::runtime_error("BUFR file not open");
    }

    // Allocate buffers for Fortran
    char subsetBuffer[9] = {0};  // Typically 8 chars plus null terminator
    int idate = 0;
    int iret = 0;

    // Call the Fortran function
    readpb_(&unitNumber_, subsetBuffer, &idate, &iret,
            sizeof(subsetBuffer) - 1);

    // Check for EOF
    if (iret != 0) {
      return false;
    }

    // Copy results to output parameters
    subset = std::string(subsetBuffer);
    date = idate;

    return true;
  }

  /**
   * @brief Read and process the next station report from a PREPBUFR file
   *
   * @param subset Output parameter that will contain the subset name
   * @param date Output parameter that will contain the date
   * @return int 0 if OK, 1 if last subset, -1 if EOF
   */
  int readPrepbufr(std::string& subset, int& date) {
    if (!isOpen_) {
      throw std::runtime_error("BUFR file not open");
    }

    // Allocate buffers for Fortran
    char subsetBuffer[9] = {0};  // Typically 8 chars plus null terminator
    int idate = 0;
    int iret = 0;

    // Call the Fortran function
    readpb_(&unitNumber_, subsetBuffer, &idate, &iret,
            sizeof(subsetBuffer) - 1);

    // Copy results to output parameters
    subset = std::string(subsetBuffer);
    date = idate;

    return iret;
  }

  /**
   * @brief Close the BUFR file
   */
  void close() {
    if (isOpen_) {
      // Use the Fortran close function
      close_bufr_file_(unitNumber_);

      // Release the unit number back to the pool
      UnitNumberManager::getInstance().releaseUnit(unitNumber_);
      unitNumber_ = -1;
      tableUnit_ = -1;

      isOpen_ = false;
      filename_.clear();
    }
  }

  /**
   * @brief Get the unit number currently in use
   *
   * @return int The current unit number, or -1 if no file is open
   */
  int getUnitNumber() const { return unitNumber_; }

 private:
  int unitNumber_ = -1;
  int tableUnit_ = -1;
  bool isOpen_ = false;
  std::string filename_;
};

}  // namespace metada::backends::io