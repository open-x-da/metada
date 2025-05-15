#pragma once

#include <cstring>
#include <fstream>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>

#include "BufrFortranAPI.h"

namespace metada::backends::io {

// Constants for BUFR array dimensions from readpb.prm
constexpr int MXR8PM = 10;   // Number of event data types
constexpr int MXR8LV = 400;  // Maximum number of levels
constexpr int MXR8VN = 10;   // Maximum number of event stacks
constexpr int MXR8VT = 6;    // Number of variable types (P,Q,T,Z,U,V)
constexpr int NHR8PM = 8;    // Number of header elements

constexpr double R8BFMS = 10.0E10;  // Missing value for real*8

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
  BufrFortranWrapper() {
    // Initialize the subset buffer with spaces (Fortran-friendly)
    memset(subsetBuffer_, ' ', 8);
    subsetBuffer_[8] = '\0';

    // Initialize the temporary header buffer to zeros
    memset(tempHeader_, 0, sizeof(tempHeader_));

    // Initialize the events buffer and level count
    evnsSize_ = MXR8PM * MXR8LV * MXR8VN * MXR8VT;
    tempEvns_ = new double[evnsSize_]();
    tempNlev_ = 0;
  }

  /**
   * @brief Destructor - ensures all files are closed
   */
  ~BufrFortranWrapper() {
    if (isOpen_) {
      close_bufr_file_(unitNumber_);
    }

    // Clean up the events buffer
    delete[] tempEvns_;
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
   * @brief Read and process the next station report from a PREPBUFR file
   *
   * @param subset Output parameter that will contain the subset name
   * @param date Output parameter that will contain the date
   * @param header Output parameter that will contain the header data
   * @param events Output parameter that will contain event data (can be
   * nullptr)
   * @param nlev Output parameter that will contain the number of levels (can be
   * nullptr)
   * @return int 0 if OK, 1 if last subset, -1 if EOF
   */
  int readPrepbufr(std::string& subset, int& date, double* header = nullptr,
                   double* events = nullptr, int* nlev = nullptr) {
    if (!isOpen_) {
      throw std::runtime_error("BUFR file not open");
    }

    int idate = 0;
    int iret = 0;
    int inlev = 0;

    // Use the provided header or the class member
    double* headerPtr = header ? header : tempHeader_;

    // Use the provided events array or the class member
    double* eventsPtr = events ? events : tempEvns_;

    // Use the provided nlev or the class member
    int* nlevPtr = nlev ? nlev : &tempNlev_;

    // Call the Fortran function
    readpb_(&unitNumber_, subsetBuffer_, &idate, headerPtr, eventsPtr, nlevPtr,
            &iret, 8);

    // Copy results to output parameters
    subset = std::string(subsetBuffer_);
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

  /**
   * @brief Get the event data for a specific variable, level, and event type
   *
   * @param events The events array returned from readPrepbufr
   * @param ii Event type index (1-8: OBS, QM, PGM, RSN, FCT, ANL, ERR, CAT)
   * @param lv Level index (1-based)
   * @param jj Event stack index (1-based)
   * @param kk Variable type index (1-6: P, Q, T, Z, U, V)
   * @return The event data value
   */
  static double getEvnsValue(const double* events, int ii, int lv, int jj,
                             int kk) {
    int idx = (ii - 1) + (lv - 1) * MXR8PM + (jj - 1) * MXR8PM * MXR8LV +
              (kk - 1) * MXR8PM * MXR8LV * MXR8VN;
    return events[idx];
  }

 private:
  int unitNumber_ = -1;
  int tableUnit_ = -1;
  bool isOpen_ = false;
  std::string filename_;
  char subsetBuffer_[9];        // Persistent buffer for Fortran string exchange
  double tempHeader_[NHR8PM];   // Persistent buffer for header data
  double* tempEvns_ = nullptr;  // Persistent buffer for event data
  int evnsSize_ = 0;            // Size of the events array
  int tempNlev_ = 0;            // Number of levels in the current report
};

}  // namespace metada::backends::io