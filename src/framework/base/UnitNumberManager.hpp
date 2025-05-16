#pragma once

#include <mutex>
#include <set>
#include <stdexcept>
#include <string>

namespace metada::framework::base {

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
  int allocate() {
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
  void release(int unit) {
    std::lock_guard<std::mutex> lock(mutex_);
    usedUnits_.erase(unit);
  }

  /**
   * @brief Check if a unit number is in use
   *
   * @param unit The unit number to check
   * @return true if the unit is in use, false otherwise
   */
  bool isInUse(int unit) {
    std::lock_guard<std::mutex> lock(mutex_);
    return usedUnits_.find(unit) != usedUnits_.end();
  }

  /**
   * @brief Manually register a unit number
   *
   * @param unit The unit number to register
   */
  void manuallyRegister(int unit) {
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

}  // namespace metada::framework::base