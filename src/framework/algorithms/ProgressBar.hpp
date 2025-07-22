#pragma once

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace metada::framework {

/**
 * @brief Configuration struct for progress bar appearance
 */
struct ProgressBarConfig {
  int bar_width = 40;                 ///< Width of the progress bar
  std::string fill_char = "█";        ///< Character for filled portion
  std::string empty_char = " ";       ///< Character for empty portion
  std::string prefix = "Progress: ";  ///< Prefix text
  bool show_percentage = true;        ///< Whether to show percentage
  bool show_count = true;             ///< Whether to show current/total count
  int precision = 1;                  ///< Decimal precision for percentage
};

/**
 * @brief Class for creating and managing progress bars
 *
 * This class provides static methods for creating visual progress bars
 * that can be used in any part of the application. The progress bars
 * are similar to those found in download/installation tools.
 *
 * Example usage:
 * ```cpp
 * ProgressBar::Config config;
 * config.prefix = "LETKF progress: ";
 *
 * for (size_t i = 0; i <= total; ++i) {
 *   std::string bar = ProgressBar::Create(i, total, config);
 *   logger.Info() << bar;
 * }
 * ```
 */
class ProgressBar {
 public:
  using Config = ProgressBarConfig;

  /**
   * @brief Create a progress bar string
   * @param current Current progress value
   * @param total Total value (100% completion)
   * @param config Configuration for the progress bar appearance
   * @return Formatted progress bar string
   */
  static std::string Create(size_t current, size_t total,
                            const Config& config = Config{}) {
    if (total == 0) return config.prefix + "[ERROR: total is zero]";

    double percentage = (static_cast<double>(current) / total) * 100.0;

    // Ensure percentage doesn't exceed 100%
    percentage = std::min(percentage, 100.0);

    // Create visual progress bar
    int filled_width =
        static_cast<int>((percentage / 100.0) * config.bar_width);
    std::string progress_bar = "[";

    for (int i = 0; i < config.bar_width; ++i) {
      if (i < filled_width) {
        progress_bar += config.fill_char;
      } else {
        progress_bar += config.empty_char;
      }
    }
    progress_bar += "]";

    // Build final string
    std::stringstream ss;
    ss << config.prefix << progress_bar;

    if (config.show_percentage) {
      ss << " " << std::fixed << std::setprecision(config.precision)
         << percentage << "%";
    }

    if (config.show_count) {
      ss << " (" << current << "/" << total << ")";
    }

    return ss.str();
  }

  /**
   * @brief Create a progress bar with percentage only
   * @param percentage Percentage value (0.0 to 100.0)
   * @param config Configuration for the progress bar appearance
   * @return Formatted progress bar string
   */
  static std::string CreateFromPercentage(double percentage,
                                          const Config& config = Config{}) {
    // Clamp percentage to valid range
    percentage = std::max(0.0, std::min(percentage, 100.0));

    // Create visual progress bar
    int filled_width =
        static_cast<int>((percentage / 100.0) * config.bar_width);
    std::string progress_bar = "[";

    for (int i = 0; i < config.bar_width; ++i) {
      if (i < filled_width) {
        progress_bar += config.fill_char;
      } else {
        progress_bar += config.empty_char;
      }
    }
    progress_bar += "]";

    // Build final string
    std::stringstream ss;
    ss << config.prefix << progress_bar;

    if (config.show_percentage) {
      ss << " " << std::fixed << std::setprecision(config.precision)
         << percentage << "%";
    }

    return ss.str();
  }

  /**
   * @brief Get a simple ASCII-compatible configuration
   * @return Configuration using ASCII characters only
   */
  static Config GetASCIIConfig() {
    Config config;
    config.fill_char = "=";
    config.empty_char = " ";
    return config;
  }

  /**
   * @brief Get a configuration with arrow indicator
   * @return Configuration with arrow-style progress bar
   */
  static Config GetArrowConfig() {
    Config config;
    config.fill_char = "=";
    config.empty_char = " ";
    // Note: Arrow indicator would be added in a more complex implementation
    return config;
  }

  /**
   * @brief Calculate progress interval for logging
   * @param total Total number of items
   * @param target_updates Target number of progress updates (default: 20 for 5%
   * intervals)
   * @return Interval for progress updates
   */
  static size_t CalculateInterval(size_t total, size_t target_updates = 20) {
    return std::max(static_cast<size_t>(1), total / target_updates);
  }

  /**
   * @brief Check if progress should be logged at current point
   * @param current Current progress value
   * @param total Total value
   * @param interval Logging interval
   * @return True if progress should be logged
   */
  static bool ShouldLog(size_t current, size_t total, size_t interval) {
    return (current % interval == 0) || (current == total);
  }
};

}  // namespace metada::framework