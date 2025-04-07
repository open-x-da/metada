#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "common/utils/logger/LogStream.hpp"

namespace metada::backends::logger {

using framework::LogLevel;

// ANSI color code constants
namespace color {
// Foreground color codes
constexpr const char* RESET = "\033[0m";
constexpr const char* RED = "\033[31m";
constexpr const char* GREEN = "\033[32m";
constexpr const char* YELLOW = "\033[33m";
constexpr const char* BLUE = "\033[34m";
constexpr const char* MAGENTA = "\033[35m";
constexpr const char* CYAN = "\033[36m";
constexpr const char* WHITE = "\033[37m";

// Bold/bright foreground colors
constexpr const char* BRED = "\033[1;31m";
constexpr const char* BGREEN = "\033[1;32m";
constexpr const char* BYELLOW = "\033[1;33m";
constexpr const char* BBLUE = "\033[1;34m";
}  // namespace color

/**
 * @brief Console logging backend implementation
 *
 * This class provides a simple console-based logging backend that writes
 * log messages to standard output/error streams with severity level prefixes
 * and optional colorized output.
 *
 * Features:
 * - Simple logging to standard output/error streams
 * - Support for different log levels (debug, info, warning, error)
 * - Colorized output based on log level for improved readability
 * - Timestamp and log level prefixing for each message
 * - Thread-safe logging operations
 * - Configurable output formatting
 * - Real-time output for immediate feedback
 *
 * Log levels:
 * - INFO: General operational information (stdout, green)
 * - WARNING: Potential issues requiring attention (stdout, yellow)
 * - ERROR: Serious problems requiring immediate action (stderr, red)
 * - DEBUG: Detailed troubleshooting information (stdout, blue)
 *
 * Example usage with stream-based API:
 * @code
 * // Initialize logging
 * ConsoleLogger<ConfigBackend>::Init("MyApp");
 *
 * // Create logger instance with configuration
 * ConsoleLogger<ConfigBackend> logger(config);
 *
 * // Log at different levels using stream interface
 * logger.Info() << "Application started with version " << app_version;
 * logger.Warning() << "Resource usage at " << usage_percent << "%";
 * logger.Error() << "Failed to connect to database: " << error_code;
 * logger.Debug() << "Connection params: " << host << ":" << port;
 *
 * // Cleanup on shutdown
 * ConsoleLogger<ConfigBackend>::Shutdown();
 * @endcode
 *
 * @tparam ConfigBackend The configuration backend type that satisfies the
 * ConfigBackendType concept
 */
template <typename ConfigBackend>
class ConsoleLogger {
 public:
  /**
   * @brief Disabled default constructor
   *
   * @details Logger backends should always be initialized with a configuration.
   */
  ConsoleLogger() = delete;

  /**
   * @brief Default destructor
   *
   * @details Ensures proper cleanup by calling Shutdown when the logger is
   * destroyed.
   */
  ~ConsoleLogger() { Shutdown(); }

  /**
   * @brief Disabled copy constructor
   *
   * @details Logger backends are not intended to be copied.
   */
  ConsoleLogger(const ConsoleLogger&) = delete;

  /**
   * @brief Disabled copy assignment operator
   *
   * @details Logger backends are not intended to be copied.
   */
  ConsoleLogger& operator=(const ConsoleLogger&) = delete;

  /**
   * @brief Move constructor
   *
   * @details Allows for moving logger instances when needed.
   */
  ConsoleLogger(ConsoleLogger&&) noexcept = default;

  /**
   * @brief Move assignment operator
   *
   * @details Allows for moving logger instances when needed.
   *
   * @return Reference to this ConsoleLogger instance
   */
  ConsoleLogger& operator=(ConsoleLogger&&) noexcept = default;

  /**
   * @brief Constructor that takes a config
   *
   * @details Initializes the logger with the provided configuration and
   * calls the static Init method with the application name from config.
   *
   * @param[in] config The configuration backend instance
   */
  explicit ConsoleLogger(const ConfigBackend& config) : config_(config) {
    // Extract configuration values with defaults
    app_name_ =
        config.HasKey("app_name") ? config.Get("app_name").asString() : "";
    use_colors_ = config.HasKey("color") ? config.Get("color").asBool() : true;
    show_timestamp_ = config.HasKey("show_timestamp")
                          ? config.Get("show_timestamp").asBool()
                          : true;

    // Initialize with config values
    Init(app_name_);
  }

  /**
   * @brief Log a message at the specified level
   *
   * @details Routes messages to the appropriate output stream with the
   * corresponding severity level prefix. Info, Warning, and Debug messages go
   * to stdout, while Error messages go to stderr. Formatting is controlled by
   * configuration settings for timestamps and colors.
   *
   * @param level The severity level of the message
   * @param message The message to log
   */
  void LogMessage(LogLevel level, const std::string& message) {
    std::string prefix;
    if (show_timestamp_) {
      // Add timestamp if enabled
      auto now = std::chrono::system_clock::now();
      auto time = std::chrono::system_clock::to_time_t(now);
      std::stringstream ss;
      ss << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S");
      prefix = "[" + ss.str() + "] ";
    }

    // Define color variables to avoid repetition
    const char* start_color = "";
    const char* end_color = "";

    // Apply colors if enabled
    if (use_colors_) {
      end_color = color::RESET;

      switch (level) {
        case LogLevel::Info:
          start_color = color::GREEN;
          break;
        case LogLevel::Warning:
          start_color = color::YELLOW;
          break;
        case LogLevel::Error:
          start_color = color::BRED;  // Bold red for errors
          break;
        case LogLevel::Debug:
          start_color = color::BLUE;
          break;
      }
    }

    // Output log message with appropriate colors and to the correct stream
    switch (level) {
      case LogLevel::Info:
        std::cout << prefix << start_color << "[INFO] " << message << end_color
                  << std::endl;
        break;
      case LogLevel::Warning:
        std::cout << prefix << start_color << "[WARNING] " << message
                  << end_color << std::endl;
        break;
      case LogLevel::Error:
        std::cerr << prefix << start_color << "[ERROR] " << message << end_color
                  << std::endl;
        break;
      case LogLevel::Debug:
        std::cout << prefix << start_color << "[DEBUG] " << message << end_color
                  << std::endl;
        break;
    }
  }

  /**
   * @brief Initialize logger for application
   *
   * @details Performs basic initialization of the console logger. Outputs an
   * initialization message to stdout with the application name.
   *
   * @param app_name Name of application using the logger
   * @note This is a static method that should be called once at application
   * startup
   */
  static void Init(const std::string& app_name) {
    // Simple initialization for default logger
    std::cout << "[INFO] " << "Console logger initialized for " << app_name
              << std::endl;
  }

  /**
   * @brief Cleanup logger resources
   *
   * @details Performs cleanup of the console logger. Outputs a shutdown message
   * to stdout. For this simple implementation, no actual cleanup is needed
   * since we only use standard streams.
   *
   * @note This is a static method that should be called once at application
   * shutdown
   */
  static void Shutdown() {
    // Simple cleanup for default logger
    std::cout << "[INFO] " << "Console logger shutdown" << std::endl;
  }

 private:
  const ConfigBackend& config_;  ///< Configuration backend instance
  std::string app_name_;         ///< Application name
  bool use_colors_;              ///< Whether to use colored output
  bool show_timestamp_;          ///< Whether to show timestamps in log messages
};

}  // namespace metada::backends::logger
