#pragma once

#include <glog/logging.h>

#include "common/utils/logger/LogStream.hpp"

namespace metada::backends::logger {

using framework::LogLevel;

/**
 * @brief Google logging backend implementation
 *
 * This class provides a backend implementation of the logger interface
 * using Google's glog library. It maps the framework's logging levels to
 * glog severity levels and provides a clean integration with the logging
 * system.
 *
 * Features:
 * - Full integration with Google's logging library
 * - Support for all framework logging levels
 * - Automatic log file creation and rotation
 * - Log levels mapped to glog severity levels
 * - Thread-safe logging operations
 * - Configurable through glog flags and command line
 *
 * Example usage with stream interface:
 * @code
 * // Initialize glog with application name
 * GoogleLogger<ConfigBackend>::Init("MyApp");
 *
 * // Create logger instance with configuration
 * GoogleLogger<ConfigBackend> logger(config);
 *
 * // Log at different levels
 * logger.Info() << "Application started with version " << version;
 * logger.Warning() << "Resource usage high: " << usage_percent << "%";
 * logger.Error() << "Failed to connect to " << host << ":" << port;
 * logger.Debug() << "Connection details: " << details;
 *
 * // Clean up at shutdown
 * GoogleLogger<ConfigBackend>::Shutdown();
 * @endcode
 *
 * @note For debug logging to work in release builds, you may need to set
 * appropriate glog verbosity flags.
 *
 * @tparam ConfigBackend The configuration backend type that satisfies the
 * ConfigBackendType concept
 */
template <typename ConfigBackend>
class GoogleLogger {
 public:
  /**
   * @brief Disabled default constructor
   *
   * @details Logger backends should always be initialized with a configuration.
   */
  GoogleLogger() = delete;

  /**
   * @brief Default destructor
   *
   * @details Ensures proper cleanup by calling Shutdown when the logger is
   * destroyed.
   */
  ~GoogleLogger() { Shutdown(); }

  /**
   * @brief Disabled copy constructor
   *
   * @details Logger backends are not intended to be copied.
   */
  GoogleLogger(const GoogleLogger&) = delete;

  /**
   * @brief Disabled copy assignment operator
   *
   * @details Logger backends are not intended to be copied.
   */
  GoogleLogger& operator=(const GoogleLogger&) = delete;

  /**
   * @brief Move constructor
   *
   * @details Allows for moving logger instances when needed.
   */
  GoogleLogger(GoogleLogger&&) noexcept = default;

  /**
   * @brief Move assignment operator
   *
   * @details Allows for moving logger instances when needed.
   *
   * @return Reference to this GoogleLogger instance
   */
  GoogleLogger& operator=(GoogleLogger&&) noexcept = default;

  /**
   * @brief Constructor that takes a config
   *
   * @details Initializes the logger with the provided configuration and
   * calls the static Init method with a default application name.
   *
   * @param[in] config The configuration backend instance
   */
  explicit GoogleLogger(const ConfigBackend& config) : config_(config) {
    Init("GoogleLogger");
  }

  /**
   * @brief Log a message at the specified level
   *
   * Maps the framework's log levels to glog severity levels and routes
   * the message to the appropriate glog logging macro.
   *
   * @param level The severity level of the message
   * @param message The message to log
   */
  void LogMessage(LogLevel level, const std::string& message) {
    switch (level) {
      case LogLevel::Info:
        LOG(INFO) << message;
        break;
      case LogLevel::Warning:
        LOG(WARNING) << message;
        break;
      case LogLevel::Error:
        LOG(ERROR) << message;
        break;
      case LogLevel::Debug:
        DLOG(INFO) << message;
        break;
    }
  }

  /**
   * @brief Initialize glog for the application
   *
   * Sets up Google's logging library with the application name,
   * configuring log file names and initial settings.
   *
   * @param app_name Name of the application for log file identification
   */
  static void Init(const std::string& app_name) {
    // Initialize glog with application name
    google::InitGoogleLogging(app_name.c_str());

    // Configure standard settings
    FLAGS_logtostderr = 1;  // Log to stderr by default

    // Log initialization message
    LOG(INFO) << "Google logger initialized for " << app_name;
  }

  /**
   * @brief Shutdown glog
   *
   * Performs proper cleanup of Google's logging library
   */
  static void Shutdown() { google::ShutdownGoogleLogging(); }

 private:
  const ConfigBackend& config_;  ///< Configuration backend instance
};

}  // namespace metada::backends::logger