#pragma once

#include <glog/logging.h>

#include "common/utils/logger/LogStream.hpp"

namespace metada::backends::logger {

using framework::LogLevel;

/**
 * @brief Google logging backend implementation
 *
 * @details
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
 * - Configurable through configuration backend
 *
 * @par Example usage with framework Logger:
 * @code
 * // Create application context with configuration
 * auto context = ApplicationContext<BackendTag>(argv[0], argv[1]);
 * auto& logger = context.getLogger();
 *
 * // Log at different levels
 * logger.Info() << "Application started with version " << version;
 * logger.Warning() << "Resource usage high: " << usage_percent << "%";
 * logger.Error() << "Failed to connect to " << host << ":" << port;
 * logger.Debug() << "Connection details: " << details;
 * @endcode
 *
 * @note For debug logging to work in release builds, you may need to set
 * appropriate glog verbosity flags through the configuration.
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
    Init(config_.Get("app_name").asString());
    FLAGS_colorlogtostderr = config_.Get("color").asBool();
    FLAGS_logtostderr = config_.Get("console").asBool();
    auto level_str = config_.Get("level").asString();
    std::transform(level_str.begin(), level_str.end(), level_str.begin(),
                   ::tolower);
    if (level_str == "info") {
      FLAGS_minloglevel = 0;  // INFO
    } else if (level_str == "warning") {
      FLAGS_minloglevel = 1;  // WARNING
    } else if (level_str == "error") {
      FLAGS_minloglevel = 2;  // ERROR
    } else if (level_str == "debug") {
      FLAGS_v = 1;            // Enable VLOG level 1
      FLAGS_minloglevel = 0;  // INFO
    } else {
      FLAGS_minloglevel = 0;  // INFO
    }
  }

  /**
   * @brief Log a message at the specified level
   *
   * @details
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
        VLOG(1) << message;  // Using VLOG for debug level
        break;
    }
  }

  /**
   * @brief Initialize glog for the application
   *
   * @details
   * Sets up Google's logging library with the application name,
   * configuring log file names and initial settings.
   *
   * @param app_name Name of the application for log file identification
   */
  static void Init(const std::string& app_name) {
    // Initialize glog with application name
    google::InitGoogleLogging(app_name.c_str());

    // Log initialization message
    LOG(INFO) << "Google logger initialized for " << app_name;
  }

  /**
   * @brief Shutdown glog
   *
   * @details
   * Performs proper cleanup of Google's logging library
   */
  static void Shutdown() { google::ShutdownGoogleLogging(); }

 private:
  const ConfigBackend& config_;  ///< Configuration backend instance
};

}  // namespace metada::backends::logger