#pragma once

#include <ng-log/logging.h>

#include "common/utils/logger/LogStream.hpp"

namespace metada::backends::logger {

using framework::LogLevel;

/**
 * @brief ng-log logging backend implementation
 *
 * @details
 * This class provides a backend implementation of the logger interface
 * using ng-log. It maps the framework's logging levels to ng-log severity
 * levels and provides a clean integration with the logging system.
 *
 * Features:
 * - Full integration with ng-log
 * - Support for all framework logging levels
 * - Automatic log file creation and rotation (if supported by ng-log)
 * - Log levels mapped to ng-log severity levels
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
 * appropriate ng-log verbosity flags through the configuration.
 *
 * @tparam ConfigBackend The configuration backend type that satisfies the
 * ConfigBackendType concept
 */
template <typename ConfigBackend>
class NgLogger {
 public:
  NgLogger() = delete;
  ~NgLogger() { Shutdown(); }
  NgLogger(const NgLogger&) = delete;
  NgLogger& operator=(const NgLogger&) = delete;
  NgLogger(NgLogger&&) noexcept = default;
  NgLogger& operator=(NgLogger&&) noexcept = default;

  explicit NgLogger(const ConfigBackend& config) : config_(config) {
    Init(config_.Get("app_name").asString());
    // TODO: Set ng-log specific flags from config if needed
  }

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
        VLOG(1) << message;
        break;
    }
  }

  static void Init(const std::string& app_name) {
    nglog::InitializeLogging(app_name.c_str());
    LOG(INFO) << "ng-log logger initialized for " << app_name;
    // TODO: Set ng-log specific flags from config if needed
  }

  static void Shutdown() { nglog::ShutdownLogging(); }

 private:
  const ConfigBackend& config_;
};

}  // namespace metada::backends::logger