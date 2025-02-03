#ifndef METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGER_H_
#define METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGER_H_

#include <glog/logging.h>

#include <string>

#include "ILogger.h"

namespace metada {
namespace backends {
namespace tools {
namespace logger {
namespace google_logger {

/**
 * @brief Google logging implementation of ILogger interface
 *
 * Provides logging functionality by delegating to Google's glog library.
 */
class GoogleLogger : public metada::framework::tools::logger::ILogger {
 public:
  /**
   * @brief Log info message using glog
   * @param message Message to log
   */
  void Info(const std::string& message) override { LOG(INFO) << message; }

  /**
   * @brief Log warning message using glog
   * @param message Message to log
   */
  void Warning(const std::string& message) override { LOG(WARNING) << message; }

  /**
   * @brief Log error message using glog
   * @param message Message to log
   */
  void Error(const std::string& message) override { LOG(ERROR) << message; }

  /**
   * @brief Log debug message using glog
   * @param message Message to log
   */
  void Debug(const std::string& message) override { DLOG(INFO) << message; }

  /**
   * @brief Initialize glog for application
   * @param app_name Name of application using logger
   */
  static void Init(const std::string& app_name) {
    google::InitGoogleLogging(app_name.c_str());
  }

  /**
   * @brief Shutdown glog logging
   */
  static void Shutdown() { google::ShutdownGoogleLogging(); }
};

}  // namespace google_logger
}  // namespace logger
}  // namespace tools
}  // namespace backends
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGER_H_
