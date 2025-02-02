#ifndef METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGER_H_
#define METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGER_H_

#include <glog/logging.h>

#include <string>

namespace metada {
namespace logger {

class GoogleLogger {
 public:
  static void Init(const std::string& app_name) {
    google::InitGoogleLogging(app_name.c_str());
    google::SetLogDestination(google::GLOG_INFO, "logs/METADA_INFO_");
    google::SetLogDestination(google::GLOG_WARNING, "logs/METADA_WARNING_");
    google::SetLogDestination(google::GLOG_ERROR, "logs/METADA_ERROR_");
  }

  static void Shutdown() { google::ShutdownGoogleLogging(); }

  static void LogInfo(const std::string& message) { LOG(INFO) << message; }

  static void LogWarning(const std::string& message) {
    LOG(WARNING) << message;
  }

  static void LogError(const std::string& message) { LOG(ERROR) << message; }

  static void LogDebug(const std::string& message) {
    LOG(INFO) << "[DEBUG] " << message;
  }
};

}  // namespace logger
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGER_H_
