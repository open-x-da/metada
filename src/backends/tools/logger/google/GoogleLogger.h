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

class GoogleLogger : public metada::framework::tools::logger::ILogger {
 public:
  void Info(const std::string& message) override { LOG(INFO) << message; }
  void Warning(const std::string& message) override { LOG(WARNING) << message; }
  void Error(const std::string& message) override { LOG(ERROR) << message; }
  void Debug(const std::string& message) override { DLOG(INFO) << message; }

  static void Init(const std::string& app_name) {
    google::InitGoogleLogging(app_name.c_str());
  }

  static void Shutdown() { google::ShutdownGoogleLogging(); }
};

}  // namespace google_logger
}  // namespace logger
}  // namespace tools
}  // namespace backends
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGER_H_
