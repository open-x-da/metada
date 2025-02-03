#ifndef METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_
#define METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_

#include <iostream>

#include "ILogger.h"

namespace metada {
namespace backends {
namespace tools {
namespace logger {
namespace default_logger {

class DefaultLogger : public metada::framework::tools::logger::ILogger {
 public:
  void Info(const std::string& message) override {
    std::cout << "[INFO] " << message << std::endl;
  }

  void Warning(const std::string& message) override {
    std::cout << "[WARNING] " << message << std::endl;
  }

  void Error(const std::string& message) override {
    std::cerr << "[ERROR] " << message << std::endl;
  }

  void Debug(const std::string& message) override {
    std::cout << "[DEBUG] " << message << std::endl;
  }

  static void Init(const std::string& app_name) {
    // Simple initialization for default logger
  }

  static void Shutdown() {
    // Simple cleanup for default logger
  }
};

}  // namespace default_logger
}  // namespace logger
}  // namespace tools
}  // namespace backends
}  // namespace metada