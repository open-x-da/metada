#ifndef METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_
#define METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_

#include <iostream>
#include <string>

namespace metada {
namespace backends {
namespace tools {
namespace logger {
namespace default_logger {

class DefaultLogger {
 public:
  static void Init(const std::string&) {}
  static void Shutdown() {}

  static void LogInfo(const std::string& message) {
    std::cout << "[INFO] " << message << std::endl;
  }

  static void LogWarning(const std::string& message) {
    std::cerr << "[WARNING] " << message << std::endl;
  }

  static void LogError(const std::string& message) {
    std::cerr << "[ERROR] " << message << std::endl;
  }

  static void LogDebug(const std::string& message) {
    std::cout << "[DEBUG] " << message << std::endl;
  }
};

}  // namespace default_logger
}  // namespace logger
}  // namespace tools
}  // namespace backends
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_