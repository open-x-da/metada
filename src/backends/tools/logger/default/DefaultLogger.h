#ifndef METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_
#define METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_

#include <iostream>
#include <string>

namespace metada {
namespace logger {

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

}  // namespace logger
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_