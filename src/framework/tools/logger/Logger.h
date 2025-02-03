#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_LOGGER_H_
#define METADA_FRAMEWORK_TOOLS_LOGGER_LOGGER_H_

#include <string>

namespace metada {
namespace framework {
namespace tools {
namespace logger {

template <typename Backend>
class Logger {
 public:
  Logger() = default;
  ~Logger() = default;

  void Info(const std::string& message) { Backend::LogInfo(message); }

  void Warning(const std::string& message) { Backend::LogWarning(message); }

  void Error(const std::string& message) { Backend::LogError(message); }

  void Debug(const std::string& message) { Backend::LogDebug(message); }
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_LOGGER_H_
