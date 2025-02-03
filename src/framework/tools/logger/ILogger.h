#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_ILOGGER_H_
#define METADA_FRAMEWORK_TOOLS_LOGGER_ILOGGER_H_

#include <string>

namespace metada {
namespace framework {
namespace tools {
namespace logger {

class ILogger {
 public:
  virtual ~ILogger() = default;

  virtual void Info(const std::string& message) = 0;
  virtual void Warning(const std::string& message) = 0;
  virtual void Error(const std::string& message) = 0;
  virtual void Debug(const std::string& message) = 0;

  // Static methods for initialization/shutdown should be implemented by
  // concrete loggers static void Init(const std::string& app_name) = 0; static
  // void Shutdown() = 0;
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_ILOGGER_H_