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

  // Note: Concrete implementations must provide the following static methods:
  //   static void Init(const std::string& app_name);
  //   static void Shutdown();
  // These methods manage process-wide logging setup and cleanup,
  // which should be done once per application rather than per logger instance.
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_ILOGGER_H_