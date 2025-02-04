#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_ILOGGER_H_
#define METADA_FRAMEWORK_TOOLS_LOGGER_ILOGGER_H_

#include <string>

namespace metada {
namespace framework {
namespace tools {
namespace logger {

/**
 * @brief Interface for logger implementations
 *
 * Defines the core logging interface that concrete logger backends must
 * implement. Also requires static Init/Shutdown methods for process-wide setup.
 */
class ILogger {
 public:
  virtual ~ILogger() = default;

  /**
   * @brief Log an informational message
   * @param message The message to log
   */
  virtual void Info(const std::string& message) = 0;

  /**
   * @brief Log a warning message
   * @param message The message to log
   */
  virtual void Warning(const std::string& message) = 0;

  /**
   * @brief Log an error message
   * @param message The message to log
   */
  virtual void Error(const std::string& message) = 0;

  /**
   * @brief Log a debug message
   * @param message The message to log
   */
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