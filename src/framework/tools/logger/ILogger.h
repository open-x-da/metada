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
 * This abstract class defines the core logging interface that concrete logger
 * backends must implement. It provides pure virtual methods for different
 * logging levels (Info, Warning, Error, Debug) that derived classes need to
 * implement.
 *
 * Additionally, concrete implementations must provide static Init() and
 * Shutdown() methods for process-wide logging setup and cleanup. These static
 * methods are not part of the interface since they need to be called directly
 * on the concrete implementation class.
 *
 * @note This interface is designed to be backend-agnostic, allowing different
 * logging implementations (console, file, third-party libraries) to be used
 * interchangeably.
 *
 * @see Logger
 * @see LoggerTraits
 * @see ConsoleLogger
 */
class ILogger {
 public:
  /**
   * @brief Virtual destructor
   *
   * Ensures proper cleanup of derived classes through base pointer
   */
  virtual ~ILogger() = default;

  /**
   * @brief Log an informational message
   *
   * Used for general operational information about the normal functioning
   * of the application. Info messages should be concise and meaningful.
   *
   * @param message The message to log at INFO level
   */
  virtual void Info(const std::string& message) = 0;

  /**
   * @brief Log a warning message
   *
   * Used for potentially harmful situations or conditions that might
   * require attention but don't prevent the application from functioning.
   *
   * @param message The message to log at WARNING level
   */
  virtual void Warning(const std::string& message) = 0;

  /**
   * @brief Log an error message
   *
   * Used for error conditions that prevent normal operation or indicate
   * serious problems that need immediate attention.
   *
   * @param message The message to log at ERROR level
   */
  virtual void Error(const std::string& message) = 0;

  /**
   * @brief Log a debug message
   *
   * Used for detailed information useful during development and
   * troubleshooting. Debug messages may include technical details not relevant
   * for normal operation.
   *
   * @param message The message to log at DEBUG level
   */
  virtual void Debug(const std::string& message) = 0;

  /**
   * @note Required static methods for concrete implementations:
   *
   * static void Init(const std::string& app_name);
   * @brief Initialize the logging system
   * @param app_name Name of the application for logging identification
   *
   * static void Shutdown();
   * @brief Cleanup and shutdown the logging system
   *
   * These methods manage process-wide logging setup and cleanup,
   * which should be done once per application rather than per logger instance.
   */
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_ILOGGER_H_