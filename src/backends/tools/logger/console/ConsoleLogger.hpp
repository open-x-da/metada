#ifndef METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGER_HPP_
#define METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGER_HPP_

#include <iostream>

#include "ILogger.hpp"

namespace metada {

/**
 * @brief Console logger backend implementation
 *
 * This class implements the ILogger interface using standard output/error
 * streams as the logging backend. It provides basic logging functionality with
 * severity level prefixes and colorized output.
 *
 * The logger writes messages to:
 * - stdout for INFO, WARNING and DEBUG levels
 * - stderr for ERROR level
 *
 * Messages are formatted with severity level prefixes:
 * - [INFO] for informational messages (green)
 * - [WARNING] for warning messages (yellow)
 * - [ERROR] for error messages (red)
 * - [DEBUG] for debug messages (blue)
 *
 * Example usage:
 * @code
 * ConsoleLogger logger;
 * logger.Info("Application started");
 * logger.Warning("Resource usage high");
 * logger.Error("Failed to connect to database");
 * logger.Debug("Connection attempt details: ...");
 * @endcode
 *
 * The logger provides static Init() and Shutdown() methods for application-wide
 * initialization and cleanup:
 * @code
 * ConsoleLogger::Init("MyApp");
 * // ... application code ...
 * ConsoleLogger::Shutdown();
 * @endcode
 *
 * @see ILogger Base interface class
 */
class ConsoleLogger : public framework::tools::logger::ILogger {
 public:
  /**
   * @brief Log info message to stdout
   *
   * Writes an informational message to standard output with [INFO] prefix.
   * Info messages are used for general operational information.
   *
   * @param message The message to log at INFO level
   */
  void Info(const std::string& message) override {
    std::cout << "[INFO] " << message << std::endl;
  }

  /**
   * @brief Log warning message to stdout
   *
   * Writes a warning message to standard output with [WARNING] prefix.
   * Warning messages indicate potential issues that may require attention.
   *
   * @param message The message to log at WARNING level
   */
  void Warning(const std::string& message) override {
    std::cout << "[WARNING] " << message << std::endl;
  }

  /**
   * @brief Log error message to stderr
   *
   * Writes an error message to standard error with [ERROR] prefix.
   * Error messages indicate serious problems requiring immediate attention.
   *
   * @param message The message to log at ERROR level
   */
  void Error(const std::string& message) override {
    std::cerr << "[ERROR] " << message << std::endl;
  }

  /**
   * @brief Log debug message to stdout
   *
   * Writes a debug message to standard output with [DEBUG] prefix.
   * Debug messages provide detailed information for troubleshooting.
   *
   * @param message The message to log at DEBUG level
   */
  void Debug(const std::string& message) override {
    std::cout << "[DEBUG] " << message << std::endl;
  }

  /**
   * @brief Initialize logger for application
   *
   * Performs basic initialization of the console logger. Outputs an
   * initialization message to stdout with the application name.
   *
   * @param app_name Name of application using the logger
   * @note This is a static method that should be called once at application
   * startup
   */
  static void Init(const std::string& app_name) {
    // Simple initialization for default logger
    std::cout << "[INFO] " << "Console logger initialized for " << app_name
              << std::endl;
  }

  /**
   * @brief Cleanup logger resources
   *
   * Performs cleanup of the console logger. Outputs a shutdown message to
   * stdout. For this simple implementation, no actual cleanup is needed since
   * we only use standard streams.
   *
   * @note This is a static method that should be called once at application
   * shutdown
   */
  static void Shutdown() {
    // Simple cleanup for default logger
    std::cout << "[INFO] " << "Console logger shutdown" << std::endl;
  }
};

}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGER_HPP_