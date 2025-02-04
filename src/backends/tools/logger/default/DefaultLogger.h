#ifndef METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_
#define METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_

#include <iostream>

#include "ILogger.h"

namespace metada {

/**
 * @brief Simple console-based logger implementation
 *
 * Provides basic logging functionality by writing messages to stdout/stderr
 * with severity level prefixes.
 */
class DefaultLogger : public ILogger {
 public:
  /**
   * @brief Log info message to stdout
   * @param message Message to log
   */
  void Info(const std::string& message) override {
    std::cout << "[INFO] " << message << std::endl;
  }

  /**
   * @brief Log warning message to stdout
   * @param message Message to log
   */
  void Warning(const std::string& message) override {
    std::cout << "[WARNING] " << message << std::endl;
  }

  /**
   * @brief Log error message to stderr
   * @param message Message to log
   */
  void Error(const std::string& message) override {
    std::cerr << "[ERROR] " << message << std::endl;
  }

  /**
   * @brief Log debug message to stdout
   * @param message Message to log
   */
  void Debug(const std::string& message) override {
    std::cout << "[DEBUG] " << message << std::endl;
  }

  /**
   * @brief Initialize logger for application
   * @param app_name Name of application using logger
   */
  static void Init(const std::string& app_name) {
    // Simple initialization for default logger
    std::cout << "[INFO] " << "Default logger initialized for " << app_name
              << std::endl;
  }

  /**
   * @brief Cleanup logger resources
   */
  static void Shutdown() {
    // Simple cleanup for default logger
    std::cout << "[INFO] " << "Default logger shutdown" << std::endl;
  }
};

template <>
struct LoggerTraits<void> {
  using LoggerBackend = DefaultLogger;
};

}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_DEFAULT_DEFAULTLOGGER_H_