#pragma once

#include <string>

namespace metada::framework {

/**
 * @brief Generic logger class that delegates to a backend implementation
 *
 * This class template provides a generic logging interface that forwards
 * logging calls to a concrete backend implementation. The backend type is
 * specified as a template parameter and must implement the ILogger interface.
 *
 * The Logger class supports different logging levels (Info, Warning, Error,
 * Debug) and delegates all logging operations to the backend implementation.
 * This allows switching between different logging backends (e.g. console, file,
 * third-party logging libraries) without changing the logging code.
 *
 * Example usage:
 * @code
 * ConsoleLogger backend;
 * Logger<ConsoleLogger> logger;
 * logger.Info("Application started");
 * logger.Error("An error occurred");
 * @endcode
 *
 * @tparam Backend The concrete logger backend type that implements the ILogger
 *                 interface
 *
 * @see ILogger
 * @see LoggerTraits
 */
template <typename Backend>
class Logger {
 public:
  /**
   * @brief Default constructor
   *
   * Creates a new Logger instance with a default-constructed backend
   */
  Logger() = default;

  /**
   * @brief Default destructor
   */
  ~Logger() = default;

  /**
   * @brief Log an informational message
   *
   * Logs a message at INFO level. Info messages are used for general
   * operational information about the normal functioning of the application.
   *
   * @param message The message to log at INFO level
   */
  void Info(const std::string& message) { backend_.Info(message); }

  /**
   * @brief Log a warning message
   *
   * Logs a message at WARNING level. Warning messages indicate potential
   * issues that are not immediately threatening but may require attention.
   *
   * @param message The message to log at WARNING level
   */
  void Warning(const std::string& message) { backend_.Warning(message); }

  /**
   * @brief Log an error message
   *
   * Logs a message at ERROR level. Error messages indicate serious problems
   * that need immediate attention and may prevent proper functioning.
   *
   * @param message The message to log at ERROR level
   */
  void Error(const std::string& message) { backend_.Error(message); }

  /**
   * @brief Log a debug message
   *
   * Logs a message at DEBUG level. Debug messages provide detailed information
   * useful during development and troubleshooting.
   *
   * @param message The message to log at DEBUG level
   */
  void Debug(const std::string& message) { backend_.Debug(message); }

  /**
   * @brief Get reference to the underlying backend
   *
   * Provides access to the concrete logger backend instance. This can be
   * useful for backend-specific configuration or testing purposes.
   *
   * @return Reference to the backend logger implementation
   */
  Backend& backend() { return backend_; }

 private:
  Backend backend_;  ///< The underlying logger backend instance that performs
                     ///< actual logging
};

}  // namespace metada::framework
