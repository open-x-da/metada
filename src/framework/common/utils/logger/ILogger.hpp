#pragma once

#include <string>

#include "utils/NonCopyable.hpp"

namespace metada::framework {

// Forward declarations
class LogStream;
enum class LogLevel;

/**
 * @brief Abstract interface for logger backend implementations
 *
 * This interface defines the contract that all logger backends must implement.
 * It provides a unified API for logging messages at different severity levels
 * in a backend-agnostic way using a stream-based interface.
 *
 * The ILogger interface is non-copyable to prevent unintended duplication of
 * logger instances, which could lead to resource contention or inconsistent
 * logging behavior. However, it supports move semantics to allow transferring
 * ownership when needed.
 *
 * Key features:
 * - Multiple severity levels (Info, Warning, Error, Debug)
 * - Process-wide initialization and cleanup
 * - Backend-agnostic interface
 * - Stream-based logging with << operator support
 *
 * Example usage:
 * @code
 * class ConsoleLogger : public ILogger {
 *   // Implement internal methods to handle messages
 *   virtual void LogMessage(LogLevel level, const std::string& message)
 * override;
 *   // Implement stream methods by delegating to base class
 *   // ... static Init() and Shutdown() methods
 * };
 * @endcode
 *
 * Concrete implementations must provide:
 * - Method to log messages at different levels
 * - Static Init() and Shutdown() methods for process-wide setup
 * - Thread-safe logging operations
 * - Proper message formatting and output
 *
 * @see Logger Main logger class template using this interface
 * @see LoggerTraits Compile-time backend selection
 * @see NonCopyable
 */
class ILogger : public NonCopyable {
 public:
  /**
   * @brief Virtual destructor
   *
   * Ensures proper cleanup of derived classes through base pointer
   */
  virtual ~ILogger() = default;

  /**
   * @brief Internal method to log a message at a specific level
   *
   * This is the core method that concrete implementations must provide.
   * It receives messages at different severity levels and handles the
   * actual logging operation.
   *
   * @param level The severity level of the message
   * @param message The formatted message to log
   */
  virtual void LogMessage(LogLevel level, const std::string& message) = 0;

  /**
   * @brief Create a LogStream for info-level logging
   *
   * @return A LogStream object for chaining with << operator
   */
  virtual LogStream InfoStream();

  /**
   * @brief Create a LogStream for warning-level logging
   *
   * @return A LogStream object for chaining with << operator
   */
  virtual LogStream WarningStream();

  /**
   * @brief Create a LogStream for error-level logging
   *
   * @return A LogStream object for chaining with << operator
   */
  virtual LogStream ErrorStream();

  /**
   * @brief Create a LogStream for debug-level logging
   *
   * @return A LogStream object for chaining with << operator
   */
  virtual LogStream DebugStream();

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

}  // namespace metada::framework