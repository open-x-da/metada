#pragma once

#include <string>

#include "utils/NonCopyable.hpp"

namespace metada::framework {

/**
 * @brief Abstract interface for logger backend implementations
 *
 * This interface defines the contract that all logger backends must implement.
 * It provides a unified API for logging messages at different severity levels
 * in a backend-agnostic way.
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
 *
 * Example usage:
 * @code
 * class ConsoleLogger : public ILogger {
 *   void Info(const std::string& message) override;
 *   void Warning(const std::string& message) override;
 *   void Error(const std::string& message) override;
 *   void Debug(const std::string& message) override;
 *   // ... static Init() and Shutdown() methods
 * };
 * @endcode
 *
 * Concrete implementations must provide:
 * - Logging methods for each severity level
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

}  // namespace metada::framework