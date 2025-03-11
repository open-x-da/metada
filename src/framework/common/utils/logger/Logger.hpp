#pragma once

#include <string>

#include "utils/NonCopyable.hpp"
#include "utils/logger/LogStream.hpp"

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
 * The Logger class is non-copyable to prevent unintended duplication of logger
 * instances, which could lead to resource contention or inconsistent logging
 * behavior. However, it supports move semantics to allow transferring ownership
 * when needed.
 *
 * Example usage with stream-based logging:
 * @code
 * ConsoleLogger backend;
 * Logger<ConsoleLogger> logger;
 * logger.Info() << "Application started with " << num_threads << " threads";
 * logger.Error() << "Failed to connect to " << server << ": " << error_msg;
 * @endcode
 *
 * @tparam Backend The concrete logger backend type that implements the ILogger
 *                 interface
 *
 * @see ILogger
 * @see LoggerTraits
 * @see NonCopyable
 * @see LogStream
 */
template <typename Backend>
class Logger : public NonCopyable {
 public:
  /**
   * @brief Default constructor
   *
   * Creates a new Logger instance with a default-constructed backend
   */
  Logger() = default;

  /**
   * @brief Move constructor - explicitly defined for compatibility with mock
   * objects
   */
  Logger(Logger&& other) noexcept : backend_(std::move(other.backend_)) {}

  /**
   * @brief Move assignment - explicitly defined for compatibility with mock
   * objects
   */
  Logger& operator=(Logger&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
    }
    return *this;
  }

  /**
   * @brief Default destructor
   */
  ~Logger() = default;

  /**
   * @brief Create a stream for info-level logging
   *
   * @return A LogStream object for stream-style logging with << operator
   */
  LogStream Info() { return LogStream(backend_, LogLevel::Info); }

  /**
   * @brief Create a stream for warning-level logging
   *
   * @return A LogStream object for stream-style logging with << operator
   */
  LogStream Warning() { return LogStream(backend_, LogLevel::Warning); }

  /**
   * @brief Create a stream for error-level logging
   *
   * @return A LogStream object for stream-style logging with << operator
   */
  LogStream Error() { return LogStream(backend_, LogLevel::Error); }

  /**
   * @brief Create a stream for debug-level logging
   *
   * @return A LogStream object for stream-style logging with << operator
   */
  LogStream Debug() { return LogStream(backend_, LogLevel::Debug); }

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
