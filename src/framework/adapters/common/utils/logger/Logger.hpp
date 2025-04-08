#pragma once

#include <string>

#include "BackendTraits.hpp"
#include "LogStream.hpp"
#include "LoggerConcepts.hpp"
#include "NonCopyable.hpp"

namespace metada::framework {

/**
 * @brief Forward declaration of Config class
 */
template <typename BackendTag>
  requires ConfigBackendType<
      typename traits::BackendTraits<BackendTag>::ConfigBackend>
class Config;

/**
 * @brief Generic logger class that delegates to a backend implementation
 *
 * @details
 * This class template provides a generic logging interface that forwards
 * logging calls to a concrete backend implementation. The backend type is
 * specified as a template parameter through a BackendTag, which is used with
 * BackendTraits to determine the actual logger backend implementation.
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
 * Logger instances are typically initialized with a Config object that provides
 * configuration settings for the logger backend.
 *
 * @par Example usage with stream-based logging:
 * @code
 * Config<MyBackendTag> config("config.yaml");
 * Logger<MyBackendTag> logger(config);
 * logger.Info() << "Application started with " << num_threads << " threads";
 * logger.Error() << "Failed to connect to " << server << ": " << error_msg;
 * @endcode
 *
 * @tparam BackendTag The backend tag type that provides the LoggerBackend type
 *                    through BackendTraits
 *
 * @see BackendTraits
 * @see NonCopyable
 * @see LogStream
 * @see Config
 */
template <LoggerBackendType BackendTag>
class Logger : public NonCopyable {
 public:
  using LoggerBackend =
      typename traits::BackendTraits<BackendTag>::LoggerBackend;
  using ConfigBackend =
      typename traits::BackendTraits<BackendTag>::ConfigBackend;

  /**
   * @brief Disabled default constructor
   */
  Logger() = delete;

  /**
   * @brief Default destructor
   */
  ~Logger() = default;

  /**
   * @brief Disabled copy constructor
   * @param[in] other The logger instance to copy from
   */
  Logger(const Logger& other) = delete;

  /**
   * @brief Disabled copy assignment operator
   * @param[in] other The logger instance to copy from
   * @return Reference to this logger instance
   */
  Logger& operator=(const Logger& other) = delete;

  /**
   * @brief Move constructor
   * @param[in] other The logger instance to move from
   */
  Logger(Logger&& other) noexcept : backend_(std::move(other.backend_)) {}

  /**
   * @brief Move assignment operator
   * @details Explicitly defined for compatibility with mock objects
   * @param[in] other The logger instance to move from
   * @return Reference to this logger instance
   */
  Logger& operator=(Logger&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
    }
    return *this;
  }

  /**
   * @brief Constructor that takes a config
   * @param[in] config The config to use for initializing the logger backend
   */
  Logger(const Config<BackendTag>& config)
      : backend_(LoggerBackend(config.backend())) {}

  /**
   * @brief Create a stream for info-level logging
   * @return A LogStream object for stream-style logging with << operator
   */
  LogStream<LoggerBackend, ConfigBackend> Info() {
    return LogStream<LoggerBackend, ConfigBackend>(backend_, LogLevel::Info);
  }

  /**
   * @brief Create a stream for warning-level logging
   * @return A LogStream object for stream-style logging with << operator
   */
  LogStream<LoggerBackend, ConfigBackend> Warning() {
    return LogStream<LoggerBackend, ConfigBackend>(backend_, LogLevel::Warning);
  }

  /**
   * @brief Create a stream for error-level logging
   * @return A LogStream object for stream-style logging with << operator
   */
  LogStream<LoggerBackend, ConfigBackend> Error() {
    return LogStream<LoggerBackend, ConfigBackend>(backend_, LogLevel::Error);
  }

  /**
   * @brief Create a stream for debug-level logging
   * @return A LogStream object for stream-style logging with << operator
   */
  LogStream<LoggerBackend, ConfigBackend> Debug() {
    return LogStream<LoggerBackend, ConfigBackend>(backend_, LogLevel::Debug);
  }

  /**
   * @brief Get reference to the underlying backend
   * @details Provides access to the concrete logger backend instance. This can
   * be useful for backend-specific configuration or testing purposes.
   * @return Reference to the backend logger implementation
   */
  LoggerBackend& backend() { return backend_; }

 private:
  LoggerBackend backend_;  ///< The underlying logger backend instance that
                           ///< performs actual logging
};

}  // namespace metada::framework
