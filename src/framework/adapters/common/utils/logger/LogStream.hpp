#pragma once

#include <sstream>
#include <string>
#include <utility>

#include "LogLevel.hpp"
#include "LoggerConcepts.hpp"

namespace metada::framework {

/**
 * @brief Stream-like class for logging that supports the << operator
 *
 * This class provides a stream interface similar to std::cout for logging
 * messages. It allows chaining of << operators to build log messages from
 * different types of data.
 *
 * The LogStream captures all data inserted via the << operator in an internal
 * string stream and forwards the complete message to the appropriate logger
 * method when the stream is destroyed or explicitly flushed.
 *
 * Example usage:
 * @code
 * logger.Info() << "User " << user_id << " logged in from " << ip_address;
 * logger.Error() << "Failed to process request: " << error_code << " - " <<
 * error_message;
 * @endcode
 *
 * @tparam Backend The logger backend type that implements the LoggerBackend
 * concept
 * @tparam ConfigBackend The configuration backend type required by the logger
 */
template <typename Backend, typename ConfigBackend>
  requires LoggerBackend<Backend, ConfigBackend>
class LogStream {
 public:
  /**
   * @brief Constructor that associates the stream with a logger and level
   *
   * @param logger Reference to the logger that will receive the message
   * @param level The severity level for this log message
   */
  LogStream(Backend& logger, LogLevel level)
      : logger_(logger), level_(level), moved_(false) {}

  /**
   * @brief Move constructor
   *
   * Transfers ownership of the stream content from another LogStream object.
   * Sets the moved_ flag on the source object to prevent it from flushing
   * in its destructor.
   *
   * @param other The LogStream to move from
   */
  LogStream(LogStream&& other) noexcept
      : logger_(other.logger_),
        level_(other.level_),
        stream_(std::move(other.stream_)),
        moved_(false) {
    other.moved_ = true;
  }

  /**
   * @brief Destructor that flushes the stream if it hasn't been moved
   *
   * If this LogStream still owns its content (hasn't been moved from),
   * the destructor will automatically flush the accumulated message to
   * the logger backend.
   */
  ~LogStream() {
    if (!moved_) {
      Flush();
    }
  }

  /**
   * @brief Move assignment is deleted to prevent potential issues with stream
   * ownership
   *
   * @return LogStream& Reference to this LogStream
   */
  LogStream& operator=(LogStream&&) = delete;

  /**
   * @brief Flush the stream content to the logger
   *
   * Sends the accumulated stream content to the logger backend using the
   * LogMessage method with the appropriate log level. After flushing,
   * the internal stream is cleared for potential reuse.
   */
  void Flush() {
    const std::string message = stream_.str();

    // Use the unified LogMessage method instead of individual severity methods
    logger_.LogMessage(level_, message);

    // Clear the stream for potential reuse
    stream_.str("");
    stream_.clear();
  }

  /**
   * @brief Generic stream insertion operator
   *
   * Allows inserting any type that can be inserted into std::ostringstream.
   *
   * @tparam T The type of value being inserted
   * @param value The value to insert into the log stream
   * @return LogStream& Reference to this LogStream for chaining
   */
  template <typename T>
  LogStream& operator<<(const T& value) {
    stream_ << value;
    return *this;
  }

 private:
  Backend& logger_;            ///< Reference to the logger backend
  LogLevel level_;             ///< Log level for this stream
  std::ostringstream stream_;  ///< Internal stream that accumulates the message
  bool moved_;  ///< Flag to indicate if this object has been moved from
};

}  // namespace metada::framework