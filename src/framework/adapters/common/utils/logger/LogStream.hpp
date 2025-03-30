#pragma once

#include <sstream>
#include <string>
#include <utility>

#include "NonCopyable.hpp"

namespace metada::framework {

/**
 * @brief Enumeration of log severity levels
 */
enum class LogLevel { Debug, Info, Warning, Error };

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
 * @tparam Backend The logger backend type that implements the LogMessage method
 */
template <typename Backend>
class LogStream : public NonCopyable {
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
   */
  ~LogStream() {
    if (!moved_) {
      Flush();
    }
  }

  /**
   * @brief Move assignment is deleted to prevent potential issues with stream
   * ownership
   */
  LogStream& operator=(LogStream&&) = delete;

  /**
   * @brief Flush the stream content to the logger
   *
   * Sends the accumulated stream content to the appropriate logger method
   * based on the log level.
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
   * @return Reference to this LogStream for chaining
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