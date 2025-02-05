#ifndef METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGER_H_
#define METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGER_H_

#include <glog/logging.h>

#include "ILogger.h"

namespace metada {

/**
 * @brief Google logging implementation of ILogger interface
 *
 * This class provides logging functionality by delegating to Google's glog
 * library. It implements the ILogger interface to provide a glog-based logging
 * backend.
 *
 * Messages are logged using glog's severity levels:
 * - INFO for informational messages
 * - WARNING for warning messages
 * - ERROR for error messages
 * - DLOG(INFO) for debug messages (only in debug builds)
 *
 * Glog provides additional features like:
 * - Automatic file naming and rotation
 * - Line numbers and timestamps
 * - Different log files per severity level
 * - Conditional logging
 *
 * @see ILogger Base interface class
 * @see https://github.com/google/glog Documentation for Google logging library
 */
class GoogleLogger : public framework::tools::logger::ILogger {
 public:
  /**
   * @brief Log info message using glog
   *
   * Writes an informational message using glog's INFO level. Info messages
   * are used for general operational information and are written to the info
   * log file.
   *
   * @param message The message to log at INFO level
   */
  void Info(const std::string& message) override { LOG(INFO) << message; }

  /**
   * @brief Log warning message using glog
   *
   * Writes a warning message using glog's WARNING level. Warning messages
   * indicate potential issues that may require attention and are written to the
   * warning log file.
   *
   * @param message The message to log at WARNING level
   */
  void Warning(const std::string& message) override { LOG(WARNING) << message; }

  /**
   * @brief Log error message using glog
   *
   * Writes an error message using glog's ERROR level. Error messages indicate
   * serious problems requiring immediate attention and are written to the error
   * log file.
   *
   * @param message The message to log at ERROR level
   */
  void Error(const std::string& message) override { LOG(ERROR) << message; }

  /**
   * @brief Log debug message using glog
   *
   * Writes a debug message using glog's debug logging (DLOG). Debug messages
   * provide detailed information for troubleshooting and are only included in
   * debug builds. In release builds, these messages are compiled out.
   *
   * @param message The message to log at DEBUG level
   */
  void Debug(const std::string& message) override { DLOG(INFO) << message; }

  /**
   * @brief Initialize glog for application
   *
   * Performs initialization of the Google logging library. This sets up log
   * file locations, naming patterns, and other glog configuration options. Must
   * be called before any logging occurs.
   *
   * @param app_name Name of application using logger, used in log file names
   * @note This is a static method that should be called once at application
   * startup
   */
  static void Init(const std::string& app_name) {
    google::InitGoogleLogging(app_name.c_str());
  }

  /**
   * @brief Shutdown glog logging
   *
   * Performs cleanup of the Google logging library. This ensures all log
   * messages are flushed and files are properly closed. Should be called before
   * application exit.
   *
   * @note This is a static method that should be called once at application
   * shutdown
   */
  static void Shutdown() { google::ShutdownGoogleLogging(); }
};

}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGER_H_
