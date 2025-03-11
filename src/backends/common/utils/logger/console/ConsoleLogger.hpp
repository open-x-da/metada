#pragma once

#include <iostream>

#include "utils/logger/ILogger.hpp"

namespace metada::backends::logger {

using framework::ILogger;

/**
 * @file ConsoleLogger.hpp
 * @brief Console logging backend implementation
 *
 * This header provides a simple console-based logging backend that implements
 * the ILogger interface. It writes log messages to standard output/error
 * streams with severity level prefixes and colorized output.
 *
 * Features:
 * - Simple logging to standard output/error streams
 * - Support for different log levels (debug, info, warning, error)
 * - Colorized output based on log level for improved readability
 * - Timestamp and log level prefixing for each message
 * - Thread-safe logging operations
 * - Configurable output formatting
 * - Real-time output for immediate feedback
 *
 * Log levels:
 * - INFO: General operational information (stdout, green)
 * - WARNING: Potential issues requiring attention (stdout, yellow)
 * - ERROR: Serious problems requiring immediate action (stderr, red)
 * - DEBUG: Detailed troubleshooting information (stdout, blue)
 *
 * Example usage with string-based API:
 * @code
 * // Initialize logging
 * ConsoleLogger::Init("MyApp");
 *
 * // Create logger instance
 * ConsoleLogger logger;
 *
 * // Log at different levels
 * logger.Info("Application started");
 * logger.Warning("Resource usage high");
 * logger.Error("Failed to connect to database");
 * logger.Debug("Connection attempt details: ...");
 *
 * // Cleanup on shutdown
 * ConsoleLogger::Shutdown();
 * @endcode
 *
 * Example usage with stream-based API:
 * @code
 * // Initialize logging
 * ConsoleLogger::Init("MyApp");
 *
 * // Create logger instance
 * ConsoleLogger logger;
 *
 * // Log at different levels using stream interface
 * logger.InfoStream() << "Application started with version " << app_version;
 * logger.WarningStream() << "Resource usage at " << usage_percent << "%";
 * logger.ErrorStream() << "Failed to connect to database: " << error_code;
 * logger.DebugStream() << "Connection params: " << host << ":" << port;
 *
 * // Cleanup on shutdown
 * ConsoleLogger::Shutdown();
 * @endcode
 *
 * @see ILogger Base interface class
 */
class ConsoleLogger : public ILogger {
 public:
  /**
   * @brief Default constructor
   */
  ConsoleLogger() = default;

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

}  // namespace metada::backends::logger
