/**
 * @file MACOMlogging.hpp
 * @brief Simple standalone logging system for MACOM components
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace metada::backends::macom {

/**
 * @brief Enum for log levels
 */
enum class MACOMLogLevel { DEBUG, INFO, WARNING, ERROR };

/**
 * @brief Simple standalone logger for MACOM components
 */
class MACOMLogger {
 public:
  /**
   * @brief Set the global minimum log level
   *
   * @param level Minimum level to display
   */
  static void setGlobalLevel(MACOMLogLevel level) { s_minLevel = level; }

  /**
   * @brief Enable or disable timestamps
   *
   * @param enable True to show timestamps, false to hide
   */
  static void enableTimestamps(bool enable) { s_showTimestamps = enable; }

  /**
   * @brief Enable or disable component tagging
   *
   * @param enable True to show component tags, false to hide
   */
  static void enableComponentTags(bool enable) { s_showComponents = enable; }

  /**
   * @brief Log a message with specified level and component
   *
   * @param level Log level
   * @param component Component name (e.g., "MACOMGeometry")
   * @param message Log message
   */
  static void log(MACOMLogLevel level, const std::string& component,
                  const std::string& message) {
    // Skip if below minimum level
    if (level < s_minLevel) {
      return;
    }

    std::stringstream ss;

    // Add timestamp if enabled
    if (s_showTimestamps) {
      auto now = std::chrono::system_clock::now();
      auto time = std::chrono::system_clock::to_time_t(now);
      ss << "[" << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S")
         << "] ";
    }

    // Add log level tag
    switch (level) {
      case MACOMLogLevel::DEBUG:
        ss << "[DEBUG] ";
        break;
      case MACOMLogLevel::INFO:
        ss << "[INFO] ";
        break;
      case MACOMLogLevel::WARNING:
        ss << "[WARNING] ";
        break;
      case MACOMLogLevel::ERROR:
        ss << "[ERROR] ";
        break;
    }

    // Add component tag if enabled
    if (s_showComponents && !component.empty()) {
      ss << "[" << component << "] ";
    }

    // Add message
    ss << message;

    // Output to appropriate stream
    if (level == MACOMLogLevel::ERROR) {
      std::cerr << ss.str() << std::endl;
    } else {
      std::cout << ss.str() << std::endl;
    }
  }

  /**
   * @brief Log a debug message
   *
   * @param component Component name
   * @param message Log message
   */
  static void debug(const std::string& component, const std::string& message) {
    log(MACOMLogLevel::DEBUG, component, message);
  }

  /**
   * @brief Log an info message
   *
   * @param component Component name
   * @param message Log message
   */
  static void info(const std::string& component, const std::string& message) {
    log(MACOMLogLevel::INFO, component, message);
  }

  /**
   * @brief Log a warning message
   *
   * @param component Component name
   * @param message Log message
   */
  static void warning(const std::string& component,
                      const std::string& message) {
    log(MACOMLogLevel::WARNING, component, message);
  }

  /**
   * @brief Log an error message
   *
   * @param component Component name
   * @param message Log message
   */
  static void error(const std::string& component, const std::string& message) {
    log(MACOMLogLevel::ERROR, component, message);
  }

 private:
  static MACOMLogLevel s_minLevel;
  static bool s_showTimestamps;
  static bool s_showComponents;
};

// Initialize static members
MACOMLogLevel MACOMLogger::s_minLevel = MACOMLogLevel::INFO;
bool MACOMLogger::s_showTimestamps = true;
bool MACOMLogger::s_showComponents = true;

// Convenience macros
#define MACOM_LOG_DEBUG(component, message) \
  metada::backends::macom::MACOMLogger::debug(component, message)

#define MACOM_LOG_INFO(component, message) \
  metada::backends::macom::MACOMLogger::info(component, message)

#define MACOM_LOG_WARNING(component, message) \
  metada::backends::macom::MACOMLogger::warning(component, message)

#define MACOM_LOG_ERROR(component, message) \
  metada::backends::macom::MACOMLogger::error(component, message)

}  // namespace metada::backends::macom
