#pragma once

namespace metada::framework {

/**
 * @brief Enumeration of log severity levels
 *
 * @details This enumeration defines the standard logging levels used throughout
 * the framework. The levels are ordered from least severe (Debug) to most
 * severe (Error). Logger implementations use these levels to filter and
 * categorize log messages.
 *
 * The levels are used by the Logger class and LogStream to specify the severity
 * of log messages. Backend implementations can use these levels to determine
 * how to format or filter messages.
 *
 * @see Logger
 * @see LogStream
 * @see LoggerBackend
 */
enum class LogLevel { Debug, Info, Warning, Error };

}  // namespace metada::framework