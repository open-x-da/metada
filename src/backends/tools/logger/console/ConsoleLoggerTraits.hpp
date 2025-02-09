#ifndef METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_HPP_
#define METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_HPP_

#include "ConsoleLogger.hpp"
#include "LoggerTraits.hpp"

/**
 * @file ConsoleLoggerTraits.hpp
 * @brief Traits specialization for console logging backend
 *
 * This header provides a specialization of the LoggerTraits template for the
 * console logging backend. It maps the generic logging interface to the
 * console-specific implementation provided by ConsoleLogger.
 *
 * The ConsoleLogger backend provides:
 * - Simple logging to standard output/error streams
 * - Support for different log levels (debug, info, warning, error)
 * - Colorized output based on log level for improved readability
 * - Timestamp and log level prefixing for each message
 * - Thread-safe logging operations
 * - Configurable output formatting
 * - Real-time output for immediate feedback
 *
 * @see ConsoleLogger Console logging backend implementation
 * @see LoggerTraits Generic logger traits template
 */
template <>
struct metada::framework::tools::logger::LoggerTraits<void> {
  /** @brief Selected logger backend type (Console implementation) */
  using LoggerBackend = metada::backends::tools::logger::console::ConsoleLogger;
};

#endif  // METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_HPP_