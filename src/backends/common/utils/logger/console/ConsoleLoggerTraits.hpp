#ifndef METADA_BACKENDS_COMMON_UTILS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_HPP_
#define METADA_BACKENDS_COMMON_UTILS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_HPP_

#include "ConsoleLogger.hpp"
#include "utils/logger/LoggerTraits.hpp"

namespace metada {
namespace framework {
namespace common {
namespace utils {
namespace logger {

/**
 * @file ConsoleLoggerTraits.hpp
 * @brief Traits specialization for console logging backend
 * @ingroup logging
 *
 * @details
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
 * Example usage:
 * @code{.cpp}
 * using Logger = LoggerTraits<void>::LoggerBackend;
 * Logger logger;
 * logger.Info("Application started");
 * @endcode
 *
 * @see ConsoleLogger Console logging backend implementation
 * @see LoggerTraits Generic logger traits template
 * @since 0.1.0
 */
template <>
struct LoggerTraits<void> {
  /**
   * @brief Selected logger backend type (Console implementation)
   * @details Defines the concrete logger implementation to be used when the
   * void specialization is selected. In this case, it maps to the
   * ConsoleLogger.
   */
  using LoggerBackend = backends::common::utils::logger::console::ConsoleLogger;
};

}  // namespace logger
}  // namespace utils
}  // namespace common
}  // namespace framework
}  // namespace metada

#endif  // METADA_BACKENDS_COMMON_UTILS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_HPP_