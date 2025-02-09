#ifndef METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_HPP_
#define METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_HPP_

#include "ConsoleLogger.hpp"
#include "LoggerTraits.hpp"

namespace metada {
namespace framework {
namespace tools {
namespace logger {

/**
 * @brief Console logger backend selector traits
 *
 * This specialization of LoggerTraits for void provides compile-time selection
 * of the console logging backend implementation.
 *
 * The ConsoleLogger backend provides:
 * - Logging to standard output/error streams
 * - Support for different log levels (debug, info, warning, error)
 * - Colorized output based on log level
 * - Timestamp and log level prefixing
 *
 * @see ConsoleLogger Console logging backend implementation
 * @see LoggerTraits Generic logger traits template
 */
template <>
struct LoggerTraits<void> {
  /** @brief Selected logger backend type (Console implementation) */
  using LoggerBackend = metada::ConsoleLogger;
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_HPP_