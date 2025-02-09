#ifndef METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGERTRAITS_HPP_
#define METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGERTRAITS_HPP_

#include "GoogleLogger.hpp"
#include "LoggerTraits.hpp"

namespace metada {
namespace framework {
namespace tools {
namespace logger {

/**
 * @brief Google logger backend selector traits
 *
 * This specialization of LoggerTraits for void provides compile-time selection
 * of the Google logging backend implementation.
 *
 * The GoogleLogger backend provides:
 * - Advanced logging with automatic file management
 * - Support for different log levels (debug, info, warning, error)
 * - Automatic log file naming and rotation
 * - Timestamps and source location information
 * - Separate log files per severity level
 * - Conditional logging and debug-only logging
 *
 * @see GoogleLogger Google logging backend implementation
 * @see LoggerTraits Generic logger traits template
 */
template <>
struct LoggerTraits<void> {
  /** @brief Selected logger backend type (Google implementation) */
  using LoggerBackend = metada::GoogleLogger;
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGERTRAITS_HPP_