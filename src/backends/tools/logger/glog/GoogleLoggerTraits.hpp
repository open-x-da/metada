#ifndef METADA_BACKENDS_TOOLS_LOGGER_GLOG_GOOGLELOGGERTRAITS_HPP_
#define METADA_BACKENDS_TOOLS_LOGGER_GLOG_GOOGLELOGGERTRAITS_HPP_

#include "GoogleLogger.hpp"
#include "LoggerTraits.hpp"

namespace metada {
namespace framework {
namespace tools {
namespace logger {

/**
 * @file GoogleLoggerTraits.hpp
 * @brief Traits specialization for Google logging backend
 * @ingroup logging
 *
 * @details
 * This header provides a specialization of the LoggerTraits template for the
 * Google logging (glog) backend. It maps the generic logging interface to
 * the glog-specific implementation provided by GoogleLogger.
 *
 * The GoogleLogger backend provides:
 * - Advanced logging with automatic file management
 * - Support for different log levels (debug, info, warning, error)
 * - Automatic log file naming and rotation
 * - Timestamps and source location information
 * - Separate log files per severity level
 * - Thread-safe logging operations
 * - Configurable log file locations
 * - Conditional logging and debug-only logging
 *
 * Example usage:
 * @code{.cpp}
 * using Logger = LoggerTraits<void>::LoggerBackend;
 * Logger logger;
 * logger.Info("Application started");
 * @endcode
 *
 * @see GoogleLogger Google logging backend implementation
 * @see LoggerTraits Generic logger traits template
 * @since 0.1.0
 */
template <>
struct LoggerTraits<void> {
  /**
   * @brief Selected logger backend type (Google implementation)
   * @details Defines the concrete logger implementation to be used when the void
   * specialization is selected. In this case, it maps to the GoogleLogger.
   */
  using LoggerBackend = backends::tools::logger::glog::GoogleLogger;
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_GLOG_GOOGLELOGGERTRAITS_HPP_