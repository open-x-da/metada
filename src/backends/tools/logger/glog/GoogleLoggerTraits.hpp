#ifndef METADA_BACKENDS_TOOLS_LOGGER_GLOG_GOOGLELOGGERTRAITS_HPP_
#define METADA_BACKENDS_TOOLS_LOGGER_GLOG_GOOGLELOGGERTRAITS_HPP_

#include "GoogleLogger.hpp"
#include "LoggerTraits.hpp"

/**
 * @file GoogleLoggerTraits.hpp
 * @brief Traits specialization for Google logging backend
 *
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
 * @see GoogleLogger Google logging backend implementation
 * @see LoggerTraits Generic logger traits template
 */
template <>
struct metada::framework::tools::logger::LoggerTraits<void> {
  /** @brief Selected logger backend type (Google implementation) */
  using LoggerBackend = metada::backends::tools::logger::glog::GoogleLogger;
};

#endif  // METADA_BACKENDS_TOOLS_LOGGER_GLOG_GOOGLELOGGERTRAITS_HPP_