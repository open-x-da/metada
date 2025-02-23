#ifndef METADA_BACKENDS_COMMON_UTILS_LOGGER_LOGGERBACKENDSELECTOR_HPP_
#define METADA_BACKENDS_COMMON_UTILS_LOGGER_LOGGERBACKENDSELECTOR_HPP_

/**
 * @file LoggerBackendSelector.hpp
 * @brief Header file for selecting the logging backend based on compile flags
 *
 * This header provides conditional inclusion of either Console or Glog logging
 * backend traits based on compile-time flags (USE_CONSOLE_LOGGER or
 * USE_GLOG_LOGGER). It allows the logging system to be built with different
 * backends without code changes.
 */

#ifdef USE_CONSOLE_LOGGER
#include "console/ConsoleLoggerTraits.hpp"
#elif defined(USE_GLOG_LOGGER)
#include "glog/GlogLoggerTraits.hpp"
#endif

#endif  // METADA_BACKENDS_COMMON_UTILS_LOGGER_LOGGERBACKENDSELECTOR_HPP_
