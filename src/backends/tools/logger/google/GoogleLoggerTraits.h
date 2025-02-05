#ifndef METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGERTRAITS_H_
#define METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGERTRAITS_H_

#include "GoogleLogger.h"
#include "LoggerTraits.h"

namespace metada {
namespace framework {
namespace tools {
namespace logger {

/**
 * @brief Template specialization of LoggerTraits for Google logging
 *
 * This specialization configures the logging system to use GoogleLogger
 * as the backend implementation. The void template parameter indicates
 * this is the default logging configuration.
 */
template <>
struct LoggerTraits<void> {
  /** @brief Type alias defining GoogleLogger as the concrete backend
   * implementation */
  using LoggerBackend = metada::GoogleLogger;
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_GOOGLE_GOOGLELOGGERTRAITS_H_