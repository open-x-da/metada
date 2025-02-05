#ifndef METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_H_
#define METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_H_

#include "ConsoleLogger.h"
#include "LoggerTraits.h"

namespace metada {
namespace framework {
namespace tools {
namespace logger {

/**
 * @brief Template specialization of LoggerTraits for console logging
 *
 * This specialization configures the logging system to use ConsoleLogger
 * as the backend implementation. The void template parameter indicates
 * this is the default logging configuration.
 */
template <>
struct LoggerTraits<void> {
  /** @brief Type alias defining ConsoleLogger as the concrete backend
   * implementation */
  using LoggerBackend = metada::ConsoleLogger;
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_LOGGER_CONSOLE_CONSOLELOGGERTRAITS_H_