#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_
#define METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_

namespace metada {
namespace framework {
namespace tools {
namespace logger {

/**
 * @brief Traits class for configuring logger backend
 *
 * This template class allows specifying which concrete logger implementation
 * to use as the backend for the logging system. The backend type is provided
 * as a template parameter and exposed through the LoggerBackend type alias.
 *
 * The LoggerTraits class is used to configure which logging backend
 * implementation will be used by the Logger class template. This allows
 * switching between different logging backends (e.g. console, file, third-party
 * libraries) by providing different specializations of LoggerTraits.
 *
 * Typical usage:
 * @code
 * // Default traits using ConsoleLogger
 * template <>
 * struct LoggerTraits<void> {
 *   using LoggerBackend = ConsoleLogger;
 * };
 * @endcode
 *
 * @tparam Backend The concrete logger implementation to use as backend. The
 * backend must implement the ILogger interface.
 *
 * @see Logger
 * @see ILogger
 * @see ConsoleLogger
 * @see GoogleLogger
 */
template <typename Backend>
struct LoggerTraits {
  /**
   * @brief Type alias for the concrete logger backend implementation
   *
   * This type will be used by Logger to instantiate the actual logging backend.
   * The backend type must implement the ILogger interface.
   */
  using LoggerBackend = Backend;
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_