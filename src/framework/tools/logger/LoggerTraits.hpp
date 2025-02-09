#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_HPP_
#define METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_HPP_

namespace metada {
namespace framework {
namespace tools {
namespace logger {

/**
 * @brief Primary template for logger backend traits
 *
 * This template class provides compile-time selection and customization of the
 * logging backend implementation used throughout the framework. It acts as a
 * traits class that maps backend types to their corresponding implementations.
 *
 * The primary template accepts any backend type and simply forwards it as the
 * LoggerBackend type. Specializations can be created to provide specific
 * backend mappings, such as mapping void to a default implementation.
 *
 * The logging backend must implement the ILogger interface and provide:
 * - Multiple severity levels (Info, Warning, Error, Debug)
 * - Process-wide initialization and cleanup
 * - Thread-safe logging operations
 * - Proper message formatting and output
 *
 * Example usage:
 * @code
 * // Default traits using ConsoleLogger
 * template <>
 * struct LoggerTraits<void> {
 *   using LoggerBackend = ConsoleLogger;
 * };
 * @endcode
 *
 * @tparam Backend The logging backend type to use
 *
 * @see Logger Main logger class using these traits
 * @see ILogger Base interface for logging backends
 */
template <typename Backend>
struct LoggerTraits {
  /** @brief Selected logging backend type */
  using LoggerBackend = Backend;
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_HPP_