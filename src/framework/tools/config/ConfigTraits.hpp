#ifndef METADA_FRAMEWORK_TOOLS_CONFIG_CONFIGTRAITS_HPP_
#define METADA_FRAMEWORK_TOOLS_CONFIG_CONFIGTRAITS_HPP_

namespace metada {
namespace framework {
namespace tools {
namespace config {

/**
 * @brief Primary template for configuration backend traits
 *
 * This template class provides compile-time selection and customization of the
 * configuration backend implementation used throughout the framework. It acts
 * as a traits class that maps backend types to their corresponding
 * implementations.
 *
 * The primary template accepts any backend type and simply forwards it as the
 * ConfigBackend type. Specializations can be created to provide specific
 * backend mappings, such as mapping void to a default implementation.
 *
 * The configuration backend must implement the IConfig interface and provide:
 * - Loading/saving from files and strings
 * - Type-safe value access with defaults
 * - Hierarchical configuration using dot notation
 * - Support for basic types and arrays
 *
 * @tparam Backend The configuration backend type to use
 *
 * @see Config Main configuration class using these traits
 * @see IConfig Base interface for configuration backends
 */
template <typename Backend>
struct ConfigTraits {
  /** @brief Selected configuration backend type */
  using ConfigBackend = Backend;
};

}  // namespace config
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_CONFIG_CONFIGTRAITS_HPP_