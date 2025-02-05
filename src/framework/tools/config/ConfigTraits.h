#ifndef METADA_FRAMEWORK_TOOLS_CONFIG_CONFIG_TRAITS_H_
#define METADA_FRAMEWORK_TOOLS_CONFIG_CONFIG_TRAITS_H_

namespace metada {
namespace framework {
namespace tools {
namespace config {

/**
 * @brief Primary template for configuration backend traits
 *
 * This template class provides a way to specify and customize the configuration
 * backend type used throughout the framework. It acts as a traits class that
 * maps backend types to their corresponding implementations.
 *
 * The primary template accepts any backend type and simply forwards it as the
 * ConfigBackend type. Specializations of this template can be created to
 * provide specific backend mappings, such as mapping void to a default
 * implementation.
 *
 * @tparam Backend The configuration backend type to use
 *
 * @see Config The main configuration class that uses these traits
 * @see YamlConfigTraits Example specialization for YAML backend
 * @see JsonConfigTraits Example specialization for JSON backend
 */
template <typename Backend>
struct ConfigTraits {
  /** @brief The configuration backend type to use */
  using ConfigBackend = Backend;
};

}  // namespace config
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_CONFIG_CONFIG_TRAITS_H_