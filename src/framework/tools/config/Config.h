#ifndef METADA_FRAMEWORK_TOOLS_CONFIG_CONFIG_H_
#define METADA_FRAMEWORK_TOOLS_CONFIG_CONFIG_H_

#include "IConfig.h"

namespace metada {
namespace framework {
namespace tools {
namespace config {

/**
 * @brief Main configuration class template providing a generic interface to
 * configuration backends
 *
 * This class template provides a static interface for loading, accessing,
 * modifying and saving configuration data using a backend specified by the
 * ConfigBackend template parameter. The backend must implement the IConfig
 * interface.
 *
 * The configuration data is stored in a hierarchical structure where keys use
 * dot notation to access nested values. For example, "database.host" would
 * access the "host" field within the "database" object.
 *
 * Supported value types include:
 * - Boolean
 * - Integer
 * - Double
 * - String
 * - Arrays of the above types
 *
 * @tparam ConfigBackend The configuration backend type that implements IConfig
 */
template <typename ConfigBackend>
class Config {
 private:
  ConfigBackend backend_;  ///< Instance of the configuration backend

 public:
  /** @brief Default constructor */
  Config() = default;

  /**
   * @brief Get direct access to the backend instance
   * @return Reference to the backend instance
   */
  ConfigBackend& backend() { return backend_; }

  /**
   * @brief Get const access to the backend instance
   * @return Const reference to the backend instance
   */
  const ConfigBackend& backend() const { return backend_; }

  /**
   * @brief Load configuration from a file
   * @param filename Path to the configuration file
   * @return true if loading was successful, false otherwise
   */
  bool LoadFromFile(const std::string& filename) {
    return backend_.LoadFromFile(filename);
  }

  /**
   * @brief Load configuration from a string
   * @param content String containing configuration data
   * @return true if loading was successful, false otherwise
   */
  bool LoadFromString(const std::string& content) {
    return backend_.LoadFromString(content);
  }

  /**
   * @brief Get a value from the configuration with a default fallback
   * @param key Dot-separated path to the configuration value (e.g.
   * "database.host")
   * @param default_value Value to return if key doesn't exist or access fails
   * @return ConfigValue containing either the requested value or the default
   * value
   */
  ConfigValue Get(const std::string& key,
                  const ConfigValue& default_value = ConfigValue()) {
    try {
      return backend_.Get(key);
    } catch (...) {
      return default_value;
    }
  }

  /**
   * @brief Set a value in the configuration
   * @param key Dot-separated path where to store the value (e.g.
   * "database.port")
   * @param value ConfigValue to store
   */
  void Set(const std::string& key, const ConfigValue& value) {
    backend_.Set(key, value);
  }

  /**
   * @brief Check if a key exists in the configuration
   * @param key Dot-separated path to check
   * @return true if the key exists, false otherwise
   */
  bool HasKey(const std::string& key) const { return backend_.HasKey(key); }

  /**
   * @brief Save configuration to a file
   * @param filename Path where to save the configuration
   * @return true if saving was successful, false otherwise
   */
  bool SaveToFile(const std::string& filename) const {
    return backend_.SaveToFile(filename);
  }

  /**
   * @brief Convert configuration to string representation
   * @return String containing the configuration data
   */
  std::string ToString() const { return backend_.ToString(); }

  /** @brief Clear all configuration data */
  void Clear() { backend_.Clear(); }
};

}  // namespace config
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_CONFIG_CONFIG_H_