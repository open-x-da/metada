#ifndef METADA_FRAMEWORK_TOOLS_CONFIG_ICONFIG_H_
#define METADA_FRAMEWORK_TOOLS_CONFIG_ICONFIG_H_

#include <string>
#include <variant>
#include <vector>

namespace metada {
namespace framework {
namespace tools {
namespace config {

/**
 * @brief Configuration value type that can hold different types of data
 *
 * This variant type supports the following value types:
 * - bool: Boolean values (true/false)
 * - int: Integer numbers
 * - double: Floating point numbers
 * - std::string: Text strings
 * - std::vector<bool>: Arrays of boolean values
 * - std::vector<int>: Arrays of integer numbers
 * - std::vector<double>: Arrays of floating point numbers
 * - std::vector<std::string>: Arrays of text strings
 */
using ConfigValue =
    std::variant<bool, int, double, std::string, std::vector<bool>,
                 std::vector<int>, std::vector<double>,
                 std::vector<std::string> >;

/**
 * @brief Interface class for configuration backends
 *
 * This abstract class defines the interface that all configuration backend
 * implementations must follow. It provides methods for loading, accessing,
 * modifying and saving configuration data in a backend-agnostic way.
 *
 * The configuration data is stored in a hierarchical structure where keys use
 * dot notation to access nested values. For example, "database.host" would
 * access the "host" field within the "database" object.
 *
 * Configuration backends implementing this interface must support all value
 * types defined in ConfigValue and handle the hierarchical key structure
 * appropriately.
 */
class IConfig {
 public:
  /**
   * @brief Virtual destructor for proper cleanup of derived classes
   */
  virtual ~IConfig() = default;

  /**
   * @brief Load configuration from a file
   * @param filename Path to the configuration file
   * @return true if loading was successful, false otherwise
   * @throws std::runtime_error if the file cannot be opened or contains invalid
   * data
   */
  virtual bool LoadFromFile(const std::string& filename) = 0;

  /**
   * @brief Load configuration from a string
   * @param content String containing the configuration data
   * @return true if loading was successful, false otherwise
   * @throws std::runtime_error if the content contains invalid data
   */
  virtual bool LoadFromString(const std::string& content) = 0;

  /**
   * @brief Get a value from the configuration
   * @param key Dot-separated path to the configuration value (e.g.
   * "database.host")
   * @return ConfigValue containing the requested value
   * @throws std::runtime_error if the key doesn't exist or value type is not
   * supported
   */
  virtual ConfigValue Get(const std::string& key) const = 0;

  /**
   * @brief Set a value in the configuration
   * @param key Dot-separated path where to store the value (e.g.
   * "database.port")
   * @param value ConfigValue to store
   * @throws std::runtime_error if the key path is invalid or value type is not
   * supported
   */
  virtual void Set(const std::string& key, const ConfigValue& value) = 0;

  /**
   * @brief Check if a key exists in the configuration
   * @param key Dot-separated path to check (e.g. "database.host")
   * @return true if the key exists and contains a valid value, false otherwise
   */
  virtual bool HasKey(const std::string& key) const = 0;

  /**
   * @brief Save configuration to a file
   * @param filename Path where to save the configuration
   * @return true if saving was successful, false otherwise
   * @throws std::runtime_error if the file cannot be created or written to
   */
  virtual bool SaveToFile(const std::string& filename) const = 0;

  /**
   * @brief Get configuration as a string
   * @return String representation of the entire configuration in the backend's
   * native format
   * @throws std::runtime_error if the configuration cannot be serialized
   */
  virtual std::string ToString() const = 0;

  /**
   * @brief Clear all configuration data
   *
   * Removes all key-value pairs from the configuration,
   * returning it to an empty state. After calling this method,
   * the configuration will be equivalent to a newly constructed instance.
   */
  virtual void Clear() = 0;
};

}  // namespace config
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_CONFIG_ICONFIG_H_