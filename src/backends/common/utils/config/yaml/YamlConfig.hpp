#pragma once

#ifdef _WIN32
#define YAML_CPP_STATIC_DEFINE
#endif
#include <yaml-cpp/yaml.h>

#include <string>
#include <vector>

#include "../utils.hpp"
#include "ConfigValue.hpp"

namespace metada::backends::config {

using framework::ConfigValue;

/**
 * @brief YAML configuration backend implementation
 *
 * This class implements the configuration backend contract using yaml-cpp as
 * the underlying storage format. It provides functionality to load, access,
 * modify and save configuration data in YAML format.
 *
 * The configuration data is stored in a hierarchical structure where keys use
 * dot notation to access nested values. For example, "database.host" would
 * access the "host" field within the "database" object.
 *
 * Features:
 * - Loading from YAML files and strings
 * - Saving to YAML files
 * - Accessing nested values with dot notation
 * - Type-safe value retrieval and storage
 * - Support for all ConfigValue types
 * - Thread-safe operations
 *
 * Example usage:
 * @code
 * YamlConfig config("config.yaml");
 *
 * // Get values
 * auto host = config.Get("database.host").AsString();
 * auto port = config.Get("database.port").AsInt();
 *
 * // Set values
 * config.Set("database.user", ConfigValue("admin"));
 * config.Set("database.timeout", ConfigValue(30));
 *
 * config.SaveToFile("updated_config.yaml");
 * @endcode
 *
 * Supported value types:
 * - Boolean
 * - Integer
 * - Double
 * - String
 * - Arrays of the above types
 *
 * This implementation satisfies the ConfigBackendType concept required by
 * the framework's Config template class.
 *
 * @see YAML::Node YAML-CPP library used for implementation
 * @see framework::ConfigValue The variant type used to store configuration
 * values
 * @see framework::ConfigBackendType The concept this class satisfies
 */
class YamlConfig {
 public:
  /** @brief Default constructor */
  YamlConfig() : root_(YAML::Node()) {}

  /** @brief Default destructor */
  ~YamlConfig() = default;

  /** @brief Disable copy constructor */
  YamlConfig(const YamlConfig&) = delete;

  /** @brief Disable copy assignment operator */
  YamlConfig& operator=(const YamlConfig&) = delete;

  /**
   * @brief Move constructor
   */
  YamlConfig(YamlConfig&& other) noexcept : root_(std::move(other.root_)) {}

  /**
   * @brief Move assignment operator
   */
  YamlConfig& operator=(YamlConfig&& other) noexcept {
    root_ = std::move(other.root_);
    return *this;
  }

  /**
   * @brief Constructor that loads configuration from a file
   * @param filename Path to the YAML configuration file
   * @throws std::runtime_error If loading fails
   */
  explicit YamlConfig(const std::string& filename) { LoadFromFile(filename); }

  /**
   * @brief Constructor that loads configuration from a ConfigMap
   * @param map ConfigMap containing the configuration data
   */
  explicit YamlConfig(const framework::ConfigMap& map);

  /**
   * @brief Load configuration from a YAML file
   * @param filename Path to the YAML configuration file
   * @return true if loading was successful, false otherwise
   * @throws std::runtime_error if file cannot be opened or contains invalid
   * YAML
   */
  bool LoadFromFile(const std::string& filename);

  /**
   * @brief Load configuration from a YAML string
   * @param content String containing YAML configuration data
   * @return true if loading was successful, false otherwise
   * @throws std::runtime_error if content contains invalid YAML
   */
  bool LoadFromString(const std::string& content);

  /**
   * @brief Get a value from the configuration
   * @param key Dot-separated path to the configuration value (e.g.
   * "database.host")
   * @return ConfigValue containing the requested value
   * @throws std::runtime_error if the key doesn't exist or value type is not
   * supported
   */
  ConfigValue Get(const std::string& key) const;

  /**
   * @brief Set a value in the configuration
   * @param key Dot-separated path where to store the value (e.g.
   * "database.port")
   * @param value ConfigValue to store
   * @throws std::runtime_error if the key path is invalid or value type is not
   * supported
   */
  void Set(const std::string& key, const ConfigValue& value);

  /**
   * @brief Check if a key exists in the configuration
   * @param key Dot-separated path to check (e.g. "database.password")
   * @return true if the key exists, false otherwise
   */
  bool HasKey(const std::string& key) const;

  /**
   * @brief Save configuration to a YAML file
   * @param filename Path where to save the configuration
   * @return true if saving was successful, false otherwise
   * @throws std::runtime_error if file cannot be created or written to
   */
  bool SaveToFile(const std::string& filename) const;

  /**
   * @brief Get configuration as a YAML string
   * @return String containing the YAML representation of the configuration
   * @throws std::runtime_error if serialization fails
   */
  std::string ToString() const;

  /**
   * @brief Clear all configuration data
   * Resets the configuration to an empty state
   */
  void Clear();

  /**
   * @brief Stream output operator for YamlConfig
   *
   * Outputs the YAML configuration as a formatted string to the given stream.
   * This is useful for debugging and logging purposes.
   *
   * @param os Output stream to write to
   * @param config YamlConfig object to output
   * @return Reference to the output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const YamlConfig& config) {
    os << config.ToString();
    return os;
  }

  /**
   * @brief Create a new YamlConfig object representing a subsection
   *
   * @param key Dot-separated path to the subsection
   * @return A new YamlConfig object representing the subsection
   */
  YamlConfig CreateSubsection(const std::string& key) const {
    YamlConfig subsection;
    subsection.root_ = GetYamlNodeConst(root_, key);
    return subsection;
  }

 private:
  /** @brief Root YAML node storing the configuration data */
  YAML::Node root_;

  /**
   * @brief Get a YAML node at the specified path
   * @param node The YAML node to traverse
   * @param key Dot-separated path to the desired value (e.g. "database.host")
   * @return YAML node at the specified path
   * @throw std::runtime_error if the path is invalid or value doesn't exist
   */
  static YAML::Node GetYamlNode(YAML::Node node, const std::string& key);

  /**
   * @brief Get a const YAML node at the specified path
   * @param node The YAML node to traverse
   * @param key Dot-separated path to the desired value (e.g. "database.host")
   * @return Const YAML node at the specified path
   * @throw std::runtime_error if the path is invalid or value doesn't exist
   */
  static const YAML::Node GetYamlNodeConst(const YAML::Node& node,
                                           const std::string& key);

  /**
   * @brief Convert a YAML value to a ConfigValue
   * @param node The YAML value to convert
   * @return ConfigValue containing the converted value
   * @throw std::runtime_error if the YAML value type is not supported
   */
  static ConfigValue YamlToConfigValue(const YAML::Node& node);

  /**
   * @brief Convert a ConfigValue to a YAML value
   * @param value The ConfigValue to convert
   * @return YAML::Node containing the converted value
   */
  static YAML::Node ConfigValueToYaml(const ConfigValue& value);
};

}  // namespace metada::backends::config