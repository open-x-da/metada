#pragma once

#include <yaml-cpp/yaml.h>

#include <memory>
#include <string>

#include "IConfig.hpp"

namespace metada::backends::common::utils::config::yaml {

/**
 * @brief YAML configuration backend implementation
 * 
 * This class implements the IConfig interface using yaml-cpp as the
 * underlying storage format. It provides functionality to load, access, modify
 * and save configuration data in YAML format.
 *
 * The configuration data is stored in a hierarchical structure where keys use
 * dot notation to access nested values. For example, "database.host" would
 * access the "host" field within the "database" object.
 *
 * Example usage:
 * @code
 * YamlConfig config;
 * config.LoadFromFile("config.yaml");
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
 * @see IConfig Base interface class
 * @see YAML::Node YAML-CPP library used for implementation
 */
class YamlConfig : public framework::tools::config::IConfig {
 public:
  /** @brief Default constructor creating an empty configuration */
  YamlConfig() = default;

  /** @brief Virtual destructor */
  ~YamlConfig() override = default;

  /**
   * @brief Load configuration from a YAML file
   * @param filename Path to the YAML configuration file
   * @return true if loading was successful, false otherwise
   * @throws std::runtime_error if file cannot be opened or contains invalid
   * YAML
   */
  bool LoadFromFile(const std::string& filename) override;

  /**
   * @brief Load configuration from a YAML string
   * @param content String containing YAML configuration data
   * @return true if loading was successful, false otherwise
   * @throws std::runtime_error if content contains invalid YAML
   */
  bool LoadFromString(const std::string& content) override;

  /**
   * @brief Get a value from the configuration
   * @param key Dot-separated path to the configuration value (e.g.
   * "database.host")
   * @return ConfigValue containing the requested value
   * @throws std::runtime_error if the key doesn't exist or value type is not
   * supported
   */
  framework::tools::config::ConfigValue Get(
      const std::string& key) const override;

  /**
   * @brief Set a value in the configuration
   * @param key Dot-separated path where to store the value (e.g.
   * "database.port")
   * @param value ConfigValue to store
   * @throws std::runtime_error if the key path is invalid or value type is not
   * supported
   */
  void Set(const std::string& key,
           const framework::tools::config::ConfigValue& value) override;

  /**
   * @brief Check if a key exists in the configuration
   * @param key Dot-separated path to check (e.g. "database.password")
   * @return true if the key exists, false otherwise
   */
  bool HasKey(const std::string& key) const override;

  /**
   * @brief Save configuration to a YAML file
   * @param filename Path where to save the configuration
   * @return true if saving was successful, false otherwise
   * @throws std::runtime_error if file cannot be created or written to
   */
  bool SaveToFile(const std::string& filename) const override;

  /**
   * @brief Get configuration as a YAML string
   * @return String containing the YAML representation of the configuration
   * @throws std::runtime_error if serialization fails
   */
  std::string ToString() const override;

  /**
   * @brief Clear all configuration data
   * Resets the configuration to an empty state
   */
  void Clear() override;

 private:
  /** @brief Root YAML node storing the configuration data */
  YAML::Node root_;

  /**
   * @brief Get a YAML node at the specified path
   * @param node YAML node to search in
   * @param key Dot-separated path to the node (e.g. "database.host")
   * @return YAML node at the specified path
   * @throws std::runtime_error if path is invalid or node doesn't exist
   */
  static YAML::Node GetNode(const YAML::Node& node, const std::string& key);

  /**
   * @brief Convert a YAML node to a ConfigValue
   * @param node YAML node to convert
   * @return ConfigValue containing the converted value
   * @throws std::runtime_error if node type is not supported
   */
  static framework::tools::config::ConfigValue NodeToConfigValue(
      const YAML::Node& node);

  /**
   * @brief Convert a ConfigValue to a YAML node
   * @param value ConfigValue to convert
   * @return YAML node containing the converted value
   * @throws std::runtime_error if value type is not supported
   */
  static YAML::Node ConfigValueToNode(
      const framework::tools::config::ConfigValue& value);
};

}  // namespace metada::backends::common::utils::config::yaml