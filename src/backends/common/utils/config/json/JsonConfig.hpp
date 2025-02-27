#pragma once

#include <memory>
#include <nlohmann/json.hpp>
#include <string>

#include "utils/config/ConfigValue.hpp"
#include "utils/config/IConfig.hpp"

namespace metada::backends::config {

using framework::ConfigValue;
using framework::IConfig;
/**
 * @brief JSON configuration backend implementation
 *
 * This class implements the IConfig interface using nlohmann::json as the
 * underlying storage format. It provides functionality to load, access, modify
 * and save configuration data in JSON format.
 *
 * The configuration data is stored in a hierarchical structure where keys use
 * dot notation to access nested values. For example, "database.host" would
 * access the "host" field within the "database" object.
 *
 * Example usage:
 * @code
 * JsonConfig config;
 * config.LoadFromFile("config.json");
 *
 * // Get values
 * auto host = config.Get("database.host").AsString();
 * auto port = config.Get("database.port").AsInt();
 *
 * // Set values
 * config.Set("database.user", ConfigValue("admin"));
 * config.Set("database.timeout", ConfigValue(30));
 *
 * config.SaveToFile("updated_config.json");
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
 * @see nlohmann::json JSON library used for implementation
 */
class JsonConfig : public IConfig {
 public:
  /** @brief Default constructor creating an empty configuration */
  JsonConfig() = default;

  /** @brief Virtual destructor */
  ~JsonConfig() override = default;

  /**
   * @brief Load configuration from a JSON file
   * @param filename Path to the JSON configuration file
   * @return true if loading was successful, false otherwise
   * @throws std::runtime_error if file cannot be opened or contains invalid
   * JSON
   */
  bool LoadFromFile(const std::string& filename) override;

  /**
   * @brief Load configuration from a JSON string
   * @param content String containing JSON configuration data
   * @return true if loading was successful, false otherwise
   * @throws std::runtime_error if content contains invalid JSON
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
  ConfigValue Get(const std::string& key) const override;

  /**
   * @brief Set a value in the configuration
   * @param key Dot-separated path where to store the value (e.g.
   * "database.port")
   * @param value ConfigValue to store
   * @throws std::runtime_error if the key path is invalid or value type is not
   * supported
   */
  void Set(const std::string& key, const ConfigValue& value) override;

  /**
   * @brief Check if a key exists in the configuration
   * @param key Dot-separated path to check (e.g. "database.password")
   * @return true if the key exists, false otherwise
   */
  bool HasKey(const std::string& key) const override;

  /**
   * @brief Save configuration to a JSON file
   * @param filename Path where to save the configuration
   * @return true if saving was successful, false otherwise
   * @throws std::runtime_error if file cannot be created or written to
   */
  bool SaveToFile(const std::string& filename) const override;

  /**
   * @brief Get configuration as a JSON string
   * @return String containing the JSON representation of the configuration
   * @throws std::runtime_error if serialization fails
   */
  std::string ToString() const override;

  /**
   * @brief Clear all configuration data
   * Resets the configuration to an empty state
   */
  void Clear() override;

 private:
  /** @brief Root JSON object storing the configuration data */
  nlohmann::json root_;

  /**
   * @brief Get a reference to a JSON value at the specified path
   * @param j JSON object to search in
   * @param key Dot-separated path to the value (e.g. "database.host")
   * @return Reference to the JSON value
   * @throws std::runtime_error if path is invalid or value doesn't exist
   */
  static nlohmann::json& GetJsonRef(nlohmann::json& j, const std::string& key);

  /**
   * @brief Get a const reference to a JSON value at the specified path
   * @param j JSON object to search in
   * @param key Dot-separated path to the value (e.g. "database.host")
   * @return Const reference to the JSON value
   * @throws std::runtime_error if path is invalid or value doesn't exist
   */
  static const nlohmann::json& GetJsonRef(const nlohmann::json& j,
                                          const std::string& key);

  /**
   * @brief Convert a JSON value to a ConfigValue
   * @param j JSON value to convert
   * @return Equivalent ConfigValue
   * @throws std::runtime_error if JSON value type is not supported
   */
  static ConfigValue JsonToConfigValue(const nlohmann::json& j);

  /**
   * @brief Convert a ConfigValue to a JSON value
   * @param value ConfigValue to convert
   * @return Equivalent JSON value
   * @throws std::runtime_error if ConfigValue type is not supported
   */
  static nlohmann::json ConfigValueToJson(const ConfigValue& value);
};

}  // namespace metada::backends::config