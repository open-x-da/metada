#pragma once

#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "../utils.hpp"
#include "common/utils/config/ConfigValue.hpp"

namespace metada::backends::config {

using framework::ConfigValue;

/**
 * @brief JSON configuration backend implementation
 *
 * This class implements the configuration backend contract using nlohmann::json
 * as the underlying storage format. It provides functionality to load, access,
 * modify and save configuration data in JSON format.
 *
 * The configuration data is stored in a hierarchical structure where keys use
 * dot notation to access nested values. For example, "database.host" would
 * access the "host" field within the "database" object.
 *
 * Features:
 * - Loading from JSON files and strings
 * - Saving to JSON files
 * - Accessing nested values with dot notation
 * - Type-safe value retrieval and storage
 * - Support for all ConfigValue types
 * - Thread-safe operations
 *
 * Example usage:
 * @code
 * JsonConfig config("config.json");
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
 * - Float
 * - String
 * - Arrays of the above types
 *
 * This implementation satisfies the ConfigBackendType concept required by
 * the framework's Config template class.
 *
 * @see nlohmann::json JSON library used for implementation
 * @see framework::ConfigValue The variant type used to store configuration
 * values
 * @see framework::ConfigBackendType The concept this class satisfies
 */
class JsonConfig {
 public:
  /** @brief Default constructor */
  JsonConfig() : root_(nlohmann::json::object()) {}

  /** @brief Default destructor */
  ~JsonConfig() = default;

  /**
   * @brief Disable copy constructor
   */
  JsonConfig(const JsonConfig&) = delete;

  /**
   * @brief Disable copy assignment operator
   */
  JsonConfig& operator=(const JsonConfig&) = delete;

  /**
   * @brief Move constructor
   */
  JsonConfig(JsonConfig&& other) noexcept : root_(std::move(other.root_)) {}

  /**
   * @brief Move assignment operator
   */
  JsonConfig& operator=(JsonConfig&& other) noexcept {
    root_ = std::move(other.root_);
    return *this;
  }

  /**
   * @brief Constructor that loads configuration from a file
   * @param filename Path to the JSON configuration file
   * @throws std::runtime_error If loading fails
   */
  explicit JsonConfig(const std::string& filename) { LoadFromFile(filename); }

  /**
   * @brief Load configuration from a JSON file
   * @param filename Path to the JSON configuration file
   * @return true if loading was successful, false otherwise
   * @throws std::runtime_error if file cannot be opened or contains invalid
   * JSON
   */
  bool LoadFromFile(const std::string& filename);

  /**
   * @brief Load configuration from a JSON string
   * @param content String containing JSON configuration data
   * @return true if loading was successful, false otherwise
   * @throws std::runtime_error if content contains invalid JSON
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
   * @brief Save configuration to a JSON file
   * @param filename Path where to save the configuration
   * @return true if saving was successful, false otherwise
   * @throws std::runtime_error if file cannot be created or written to
   */
  bool SaveToFile(const std::string& filename) const;

  /**
   * @brief Get configuration as a JSON string
   * @return String containing the JSON representation of the configuration
   * @throws std::runtime_error if serialization fails
   */
  std::string ToString() const;

  /**
   * @brief Clear all configuration data
   * Resets the configuration to an empty state
   */
  void Clear();

  /**
   * @brief Create a new JsonConfig object representing a subsection
   *
   * @param key Dot-separated path to the subsection
   * @return A new JsonConfig object representing the subsection
   */
  JsonConfig CreateSubsection(const std::string& key) const {
    JsonConfig subsection;
    subsection.root_ = GetJsonRef(root_, key);
    return subsection;
  }

 private:
  /** @brief Root JSON node storing the configuration data */
  nlohmann::json root_;

  /**
   * @brief Get a reference to a JSON value at the specified path
   * @param j The JSON object to traverse
   * @param key Dot-separated path to the desired value
   * @return Reference to the JSON value at the specified path
   * @throws std::runtime_error if the path is invalid or contains invalid
   * characters
   */
  nlohmann::json& GetJsonRef(nlohmann::json& j, const std::string& key);

  /**
   * @brief Get a const reference to a JSON value at the specified path
   * @param j The JSON object to traverse
   * @param key Dot-separated path to the desired value
   * @return Const reference to the JSON value at the specified path
   * @throws std::runtime_error if the path is invalid or value doesn't
   * exist
   */
  const nlohmann::json& GetJsonRef(const nlohmann::json& j,
                                   const std::string& key) const;

  /**
   * @brief Convert a JSON value to a ConfigValue
   * @param j The JSON value to convert
   * @return ConfigValue containing the converted value
   * @throw std::runtime_error if the JSON value type is not supported
   */
  static ConfigValue JsonToConfigValue(const nlohmann::json& j);

  /**
   * @brief Convert a ConfigValue to a JSON value
   * @param value The ConfigValue to convert
   * @return nlohmann::json containing the converted value
   */
  static nlohmann::json ConfigValueToJson(const ConfigValue& value);
};

}  // namespace metada::backends::config