#pragma once

#include <string>
#include <vector>

#include "ConfigValue.hpp"

namespace metada::framework {

/**
 * @brief Abstract interface for configuration backend implementations
 *
 * This interface defines the contract that all configuration backends must
 * implement. It provides a unified API for loading, accessing, modifying and
 * persisting configuration data in a backend-agnostic way.
 *
 * Key features:
 * - Hierarchical configuration using dot notation (e.g. "database.host")
 * - Type-safe value storage and retrieval via ConfigValue
 * - File and string-based loading/saving
 * - Key existence checking
 * - Full configuration clearing
 *
 * Example usage:
 * @code
 * class YamlConfig : public IConfig {
 *   bool LoadFromFile(const std::string& filename) override;
 *   ConfigValue Get(const std::string& key) const override;
 *   void Set(const std::string& key, const ConfigValue& value) override;
 *   // ... implement other methods
 * };
 * @endcode
 *
 * @see ConfigValue Type-safe variant for configuration values
 * @see Config Main configuration class template using this interface
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

}  // namespace metada::framework