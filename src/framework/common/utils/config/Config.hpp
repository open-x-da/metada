#pragma once

#include "../NonCopyable.hpp"
#include "ConfigValue.hpp"
#include "IConfig.hpp"

namespace metada::framework {

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
 * Example usage:
 * @code
 * Config<YamlConfig> config;
 * config.LoadFromFile("config.yaml");
 * auto host = config.Get("database.host", "localhost");
 * auto port = config.Get("database.port", 5432);
 * config.Set("database.username", "admin");
 * config.SaveToFile("config.yaml");
 * @endcode
 *
 * Key features:
 * - Hierarchical configuration structure using dot notation
 * - Type-safe value access with default fallbacks
 * - File and string-based loading/saving
 * - Backend-agnostic interface
 * - Exception safety through Get() vs GetUnsafe()
 *
 * Supported value types:
 * - Boolean
 * - Integer
 * - Double
 * - String
 * - Arrays of the above types
 *
 * @tparam Backend The configuration backend type that implements IConfig
 * @see IConfig Base interface class for configuration backends
 */
template <typename Backend>
class Config : public NonCopyable {
 private:
  Backend backend_;  ///< Instance of the configuration backend

 public:
  /** @brief Default constructor */
  Config() = default;

  /**
   * @brief Move constructor - explicitly defined for compatibility with mock
   * objects
   */
  Config(Config&& other) noexcept : backend_(std::move(other.backend_)) {}

  /**
   * @brief Move assignment - explicitly defined for compatibility with mock
   * objects
   */
  Config& operator=(Config&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
    }
    return *this;
  }

  /**
   * @brief Get direct access to the backend instance
   * @return Reference to the backend instance
   */
  Backend& backend() { return backend_; }

  /**
   * @brief Get const access to the backend instance
   * @return Const reference to the backend instance
   */
  const Backend& backend() const { return backend_; }

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
   * @brief Get a value from the configuration with a default fallback (const
   * version)
   */
  ConfigValue Get(const std::string& key,
                  const ConfigValue& default_value = ConfigValue()) const {
    try {
      return backend_.Get(key);
    } catch (...) {
      return default_value;
    }
  }

  /**
   * @brief Get a value from the configuration without catching exceptions
   */
  ConfigValue GetUnsafe(const std::string& key) { return backend_.Get(key); }

  /**
   * @brief Get a value from the configuration without catching exceptions
   * (const version)
   */
  ConfigValue GetUnsafe(const std::string& key) const {
    return backend_.Get(key);
  }

  /**
   * @brief Set a value in the configuration
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

}  // namespace metada::framework
