#pragma once

#include <concepts>

#include "BackendTraits.hpp"
#include "NonCopyable.hpp"
#include "common/utils/config/ConfigValue.hpp"

namespace metada::framework {

/**
 * @brief Concept defining requirements for a configuration backend
 */
template <typename T>
concept ConfigBackendType =
    requires(T t, const std::string& filename, const ConfigValue& value) {
      // Required methods
      { t.LoadFromFile(filename) } -> std::same_as<bool>;
      { t.LoadFromString(filename) } -> std::same_as<bool>;
      { t.Get(filename) } -> std::same_as<ConfigValue>;
      { t.Set(filename, value) } -> std::same_as<void>;
      { t.HasKey(filename) } -> std::same_as<bool>;
      { t.SaveToFile(filename) } -> std::same_as<bool>;
      { t.ToString() } -> std::same_as<std::string>;
      { t.Clear() } -> std::same_as<void>;
    };

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
template <typename BackendTag>
  requires ConfigBackendType<
      typename traits::BackendTraits<BackendTag>::ConfigBackend>
class Config : public NonCopyable {
 public:
  using ConfigBackend =
      typename traits::BackendTraits<BackendTag>::ConfigBackend;

  /** @brief Disabled default constructor */
  Config() = delete;

  /** @brief Default destructor */
  ~Config() = default;

  /** @brief Disabled copy constructor */
  Config(const Config&) = delete;

  /** @brief Disabled copy assignment */
  Config& operator=(const Config&) = delete;

  /**
   * @brief Move constructor - explicitly defined for compatibility with mock
   * objects
   */
  explicit Config(Config&& other) noexcept
      : backend_(std::move(other.backend_)) {}

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
   * @brief Constructor that loads configuration from a file
   * @param filename Path to the configuration file
   * @throws std::runtime_error If loading fails
   */
  explicit Config(const std::string& filename) {
    if (!backend_.LoadFromFile(filename)) {
      throw std::runtime_error("Failed to load configuration from file: " +
                               filename);
    }
  }

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

 private:
  ConfigBackend backend_;  ///< Instance of the configuration backend
};

}  // namespace metada::framework
