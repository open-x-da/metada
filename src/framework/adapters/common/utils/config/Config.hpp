#pragma once

#include <concepts>
#include <filesystem>

#include "BackendTraits.hpp"
#include "ConfigValue.hpp"
#include "NonCopyable.hpp"

namespace metada::framework {

/**
 * @brief Concept defining requirements for a configuration backend
 *
 * @details This concept enforces the contract that any configuration backend
 * must implement. It provides compile-time validation of the required methods
 * and their signatures, ensuring that backends can be properly used with the
 * Config class.
 *
 * The concept requires specific method signatures for loading, accessing,
 * modifying, and saving configuration data, as well as a constructor that takes
 * a filename.
 */
template <typename T>
concept ConfigBackendType =
    requires(T t, const std::string& filename, const ConfigValue& value) {
      // Required constructor
      { T(filename) } -> std::same_as<T>;
      // Required methods
      { t.LoadFromFile(filename) } -> std::same_as<bool>;
      { t.LoadFromString(filename) } -> std::same_as<bool>;
      { t.Get(filename) } -> std::same_as<ConfigValue>;
      { t.Set(filename, value) } -> std::same_as<void>;
      { t.HasKey(filename) } -> std::same_as<bool>;
      { t.SaveToFile(filename) } -> std::same_as<bool>;
      { t.ToString() } -> std::same_as<std::string>;
      { t.Clear() } -> std::same_as<void>;
      { t.CreateSubsection(filename) } -> std::same_as<T>;
    };

/**
 * @brief Main configuration class template providing a generic interface to
 * configuration backends
 *
 * @details This class template provides a static interface for loading,
 * accessing, modifying and saving configuration data using a backend specified
 * by the BackendTag template parameter. The backend must satisfy the
 * ConfigBackendType concept.
 *
 * The configuration data is stored in a hierarchical structure where keys use
 * dot notation to access nested values. For example, "database.host" would
 * access the "host" field within the "database" object.
 *
 * @par Example usage:
 * @code
 * // Create application context with configuration
 * auto context = ApplicationContext<BackendTag>(argv[0], argv[1]);
 * auto& config = context.getConfig();
 *
 * // Access configuration values with defaults
 * auto host = config.Get("database.host", "localhost");
 * auto port = config.Get("database.port", 5432);
 *
 * // Modify configuration
 * config.Set("database.username", "admin");
 * config.SaveToFile("updated_config.yaml");
 *
 * // Work with subsections
 * auto dbConfig = config.GetSubsection("database");
 * auto timeout = dbConfig.Get("timeout", 30);
 * @endcode
 *
 * @par Key features:
 * - Hierarchical configuration structure using dot notation
 * - Type-safe value access with default fallbacks
 * - File and string-based loading/saving
 * - Backend-agnostic interface using C++20 concepts
 * - Exception safety through Get() with default fallbacks
 * - Support for configuration subsections
 *
 * @par Supported value types:
 * The ConfigValue class supports various types including:
 * - Boolean values (true/false)
 * - Integer numbers
 * - Floating-point numbers
 * - Strings
 * - Arrays/sequences of values
 * - Maps/objects (nested configurations)
 * - Null values
 *
 * @tparam BackendTag The tag type that identifies the configuration backend via
 * BackendTraits
 */
template <typename BackendTag>
  requires ConfigBackendType<
      typename traits::BackendTraits<BackendTag>::ConfigBackend>
class Config : public NonCopyable {
 public:
  using ConfigBackend =
      typename traits::BackendTraits<BackendTag>::ConfigBackend;

  /**
   * @brief Disabled default constructor
   *
   * @details Configuration must be initialized with a source.
   */
  Config() = delete;

  /** @brief Default destructor */
  ~Config() = default;

  /**
   * @brief Disabled copy constructor
   *
   * @details Configuration instances are not intended to be copied.
   */
  Config(const Config&) = delete;

  /**
   * @brief Disabled copy assignment
   *
   * @details Configuration instances are not intended to be copied.
   */
  Config& operator=(const Config&) = delete;

  /**
   * @brief Move constructor
   *
   * @details Explicitly defined for compatibility with mock objects and
   * testing. Transfers ownership of the backend from another Config instance.
   *
   * @param other The Config instance to move from
   */
  explicit Config(Config&& other) noexcept
      : backend_(std::move(other.backend_)) {}

  /**
   * @brief Move assignment operator
   *
   * @details Explicitly defined for compatibility with mock objects and
   * testing. Transfers ownership of the backend from another Config instance.
   *
   * @param other The Config instance to move from
   * @return Reference to this Config instance
   */
  Config& operator=(Config&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
    }
    return *this;
  }

  /**
   * @brief Constructor that loads configuration from a file
   *
   * @details Initializes the configuration by loading from the specified file.
   * Verifies that the file exists before attempting to load it.
   *
   * @param filename Path to the configuration file
   * @throws std::runtime_error If file doesn't exist or loading fails
   */
  explicit Config(const std::string& filename)
      : backend_([&filename]() {
          if (!std::filesystem::exists(filename)) {
            throw std::runtime_error("Configuration file does not exist: " +
                                     filename);
          }
          return ConfigBackend(filename);
        }()) {}

  /**
   * @brief Constructor that takes an existing backend
   *
   * @details Creates a Config object from an existing backend instance.
   * This is used internally for creating subsections.
   *
   * @param backend The backend instance to use
   */
  explicit Config(ConfigBackend&& backend) : backend_(std::move(backend)) {}

  /**
   * @brief Get direct access to the backend instance
   *
   * @details Provides mutable access to the underlying configuration backend.
   * This can be used for backend-specific operations not exposed by the Config
   * interface.
   *
   * @return Reference to the backend instance
   */
  ConfigBackend& backend() { return backend_; }

  /**
   * @brief Get const access to the backend instance
   *
   * @details Provides read-only access to the underlying configuration backend.
   * This can be used for backend-specific operations not exposed by the Config
   * interface.
   *
   * @return Const reference to the backend instance
   */
  const ConfigBackend& backend() const { return backend_; }

  /**
   * @brief Get a value from the configuration with a default fallback
   *
   * @details Retrieves a value from the configuration using the specified key.
   * If the key doesn't exist or an error occurs, returns the provided default
   * value. This method provides exception safety when accessing configuration
   * values.
   *
   * @param key Dot-separated path to the configuration value
   * @param default_value Value to return if the key doesn't exist or an error
   * occurs
   * @return The configuration value or the default value
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
   *
   * @details Retrieves a value from the configuration using the specified key.
   * If the key doesn't exist or an error occurs, returns the provided default
   * value. This method provides exception safety when accessing configuration
   * values.
   *
   * @param key Dot-separated path to the configuration value
   * @param default_value Value to return if the key doesn't exist or an error
   * occurs
   * @return The configuration value or the default value
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
   *
   * @details Updates or creates a configuration value at the specified key
   * path. If intermediate nodes in the path don't exist, they will be created.
   *
   * @param key Dot-separated path where to set the value
   * @param value The value to set
   */
  void Set(const std::string& key, const ConfigValue& value) {
    backend_.Set(key, value);
  }

  /**
   * @brief Check if a key exists in the configuration
   *
   * @details Verifies whether the specified key path exists in the
   * configuration.
   *
   * @param key Dot-separated path to check
   * @return true if the key exists, false otherwise
   */
  bool HasKey(const std::string& key) const { return backend_.HasKey(key); }

  /**
   * @brief Save configuration to a file
   *
   * @details Writes the current configuration state to the specified file.
   * The format of the file depends on the backend implementation.
   *
   * @param filename Path where to save the configuration
   * @return true if saving was successful, false otherwise
   */
  bool SaveToFile(const std::string& filename) const {
    return backend_.SaveToFile(filename);
  }

  /**
   * @brief Convert configuration to string representation
   *
   * @details Serializes the entire configuration to a string.
   * The format of the string depends on the backend implementation.
   *
   * @return String containing the configuration data
   */
  std::string ToString() const { return backend_.ToString(); }

  /**
   * @brief Clear all configuration data
   *
   * @details Removes all key-value pairs from the configuration,
   * resulting in an empty configuration state.
   */
  void Clear() { backend_.Clear(); }

  /**
   * @brief Get a subsection of the configuration as a new Config object
   *
   * @details Creates a new Config object that represents a subsection of the
   * current configuration. The new Config object uses the same backend as the
   * current one, but its root is set to the specified subsection.
   *
   * @param key Dot-separated path to the subsection
   * @return A new Config object representing the subsection
   * @throws std::runtime_error if the subsection doesn't exist
   */
  Config GetSubsection(const std::string& key) const {
    // Get the subsection as a ConfigValue
    ConfigValue subsectionValue = Get(key);

    // Check if the subsection exists and is a map
    if (!subsectionValue.isMap()) {
      throw std::runtime_error("Subsection does not exist or is not a map: " +
                               key);
    }

    // Create a new Config object by moving the backend subsection
    return Config(std::move(backend_.CreateSubsection(key)));
  }

  /**
   * @brief Check if a subsection exists in the configuration
   *
   * @details Verifies whether the specified key path exists in the
   * configuration and contains a map (object).
   *
   * @param key Dot-separated path to check
   * @return true if the subsection exists and contains a map, false otherwise
   */
  bool HasSubsection(const std::string& key) const {
    auto value = Get(key);
    return value.isMap();
  }

 private:
  ConfigBackend backend_;  ///< Instance of the configuration backend
};

}  // namespace metada::framework
