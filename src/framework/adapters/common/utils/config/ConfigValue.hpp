#pragma once

#include <map>
#include <string>
#include <variant>
#include <vector>

namespace metada::framework {

// Forward declaration
class ConfigValue;

// Define the map type
using ConfigMap = std::map<std::string, ConfigValue>;

/**
 * @brief Configuration value type supporting multiple data types including
 * nested objects
 *
 * @details This class provides a flexible container that can hold
 * different types of configuration values, including nested objects represented
 * as maps. It enables type-safe storage and access of configuration data while
 * maintaining a consistent interface.
 *
 * The ConfigValue is used throughout the configuration system to represent
 * values retrieved from or stored in configuration backends. It's designed to
 * work with the Config class template and various backend implementations.
 *
 * Supported value types:
 * - bool: Boolean values (true/false)
 * - int: Integer numbers
 * - float: Floating point numbers (note: using float, not double)
 * - std::string: Text strings
 * - std::vector<bool>: Arrays of boolean values
 * - std::vector<int>: Arrays of integer numbers
 * - std::vector<float>: Arrays of floating point numbers
 * - std::vector<std::string>: Arrays of text strings
 * - ConfigMap: Maps of string keys to ConfigValue values (for nested objects)
 *
 * Example usage:
 * @code
 * ConfigValue value = config.Get("server.port");
 * if (value.isInt()) {
 *   int port = value.asInt();
 *   // Use port value
 * }
 * @endcode
 */
class ConfigValue {
 public:
  // Default constructor
  ConfigValue() : value_(std::monostate{}) {}

  // Constructors for primitive types
  ConfigValue(bool value) : value_(value) {}
  ConfigValue(int value) : value_(value) {}
  ConfigValue(float value) : value_(value) {}
  ConfigValue(const std::string& value) : value_(value) {}
  ConfigValue(const char* value) : value_(std::string(value)) {}

  // Constructors for vectors
  ConfigValue(const std::vector<bool>& value) : value_(value) {}
  ConfigValue(const std::vector<int>& value) : value_(value) {}
  ConfigValue(const std::vector<float>& value) : value_(value) {}
  ConfigValue(const std::vector<std::string>& value) : value_(value) {}

  // Constructor for maps
  ConfigValue(const ConfigMap& value) : value_(value) {}

  // Type checking methods
  bool isNull() const { return std::holds_alternative<std::monostate>(value_); }
  bool isBool() const { return std::holds_alternative<bool>(value_); }
  bool isInt() const { return std::holds_alternative<int>(value_); }
  bool isFloat() const { return std::holds_alternative<float>(value_); }
  bool isString() const { return std::holds_alternative<std::string>(value_); }
  bool isVectorBool() const {
    return std::holds_alternative<std::vector<bool>>(value_);
  }
  bool isVectorInt() const {
    return std::holds_alternative<std::vector<int>>(value_);
  }
  bool isVectorFloat() const {
    return std::holds_alternative<std::vector<float>>(value_);
  }
  bool isVectorString() const {
    return std::holds_alternative<std::vector<std::string>>(value_);
  }
  bool isMap() const { return std::holds_alternative<ConfigMap>(value_); }

  // Value access methods
  bool asBool() const { return std::get<bool>(value_); }
  int asInt() const { return std::get<int>(value_); }
  float asFloat() const { return std::get<float>(value_); }
  const std::string& asString() const { return std::get<std::string>(value_); }
  const std::vector<bool>& asVectorBool() const {
    return std::get<std::vector<bool>>(value_);
  }
  const std::vector<int>& asVectorInt() const {
    return std::get<std::vector<int>>(value_);
  }
  const std::vector<float>& asVectorFloat() const {
    return std::get<std::vector<float>>(value_);
  }
  const std::vector<std::string>& asVectorString() const {
    return std::get<std::vector<std::string>>(value_);
  }
  const ConfigMap& asMap() const { return std::get<ConfigMap>(value_); }

  // Map access methods
  bool hasKey(const std::string& key) const {
    if (!isMap()) return false;
    return asMap().find(key) != asMap().end();
  }

  ConfigValue get(const std::string& key) const {
    if (!isMap()) return ConfigValue();
    auto it = asMap().find(key);
    if (it == asMap().end()) return ConfigValue();
    return it->second;
  }

  // Get a value with a default if the key doesn't exist
  ConfigValue get(const std::string& key,
                  const ConfigValue& default_value) const {
    if (!isMap()) return default_value;
    auto it = asMap().find(key);
    if (it == asMap().end()) return default_value;
    return it->second;
  }

  // Set a value in the map
  void set(const std::string& key, const ConfigValue& value) {
    if (!isMap()) {
      value_ = ConfigMap();
    }
    std::get<ConfigMap>(value_)[key] = value;
  }

  bool operator==(const ConfigValue& other) const {
    return value_ == other.value_;
  }

  bool operator!=(const ConfigValue& other) const { return !(*this == other); }

 private:
  std::variant<std::monostate, bool, int, float, std::string, std::vector<bool>,
               std::vector<int>, std::vector<float>, std::vector<std::string>,
               ConfigMap>
      value_;
};

}  // namespace metada::framework
