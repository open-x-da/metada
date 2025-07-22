#pragma once

#include <algorithm>
#include <iomanip>
#include <map>
#include <ostream>
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

  // Constructor for vectors of ConfigValue
  ConfigValue(const std::vector<ConfigValue>& value) : value_(value) {}

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
  bool isVectorConfigValue() const {
    return std::holds_alternative<std::vector<ConfigValue>>(value_);
  }
  bool isMap() const { return std::holds_alternative<ConfigMap>(value_); }
  bool isVectorMap() const {
    if (!isVectorConfigValue()) return false;
    const auto& vec = asVectorConfigValue();
    return std::all_of(vec.begin(), vec.end(),
                       [](const ConfigValue& v) { return v.isMap(); });
  }

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
  const std::vector<ConfigValue>& asVectorConfigValue() const {
    return std::get<std::vector<ConfigValue>>(value_);
  }
  const ConfigMap& asMap() const { return std::get<ConfigMap>(value_); }
  std::vector<ConfigMap> asVectorMap() const {
    const auto& vec = asVectorConfigValue();
    std::vector<ConfigMap> result;
    result.reserve(vec.size());
    for (const auto& v : vec) result.push_back(v.asMap());
    return result;
  }

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

  friend std::ostream& operator<<(std::ostream& os, const ConfigValue& value) {
    if (value.isNull()) {
      os << "null";
    } else if (value.isBool()) {
      os << (value.asBool() ? "true" : "false");
    } else if (value.isInt()) {
      os << value.asInt();
    } else if (value.isFloat()) {
      os << std::fixed << std::setprecision(1) << value.asFloat();
    } else if (value.isString()) {
      os << value.asString();
    } else if (value.isVectorBool()) {
      os << "[";
      const auto& vec = value.asVectorBool();
      for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) os << ", ";
        os << (vec[i] ? "true" : "false");
      }
      os << "]";
    } else if (value.isVectorInt()) {
      os << "[";
      const auto& vec = value.asVectorInt();
      for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) os << ", ";
        os << vec[i];
      }
      os << "]";
    } else if (value.isVectorFloat()) {
      os << "[";
      const auto& vec = value.asVectorFloat();
      for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) os << ", ";
        os << std::fixed << std::setprecision(1) << vec[i];
      }
      os << "]";
    } else if (value.isVectorString()) {
      os << "[";
      const auto& vec = value.asVectorString();
      for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) os << ", ";
        os << '"' << vec[i] << '"';
      }
      os << "]";
    } else if (value.isMap()) {
      os << "{";
      const auto& map = value.asMap();
      bool first = true;
      for (const auto& [k, v] : map) {
        if (!first) os << ", ";
        os << '"' << k << ": " << v;
        first = false;
      }
      os << "}";
    } else if (value.isVectorConfigValue()) {
      os << "[";
      const auto& vec = value.asVectorConfigValue();
      for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) os << ", ";
        os << vec[i];
      }
      os << "]";
    }
    return os;
  }

  // Public getter for the variant type index
  size_t typeIndex() const { return value_.index(); }

  // Public getter for the underlying variant (const)
  const auto& constVariant() const { return value_; }

 private:
  std::variant<std::monostate, bool, int, float, std::string, std::vector<bool>,
               std::vector<int>, std::vector<float>, std::vector<std::string>,
               std::vector<ConfigValue>, ConfigMap>
      value_;
};

}  // namespace metada::framework
