#include "JsonConfig.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace metada::backends::config {

using framework::ConfigMap;
using framework::ConfigValue;

/**
 * @brief Load configuration from a JSON file
 *
 * Attempts to load and parse a JSON configuration file. The file must contain
 * valid JSON syntax.
 *
 * @param filename Path to the JSON configuration file
 * @return true if loading was successful, false otherwise
 * @throws std::runtime_error if file cannot be opened or contains invalid JSON
 */
bool JsonConfig::LoadFromFile(const std::string& filename) {
  try {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Failed to open file: " + filename);
    }
    file >> root_;
    return true;
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to load JSON from file: " +
                             std::string(e.what()));
  }
}

/**
 * @brief Load configuration from a JSON string
 *
 * Attempts to parse a string containing JSON configuration data.
 *
 * @param content String containing JSON configuration data
 * @return true if parsing was successful, false otherwise
 * @throws std::runtime_error if content contains invalid JSON
 */
bool JsonConfig::LoadFromString(const std::string& content) {
  try {
    root_ = nlohmann::json::parse(content);
    return true;
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to parse JSON string: " +
                             std::string(e.what()));
  }
}

/**
 * @brief Get a value from the configuration
 *
 * Retrieves a value from the configuration using a dot-separated path key.
 * The value is returned as a ConfigValue variant that can contain:
 * - bool
 * - int
 * - float
 * - string
 * - vector<bool>
 * - vector<int>
 * - vector<float>
 * - vector<string>
 *
 * @param key Dot-separated path to the configuration value (e.g.
 * "database.host")
 * @return ConfigValue containing the requested value
 * @throws std::runtime_error if the key doesn't exist or value type is not
 * supported
 */
ConfigValue JsonConfig::Get(const std::string& key) const {
  try {
    const auto& value = GetJsonRef(root_, key);
    return JsonToConfigValue(value);
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to get value for key '" + key +
                             "': " + std::string(e.what()));
  }
}

/**
 * @brief Set a value in the configuration
 *
 * Stores a value in the configuration at the specified dot-separated path key.
 * Creates intermediate objects as needed. Supported value types are:
 * - bool
 * - int
 * - float
 * - string
 * - vector<bool>
 * - vector<int>
 * - vector<float>
 * - vector<string>
 *
 * @param key Dot-separated path where to store the value (e.g. "database.host")
 * @param value ConfigValue containing the value to store
 * @throws std::runtime_error if the key path is invalid or contains invalid
 * characters
 */
void JsonConfig::Set(const std::string& key, const ConfigValue& value) {
  try {
    auto& target = GetJsonRef(root_, key);
    target = ConfigValueToJson(value);
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to set value for key '" + key +
                             "': " + std::string(e.what()));
  }
}

/**
 * @brief Check if a key exists in the configuration
 *
 * Tests whether a value exists at the specified dot-separated path key.
 *
 * @param key Dot-separated path to check (e.g. "database.host")
 * @return true if the key exists and contains a value, false otherwise
 */
bool JsonConfig::HasKey(const std::string& key) const {
  try {
    GetJsonRef(root_, key);
    return true;
  } catch (const std::exception&) {
    return false;
  }
}

/**
 * @brief Save configuration to a JSON file
 *
 * Writes the current configuration state to a JSON file. The output is
 * pretty-printed with 2-space indentation for readability.
 *
 * @param filename Path where to save the configuration file
 * @return true if saving was successful, false otherwise
 */
bool JsonConfig::SaveToFile(const std::string& filename) const {
  try {
    std::ofstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Failed to open file for writing: " + filename);
    }
    file << std::setw(2) << root_ << std::endl;
    return true;
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to save JSON to file: " +
                             std::string(e.what()));
  }
}

/**
 * @brief Convert configuration to a JSON string
 *
 * Returns a string representation of the current configuration state.
 * The output is pretty-printed with 2-space indentation for readability.
 *
 * @return Pretty-printed JSON string representation of the configuration
 */
std::string JsonConfig::ToString() const {
  try {
    return root_.dump(2);
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to convert JSON to string: " +
                             std::string(e.what()));
  }
}

/**
 * @brief Clear all configuration data
 *
 * Removes all key-value pairs from the configuration, resetting it to
 * an empty state.
 */
void JsonConfig::Clear() {
  root_ = nlohmann::json::object();
}

/**
 * @brief Get a reference to a JSON value at the specified path
 *
 * Helper function that traverses the JSON object using a dot-separated path
 * and returns a reference to the JSON value at that location. Creates
 * intermediate objects as needed when used with Set().
 *
 * @param j The JSON object to traverse
 * @param key Dot-separated path to the desired value (e.g. "database.host")
 * @return Reference to the JSON value at the specified path
 * @throws nlohmann::json::exception if the path is invalid or contains invalid
 * characters
 */
nlohmann::json& JsonConfig::GetJsonRef(nlohmann::json& j,
                                       const std::string& key) {
  auto parts = SplitKey(key);
  auto* current = &j;

  for (size_t i = 0; i < parts.size() - 1; ++i) {
    if (!current->contains(parts[i])) {
      (*current)[parts[i]] = nlohmann::json::object();
    }
    current = &(*current)[parts[i]];
    if (!current->is_object()) {
      throw std::runtime_error(
          "Invalid path: intermediate node is not an object");
    }
  }

  return (*current)[parts.back()];
}

/**
 * @brief Get a const reference to a JSON value at the specified path
 *
 * Const version of GetJsonRef() that returns a const reference. Used internally
 * by Get() and HasKey() to safely access values without modification.
 *
 * @param j The JSON object to traverse
 * @param key Dot-separated path to the desired value (e.g. "database.host")
 * @return Const reference to the JSON value at the specified path
 * @throws nlohmann::json::exception if the path is invalid or value doesn't
 * exist
 */
const nlohmann::json& JsonConfig::GetJsonRef(const nlohmann::json& j,
                                             const std::string& key) const {
  auto parts = SplitKey(key);
  const nlohmann::json* current = &j;

  for (size_t i = 0; i < parts.size() - 1; ++i) {
    if (!current->contains(parts[i])) {
      throw std::runtime_error("Key not found: " + key);
    }
    current = &(*current)[parts[i]];
    if (!current->is_object()) {
      throw std::runtime_error(
          "Invalid path: intermediate node is not an object");
    }
  }

  if (!current->contains(parts.back())) {
    throw std::runtime_error("Key not found: " + key);
  }

  return (*current)[parts.back()];
}

/**
 * @brief Convert a JSON value to a ConfigValue
 *
 * Helper function that converts nlohmann::json values to the appropriate
 * ConfigValue variant type. Supports conversion of:
 * - Boolean values to bool
 * - Integer numbers to int
 * - Floating point numbers to float
 * - Strings to string
 * - Arrays of the above types to corresponding vector types
 *
 * For arrays, the element type is determined by the first element. Empty arrays
 * are converted to vector<string> by default.
 *
 * @param j The JSON value to convert
 * @return ConfigValue containing the converted value
 * @throws std::runtime_error if the JSON value type is not supported (e.g.
 * objects or null)
 */
ConfigValue JsonConfig::JsonToConfigValue(const nlohmann::json& j) {
  if (j.is_null()) {
    return ConfigValue();
  } else if (j.is_boolean()) {
    return ConfigValue(j.get<bool>());
  } else if (j.is_number_integer()) {
    return ConfigValue(j.get<int>());
  } else if (j.is_number_float()) {
    return ConfigValue(j.get<float>());
  } else if (j.is_string()) {
    return ConfigValue(j.get<std::string>());
  } else if (j.is_array()) {
    if (j.empty()) {
      return ConfigValue(std::vector<std::string>());
    }

    // Determine array type based on first element
    if (j[0].is_boolean()) {
      std::vector<bool> vec;
      for (const auto& item : j) {
        vec.push_back(item.get<bool>());
      }
      return ConfigValue(vec);
    } else if (j[0].is_number_integer()) {
      std::vector<int> vec;
      for (const auto& item : j) {
        vec.push_back(item.get<int>());
      }
      return ConfigValue(vec);
    } else if (j[0].is_number_float()) {
      std::vector<float> vec;
      for (const auto& item : j) {
        vec.push_back(item.get<float>());
      }
      return ConfigValue(vec);
    } else if (j[0].is_string()) {
      std::vector<std::string> vec;
      for (const auto& item : j) {
        vec.push_back(item.get<std::string>());
      }
      return ConfigValue(vec);
    } else if (j[0].is_object()) {
      // Support for array of objects (nested maps)
      std::vector<framework::ConfigValue> vec;
      for (const auto& item : j) {
        vec.push_back(JsonToConfigValue(item));
      }
      return ConfigValue(vec);
    }
  } else if (j.is_object()) {
    ConfigMap map;
    for (const auto& [key, value] : j.items()) {
      map[key] = JsonToConfigValue(value);
    }
    return ConfigValue(map);
  }

  throw std::runtime_error("Unsupported JSON value type");
}

/**
 * @brief Convert a ConfigValue to a JSON value
 *
 * Helper function that converts ConfigValue variants to nlohmann::json values.
 * Uses std::visit to handle all possible types stored in the ConfigValue:
 * - bool -> JSON boolean
 * - int -> JSON integer
 * - float -> JSON number
 * - string -> JSON string
 * - vector<T> -> JSON array of corresponding type
 *
 * @param value The ConfigValue to convert
 * @return nlohmann::json containing the converted value
 */
nlohmann::json JsonConfig::ConfigValueToJson(const ConfigValue& value) {
  if (value.isNull()) {
    return nlohmann::json();
  } else if (value.isBool()) {
    return nlohmann::json(value.asBool());
  } else if (value.isInt()) {
    return nlohmann::json(value.asInt());
  } else if (value.isFloat()) {
    return nlohmann::json(value.asFloat());
  } else if (value.isString()) {
    return nlohmann::json(value.asString());
  } else if (value.isVectorBool()) {
    return nlohmann::json(value.asVectorBool());
  } else if (value.isVectorInt()) {
    return nlohmann::json(value.asVectorInt());
  } else if (value.isVectorFloat()) {
    return nlohmann::json(value.asVectorFloat());
  } else if (value.isVectorString()) {
    return nlohmann::json(value.asVectorString());
  } else if (value.isVectorConfigValue()) {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& item : value.asVectorConfigValue()) {
      arr.push_back(ConfigValueToJson(item));
    }
    return arr;
  } else if (value.isMap()) {
    nlohmann::json result = nlohmann::json::object();
    for (const auto& [key, val] : value.asMap()) {
      result[key] = ConfigValueToJson(val);
    }
    return result;
  }
  throw std::runtime_error("Unsupported ConfigValue type");
}

JsonConfig::JsonConfig(const framework::ConfigMap& map) {
  // Convert ConfigMap to nlohmann::json object
  root_ = nlohmann::json::object();
  for (const auto& [key, value] : map) {
    root_[key] = ConfigValueToJson(value);
  }
}

}  // namespace metada::backends::config