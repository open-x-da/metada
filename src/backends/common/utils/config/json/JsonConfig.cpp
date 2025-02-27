#include "JsonConfig.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace metada::backends::common::utils::config::json {

/**
 * @brief Load configuration from a JSON file
 *
 * Attempts to load and parse a JSON configuration file. The file must contain
 * valid JSON syntax.
 *
 * @param filename Path to the JSON configuration file
 * @return true if loading was successful, false if file cannot be opened or
 * contains invalid JSON
 */
bool JsonConfig::LoadFromFile(const std::string& filename) {
  try {
    std::ifstream fin(filename);
    if (!fin) return false;
    fin >> root_;
    return true;
  } catch (const nlohmann::json::exception& e) {
    return false;
  }
}

/**
 * @brief Load configuration from a JSON string
 *
 * Attempts to parse a string containing JSON configuration data.
 *
 * @param content String containing JSON configuration data
 * @return true if parsing was successful, false if content contains invalid
 * JSON
 */
bool JsonConfig::LoadFromString(const std::string& content) {
  try {
    root_ = nlohmann::json::parse(content);
    return true;
  } catch (const nlohmann::json::exception& e) {
    return false;
  }
}

/**
 * @brief Get a value from the configuration
 *
 * Retrieves a value from the configuration using a dot-separated path key.
 * The value is returned as a ConfigValue variant that can contain:
 * - bool
 * - int
 * - double
 * - string
 * - vector<bool>
 * - vector<int>
 * - vector<double>
 * - vector<string>
 *
 * @param key Dot-separated path to the configuration value (e.g.
 * "database.host")
 * @return ConfigValue containing the requested value
 * @throw std::runtime_error if the key doesn't exist or value type is not
 * supported
 */
framework::ConfigValue JsonConfig::Get(const std::string& key) const {
  try {
    const auto& j = GetJsonRef(root_, key);
    return JsonToConfigValue(j);
  } catch (const std::exception& e) {
    throw std::runtime_error("Key not found: " + key);
  }
}

/**
 * @brief Set a value in the configuration
 *
 * Stores a value in the configuration at the specified dot-separated path key.
 * Creates intermediate objects as needed. Supported value types are:
 * - bool
 * - int
 * - double
 * - string
 * - vector<bool>
 * - vector<int>
 * - vector<double>
 * - vector<string>
 *
 * @param key Dot-separated path where to store the value (e.g. "database.host")
 * @param value ConfigValue containing the value to store
 * @throw std::runtime_error if the key path is invalid or contains invalid
 * characters
 */
void JsonConfig::Set(const std::string& key, const ConfigValue& value) {
  auto& target = GetJsonRef(root_, key);
  target = ConfigValueToJson(value);
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
  } catch (...) {
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
 * @return true if saving was successful, false if file cannot be written
 */
bool JsonConfig::SaveToFile(const std::string& filename) const {
  try {
    std::ofstream fout(filename);
    if (!fout) return false;
    fout << root_.dump(2);  // Pretty print with indent=2
    return true;
  } catch (const std::exception& e) {
    return false;
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
  return root_.dump(2);  // Pretty print with indent=2
}

/**
 * @brief Clear all configuration data
 *
 * Removes all key-value pairs from the configuration, resetting it to
 * an empty state.
 */
void JsonConfig::Clear() {
  root_.clear();
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
 * @throw nlohmann::json::exception if the path is invalid or contains invalid
 * characters
 */
nlohmann::json& JsonConfig::GetJsonRef(nlohmann::json& j,
                                       const std::string& key) {
  std::istringstream key_stream(key);
  std::string part;
  nlohmann::json* current = &j;

  while (std::getline(key_stream, part, '.')) {
    current = &((*current)[part]);
  }
  return *current;
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
 * @throw nlohmann::json::exception if the path is invalid or value doesn't
 * exist
 */
const nlohmann::json& JsonConfig::GetJsonRef(const nlohmann::json& j,
                                             const std::string& key) {
  std::istringstream key_stream(key);
  std::string part;
  const nlohmann::json* current = &j;

  while (std::getline(key_stream, part, '.')) {
    current = &((*current)[part]);
  }
  return *current;
}

/**
 * @brief Convert a JSON value to a ConfigValue
 *
 * Helper function that converts nlohmann::json values to the appropriate
 * ConfigValue variant type. Supports conversion of:
 * - Boolean values to bool
 * - Integer numbers to int
 * - Floating point numbers to double
 * - Strings to string
 * - Arrays of the above types to corresponding vector types
 *
 * For arrays, the element type is determined by the first element. Empty arrays
 * are converted to vector<string> by default.
 *
 * @param j The JSON value to convert
 * @return ConfigValue containing the converted value
 * @throw std::runtime_error if the JSON value type is not supported (e.g.
 * objects or null)
 */
framework::ConfigValue JsonConfig::JsonToConfigValue(const nlohmann::json& j) {
  using framework::ConfigValue;

  if (j.is_boolean()) {
    return j.get<bool>();
  } else if (j.is_number_integer()) {
    return j.get<int>();
  } else if (j.is_number_float()) {
    return j.get<double>();
  } else if (j.is_string()) {
    return j.get<std::string>();
  } else if (j.is_array()) {
    if (j.empty()) return std::vector<std::string>();

    // Try to determine array type from first element
    if (j[0].is_boolean()) {
      return j.get<std::vector<bool>>();
    } else if (j[0].is_number_integer()) {
      return j.get<std::vector<int>>();
    } else if (j[0].is_number_float()) {
      return j.get<std::vector<double>>();
    } else {
      return j.get<std::vector<std::string>>();
    }
  }

  throw std::runtime_error("Unsupported JSON type");
}

/**
 * @brief Convert a ConfigValue to a JSON value
 *
 * Helper function that converts ConfigValue variants to nlohmann::json values.
 * Uses std::visit to handle all possible types stored in the ConfigValue:
 * - bool -> JSON boolean
 * - int -> JSON integer
 * - double -> JSON number
 * - string -> JSON string
 * - vector<T> -> JSON array of corresponding type
 *
 * @param value The ConfigValue to convert
 * @return nlohmann::json containing the converted value
 */
nlohmann::json JsonConfig::ConfigValueToJson(const ConfigValue& value) {
  return std::visit(
      [](const auto& v) -> nlohmann::json { return nlohmann::json(v); }, value);
}

}  // namespace metada::backends::common::utils::config::json