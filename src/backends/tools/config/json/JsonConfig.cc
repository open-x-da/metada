#include "JsonConfig.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace metada {
namespace backends {
namespace tools {
namespace config {
namespace json {

/**
 * @brief Load configuration from a JSON file
 * @param filename Path to the JSON configuration file
 * @return true if loading was successful, false otherwise
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
 * @param content String containing the JSON configuration data
 * @return true if loading was successful, false otherwise
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
 * @param key Dot-separated path to the configuration value (e.g.
 * "parent.child.value")
 * @return ConfigValue containing the requested value
 * @throw std::runtime_error if the key doesn't exist or value type is not
 * supported
 */
framework::tools::config::ConfigValue JsonConfig::Get(
    const std::string& key) const {
  try {
    const auto& j = GetJsonRef(root_, key);
    return JsonToConfigValue(j);
  } catch (const std::exception& e) {
    throw std::runtime_error("Key not found: " + key);
  }
}

/**
 * @brief Set a value in the configuration
 * @param key Dot-separated path where to store the value (e.g.
 * "parent.child.value")
 * @param value ConfigValue containing the value to store
 * @throw std::runtime_error if the key path is invalid
 */
void JsonConfig::Set(const std::string& key,
                     const framework::tools::config::ConfigValue& value) {
  auto& target = GetJsonRef(root_, key);
  target = ConfigValueToJson(value);
}

/**
 * @brief Check if a key exists in the configuration
 * @param key Dot-separated path to check (e.g. "parent.child.value")
 * @return true if the key exists, false otherwise
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
 * @param filename Path where to save the configuration file
 * @return true if saving was successful, false otherwise
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
 * @return Pretty-printed JSON string representation of the configuration
 */
std::string JsonConfig::ToString() const {
  return root_.dump(2);  // Pretty print with indent=2
}

/**
 * @brief Clear all configuration data
 */
void JsonConfig::Clear() {
  root_.clear();
}

/**
 * @brief Get a reference to a JSON value at the specified path
 *
 * This helper function traverses the JSON object using a dot-separated path
 * and returns a reference to the JSON value at that path.
 *
 * @param j The JSON object to traverse
 * @param key Dot-separated path to the desired value (e.g.
 * "parent.child.value")
 * @return Reference to the JSON value at the specified path
 * @throw nlohmann::json::exception if the path is invalid or value doesn't
 * exist
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
 * Const version of GetJsonRef that returns a const reference.
 *
 * @param j The JSON object to traverse
 * @param key Dot-separated path to the desired value (e.g.
 * "parent.child.value")
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
 * This helper function converts nlohmann::json values to the appropriate
 * ConfigValue type. Supports conversion of:
 * - Boolean values
 * - Integer numbers
 * - Floating point numbers
 * - Strings
 * - Arrays of the above types
 *
 * For arrays, the type is determined by the first element. Empty arrays default
 * to string arrays.
 *
 * @param j The JSON value to convert
 * @return ConfigValue containing the converted value
 * @throw std::runtime_error if the JSON value type is not supported (e.g.
 * objects or null)
 */
framework::tools::config::ConfigValue JsonConfig::JsonToConfigValue(
    const nlohmann::json& j) {
  using framework::tools::config::ConfigValue;

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
 * This helper function converts ConfigValue variants to nlohmann::json values.
 * Uses std::visit to handle all possible types stored in the ConfigValue
 * variant:
 * - bool
 * - int
 * - double
 * - string
 * - vector<bool>
 * - vector<int>
 * - vector<double>
 * - vector<string>
 *
 * @param value The ConfigValue to convert
 * @return nlohmann::json containing the converted value
 */
nlohmann::json JsonConfig::ConfigValueToJson(
    const framework::tools::config::ConfigValue& value) {
  return std::visit(
      [](const auto& v) -> nlohmann::json { return nlohmann::json(v); }, value);
}

}  // namespace json
}  // namespace config
}  // namespace tools
}  // namespace backends
}  // namespace metada