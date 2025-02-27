#include "YamlConfig.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace metada::backends::config {

/**
 * @brief Load configuration from a YAML file
 * 
 * Attempts to load and parse a YAML configuration file. The file must contain
 * valid YAML syntax.
 *
 * @param filename Path to the YAML configuration file
 * @return true if loading was successful, false if file cannot be opened or contains invalid YAML
 */
bool YamlConfig::LoadFromFile(const std::string& filename) {
  try {
    root_ = YAML::LoadFile(filename);
    return true;
  } catch (const YAML::Exception& e) {
    return false;
  }
}

/**
 * @brief Load configuration from a YAML string
 * 
 * Attempts to parse a string containing YAML configuration data.
 *
 * @param content String containing YAML configuration data
 * @return true if parsing was successful, false if content contains invalid YAML
 */
bool YamlConfig::LoadFromString(const std::string& content) {
  try {
    root_ = YAML::Load(content);
    return true;
  } catch (const YAML::Exception& e) {
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
 * @param key Dot-separated path to the configuration value (e.g. "database.host")
 * @return ConfigValue containing the requested value
 * @throws std::runtime_error if the key doesn't exist or value type is not supported
 */
ConfigValue YamlConfig::Get(const std::string& key) const {
  auto node = GetNode(root_, key);
  if (!node.IsDefined()) {
    throw std::runtime_error("Key not found: " + key);
  }
  return NodeToConfigValue(node);
}

// Helper function to split string by delimiter
std::vector<std::string> SplitString(const std::string& str, char delim) {
  std::vector<std::string> tokens;
  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, delim)) {
    tokens.push_back(token);
  }
  return tokens;
}

/**
 * @brief Set a value in the configuration
 * 
 * Sets a value in the configuration using a dot-separated path key.
 * Creates intermediate nodes in the path if they don't exist.
 * Accepts ConfigValue variants containing:
 * - bool
 * - int
 * - double
 * - string
 * - vector<bool>
 * - vector<int>
 * - vector<double>
 * - vector<string>
 *
 * @param key Dot-separated path where to set the value (e.g. "database.host")
 * @param value ConfigValue containing the value to set
 */
void YamlConfig::Set(const std::string& key,
                     const ConfigValue& value) {
  auto keys = SplitString(key, '.');
  YAML::Node current = root_;

  // Build the path, creating nodes as needed
  for (size_t i = 0; i < keys.size() - 1; ++i) {
    if (!current[keys[i]]) {
      current[keys[i]] = YAML::Node(YAML::NodeType::Map);
    }
    current = current[keys[i]];
  }

  // Set the final value
  current[keys.back()] = ConfigValueToNode(value);
  root_ = current;
}

/**
 * @brief Check if a key exists in the configuration
 * 
 * Verifies if a value exists at the specified dot-separated path key.
 *
 * @param key Dot-separated path to check (e.g. "database.password")
 * @return true if the key exists and has a value, false otherwise
 */
bool YamlConfig::HasKey(const std::string& key) const {
  return GetNode(root_, key).IsDefined();
}

/**
 * @brief Save configuration to a YAML file
 * 
 * Writes the current configuration state to a YAML file.
 * Creates the file if it doesn't exist, overwrites if it does.
 *
 * @param filename Path where to save the configuration file
 * @return true if saving was successful, false if file cannot be created or written
 */
bool YamlConfig::SaveToFile(const std::string& filename) const {
  try {
    std::ofstream fout(filename);
    if (!fout) return false;
    fout << ToString();
    return true;
  } catch (const std::exception& e) {
    return false;
  }
}

/**
 * @brief Get configuration as a YAML string
 * 
 * Serializes the current configuration state to a YAML-formatted string.
 *
 * @return String containing the YAML representation of the configuration
 */
std::string YamlConfig::ToString() const {
  std::stringstream ss;
  ss << YAML::Dump(root_);
  return ss.str();
}

/**
 * @brief Clear all configuration data
 * 
 * Resets the configuration to an empty state by creating a new empty YAML node.
 */
void YamlConfig::Clear() {
  root_ = YAML::Node();
}

/**
 * @brief Helper function to get a YAML node at the specified path
 * 
 * Traverses the YAML node tree following a dot-separated path key.
 * Returns an undefined node if any part of the path doesn't exist.
 *
 * @param node Root YAML node to search in
 * @param key Dot-separated path to the node (e.g. "database.host")
 * @return YAML::Node at the specified path, or undefined node if path doesn't exist
 */
YAML::Node YamlConfig::GetNode(const YAML::Node& node, const std::string& key) {
  auto keys = SplitString(key, '.');
  YAML::Node current = node;

  for (const auto& k : keys) {
    if (!current[k].IsDefined()) {
      return YAML::Node();
    }
    current = current[k];
  }

  return current;
}

/**
 * @brief Convert a YAML node to a ConfigValue
 * 
 * Converts YAML nodes to ConfigValue variants by attempting to parse the node
 * as different supported types in order of precedence.
 *
 * Supported scalar types (in order of precedence):
 * - null (converted to empty string)
 * - boolean
 * - integer
 * - double
 * - string
 *
 * Supported sequence types (in order of precedence):
 * - empty sequence (converted to empty string vector)
 * - vector<bool>
 * - vector<int>
 * - vector<double>
 * - vector<string>
 *
 * @param node The YAML node to convert
 * @return ConfigValue containing the converted value
 * @throws std::runtime_error if node type is not supported
 */
ConfigValue YamlConfig::NodeToConfigValue(const YAML::Node& node) {

  if (node.IsScalar()) {
    // Try to parse as different types
    try {
      if (node.IsNull()) return std::string("");
      if (node.as<bool>()) return node.as<bool>();
      return node.as<int>();
    } catch (...) {
      try {
        return node.as<double>();
      } catch (...) {
        return node.as<std::string>();
      }
    }
  } else if (node.IsSequence()) {
    // Try to parse as vector of different types
    if (node.size() == 0) return std::vector<std::string>();

    try {
      return node.as<std::vector<bool>>();
    } catch (...) {
      try {
        return node.as<std::vector<int>>();
      } catch (...) {
        try {
          return node.as<std::vector<double>>();
        } catch (...) {
          return node.as<std::vector<std::string>>();
        }
      }
    }
  }

  throw std::runtime_error("Unsupported YAML node type");
}

/**
 * @brief Convert a ConfigValue to a YAML node
 * 
 * Converts ConfigValue variants to YAML nodes using std::visit to handle
 * all possible stored types.
 *
 * Supported ConfigValue types:
 * - bool -> YAML boolean
 * - int -> YAML integer
 * - double -> YAML float
 * - string -> YAML string
 * - vector<bool> -> YAML sequence of booleans
 * - vector<int> -> YAML sequence of integers
 * - vector<double> -> YAML sequence of floats
 * - vector<string> -> YAML sequence of strings
 *
 * @param value The ConfigValue to convert
 * @return YAML::Node containing the converted value
 */
YAML::Node YamlConfig::ConfigValueToNode(
      const ConfigValue& value) {
  return std::visit([](const auto& v) -> YAML::Node { return YAML::Node(v); },
                    value);
}

}  // namespace metada::backends::common::utils::config::yaml