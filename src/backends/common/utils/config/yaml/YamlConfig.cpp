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
 * @return true if loading was successful, false otherwise
 */
bool YamlConfig::LoadFromFile(const std::string& filename) {
  try {
    root_ = YAML::LoadFile(filename);
    return true;
  } catch (const YAML::Exception& e) {
    throw std::runtime_error("Failed to load YAML from file: " + std::string(e.what()));
  }
}

/**
 * @brief Load configuration from a YAML string
 * 
 * Attempts to parse a string containing YAML configuration data.
 *
 * @param content String containing YAML configuration data
 * @return true if parsing was successful, false otherwise
 */
bool YamlConfig::LoadFromString(const std::string& content) {
  try {
    root_ = YAML::Load(content);
    return true;
  } catch (const YAML::Exception& e) {
    throw std::runtime_error("Failed to parse YAML string: " + std::string(e.what()));
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
  try {
    const auto& node = GetYamlRef(root_, key);
    return YamlToConfigValue(node);
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to get value for key '" + key + "': " + std::string(e.what()));
  }
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
  try {
    auto& target = GetYamlRef(root_, key);
    target = ConfigValueToYaml(value);
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to set value for key '" + key + "': " + std::string(e.what()));
  }
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
  try {
    GetYamlRef(root_, key);
    return true;
  } catch (const std::exception&) {
    return false;
  }
}

/**
 * @brief Save configuration to a YAML file
 * 
 * Writes the current configuration state to a YAML file.
 * Creates the file if it doesn't exist, overwrites if it does.
 *
 * @param filename Path where to save the configuration file
 * @return true if saving was successful, false otherwise
 */
bool YamlConfig::SaveToFile(const std::string& filename) const {
  try {
    std::ofstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Failed to open file for writing: " + filename);
    }
    file << root_;
    return true;
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to save YAML to file: " + std::string(e.what()));
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
  try {
    std::stringstream ss;
    ss << root_;
    return ss.str();
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to convert YAML to string: " + std::string(e.what()));
  }
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
YAML::Node& YamlConfig::GetYamlRef(YAML::Node& node, const std::string& key) {
  auto parts = SplitKey(key);
  YAML::Node* current = &node;
  
  for (size_t i = 0; i < parts.size() - 1; ++i) {
    if (!(*current)[parts[i]]) {
      (*current)[parts[i]] = YAML::Node();
    }
    current = &(*current)[parts[i]];
    if (!current->IsMap()) {
      throw std::runtime_error("Invalid path: intermediate node is not a map");
    }
  }
  
  return (*current)[parts.back()];
}

const YAML::Node& YamlConfig::GetYamlRef(const YAML::Node& node, const std::string& key) {
  auto parts = SplitKey(key);
  const YAML::Node* current = &node;
  
  for (size_t i = 0; i < parts.size() - 1; ++i) {
    if (!(*current)[parts[i]]) {
      throw std::runtime_error("Key not found: " + key);
    }
    current = &(*current)[parts[i]];
    if (!current->IsMap()) {
      throw std::runtime_error("Invalid path: intermediate node is not a map");
    }
  }
  
  if (!(*current)[parts.back()]) {
    throw std::runtime_error("Key not found: " + key);
  }
  
  return (*current)[parts.back()];
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
ConfigValue YamlConfig::YamlToConfigValue(const YAML::Node& node) {
  if (!node) {
    return ConfigValue();
  }
  
  if (node.IsNull()) {
    return ConfigValue();
  } else if (node.IsScalar()) {
    try {
      return ConfigValue(node.as<bool>());
    } catch (const YAML::Exception&) {
      try {
        return ConfigValue(node.as<int>());
      } catch (const YAML::Exception&) {
        try {
          return ConfigValue(node.as<float>());
        } catch (const YAML::Exception&) {
          return ConfigValue(node.as<std::string>());
        }
      }
    }
  } else if (node.IsSequence()) {
    if (node.size() == 0) {
      return ConfigValue(std::vector<std::string>());
    }
    
    // Determine sequence type based on first element
    try {
      std::vector<bool> vec;
      for (const auto& item : node) {
        vec.push_back(item.as<bool>());
      }
      return ConfigValue(vec);
    } catch (const YAML::Exception&) {
      try {
        std::vector<int> vec;
        for (const auto& item : node) {
          vec.push_back(item.as<int>());
        }
        return ConfigValue(vec);
      } catch (const YAML::Exception&) {
        try {
          std::vector<float> vec;
          for (const auto& item : node) {
            vec.push_back(item.as<float>());
          }
          return ConfigValue(vec);
        } catch (const YAML::Exception&) {
          std::vector<std::string> vec;
          for (const auto& item : node) {
            vec.push_back(item.as<std::string>());
          }
          return ConfigValue(vec);
        }
      }
    }
  } else if (node.IsMap()) {
    ConfigMap map;
    for (const auto& it : node) {
      map[it.first.as<std::string>()] = YamlToConfigValue(it.second);
    }
    return ConfigValue(map);
  }
  
  throw std::runtime_error("Unsupported YAML value type");
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
YAML::Node YamlConfig::ConfigValueToYaml(const ConfigValue& value) {
  return std::visit(
      [](const auto& v) -> YAML::Node {
        using T = std::decay_t<decltype(v)>;
        if constexpr (std::is_same_v<T, std::monostate>) {
          return YAML::Node();
        } else if constexpr (std::is_same_v<T, bool>) {
          return YAML::Node(v);
        } else if constexpr (std::is_same_v<T, int>) {
          return YAML::Node(v);
        } else if constexpr (std::is_same_v<T, float>) {
          return YAML::Node(v);
        } else if constexpr (std::is_same_v<T, std::string>) {
          return YAML::Node(v);
        } else if constexpr (std::is_same_v<T, std::vector<bool>>) {
          YAML::Node node;
          for (const auto& item : v) {
            node.push_back(item);
          }
          return node;
        } else if constexpr (std::is_same_v<T, std::vector<int>>) {
          YAML::Node node;
          for (const auto& item : v) {
            node.push_back(item);
          }
          return node;
        } else if constexpr (std::is_same_v<T, std::vector<float>>) {
          YAML::Node node;
          for (const auto& item : v) {
            node.push_back(item);
          }
          return node;
        } else if constexpr (std::is_same_v<T, std::vector<std::string>>) {
          YAML::Node node;
          for (const auto& item : v) {
            node.push_back(item);
          }
          return node;
        } else if constexpr (std::is_same_v<T, ConfigMap>) {
          YAML::Node node;
          for (const auto& [key, val] : v) {
            node[key] = ConfigValueToYaml(val);
          }
          return node;
        }
      },
      value);
}

std::vector<std::string> YamlConfig::SplitKey(const std::string& key) {
  std::vector<std::string> parts;
  std::string current;
  
  for (char c : key) {
    if (c == '.') {
      if (!current.empty()) {
        parts.push_back(current);
        current.clear();
      }
    } else {
      current += c;
    }
  }
  
  if (!current.empty()) {
    parts.push_back(current);
  }
  
  return parts;
}

}  // namespace metada::backends::config
}  // namespace metada::backends::config
}  // namespace metada::backends::config