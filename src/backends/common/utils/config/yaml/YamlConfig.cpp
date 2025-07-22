#include "YamlConfig.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>

namespace metada::backends::config {

/**
 * @brief Load configuration from a YAML file
 *
 * Attempts to load and parse a YAML configuration file. The file must contain
 * valid YAML syntax.
 *
 * @param filename Path to the YAML configuration file
 * @return true if loading was successful, false otherwise
 * @throws std::runtime_error if file cannot be opened or contains invalid YAML
 */
bool YamlConfig::LoadFromFile(const std::string& filename) {
  try {
    root_ = YAML::LoadFile(filename);
    return true;
  } catch (const YAML::Exception& e) {
    throw std::runtime_error("Failed to load YAML from file: " +
                             std::string(e.what()));
  }
}

/**
 * @brief Load configuration from a YAML string
 *
 * Attempts to parse a string containing YAML configuration data.
 *
 * @param content String containing YAML configuration data
 * @return true if parsing was successful, false otherwise
 * @throws std::runtime_error if content contains invalid YAML
 */
bool YamlConfig::LoadFromString(const std::string& content) {
  try {
    root_ = YAML::Load(content);
    return true;
  } catch (const YAML::Exception& e) {
    throw std::runtime_error("Failed to parse YAML string: " +
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
 * - ConfigMap (nested configuration)
 *
 * @param key Dot-separated path to the configuration value (e.g.
 * "database.host")
 * @return ConfigValue containing the requested value
 * @throws std::runtime_error if the key doesn't exist or value type is not
 * supported
 */
ConfigValue YamlConfig::Get(const std::string& key) const {
  try {
    const auto& node = GetYamlNodeConst(root_, key);
    return YamlToConfigValue(node);
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to get value for key '" + key +
                             "': " + std::string(e.what()));
  }
}

/**
 * @brief Set a value in the configuration
 *
 * Sets a value in the configuration using a dot-separated path key.
 * Creates intermediate nodes in the path if they don't exist.
 * Accepts ConfigValue variants containing:
 * - bool
 * - int
 * - float
 * - string
 * - vector<bool>
 * - vector<int>
 * - vector<float>
 * - vector<string>
 * - ConfigMap (nested configuration)
 *
 * @param key Dot-separated path where to set the value (e.g. "database.host")
 * @param value ConfigValue containing the value to set
 * @throws std::runtime_error if the key path is invalid or value type is not
 * supported
 */
void YamlConfig::Set(const std::string& key, const ConfigValue& value) {
  try {
    // Use our path resolver to find the right node to modify
    auto parts = SplitKey(key);
    if (parts.empty()) {
      // Setting the root node
      root_ = ConfigValueToYaml(value);
      return;
    }

    // Navigate to the parent node
    YAML::Node* current = &root_;
    for (size_t i = 0; i < parts.size() - 1; ++i) {
      if (!current->IsMap()) {
        *current = YAML::Node(YAML::NodeType::Map);
      }
      // Ensure the child exists as a map
      if (!(*current)[parts[i]]) {
        (*current)[parts[i]] = YAML::Node(YAML::NodeType::Map);
      }
      // Use a temporary variable to store the result
      YAML::Node temp = (*current)[parts[i]];
      current = &temp;
    }

    // Set the value at the leaf node
    (*current)[parts.back()] = ConfigValueToYaml(value);
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to set value for key '" + key +
                             "': " + std::string(e.what()));
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
    GetYamlNodeConst(root_, key);
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
 * @throws std::runtime_error if file cannot be created or written to
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
    throw std::runtime_error("Failed to save YAML to file: " +
                             std::string(e.what()));
  }
}

/**
 * @brief Get configuration as a YAML string
 *
 * Serializes the current configuration state to a YAML-formatted string.
 *
 * @return String containing the YAML representation of the configuration
 * @throws std::runtime_error if serialization fails
 */
std::string YamlConfig::ToString() const {
  try {
    std::stringstream ss;
    ss << root_;
    return ss.str();
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to convert YAML to string: " +
                             std::string(e.what()));
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
 * @brief Convert a YAML node to a ConfigValue
 *
 * Converts YAML nodes to ConfigValue variants by attempting to parse the node
 * as different supported types in order of precedence.
 *
 * Supported scalar types (in order of precedence):
 * - null (converted to empty ConfigValue)
 * - boolean
 * - integer
 * - float
 * - string
 *
 * Supported sequence types (in order of precedence):
 * - empty sequence (converted to empty string vector)
 * - vector<bool>
 * - vector<int>
 * - vector<float>
 * - vector<string>
 *
 * Supported map type:
 * - ConfigMap (nested configuration)
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
    // Use YAML-cpp's built-in type detection to respect original YAML types
    // This properly handles quoted strings vs unquoted values

    // Try bool first
    try {
      auto b = node.as<bool>();
      return ConfigValue(b);
    } catch (const YAML::Exception&) {
    }

    // Try int - but only if the YAML library thinks it's actually an integer
    try {
      // This will only succeed if the node is genuinely an integer in YAML
      // It won't convert "3DVAR" to 3 because YAML-cpp respects the quotes
      auto i = node.as<int>();
      return ConfigValue(i);
    } catch (const YAML::Exception&) {
    }

    // Try float - but only if the YAML library thinks it's actually a float
    try {
      auto f = node.as<float>();
      return ConfigValue(f);
    } catch (const YAML::Exception&) {
    }

    // Fallback to string - this will handle quoted strings like "3DVAR"
    // properly
    return ConfigValue(node.as<std::string>());
  } else if (node.IsSequence()) {
    if (node.size() == 0) {
      return ConfigValue(std::vector<std::string>());
    }

    // Try to detect if this is a sequence of maps
    if (node[0].IsMap()) {
      std::vector<ConfigValue> vec;
      for (const auto& item : node) {
        vec.push_back(YamlToConfigValue(item));
      }
      return ConfigValue(vec);
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
    framework::ConfigMap map;
    for (const auto& it : node) {
      map[it.first.as<std::string>()] = YamlToConfigValue(it.second);
    }
    return ConfigValue(map);
  }

  throw std::runtime_error("Unsupported YAML value type");
}

// Custom visitor for converting ConfigValue to YAML
class ConfigValueToYamlVisitor {
 public:
  YAML::Node operator()(std::monostate) const { return YAML::Node(); }

  YAML::Node operator()(bool v) const { return YAML::Node(v); }

  YAML::Node operator()(int v) const { return YAML::Node(v); }

  YAML::Node operator()(float v) const {
    // Always emit as string with decimal to preserve float-ness
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << v;
    return YAML::Node(oss.str());
  }

  YAML::Node operator()(const std::string& v) const { return YAML::Node(v); }

  YAML::Node operator()(const std::vector<bool>& v) const {
    YAML::Node node;
    for (bool item : v) {
      node.push_back(item);
    }
    return node;
  }

  YAML::Node operator()(const std::vector<int>& v) const {
    YAML::Node node;
    for (int item : v) {
      node.push_back(item);
    }
    return node;
  }

  YAML::Node operator()(const std::vector<float>& v) const {
    YAML::Node node;
    for (float item : v) {
      node.push_back(item);
    }
    return node;
  }

  YAML::Node operator()(const std::vector<std::string>& v) const {
    YAML::Node node;
    for (const std::string& item : v) {
      node.push_back(item);
    }
    return node;
  }

  YAML::Node operator()(
      const std::vector<metada::framework::ConfigValue>& v) const {
    YAML::Node node;
    ConfigValueToYamlVisitor tempVisitor;
    for (const auto& item : v) {
      node.push_back(std::visit(tempVisitor, item.constVariant()));
    }
    return node;
  }

  YAML::Node operator()(const framework::ConfigMap& v) const {
    YAML::Node node;
    for (const auto& [key, val] : v) {
      // Create a temporary visitor and use it directly
      ConfigValueToYamlVisitor tempVisitor;
      if (val.isNull())
        node[key] = tempVisitor(std::monostate{});
      else if (val.isBool())
        node[key] = tempVisitor(val.asBool());
      else if (val.isInt())
        node[key] = tempVisitor(val.asInt());
      else if (val.isFloat())
        node[key] = tempVisitor(val.asFloat());
      else if (val.isString())
        node[key] = tempVisitor(val.asString());
      else if (val.isVectorBool())
        node[key] = tempVisitor(val.asVectorBool());
      else if (val.isVectorInt())
        node[key] = tempVisitor(val.asVectorInt());
      else if (val.isVectorFloat())
        node[key] = tempVisitor(val.asVectorFloat());
      else if (val.isVectorString())
        node[key] = tempVisitor(val.asVectorString());
      else if (val.isVectorConfigValue())
        node[key] = tempVisitor(val.asVectorConfigValue());
      else if (val.isMap())
        node[key] = tempVisitor(val.asMap());
      else
        node[key] = YAML::Node();
    }
    return node;
  }
};

YAML::Node YamlConfig::ConfigValueToYaml(const ConfigValue& value) {
  ConfigValueToYamlVisitor visitor;

  if (value.isNull()) return visitor(std::monostate{});
  if (value.isBool()) return visitor(value.asBool());
  if (value.isInt()) return visitor(value.asInt());
  if (value.isFloat()) return visitor(value.asFloat());
  if (value.isString()) return visitor(value.asString());
  if (value.isVectorBool()) return visitor(value.asVectorBool());
  if (value.isVectorInt()) return visitor(value.asVectorInt());
  if (value.isVectorFloat()) return visitor(value.asVectorFloat());
  if (value.isVectorString()) return visitor(value.asVectorString());
  if (value.isVectorConfigValue()) return visitor(value.asVectorConfigValue());
  if (value.isMap()) return visitor(value.asMap());

  // Default case
  return YAML::Node();
}

// Helper to get a YAML node at a path
YAML::Node YamlConfig::GetYamlNode(YAML::Node node, const std::string& key) {
  auto parts = SplitKey(key);

  for (const auto& part : parts) {
    if (!node.IsMap()) {
      throw std::runtime_error("Invalid path: intermediate node is not a map");
    }

    node = node[part];
    if (!node) {
      throw std::runtime_error("Key not found: " + part);
    }
  }

  return node;
}

// Helper to get a const YAML node at a path
const YAML::Node YamlConfig::GetYamlNodeConst(const YAML::Node& node,
                                              const std::string& key) {
  auto parts = SplitKey(key);

  if (parts.empty()) {
    return node;
  }

  // Handle the first part
  const auto& part = parts.front();

  // Check if node is a map
  if (!node.IsMap()) {
    throw std::runtime_error("Invalid path: node is not a map");
  }

  // Get the child node (creates a copy)
  YAML::Node child = node[part];

  // Check if child exists
  if (!child) {
    throw std::runtime_error("Key not found: " + part);
  }

  // If we have more parts, recurse with the remaining parts
  if (parts.size() > 1) {
    // Create a new path with the remaining parts
    std::string remaining;
    for (size_t i = 1; i < parts.size(); ++i) {
      if (i > 1) remaining += ".";
      remaining += parts[i];
    }
    return GetYamlNodeConst(child, remaining);
  }

  // Return the child node (last part of the path)
  return child;
}

YamlConfig::YamlConfig(const framework::ConfigMap& map) {
  root_ = ConfigValueToYaml(ConfigValue(map));
}

}  // namespace metada::backends::config