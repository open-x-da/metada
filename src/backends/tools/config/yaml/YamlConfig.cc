#include "YamlConfig.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace metada {
namespace backends {
namespace tools {
namespace config {
namespace yaml {

/**
 * @brief Load configuration from a YAML file
 * @param filename Path to the YAML configuration file
 * @return true if loading was successful, false otherwise
 * @throws std::runtime_error if file cannot be opened or contains invalid YAML
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
 * @param content String containing YAML configuration data
 * @return true if loading was successful, false otherwise
 * @throws std::runtime_error if content contains invalid YAML
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
 * @param key Dot-separated path to the configuration value (e.g.
 * "database.host")
 * @return ConfigValue containing the requested value
 * @throws std::runtime_error if the key doesn't exist or value type is not
 * supported
 */
framework::tools::config::ConfigValue YamlConfig::Get(
    const std::string& key) const {
  auto node = GetNode(root_, key);
  if (!node.IsDefined()) {
    throw std::runtime_error("Key not found: " + key);
  }
  return NodeToConfigValue(node);
}

/**
 * @brief Set a value in the configuration
 * @param key Dot-separated path where to store the value (e.g. "database.port")
 * @param value ConfigValue to store
 * @throws std::runtime_error if the key path is invalid or value type is not
 * supported
 */
void YamlConfig::Set(const std::string& key,
                     const framework::tools::config::ConfigValue& value) {
  auto keys = YAML::Split(key, '.');
  YAML::Node* current = &root_;

  for (size_t i = 0; i < keys.size() - 1; ++i) {
    current = &(*current)[keys[i]];
  }

  (*current)[keys.back()] = ConfigValueToNode(value);
}

/**
 * @brief Check if a key exists in the configuration
 * @param key Dot-separated path to check (e.g. "database.password")
 * @return true if the key exists, false otherwise
 */
bool YamlConfig::HasKey(const std::string& key) const {
  return GetNode(root_, key).IsDefined();
}

/**
 * @brief Save configuration to a YAML file
 * @param filename Path where to save the configuration
 * @return true if saving was successful, false otherwise
 * @throws std::runtime_error if file cannot be created or written to
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
 * @return String containing the YAML representation of the configuration
 * @throws std::runtime_error if serialization fails
 */
std::string YamlConfig::ToString() const {
  std::stringstream ss;
  ss << root_;
  return ss.str();
}

/**
 * @brief Clear all configuration data
 * Resets the configuration to an empty state
 */
void YamlConfig::Clear() {
  root_ = YAML::Node();
}

/**
 * @brief Helper function to get a YAML node at the specified path
 * @param node Root YAML node to search in
 * @param key Dot-separated path to the node (e.g. "database.host")
 * @return YAML::Node at the specified path, or undefined node if path doesn't
 * exist
 */
YAML::Node YamlConfig::GetNode(const YAML::Node& node, const std::string& key) {
  auto keys = YAML::Split(key, '.');
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
 * This helper function converts YAML nodes to ConfigValue variants.
 * Supports the following YAML types:
 * - Scalar: null, boolean, integer, double, string
 * - Sequence: vector of boolean, integer, double, or string
 *
 * @param node The YAML node to convert
 * @return ConfigValue containing the converted value
 * @throws std::runtime_error if node type is not supported
 */
framework::tools::config::ConfigValue YamlConfig::NodeToConfigValue(
    const YAML::Node& node) {
  using framework::tools::config::ConfigValue;

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
 * This helper function converts ConfigValue variants to YAML nodes.
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
 * @return YAML::Node containing the converted value
 */
YAML::Node YamlConfig::ConfigValueToNode(
    const framework::tools::config::ConfigValue& value) {
  return std::visit([](const auto& v) -> YAML::Node { return YAML::Node(v); },
                    value);
}

}  // namespace yaml
}  // namespace config
}  // namespace tools
}  // namespace backends
}  // namespace metada