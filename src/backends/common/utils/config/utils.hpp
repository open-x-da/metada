#pragma once

#include <string>
#include <vector>

namespace metada::backends::config {

/**
 * @brief Split a dot-separated key into its components
 *
 * Parses a string with dot notation into individual path components.
 * For example, "database.host" would be split into ["database", "host"].
 * 
 * @param key The dot-separated key to split
 * @return Vector of key components
 */
inline std::vector<std::string> SplitKey(const std::string& key) {
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