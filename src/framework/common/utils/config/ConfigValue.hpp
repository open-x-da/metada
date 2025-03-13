#pragma once

#include <string>
#include <variant>
#include <vector>

namespace metada::framework {

/**
 * @brief Configuration value type supporting multiple data types
 *
 * This variant type provides a flexible container that can hold different types
 * of configuration values. It enables type-safe storage and access of
 * configuration data while maintaining a consistent interface.
 *
 * Supported value types:
 * - bool: Boolean values (true/false)
 * - int: Integer numbers
 * - double: Floating point numbers
 * - std::string: Text strings
 * - std::vector<bool>: Arrays of boolean values
 * - std::vector<int>: Arrays of integer numbers
 * - std::vector<double>: Arrays of floating point numbers
 * - std::vector<std::string>: Arrays of text strings
 */
using ConfigValue = std::variant<bool, int, float, std::string,
                                 std::vector<bool>, std::vector<int>,
                                 std::vector<float>, std::vector<std::string>>;

}  // namespace metada::framework
