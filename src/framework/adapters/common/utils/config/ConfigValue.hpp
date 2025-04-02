#pragma once

#include <string>
#include <variant>
#include <vector>

namespace metada::framework {

/**
 * @brief Configuration value type supporting multiple data types
 *
 * @details This variant type provides a flexible container that can hold
 * different types of configuration values. It enables type-safe storage and
 * access of configuration data while maintaining a consistent interface.
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
 *
 * Example usage:
 * @code
 * ConfigValue value = config.Get("server.port");
 * if (std::holds_alternative<int>(value)) {
 *   int port = std::get<int>(value);
 *   // Use port value
 * }
 * @endcode
 */
using ConfigValue = std::variant<bool, int, float, std::string,
                                 std::vector<bool>, std::vector<int>,
                                 std::vector<float>, std::vector<std::string>>;

}  // namespace metada::framework
