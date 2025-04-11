#pragma once

#include <concepts>
#include <string>

#include "../../../CommonConcepts.hpp"
#include "ConfigValue.hpp"

namespace metada::framework {

/**
 * @brief Concept defining requirements for a configuration backend
 *
 * @details This concept enforces the contract that any configuration backend
 * must implement. It provides compile-time validation of the required methods
 * and their signatures, ensuring that backends can be properly used with the
 * Config class.
 *
 * The concept requires specific method signatures for loading, accessing,
 * modifying, and saving configuration data, as well as a constructor that takes
 * a filename.
 */
template <typename T>
concept ConfigBackendType =
    requires(T t, const std::string& filename, const ConfigValue& value) {
      // Required constructor
      { T(filename) } -> std::same_as<T>;
      // Required methods
      { t.LoadFromFile(filename) } -> std::same_as<bool>;
      { t.LoadFromString(filename) } -> std::same_as<bool>;
      { t.Get(filename) } -> std::same_as<ConfigValue>;
      { t.Set(filename, value) } -> std::same_as<void>;
      { t.HasKey(filename) } -> std::same_as<bool>;
      { t.SaveToFile(filename) } -> std::same_as<bool>;
      { t.ToString() } -> std::same_as<std::string>;
      { t.Clear() } -> std::same_as<void>;
      { t.CreateSubsection(filename) } -> std::same_as<T>;
    };

}  // namespace metada::framework