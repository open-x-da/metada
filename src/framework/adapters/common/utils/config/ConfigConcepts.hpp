#pragma once

#include <concepts>
#include <string>

#include "../../../CommonConcepts.hpp"
#include "ConfigValue.hpp"

namespace metada::framework {

//-----------------------------------------------------------------------------
// Basic capability concepts
//-----------------------------------------------------------------------------

/**
 * @brief Checks if a type can be initialized with a filename
 */
template <typename T>
concept HasFileConstructor = requires(const std::string& filename) {
  { T(filename) } -> std::same_as<T>;
};

/**
 * @brief Checks if a type supports file-based loading operations
 */
template <typename T>
concept HasFileLoading = requires(T t, const std::string& filename) {
  { t.LoadFromFile(filename) } -> std::same_as<bool>;
  { t.LoadFromString(filename) } -> std::same_as<bool>;
};

/**
 * @brief Checks if a type supports value access operations
 */
template <typename T>
concept HasValueAccess =
    requires(T t, const std::string& key, const ConfigValue& value) {
      { t.Get(key) } -> std::same_as<ConfigValue>;
      { t.Set(key, value) } -> std::same_as<void>;
      { t.HasKey(key) } -> std::same_as<bool>;
    };

/**
 * @brief Checks if a type supports persistence operations
 */
template <typename T>
concept HasPersistence = requires(T t, const std::string& filename) {
  { t.SaveToFile(filename) } -> std::same_as<bool>;
  { t.ToString() } -> std::same_as<std::string>;
};

/**
 * @brief Checks if a type supports structure management operations
 */
template <typename T>
concept HasStructureManagement = requires(T t, const std::string& key) {
  { t.Clear() } -> std::same_as<void>;
  { t.CreateSubsection(key) } -> std::same_as<T>;
};

//-----------------------------------------------------------------------------
// Component concepts
//-----------------------------------------------------------------------------

/**
 * @brief Defines requirements for a configuration backend implementation
 *
 * @details A valid configuration backend must:
 * - Be constructible with a filename
 * - Support loading configuration from files and strings
 * - Provide methods to get, set, and check for configuration values
 * - Support saving configurations to files and string representation
 * - Provide structure management capabilities (clearing and creating
 * subsections)
 *
 * This concept constrains the configuration backends that can be used with
 * the Config class, ensuring they provide all required functionality.
 *
 * @tparam T The configuration backend implementation type
 */
template <typename T>
concept ConfigBackendType =
    HasFileConstructor<T> && HasFileLoading<T> && HasValueAccess<T> &&
    HasPersistence<T> && HasStructureManagement<T>;

}  // namespace metada::framework