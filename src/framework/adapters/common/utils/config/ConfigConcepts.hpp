#pragma once

#include <concepts>
#include <string>

#include "../../../CommonConcepts.hpp"
#include "ConfigValue.hpp"

namespace metada::framework {

//-----------------------------------------------------------------------------
// Backend capability concepts
//-----------------------------------------------------------------------------

/**
 * @brief Concept that checks if a type can be initialized with a filename
 *
 * @details Verifies that the type has a constructor that accepts a string
 * filename parameter and returns an instance of the type.
 *
 * @tparam T The type to check for file constructor capability
 */
template <typename T>
concept HasFileConstructor = requires(const std::string& filename) {
  { T(filename) } -> std::same_as<T>;
};

/**
 * @brief Concept that checks if a type supports file-based loading operations
 *
 * @details Verifies that the type provides methods to load configuration
 * from both files and strings, with both methods returning a boolean
 * success indicator.
 *
 * @tparam T The type to check for file loading capability
 */
template <typename T>
concept HasFileLoading = requires(T t, const std::string& filename) {
  { t.LoadFromFile(filename) } -> std::same_as<bool>;
  { t.LoadFromString(filename) } -> std::same_as<bool>;
};

/**
 * @brief Concept that checks if a type supports value access operations
 *
 * @details Verifies that the type provides methods to get, set, and check
 * for configuration values using string keys. The Get method must return
 * a ConfigValue, Set must accept a key and value, and HasKey must return
 * a boolean indicating if the key exists.
 *
 * @tparam T The type to check for value access capability
 */
template <typename T>
concept HasValueAccess =
    requires(T t, const std::string& key, const ConfigValue& value) {
      { t.Get(key) } -> std::same_as<ConfigValue>;
      { t.Set(key, value) } -> std::same_as<void>;
      { t.HasKey(key) } -> std::same_as<bool>;
    };

/**
 * @brief Concept that checks if a type supports persistence operations
 *
 * @details Verifies that the type provides methods to save configuration
 * to a file and convert it to a string representation. SaveToFile must
 * return a boolean success indicator, and ToString must return a string.
 *
 * @tparam T The type to check for persistence capability
 */
template <typename T>
concept HasPersistence = requires(T t, const std::string& filename) {
  { t.SaveToFile(filename) } -> std::same_as<bool>;
  { t.ToString() } -> std::same_as<std::string>;
};

/**
 * @brief Concept that checks if a type supports structure management operations
 *
 * @details Verifies that the type provides methods to clear the configuration
 * and create subsections. Clear must return void, and CreateSubsection must
 * return a new instance of the same type representing the subsection.
 *
 * @tparam T The type to check for structure management capability
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
 * @brief Concept that defines requirements for a configuration backend implementation
 *
 * @details A valid configuration backend must:
 * - Be constructible with a filename
 * - Support loading configuration from files and strings
 * - Provide methods to get, set, and check for configuration values
 * - Support saving configurations to files and string representation
 * - Provide structure management capabilities (clearing and creating
 *   subsections)
 * - Have deleted copy constructor and assignment operator to prevent unintended copies
 *
 * This concept is used to ensure that backend implementations provide
 * all the necessary functionality required by the Config class.
 *
 * @tparam T The configuration backend implementation type
 *
 * @see HasFileConstructor
 * @see HasFileLoading
 * @see HasValueAccess
 * @see HasPersistence
 * @see HasStructureManagement
 * @see HasDeletedCopyConstructor
 * @see HasDeletedCopyAssignment
 */
template <typename T>
concept ConfigBackendImpl =
    HasFileConstructor<T> && HasFileLoading<T> && HasValueAccess<T> &&
    HasPersistence<T> && HasStructureManagement<T> &&
    HasDeletedCopyConstructor<T> && HasDeletedCopyAssignment<T>;

/**
 * @brief Concept that defines requirements for a configuration backend tag type
 *
 * @details A valid backend tag must:
 * - Provide a ConfigBackend type through BackendTraits
 * - Ensure the ConfigBackend type satisfies the ConfigBackendImpl concept
 *
 * This concept constrains the template parameter of the Config class,
 * ensuring that only valid backend configurations can be used. It is used
 * in template constraints for classes like Config, Logger, Model, and State
 * that depend on configuration backends.
 *
 * @tparam T The backend tag type to check
 *
 * @see HasConfigBackend
 * @see ConfigBackendImpl
 * @see Config
 */
template <typename T>
concept ConfigBackendType =
    HasConfigBackend<T> &&
    ConfigBackendImpl<typename traits::BackendTraits<T>::ConfigBackend>;

}  // namespace metada::framework