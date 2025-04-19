/**
 * @file ConfigConcepts.hpp
 * @brief Concept definitions for configuration classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the Config adapter class. These concepts ensure that backend
 * implementations provide all the necessary functionality required by the
 * configuration operations.
 */

#pragma once

#include <concepts>
#include <string>

#include "CommonConcepts.hpp"
#include "ConfigValue.hpp"

namespace metada::framework {

//-----------------------------------------------------------------------------
// Component concepts
//-----------------------------------------------------------------------------

/**
 * @brief Concept that defines requirements for a configuration backend
 * implementation
 *
 * @details A valid configuration backend implementation must provide:
 * - File-based construction and loading operations
 * - Value access operations (get, set, check existence)
 * - Persistence operations (save, string conversion)
 * - Structure management (clear, subsection creation)
 * - Proper resource management with deleted copy constructor and copy
 * assignment operator
 *
 * This concept is used to ensure that backend implementations provide
 * all the necessary functionality required by the Config class.
 *
 * @tparam T The configuration backend implementation type
 *
 * @see HasDeletedCopyConstructor
 * @see HasDeletedCopyAssignment
 */
template <typename T>
concept ConfigBackendImpl =
    requires(T& t, const T& ct, const std::string& key,
             const std::string& filename, const ConfigValue& value) {
      // File construction and loading
      { T(filename) } -> std::same_as<T>;
      { t.LoadFromFile(filename) } -> std::same_as<bool>;
      { t.LoadFromString(filename) } -> std::same_as<bool>;

      // Value access
      { t.Get(key) } -> std::same_as<ConfigValue>;
      { t.Set(key, value) } -> std::same_as<void>;
      { t.HasKey(key) } -> std::same_as<bool>;

      // Persistence
      { t.SaveToFile(filename) } -> std::same_as<bool>;
      { t.ToString() } -> std::same_as<std::string>;

      // Structure management
      { t.Clear() } -> std::same_as<void>;
      { t.CreateSubsection(key) } -> std::same_as<T>;

      // Resource management constraints
      requires HasDeletedCopyConstructor<T>;
      requires HasDeletedCopyAssignment<T>;
    };

/**
 * @brief Concept that defines requirements for a configuration backend tag type
 *
 * @details A valid backend tag must:
 * - Provide a ConfigBackend type through BackendTraits
 * - Ensure the ConfigBackend type satisfies the ConfigBackendImpl concept
 *
 * This concept constrains the template parameter of the Config class,
 * ensuring that only valid backend configurations can be used. It provides
 * compile-time validation of backend compatibility.
 *
 * @tparam T The backend tag type to check
 *
 * @see HasConfigBackend
 * @see ConfigBackendImpl
 */
template <typename T>
concept ConfigBackendType =
    HasConfigBackend<T> &&
    ConfigBackendImpl<typename traits::BackendTraits<T>::ConfigBackend>;

}  // namespace metada::framework