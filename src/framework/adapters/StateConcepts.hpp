#pragma once

#include <concepts>
#include <memory>
#include <string>
#include <vector>

#include "CommonConcepts.hpp"

namespace metada::framework {

//-----------------------------------------------------------------------------
// Core state capability concepts
//-----------------------------------------------------------------------------

/**
 * @brief Concept that checks if a type provides data access methods
 *
 * @details Verifies that a type T provides methods to access the underlying
 * state data as both mutable and const void pointers, allowing for type-safe
 * casting at the interface level.
 *
 * @tparam T The type to check for data access capability
 */
template <typename T>
concept HasGetData = requires(T& t, const T& ct) {
  { t.getData() } -> std::same_as<void*>;
  { ct.getData() } -> std::same_as<const void*>;
};

/**
 * @brief Concept that checks if a type provides variable name access
 *
 * @details Verifies that a type T provides a method to retrieve the names
 * of state variables as a vector of strings, supporting discovery and
 * introspection of state components.
 *
 * @tparam T The type to check for variable name access capability
 */
template <typename T>
concept HasGetVariableNames = requires(const T& t) {
  { t.getVariableNames() } -> std::same_as<const std::vector<std::string>&>;
};

/**
 * @brief Concept that checks if a type provides dimension information
 *
 * @details Verifies that a type T provides a method to retrieve the dimensions
 * of each state variable, supporting multi-dimensional state representations.
 *
 * @tparam T The type to check for dimension access capability
 */
template <typename T>
concept HasGetDimensions = requires(const T& t) {
  {
    t.getDimensions(std::string{})
  } -> std::same_as<const std::vector<size_t>&>;
};

//-----------------------------------------------------------------------------
// Component concepts
//-----------------------------------------------------------------------------

/**
 * @brief Concept that defines requirements for a state backend implementation
 *
 * @details A valid state backend must:
 * - Be constructible from a ConfigBackend instance
 * - Provide a clone method for deep copying
 * - Provide methods to access underlying data
 * - Provide methods to retrieve variable names and dimensions
 * - Have deleted copy constructor and assignment operator
 *
 * This concept is used to ensure that backend implementations provide
 * all the necessary functionality required by the State class.
 *
 * @tparam T The state backend implementation type
 * @tparam ConfigBackend The configuration backend type
 *
 * @see HasConfigConstructor
 * @see HasClone
 * @see HasGetData
 * @see HasGetVariableNames
 * @see HasGetDimensions
 */
template <typename T, typename ConfigBackend>
concept StateBackendImpl =
    HasConfigConstructor<T, ConfigBackend> && HasClone<T> && HasGetData<T> &&
    HasGetVariableNames<T> && HasGetDimensions<T> &&
    HasDeletedCopyConstructor<T> && HasDeletedCopyAssignment<T>;

/**
 * @brief Concept that defines requirements for a state backend tag type
 *
 * @details A valid backend tag must:
 * - Provide both StateBackend and ConfigBackend types via BackendTraits
 * - Ensure the StateBackend type satisfies the StateBackendImpl concept
 *   when paired with the ConfigBackend type
 *
 * This concept constrains the template parameter of the State class,
 * ensuring that only valid backend configurations can be used. It provides
 * compile-time validation of backend compatibility.
 *
 * @tparam T The backend tag type to check
 *
 * @see HasStateBackend
 * @see HasConfigBackend
 * @see StateBackendImpl
 */
template <typename T>
concept StateBackendType =
    HasStateBackend<T> && HasConfigBackend<T> &&
    StateBackendImpl<typename traits::BackendTraits<T>::StateBackend,
                     typename traits::BackendTraits<T>::ConfigBackend>;

}  // namespace metada::framework