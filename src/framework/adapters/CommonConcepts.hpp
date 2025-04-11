/**
 * @file CommonConcepts.hpp
 * @brief Common concepts shared across different adapter types
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that are used by multiple adapter
 * types across the framework, ensuring consistent definitions and avoiding ODR
 * violations.
 */

#pragma once

#include <concepts>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace metada::framework {

/**
 * @brief Checks if a type provides a ConfigBackend type through BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines a ConfigBackend type through the BackendTraits specialization.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasConfigBackend =
    requires { typename traits::BackendTraits<T>::ConfigBackend; };

/**
 * @brief Checks if a type provides a LoggerBackend type through BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines a LoggerBackend type through the BackendTraits specialization.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasLoggerBackend =
    requires { typename traits::BackendTraits<T>::LoggerBackend; };

/**
 * @brief Checks if a type provides a StateBackend type through BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines a StateBackend type through the BackendTraits specialization.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasStateBackend =
    requires { typename traits::BackendTraits<T>::StateBackend; };

/**
 * @brief Checks if a type provides a ModelBackend type through BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines a ModelBackend type through the BackendTraits specialization.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasModelBackend =
    requires { typename traits::BackendTraits<T>::ModelBackend; };

/**
 * @brief Concept that checks if a type has a constructor from a ConfigBackend
 *
 * @details This concept verifies that a type T has a constructor that takes a
 * ConfigBackend parameter, allowing components to be constructed from
 * configuration.
 */
template <typename T, typename ConfigBackend>
concept HasConfigConstructor = requires(const ConfigBackend& config) {
  {
    T(config)
  } -> std::same_as<T>;  // Check if T can be constructed from config
};

/**
 * @brief Concept that checks if a type provides a clone method
 *
 * @details Verifies that a type T provides a clone method that returns
 * a unique_ptr to a new instance of the same type, enabling deep copying
 * of state objects.
 *
 * @tparam T The type to check for cloning capability
 */
template <typename T>
concept HasClone = requires(const T& t) {
  { t.clone() } -> std::convertible_to<std::unique_ptr<T>>;
};

/**
 * @brief Concept to check if default construction is deleted
 *
 * @details This concept verifies that a type cannot be default-constructed,
 * which is used to ensure classes are always properly initialized with
 * parameters.
 */
template <typename T>
concept HasDeletedDefaultConstructor = !std::is_default_constructible_v<T>;

/**
 * @brief Concept to check if copy construction is deleted
 *
 * @details This concept verifies that a type cannot be copy-constructed,
 * which is a requirement for certain backends to prevent unintended copying.
 */
template <typename T>
concept HasDeletedCopyConstructor = !std::is_copy_constructible_v<T>;

/**
 * @brief Concept to check if copy assignment is deleted
 *
 * @details This concept verifies that a type cannot be copy-assigned,
 * which is a requirement for certain backends to prevent unintended copying.
 */
template <typename T>
concept HasDeletedCopyAssignment = !std::is_copy_assignable_v<T>;

/**
 * @brief Concept requiring initialization capability
 *
 * @details Verifies that a component can be initialized with a configuration
 * object, setting up its internal state before use.
 */
template <typename T, typename ConfigBackend>
concept HasInitialize = requires(T& component, const ConfigBackend& config) {
  { component.initialize(config) } -> std::same_as<void>;
};

}  // namespace metada::framework