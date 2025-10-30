/**
 * @file CommonConcepts.hpp
 * @brief Common concepts shared across different adapter types
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that are used by multiple adapter
 * types across the framework, ensuring consistent definitions and avoiding ODR
 * violations. These concepts define requirements for backend implementations
 * and provide compile-time validation of adapter-backend compatibility.
 */

#pragma once

#include <concepts>
#include <memory>
#include <type_traits>

namespace metada::framework {

/**
 * @brief Checks if a type provides a ConfigBackend type through BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines a ConfigBackend type through the BackendTraits specialization.
 * The ConfigBackend is responsible for handling configuration parameters.
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
 * The LoggerBackend is responsible for handling logging operations.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasLoggerBackend =
    requires { typename traits::BackendTraits<T>::LoggerBackend; };

/**
 * @brief Checks if a type provides a GeometryBackend type through BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines a GeometryBackend type through the BackendTraits specialization.
 * The GeometryBackend is responsible for handling spatial domain
 * representations.
 *
 * Note: GeometryBackend does not need to be a template class; it can be a
 * concrete type or a class with a templated constructor.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasGeometryBackend =
    requires { typename traits::BackendTraits<T>::GeometryBackend; };

/**
 * @brief Checks if a type provides a GeometryIteratorBackend type through
 * BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines a GeometryIteratorBackend type through the BackendTraits
 * specialization. The GeometryIteratorBackend enables traversal of geometry
 * elements.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasGeometryIteratorBackend =
    requires { typename traits::BackendTraits<T>::GeometryIteratorBackend; };

/**
 * @brief Checks if a type provides a StateBackend type through BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines a StateBackend type through the BackendTraits specialization.
 * The StateBackend is responsible for handling model state representations
 * and operations such as vector arithmetic.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasStateBackend =
    requires { typename traits::BackendTraits<T>::StateBackend; };

/**
 * @brief Checks if a type provides an IncrementBackend type through
 * BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines an IncrementBackend type through the BackendTraits specialization.
 * The IncrementBackend is responsible for handling analysis increments
 * in variational data assimilation, independent of state representations.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasIncrementBackend =
    requires { typename traits::BackendTraits<T>::IncrementBackend; };

/**
 * @brief Checks if a type provides a ModelBackend type through BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines a ModelBackend type through the BackendTraits specialization.
 * The ModelBackend is responsible for handling model integrations and related
 * operations.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasModelBackend =
    requires { typename traits::BackendTraits<T>::ModelBackend; };

/**
 * @brief Checks if a type provides an ObservationBackend type through
 * BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines an ObservationBackend type through the BackendTraits specialization.
 * The ObservationBackend is responsible for handling observational data and
 * related operations in data assimilation systems.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasObservationBackend =
    requires { typename traits::BackendTraits<T>::ObservationBackend; };

/**
 * @brief Checks if a type provides an ObsOperatorBackend type through
 * BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines an ObsOperatorBackend type through the BackendTraits specialization.
 * The ObsOperatorBackend is responsible for handling the observation operator
 * and related operations in data assimilation systems.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasObsOperatorBackend =
    requires { typename traits::BackendTraits<T>::ObsOperatorBackend; };

/**
 * @brief Checks if a type provides an ObservationIteratorBackend type through
 * BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines an ObservationIteratorBackend type through the BackendTraits
 * specialization. The ObservationIteratorBackend enables traversal of
 * observation data points.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept HasObservationIteratorBackend =
    requires { typename traits::BackendTraits<T>::ObservationIteratorBackend; };

/**
 * @brief Checks if a type provides an ObsIOBackend type through BackendTraits
 *
 * @details This concept is used to verify that a backend tag type correctly
 * defines an ObsIOBackend type through the BackendTraits specialization.
 * The ObsIOBackend is responsible for handling the observation I/O operations.
 */
template <typename T>
concept HasObsIOBackend =
    requires { typename traits::BackendTraits<T>::ObsIOBackend; };

/**
 * @brief Concept that checks if a type has a constructor from a ConfigBackend
 *
 * @details This concept verifies that a type T has a constructor that takes a
 * ConfigBackend parameter, allowing components to be constructed from
 * configuration. This enables consistent initialization patterns across the
 * framework.
 *
 * @tparam T The type to check for configuration constructor
 * @tparam ConfigBackend The configuration backend type
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
 * of objects. This is essential for operations like ensemble generation,
 * state perturbation, and various data assimilation algorithms.
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
 * parameters. This prevents the creation of objects in an invalid state.
 *
 * @tparam T The type to check for deleted default constructor
 */
template <typename T>
concept HasDeletedDefaultConstructor = !std::is_default_constructible_v<T>;

/**
 * @brief Concept to check if copy construction is deleted
 *
 * @details This concept verifies that a type cannot be copy-constructed,
 * which is a requirement for certain backends to prevent unintended copying.
 * This is particularly important for resource-heavy objects where copying
 * should be explicit through clone methods.
 *
 * @tparam T The type to check for deleted copy constructor
 */
template <typename T>
concept HasDeletedCopyConstructor = !std::is_copy_constructible_v<T>;

/**
 * @brief Concept to check if copy assignment is deleted
 *
 * @details This concept verifies that a type cannot be copy-assigned,
 * which is a requirement for certain backends to prevent unintended copying.
 * This complements the HasDeletedCopyConstructor concept to ensure complete
 * control over object copying.
 *
 * @tparam T The type to check for deleted copy assignment
 */
template <typename T>
concept HasDeletedCopyAssignment = !std::is_copy_assignable_v<T>;

/**
 * @brief Concept requiring initialization capability
 *
 * @details Verifies that a component can be initialized with a configuration
 * object, setting up its internal state before use. This allows for two-phase
 * construction patterns where an object is created and then initialized
 * separately.
 *
 * @tparam T The component type to check
 * @tparam ConfigBackend The configuration backend type
 */
template <typename T, typename ConfigBackend>
concept HasInitialize = requires(T& component, const ConfigBackend& config) {
  { component.initialize(config) } -> std::same_as<void>;
};

}  // namespace metada::framework