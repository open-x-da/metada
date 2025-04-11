/**
 * @file ModelConcepts.hpp
 * @brief Concepts defining requirements for model backend types
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that specify the requirements
 * for types to be used as model backends in the framework. These concepts
 * ensure that model backends implement the necessary methods with the
 * correct signatures to be compatible with the framework's Model class.
 */

#pragma once

#include <concepts>
#include <string>

#include "CommonConcepts.hpp"

namespace metada::framework {

//-----------------------------------------------------------------------------
// Core model capability concepts
//-----------------------------------------------------------------------------

/**
 * @brief Concept requiring reset capability
 *
 * @details Ensures that a model can be reset to its initial state,
 * allowing it to be reused without reconstruction.
 *
 * @tparam T The type to check for reset capability
 */
template <typename T>
concept HasReset = requires(T& model) {
  { model.reset() } -> std::same_as<void>;
};

/**
 * @brief Concept requiring finalization capability
 *
 * @details Verifies that a model can properly clean up resources
 * when it is no longer needed.
 *
 * @tparam T The type to check for finalization capability
 */
template <typename T>
concept HasFinalize = requires(T& model) {
  { model.finalize() } -> std::same_as<void>;
};

/**
 * @brief Concept requiring parameter get capability
 *
 * @details Ensures that a model can retrieve parameter values by name,
 * returning them as strings for generic parameter access.
 *
 * @tparam T The type to check for parameter retrieval capability
 */
template <typename T>
concept HasGetParameter = requires(const T& model, const std::string& name) {
  { model.getParameter(name) } -> std::convertible_to<std::string>;
};

/**
 * @brief Concept requiring parameter set capability
 *
 * @details Verifies that a model can set parameter values by name,
 * accepting string values for generic parameter modification.
 *
 * @tparam T The type to check for parameter modification capability
 */
template <typename T>
concept HasSetParameter =
    requires(T& model, const std::string& name, const std::string& value) {
      { model.setParameter(name, value) } -> std::same_as<void>;
    };

/**
 * @brief Concept requiring model run capability
 *
 * @details Ensures that a model can execute a simulation from a given
 * initial state to produce a final state over a specified time interval.
 *
 * @tparam T The type to check for model execution capability
 * @tparam StateBackend The state backend type used for model execution
 */
template <typename T, typename StateBackend>
concept HasRun =
    requires(T& model, const StateBackend& initialState,
             StateBackend& finalState, double startTime, double endTime) {
      {
        model.run(initialState, finalState, startTime, endTime)
      } -> std::same_as<void>;
    };

//-----------------------------------------------------------------------------
// Component concepts
//-----------------------------------------------------------------------------

/**
 * @brief Concept that defines requirements for a model backend implementation
 *
 * @details A valid model backend implementation must:
 * - Be constructible from a ConfigBackend instance
 * - Support initialization from a ConfigBackend
 * - Provide reset and finalization methods
 * - Support parameter get/set operations
 * - Implement the run method for model execution
 * - Have deleted copy constructor, copy assignment, and default constructor
 *
 * This concept is used to ensure that backend implementations provide
 * all the necessary functionality required by the Model class.
 *
 * @tparam T The model backend implementation type
 * @tparam ConfigBackend The configuration backend type
 * @tparam StateBackend The state backend type
 *
 * @see HasConfigConstructor
 * @see HasInitialize
 * @see HasReset
 * @see HasFinalize
 * @see HasGetParameter
 * @see HasSetParameter
 * @see HasRun
 */
template <typename T, typename ConfigBackend, typename StateBackend>
concept ModelBackendImpl =
    HasConfigConstructor<T, ConfigBackend> && HasInitialize<T, ConfigBackend> &&
    HasReset<T> && HasFinalize<T> && HasGetParameter<T> && HasSetParameter<T> &&
    HasRun<T, StateBackend> && HasDeletedDefaultConstructor<T> &&
    HasDeletedCopyConstructor<T> && HasDeletedCopyAssignment<T>;

/**
 * @brief Concept that defines requirements for a model backend tag type
 *
 * @details A valid model backend tag must:
 * - Provide ModelBackend, ConfigBackend, and StateBackend types via
 * BackendTraits
 * - Ensure the ModelBackend type satisfies the ModelBackendImpl concept
 *   when paired with the ConfigBackend and StateBackend types
 *
 * This concept constrains the template parameter of the Model class,
 * ensuring that only valid backend configurations can be used. It provides
 * compile-time validation of backend compatibility.
 *
 * @tparam T The backend tag type to check
 *
 * @see HasModelBackend
 * @see HasConfigBackend
 * @see HasStateBackend
 * @see ModelBackendImpl
 */
template <typename T>
concept ModelBackendType =
    HasModelBackend<T> && HasConfigBackend<T> && HasStateBackend<T> &&
    ModelBackendImpl<typename traits::BackendTraits<T>::ModelBackend,
                     typename traits::BackendTraits<T>::ConfigBackend,
                     typename traits::BackendTraits<T>::StateBackend>;

}  // namespace metada::framework