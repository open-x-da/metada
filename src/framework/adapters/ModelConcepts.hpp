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

// Individual concepts for model backend requirements

/**
 * @brief Concept requiring reset capability
 *
 * @details Ensures that a model can be reset to its initial state,
 * allowing it to be reused without reconstruction.
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
 */
template <typename T, typename StateBackend>
concept HasRun =
    requires(T& model, const StateBackend& initialState,
             StateBackend& finalState, double startTime, double endTime) {
      {
        model.run(initialState, finalState, startTime, endTime)
      } -> std::same_as<void>;
    };

/**
 * @brief Combined concept for all model backend requirements
 *
 * @details
 * This concept specifies the complete set of operations that must be supported
 * by any type that will be used as a model backend. It includes initialization,
 * parameter access, execution, and resource management operations.
 *
 * The concept provides compile-time validation of the required methods
 * and their signatures, ensuring that backends can be properly used with the
 * Model class.
 *
 * @tparam T The type to check against the ModelBackendType requirements
 * @tparam ConfigBackend The configuration backend type used for initialization
 * @tparam StateBackend The state backend type used for model execution
 */
template <typename T, typename ConfigBackend, typename StateBackend>
concept ModelBackendType =
    HasConfigConstructor<T, ConfigBackend> && HasInitialize<T, ConfigBackend> &&
    HasReset<T> && HasFinalize<T> && HasGetParameter<T> && HasSetParameter<T> &&
    HasRun<T, StateBackend> && HasDeletedCopyConstructor<T> &&
    HasDeletedCopyAssignment<T>;

}  // namespace metada::framework