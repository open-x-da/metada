/**
 * @file ModelConcepts.hpp
 * @brief Concepts defining requirements for model backend types
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that specify the requirements
 * for types to be used as model backends in the framework.
 */

#pragma once

#include <concepts>
#include <string>

namespace metada::framework {

// Individual concepts for model backend requirements

/**
 * @brief Concept that checks if a type has a ModelBackend type member
 *
 * @details This concept verifies that a type T provides a ModelBackend type
 * through the BackendTraits specialization.
 */
template <typename T>
concept HasModelBackend =
    requires { typename traits::BackendTraits<T>::ModelBackend; };

/**
 * @brief Concept that checks if a type has a StateBackend type member
 *
 * @details This concept verifies that a type T provides a StateBackend type
 * through the BackendTraits specialization.
 */
template <typename T>
concept HasStateBackend =
    requires { typename traits::BackendTraits<T>::StateBackend; };

/**
 * @brief Concept requiring initialization capability
 */
template <typename T, typename ConfigBackend>
concept HasInitialize = requires(T& model, const ConfigBackend& config) {
  { model.initialize(config) } -> std::same_as<void>;
};

/**
 * @brief Concept requiring reset capability
 */
template <typename T, typename ConfigBackend>
concept HasReset = requires(T& model, const ConfigBackend& config) {
  { model.reset() } -> std::same_as<void>;
};

/**
 * @brief Concept requiring finalization capability
 */
template <typename T, typename ConfigBackend>
concept HasFinalize = requires(T& model, const ConfigBackend& config) {
  { model.finalize() } -> std::same_as<void>;
};

/**
 * @brief Concept requiring parameter get capability
 */
template <typename T, typename ConfigBackend>
concept HasGetParameter = requires(const T& model, const std::string& name) {
  { model.getParameter(name) } -> std::convertible_to<std::string>;
};

/**
 * @brief Concept requiring parameter set capability
 */
template <typename T, typename ConfigBackend>
concept HasSetParameter =
    requires(T& model, const std::string& name, const std::string& value) {
      { model.setParameter(name, value) } -> std::same_as<void>;
    };

/**
 * @brief Concept requiring model run capability
 */
template <typename T, typename ConfigBackend, typename StateBackend>
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
 * @tparam T The type to check against the ModelBackendType requirements
 */
template <typename T, typename ConfigBackend, typename StateBackend>
concept ModelBackendType =
    HasInitialize<T, ConfigBackend> && HasReset<T, ConfigBackend> &&
    HasFinalize<T, ConfigBackend> && HasGetParameter<T, ConfigBackend> &&
    HasSetParameter<T, ConfigBackend> &&
    HasRun<T, ConfigBackend, StateBackend> && HasStateBackend<T> &&
    HasModelBackend<T>;

}  // namespace metada::framework