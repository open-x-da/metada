#pragma once

#include <concepts>
#include <memory>
#include <string>
#include <vector>

#include "CommonConcepts.hpp"

namespace metada::framework {

/**
 * @brief Concept that defines requirements for an observation operator backend
 * implementation
 *
 * @details A valid observation operator backend implementation must provide:
 * - Constructor that accepts a configuration backend
 * - Initialization method that configures the operator
 * - Method to check initialization status
 * - Apply method for mapping state to observation space
 * - Methods to retrieve required state and observation variables
 * - Proper resource management (deleted default constructor, copy constructor,
 * and copy assignment)
 *
 * This concept ensures that any observation operator backend can be properly
 * initialized, can transform state variables into observation space, and
 * declares which variables it requires from both state and observation spaces.
 *
 * @tparam T The observation operator backend implementation type
 * @tparam ConfigBackend The configuration backend type
 * @tparam StateBackend The state backend type
 * @tparam ObsBackend The observation backend type
 */
template <typename T, typename ConfigBackend, typename StateBackend,
          typename ObsBackend>
concept ObsOperatorBackendImpl =
    requires(T& t, const T& ct, const ConfigBackend& config,
             const StateBackend& state, ObsBackend& obs) {
      // Constructor
      { T(config) } -> std::same_as<T>;

      // Initialization
      { t.initialize(config) } -> std::same_as<void>;
      { t.isInitialized() } -> std::same_as<bool>;

      // Apply method
      { ct.apply(state, obs) } -> std::same_as<void>;

      // Required variables
      {
        ct.getRequiredStateVars()
      } -> std::same_as<const std::vector<std::string>&>;
      {
        ct.getRequiredObsVars()
      } -> std::same_as<const std::vector<std::string>&>;

      // Resource management constraints
      requires HasDeletedDefaultConstructor<T>;
      requires HasDeletedCopyConstructor<T>;
      requires HasDeletedCopyAssignment<T>;
    };

/**
 * @brief Concept that defines requirements for an observation operator backend
 * tag type
 *
 * @details A valid backend tag must:
 * - Provide an ObsOperatorBackend type via BackendTraits
 * - Provide ConfigBackend, StateBackend, and ObservationBackend types
 * - Ensure the ObsOperatorBackend type satisfies the ObsOperatorBackendImpl
 * concept with the corresponding backend types
 *
 * This concept is used to validate that a backend tag provides all necessary
 * type information and that the associated observation operator backend
 * implementation meets the requirements defined in ObsOperatorBackendImpl.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept ObsOperatorBackendType =
    HasObsOperatorBackend<T> && HasConfigBackend<T> && HasStateBackend<T> &&
    HasObservationBackend<T> &&
    ObsOperatorBackendImpl<
        typename traits::BackendTraits<T>::ObsOperatorBackend,
        typename traits::BackendTraits<T>::ConfigBackend,
        typename traits::BackendTraits<T>::StateBackend,
        typename traits::BackendTraits<T>::ObservationBackend>;

}  // namespace metada::framework