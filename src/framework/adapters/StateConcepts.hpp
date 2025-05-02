/**
 * @file StateConcepts.hpp
 * @brief Concept definitions for state classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the State adapter class. These concepts ensure that backend
 * implementations provide all the necessary functionality required by the
 * state operations.
 */

#pragma once

#include <concepts>
#include <memory>
#include <string>
#include <vector>

#include "CommonConcepts.hpp"

namespace metada::framework {

/**
 * @brief Concept that defines requirements for a state backend implementation
 *
 * @details A valid state backend implementation must provide:
 * - Data access methods (getData for both mutable and const access)
 * - Variable name access method
 * - Dimension information for each variable
 * - Support for cloning the state
 * - Configuration constructor and proper resource management
 *
 * @tparam T The state backend implementation type
 * @tparam ConfigBackend The configuration backend type
 */
template <typename T, typename ConfigBackend>
concept StateBackendImpl =
    requires(T& t, const T& ct, const std::string& varName,
             const ConfigBackend& config) {
      // Data access
      { t.getData() } -> std::same_as<void*>;
      // TODO: Add const data access
      //{ ct.getData() } -> std::same_as<const void*>;

      // Variable information
      {
        ct.getVariableNames()
      } -> std::same_as<const std::vector<std::string>&>;
      { ct.getDimensions(varName) } -> std::same_as<const std::vector<size_t>&>;

      // Construction and cloning
      { T(config) } -> std::same_as<T>;
      { ct.clone() } -> std::convertible_to<std::unique_ptr<T>>;

      // Resource management constraints
      requires HasDeletedDefaultConstructor<T>;
      requires HasDeletedCopyConstructor<T>;
      requires HasDeletedCopyAssignment<T>;
    };

/**
 * @brief Concept that defines requirements for a state backend tag type
 *
 * @details A valid backend tag must:
 * - Provide both StateBackend and ConfigBackend types via BackendTraits
 * - Ensure the StateBackend type satisfies the StateBackendImpl concept
 *   when paired with the ConfigBackend type
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept StateBackendType =
    HasStateBackend<T> && HasConfigBackend<T> &&
    StateBackendImpl<typename traits::BackendTraits<T>::StateBackend,
                     typename traits::BackendTraits<T>::ConfigBackend>;

}  // namespace metada::framework