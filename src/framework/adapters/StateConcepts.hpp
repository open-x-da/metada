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
 * state operations in data assimilation systems, including data access,
 * variable information, arithmetic operations, and resource management.
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
 * - Vector arithmetic operations (add, dot product, norm)
 * - Equality comparison
 * - Zero initialization
 * - Construction from configuration and geometry backends
 * - Proper resource management
 *
 * These requirements enable the State adapter to perform common
 * state-space operations needed in data assimilation algorithms.
 *
 * @tparam T The state backend implementation type
 * @tparam ConfigBackend The configuration backend type
 * @tparam GeometryBackend The geometry backend type
 */
template <typename T, typename ConfigBackend, typename GeometryBackend>
concept StateBackendImpl =
    requires(T& t, const T& ct, const T& other, const std::string& varName,
             const ConfigBackend& config, const GeometryBackend& geometry) {
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
      { T(config, geometry) } -> std::same_as<T>;
      { ct.clone() } -> std::convertible_to<std::unique_ptr<T>>;

      // Vector arithmetic
      { t.zero() } -> std::same_as<void>;
      { t.add(other) } -> std::same_as<void>;
      { ct.dot(other) } -> std::convertible_to<double>;
      { ct.norm() } -> std::convertible_to<double>;

      // Comparison
      { ct.equals(other) } -> std::convertible_to<bool>;

      // Resource management constraints
      requires HasDeletedDefaultConstructor<T>;
      requires HasDeletedCopyConstructor<T>;
      requires HasDeletedCopyAssignment<T>;
    };

/**
 * @brief Concept that defines requirements for a state backend tag type
 *
 * @details A valid backend tag must:
 * - Provide StateBackend, ConfigBackend, and GeometryBackend types via
 * BackendTraits
 * - Ensure the StateBackend type satisfies the StateBackendImpl concept
 *   when paired with the ConfigBackend and GeometryBackend types
 *
 * This concept is used by the State adapter to validate that a backend
 * implementation provides all necessary functionality at compile time.
 *
 * @tparam T The backend tag type to check
 *
 * @see StateBackendImpl
 * @see HasStateBackend
 * @see HasConfigBackend
 * @see HasGeometryBackend
 */
template <typename T>
concept StateBackendType =
    HasStateBackend<T> && HasConfigBackend<T> && HasGeometryBackend<T> &&
    StateBackendImpl<typename traits::BackendTraits<T>::StateBackend,
                     typename traits::BackendTraits<T>::ConfigBackend,
                     typename traits::BackendTraits<T>::GeometryBackend>;

}  // namespace metada::framework