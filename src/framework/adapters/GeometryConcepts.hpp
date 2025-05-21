/**
 * @file GeometryConcepts.hpp
 * @brief Concept definitions for geometry classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the Geometry adapter class. These concepts ensure that backend
 * implementations provide all the necessary functionality required by the
 * geometry operations.
 */

#pragma once

#include <concepts>
#include <cstddef>
#include <iterator>
#include <memory>

#include "CommonConcepts.hpp"

namespace metada::framework {

/**
 * @brief Concept that defines requirements for a geometry backend
 * implementation
 *
 * @details A valid geometry backend implementation must provide:
 * - Iteration capabilities (begin/end)
 * - Size information
 * - Periodicity queries
 * - Cloning capability
 * - Halo exchange functionality
 * - Proper resource management with deleted default constructor, copy
 * constructor, and copy assignment operator
 *
 * This concept is used to ensure that backend implementations provide
 * all the necessary functionality required by the Geometry class.
 *
 * @tparam T The geometry backend implementation type
 * @tparam SB The state backend type to use for halo exchange
 * @tparam ConfigBackend The configuration backend type
 *
 * @see HasDeletedDefaultConstructor
 * @see HasDeletedCopyConstructor
 * @see HasDeletedCopyAssignment
 */
template <typename T, typename SB, typename ConfigBackend>
concept GeometryBackendImpl = requires(T t, T& t_ref, const T& t_const,
                                       SB& state, const ConfigBackend& config) {
  // Iteration
  { t.begin() } -> std::same_as<typename T::iterator>;
  { t.end() } -> std::same_as<typename T::iterator>;
  { t_const.begin() } -> std::same_as<typename T::const_iterator>;
  { t_const.end() } -> std::same_as<typename T::const_iterator>;

  // Size information
  { t.totalGridSize() } -> std::convertible_to<std::size_t>;

  // Periodicity checks
  { t.isPeriodicX() } -> std::same_as<bool>;
  { t.isPeriodicY() } -> std::same_as<bool>;
  { t.isPeriodicZ() } -> std::same_as<bool>;

  // Cloning
  { t_const.clone() } -> std::same_as<T>;

  // Halo exchange functionality (required)
  { t.haloExchange(state) } -> std::same_as<void>;

  // Construction and resource management
  { T(config) } -> std::same_as<T>;
  requires HasDeletedDefaultConstructor<T>;
  requires HasDeletedCopyConstructor<T>;
  requires HasDeletedCopyAssignment<T>;
};

/**
 * @brief Concept that defines requirements for a geometry backend tag type
 *
 * @details A valid backend tag must:
 * - Provide a GeometryBackend type through BackendTraits
 * - Provide a StateBackend type through BackendTraits
 * - Provide a ConfigBackend type through BackendTraits
 * - Ensure the GeometryBackend type satisfies the GeometryBackendImpl concept
 *   when paired with the StateBackend and ConfigBackend types
 *
 * This concept constrains the template parameter of the Geometry class,
 * ensuring that only valid backend configurations can be used. It provides
 * compile-time validation of backend compatibility.
 *
 * @tparam T The backend tag type to check
 *
 * @see HasGeometryBackend
 * @see HasStateBackend
 * @see HasConfigBackend
 * @see GeometryBackendImpl
 */
template <typename T>
concept GeometryBackendType =
    HasGeometryBackend<T> && HasStateBackend<T> && HasConfigBackend<T> &&
    GeometryBackendImpl<typename traits::BackendTraits<T>::GeometryBackend,
                        typename traits::BackendTraits<T>::StateBackend,
                        typename traits::BackendTraits<T>::ConfigBackend>;

}  // namespace metada::framework