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
 *
 * @tparam T The geometry backend implementation type
 * @tparam SB The state backend type to use for halo exchange
 */
template <typename T, typename SB>
concept GeometryBackendImpl =
    requires(T t, T& t_ref, const T& t_const, SB& state) {
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
    };

/**
 * @brief Concept that defines requirements for a geometry backend tag type
 *
 * @details A valid backend tag must:
 * - Provide a GeometryBackend type through BackendTraits
 * - Provide a StateBackend type through BackendTraits
 * - Ensure the GeometryBackend type satisfies the GeometryBackendImpl concept
 *   when paired with the StateBackend type
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept GeometryBackendType =
    HasGeometryBackend<T> && HasStateBackend<T> &&
    GeometryBackendImpl<typename traits::BackendTraits<T>::GeometryBackend,
                        typename traits::BackendTraits<T>::StateBackend>;

}  // namespace metada::framework