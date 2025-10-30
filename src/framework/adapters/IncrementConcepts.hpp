/**
 * @file IncrementConcepts.hpp
 * @brief Concept definitions for increment classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the Increment adapter class. These concepts ensure that backend
 * implementations provide all the necessary functionality required for
 * increment operations in variational data assimilation systems.
 */

#pragma once

#include <concepts>
#include <vector>

#include "CommonConcepts.hpp"

namespace metada::framework {

/**
 * @brief Concept that defines requirements for an increment backend
 * implementation
 *
 * @details A valid increment backend implementation must provide:
 * - Vector space operations (zero, scale, axpy)
 * - Inner product operations (dot product, norm)
 * - Arithmetic operators (+=, -=, *=, /=)
 * - Geometry access
 * - Generic data access (getData) for testing and gradient checks
 * - Randomization for testing
 * - Construction from geometry backend
 *
 * These requirements enable the Increment adapter to perform vector space
 * operations needed in variational data assimilation algorithms without
 * exposing backend-specific implementation details such as grid structure,
 * field layout, or data representation.
 *
 * @tparam T The increment backend implementation type
 * @tparam GeometryBackend The geometry backend type
 */
template <typename T, typename GeometryBackend>
concept IncrementBackendImpl =
    requires(T& t, const T& ct, const T& other, double alpha,
             const GeometryBackend& geometry) {
      // Construction from geometry
      { T(geometry) } -> std::same_as<T>;

      // Vector space operations
      { t.zero() } -> std::same_as<void>;
      { t.scale(alpha) } -> std::same_as<void>;
      { t.axpy(alpha, other) } -> std::same_as<void>;

      // Inner product operations
      { ct.dot(other) } -> std::convertible_to<double>;
      { ct.norm() } -> std::convertible_to<double>;

      // Arithmetic operators
      { t += other } -> std::same_as<T&>;
      { t -= other } -> std::same_as<T&>;
      { t *= alpha } -> std::same_as<T&>;
      { t /= alpha } -> std::same_as<T&>;

      // Geometry access
      { ct.geometry() } -> std::convertible_to<const GeometryBackend&>;

      // Generic data access for testing and gradient checks
      { ct.getData() } -> std::convertible_to<std::vector<double>>;

      // Randomization for testing
      { t.randomize() } -> std::same_as<void>;
    };

/**
 * @brief Concept that defines requirements for an increment backend tag type
 *
 * @details A valid backend tag must:
 * - Provide IncrementBackend and GeometryBackend types via BackendTraits
 * - Ensure the IncrementBackend type satisfies the IncrementBackendImpl concept
 *   when paired with the GeometryBackend type
 *
 * This concept is used by the Increment adapter to validate that a backend
 * implementation provides all necessary functionality at compile time.
 *
 * @tparam T The backend tag type to check
 *
 * @see IncrementBackendImpl
 * @see HasIncrementBackend
 * @see HasGeometryBackend
 */
template <typename T>
concept IncrementBackendType =
    HasIncrementBackend<T> && HasGeometryBackend<T> &&
    IncrementBackendImpl<typename traits::BackendTraits<T>::IncrementBackend,
                         typename traits::BackendTraits<T>::GeometryBackend>;

}  // namespace metada::framework
