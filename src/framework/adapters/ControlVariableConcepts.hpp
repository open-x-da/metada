/**
 * @file ControlVariableConcepts.hpp
 * @brief Concept definitions for control variable classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the ControlVariable adapter class. These concepts ensure that
 * backend implementations provide all the necessary functionality required for
 * control variable operations in incremental variational data assimilation.
 *
 * In incremental variational DA, the control variable v is the actual
 * optimization variable, which is transformed to state-space increments δx
 * through a transformation operator (e.g., B^(1/2) in WRFDA).
 */

#pragma once

#include <concepts>
#include <vector>

#include "BackendTraits.hpp"

namespace metada::framework {

/**
 * @brief Concept that defines requirements for a control variable backend
 * implementation
 *
 * @details A valid control variable backend implementation must provide:
 * - Vector space operations (zero, scale, axpy) in control space
 * - Inner product operations (dot product, norm) in control space
 * - Arithmetic operators (+=, -=, *=, /=) in control space
 * - Geometry access (control variables are still tied to model geometry)
 * - Generic data access (getData) for optimization algorithms
 * - Randomization for testing
 * - Construction from geometry backend
 *
 * These requirements enable the ControlVariable adapter to perform vector
 * space operations needed by optimization algorithms (L-BFGS, CG, etc.)
 * without exposing backend-specific implementation details.
 *
 * Note: Control variables live in control space, which may have different
 * dimensions and structure than state space. The transformation between
 * control space and state space is handled by ControlVariableBackend.
 *
 * @tparam T The control variable backend implementation type
 * @tparam GeometryBackend The geometry backend type
 */
template <typename T, typename GeometryBackend>
concept ControlVariableBackendImpl =
    requires(T& t, const T& ct, const T& other, double alpha,
             const GeometryBackend& geometry) {
      // Construction from geometry
      { T(geometry) } -> std::same_as<T>;

      // Vector space operations in control space
      { t.zero() } -> std::same_as<void>;
      { t.scale(alpha) } -> std::same_as<void>;
      { t.axpy(alpha, other) } -> std::same_as<void>;

      // Inner product operations in control space
      { ct.dot(other) } -> std::convertible_to<double>;
      { ct.norm() } -> std::convertible_to<double>;

      // Arithmetic operators in control space
      { t += other } -> std::same_as<T&>;
      { t -= other } -> std::same_as<T&>;
      { t *= alpha } -> std::same_as<T&>;
      { t /= alpha } -> std::same_as<T&>;

      // Geometry access
      { ct.geometry() } -> std::convertible_to<const GeometryBackend&>;

      // Generic data access for optimization algorithms
      { ct.getData() } -> std::convertible_to<std::vector<double>>;

      // Randomization for testing
      { t.randomize() } -> std::same_as<void>;

      // Set from vector (for initialization from optimization algorithm)
      { t.setFromVector(std::vector<double>{}) } -> std::same_as<void>;
    };

/**
 * @brief Concept that defines requirements for a control variable backend tag
 * type
 *
 * @details A valid backend tag must:
 * - Provide ControlVariableBackend and GeometryBackend types via BackendTraits
 * - Ensure the ControlVariableBackend type satisfies the
 *   ControlVariableBackendImpl concept when paired with the GeometryBackend
 *   type
 *
 * This concept is used by the ControlVariable adapter to validate that a
 * backend implementation provides all necessary functionality at compile time.
 *
 * @tparam T The backend tag type to check
 *
 * @see ControlVariableBackendImpl
 * @see HasControlVariableBackend
 * @see HasGeometryBackend
 */
template <typename T>
concept ControlVariableBackendType =
    HasControlVariableBackend<T> && HasGeometryBackend<T> &&
    ControlVariableBackendImpl<
        typename traits::BackendTraits<T>::ControlVariableBackend,
        typename traits::BackendTraits<T>::GeometryBackend>;

}  // namespace metada::framework
