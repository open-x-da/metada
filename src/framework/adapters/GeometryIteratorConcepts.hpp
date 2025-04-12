/**
 * @file GeometryIteratorConcepts.hpp
 * @brief Concept definitions for geometry iterator classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the GeometryIterator adapter class. These concepts ensure that
 * backend iterator implementations provide all the necessary functionality
 * required by the geometry iteration operations.
 */

#pragma once

#include <concepts>
#include <iterator>
#include <type_traits>

namespace metada::framework {

/**
 * @brief Concept that defines requirements for a geometry iterator backend
 * implementation
 *
 * @details A valid geometry iterator backend implementation must provide:
 * - Dereference operator
 * - Pre/post increment operators
 * - Equality/inequality comparison operators
 * - Forward iterator behavior
 *
 * This concept ensures that the backend iterator can be used with the standard
 * iterator interface.
 *
 * @tparam T The geometry iterator backend implementation type
 */
template <typename T>
concept GeometryIteratorBackendImpl = requires(T it, const T& const_it) {
  // Dereference
  *it;
  *const_it;

  // Pre-increment
  { ++it } -> std::same_as<T&>;

  // Post-increment (result type not checked, just that it exists)
  it++;

  // Equality comparison
  { it == it } -> std::same_as<bool>;
  { const_it == const_it } -> std::same_as<bool>;

  // Inequality comparison
  { it != it } -> std::same_as<bool>;
  { const_it != const_it } -> std::same_as<bool>;
};

/**
 * @brief Concept that defines requirements for a geometry iterator backend tag
 * type
 *
 * @details A valid backend tag must:
 * - Provide a GeometryIteratorBackend type through BackendTraits
 * - Ensure the GeometryIteratorBackend type satisfies the
 * GeometryIteratorBackendImpl concept
 *
 * This concept constrains the template parameter of the GeometryIterator class,
 * ensuring that only valid backend iterator implementations can be used.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept GeometryIteratorBackendType =
    HasGeometryIteratorBackend<T> &&
    GeometryIteratorBackendImpl<
        typename traits::BackendTraits<T>::GeometryIteratorBackend>;

}  // namespace metada::framework