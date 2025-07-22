/**
 * @file ObservationIteratorConcepts.hpp
 * @brief Concept definitions for observation iterator classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the ObservationIterator adapter class. These concepts ensure that
 * backend iterator implementations provide all the necessary functionality
 * required by the observation iteration operations.
 */

#pragma once

#include <concepts>
#include <iterator>
#include <type_traits>

#include "CommonConcepts.hpp"

namespace metada::framework {

/**
 * @brief Concept that defines requirements for an observation iterator backend
 * implementation
 *
 * @details A valid observation iterator backend implementation must provide:
 * - Dereference operator
 * - Pre/post increment operators
 * - Equality/inequality comparison operators
 * - Forward iterator behavior
 * - Proper resource management with default constructor
 *
 * This concept ensures that the backend iterator can be used with the standard
 * iterator interface and follows the framework's resource management patterns.
 *
 * @tparam T The observation iterator backend implementation type
 *
 * @see ObservationIteratorBackendType
 */
template <typename T>
concept ObservationIteratorBackendImpl = requires(T it, const T& const_it) {
  // Type aliases
  typename T::iterator_category;
  typename T::value_type;
  typename T::difference_type;
  typename T::pointer;
  typename T::reference;

  // Default constructor
  { T() };

  // Dereference
  *it;
  *const_it;
  it.operator->();
  const_it.operator->();

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
 * @brief Concept that defines requirements for an observation iterator backend
 * tag type
 *
 * @details A valid backend tag must:
 * - Provide an ObservationIteratorBackend type through BackendTraits
 * - Ensure the ObservationIteratorBackend type satisfies the
 *   ObservationIteratorBackendImpl concept
 *
 * This concept constrains the template parameter of the ObservationIterator
 * class, ensuring that only valid backend iterator implementations can be used.
 * It provides compile-time validation of backend compatibility.
 *
 * @tparam T The backend tag type to check
 *
 * @see HasObservationIteratorBackend
 * @see ObservationIteratorBackendImpl
 */
template <typename T>
concept ObservationIteratorBackendType =
    HasObservationIteratorBackend<T> &&
    ObservationIteratorBackendImpl<
        typename traits::BackendTraits<T>::ObservationIteratorBackend>;

}  // namespace metada::framework