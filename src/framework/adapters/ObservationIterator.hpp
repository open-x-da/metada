/**
 * @file ObservationIterator.hpp
 * @brief Iterator class for traversing observation data points
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This header provides an iterator class that enables traversal of observation
 * data points in an Observation object. The ObservationIterator class wraps
 * backend iterator implementations and provides a consistent interface across
 * different backends.
 *
 * The ObservationIterator class is designed to:
 * - Provide standard iterator interface (forward iterator)
 * - Delegate operations to the backend iterator
 * - Enable range-based for loop usage with Observation objects
 * - Support proper move and copy semantics
 *
 * @see Observation
 * @see ObservationBackendType
 */

#pragma once
#include <iterator>  // for std::forward_iterator_tag

#include "BackendTraits.hpp"                // For BackendTraits
#include "ConfigConcepts.hpp"               // For ConfigBackendType concept
#include "ObservationIteratorConcepts.hpp"  // For ObservationIteratorBackendType concept

namespace metada::framework {

// Forward declaration
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

/**
 * @brief Iterator class for traversing observation data points
 *
 * @details
 * This class template wraps an observation iterator backend implementation and
 * provides a standard iterator interface for traversing observation data
 * points. It delegates operations to the backend while providing a consistent
 * interface across different backends.
 *
 * The iterator satisfies the requirements of a ForwardIterator, allowing it to
 * be used in standard algorithms and range-based for loops.
 *
 * The iterator supports both copy and move operations to ensure proper resource
 * management and consistency with backend implementations.
 *
 * Example usage:
 * @code
 * Observation<Backend> observation(config);
 * for (auto& point : observation) {
 *     // Process each observation point
 * }
 * @endcode
 *
 * @tparam BackendTag The tag type that defines the observation backend through
 * BackendTraits
 *
 * @see Observation
 * @see ObservationBackendType
 */
template <typename BackendTag>
  requires ObservationIteratorBackendType<BackendTag>
class ObservationIterator {
 public:
  /** @brief Backend iterator type from traits */
  using ObservationIteratorBackend =
      typename traits::BackendTraits<BackendTag>::ObservationIteratorBackend;
  using iterator_category = std::forward_iterator_tag;
  using value_type =
      typename std::iterator_traits<ObservationIteratorBackend>::value_type;
  using difference_type = typename std::iterator_traits<
      ObservationIteratorBackend>::difference_type;
  using pointer =
      typename std::iterator_traits<ObservationIteratorBackend>::pointer;
  using reference =
      typename std::iterator_traits<ObservationIteratorBackend>::reference;

  /** @brief Default constructor (for end/sentinels) */
  ObservationIterator() = default;

  /** @brief Copy constructor */
  ObservationIterator(const ObservationIterator&) = default;

  /** @brief Copy assignment operator */
  ObservationIterator& operator=(const ObservationIterator&) = default;

  /** @brief Move constructor */
  ObservationIterator(ObservationIterator&&) noexcept = default;

  /** @brief Move assignment operator */
  ObservationIterator& operator=(ObservationIterator&&) noexcept = default;

  /** @brief Destructor */
  ~ObservationIterator() = default;

  /** @brief Constructor from a backend iterator */
  explicit ObservationIterator(ObservationIteratorBackend iter)
      : iter_(std::move(iter)) {}

  /**
   * @brief Dereference operator to access the current observation point
   * @return Reference to the current observation point
   */
  reference operator*() const { return *iter_; }

  /**
   * @brief Arrow operator to access the current observation point
   * @return Pointer to the current observation point
   */
  pointer operator->() const { return iter_.operator->(); }

  /**
   * @brief Pre-increment operator to advance to the next observation point
   * @return Reference to this iterator after advancement
   */
  ObservationIterator& operator++() {
    ++iter_;
    return *this;
  }

  /**
   * @brief Post-increment operator to advance to the next observation point
   * @return Value of this iterator before advancement
   */
  ObservationIterator operator++(int) {
    ObservationIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  /**
   * @brief Equality comparison operator
   * @param other Iterator to compare with
   * @return True if iterators point to the same position
   */
  bool operator==(const ObservationIterator& other) const {
    return iter_ == other.iter_;
  }

  /**
   * @brief Inequality comparison operator
   * @param other Iterator to compare with
   * @return True if iterators point to different positions
   */
  bool operator!=(const ObservationIterator& other) const {
    return iter_ != other.iter_;
  }

  /** @brief Get the backend iterator */
  ObservationIteratorBackend& backend() { return iter_; }

  /** @brief Get the backend iterator (const version) */
  const ObservationIteratorBackend& backend() const { return iter_; }

 private:
  /** @brief The wrapped backend iterator */
  ObservationIteratorBackend iter_;
};

}  // namespace metada::framework