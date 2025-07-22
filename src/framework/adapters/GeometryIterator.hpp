/**
 * @file GeometryIterator.hpp
 * @brief Iterator class for traversing geometry grid points
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This header provides an iterator class that enables traversal of grid points
 * in a Geometry object. The GeometryIterator class wraps backend iterator
 * implementations and provides a consistent interface across different
 * backends.
 *
 * The GeometryIterator class is designed to:
 * - Provide standard iterator interface (forward iterator)
 * - Delegate operations to the backend iterator
 * - Enable range-based for loop usage with Geometry objects
 * - Support proper move and copy semantics
 *
 * @see Geometry
 * @see GeometryBackendType
 */

#pragma once
#include <iterator>  // for std::forward_iterator_tag

#include "BackendTraits.hpp"             // For BackendTraits
#include "ConfigConcepts.hpp"            // For ConfigBackendType concept
#include "GeometryIteratorConcepts.hpp"  // For GeometryIteratorBackendType concept

namespace metada::framework {

// Forward declaration
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

/**
 * @brief Iterator class for traversing geometry grid points
 *
 * @details
 * This class template wraps a geometry iterator backend implementation and
 * provides a standard iterator interface for traversing grid points. It
 * delegates operations to the backend while providing a consistent interface
 * across different backends.
 *
 * The iterator satisfies the requirements of a ForwardIterator, allowing it to
 * be used in standard algorithms and range-based for loops.
 *
 * The iterator supports both copy and move operations to ensure proper resource
 * management and consistency with backend implementations.
 *
 * Example usage:
 * @code
 * Geometry<Backend> geometry(config);
 * for (auto& point : geometry) {
 *     // Process each grid point
 * }
 * @endcode
 *
 * @tparam BackendTag The tag type that defines the geometry backend through
 * BackendTraits
 *
 * @see Geometry
 * @see GeometryBackendType
 */
template <typename BackendTag>
  requires GeometryIteratorBackendType<BackendTag>
class GeometryIterator {
 public:
  /** @brief Backend iterator type from traits */
  using GeometryIteratorBackend =
      typename traits::BackendTraits<BackendTag>::GeometryIteratorBackend;
  using iterator_category = std::forward_iterator_tag;
  using value_type =
      typename std::iterator_traits<GeometryIteratorBackend>::value_type;
  using difference_type =
      typename std::iterator_traits<GeometryIteratorBackend>::difference_type;
  using pointer =
      typename std::iterator_traits<GeometryIteratorBackend>::pointer;
  using reference =
      typename std::iterator_traits<GeometryIteratorBackend>::reference;

  /** @brief Default constructor (for end/sentinels) */
  GeometryIterator() = default;

  /** @brief Copy constructor */
  GeometryIterator(const GeometryIterator&) = default;

  /** @brief Copy assignment operator */
  GeometryIterator& operator=(const GeometryIterator&) = default;

  /** @brief Move constructor */
  GeometryIterator(GeometryIterator&&) noexcept = default;

  /** @brief Move assignment operator */
  GeometryIterator& operator=(GeometryIterator&&) noexcept = default;

  /** @brief Destructor */
  ~GeometryIterator() = default;

  /** @brief Constructor from a backend iterator */
  explicit GeometryIterator(GeometryIteratorBackend iter)
      : iter_(std::move(iter)) {}

  /**
   * @brief Dereference operator to access the current grid point
   * @return Reference to the current grid point
   */
  reference operator*() const { return *iter_; }

  /**
   * @brief Arrow operator to access the current grid point
   * @return Pointer to the current grid point
   */
  pointer operator->() const { return iter_.operator->(); }

  /**
   * @brief Pre-increment operator to advance to the next grid point
   * @return Reference to this iterator after advancement
   */
  GeometryIterator& operator++() {
    ++iter_;
    return *this;
  }

  /**
   * @brief Post-increment operator to advance to the next grid point
   * @return Value of this iterator before advancement
   */
  GeometryIterator operator++(int) {
    GeometryIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  /**
   * @brief Equality comparison operator
   * @param other Iterator to compare with
   * @return True if iterators point to the same position
   */
  bool operator==(const GeometryIterator& other) const {
    return iter_ == other.iter_;
  }

  /**
   * @brief Inequality comparison operator
   * @param other Iterator to compare with
   * @return True if iterators point to different positions
   */
  bool operator!=(const GeometryIterator& other) const {
    return iter_ != other.iter_;
  }

  /** @brief Get the backend iterator */
  GeometryIteratorBackend& backend() { return iter_; }

  /** @brief Get the backend iterator (const version) */
  const GeometryIteratorBackend& backend() const { return iter_; }

 private:
  /** @brief The wrapped backend iterator */
  GeometryIteratorBackend iter_;
};

}  // namespace metada::framework
