/**
 * @file IGeometryIterator.hpp
 * @brief Interface for iterators traversing geometry grid points
 * @ingroup interfaces
 *
 * @details
 * This header provides the interface that geometry iterators must implement.
 * It defines the contract for traversing grid points in N-dimensional space.
 */

#pragma once

// Standard library includes
#include <iterator>
#include <memory>
#include <vector>

namespace metada::framework {

/**
 * @brief Interface for forward iterators traversing geometry grid points
 *
 * @tparam T Type of the coordinate value (typically double)
 */
template <typename T>
class IGeometryIterator {
 public:
  // Iterator traits typedefs for STL compatibility
  using iterator_category = std::forward_iterator_tag;
  using value_type = std::vector<T>;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  /**
   * @brief Virtual destructor
   */
  virtual ~IGeometryIterator() = default;

  /**
   * @brief Dereference operator returns current point coordinates
   */
  virtual const value_type& operator*() const = 0;

  /**
   * @brief Arrow operator for accessing point coordinates
   */
  virtual const value_type* operator->() const = 0;

  /**
   * @brief Pre-increment operator
   */
  virtual IGeometryIterator& operator++() = 0;

  /**
   * @brief Equality comparison
   */
  virtual bool operator==(const IGeometryIterator& other) const = 0;

  /**
   * @brief Inequality comparison
   */
  virtual bool operator!=(const IGeometryIterator& other) const = 0;

  /**
   * @brief Clone this iterator
   * @return Unique pointer to a new iterator instance
   */
  virtual std::unique_ptr<IGeometryIterator> clone() const = 0;
};

}  // namespace metada::framework