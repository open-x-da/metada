/**
 * @file WRFGeometryIterator.hpp
 * @brief WRF geometry iterator implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <cstddef>

namespace metada::backends::wrf {

// Forward declaration
class WRFGeometry;

/**
 * @brief Iterator for the WRF geometry grid
 *
 * Allows iteration over all grid points in the WRF geometry.
 * Implements the required iterator interface for use with STL algorithms.
 */
class WRFGeometry::iterator {
 public:
  // Iterator traits
  using iterator_category = std::forward_iterator_tag;
  using value_type = std::tuple<size_t, size_t, size_t>;  // i, j, k indices
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  /**
   * @brief Default constructor - creates an invalid iterator
   */
  iterator() : geometry_(nullptr), i_(0), j_(0), k_(0), index_(0) {}

  /**
   * @brief Constructor for creating a valid iterator
   *
   * @param geometry Pointer to the WRF geometry backend
   * @param index Linear index into the grid (0 for begin, grid size for end)
   */
  iterator(const WRFGeometry* geometry, size_t index);

  /**
   * @brief Dereference operator
   *
   * @return Tuple containing (i, j, k) indices of current position
   */
  value_type operator*() const;

  /**
   * @brief Pre-increment operator
   *
   * @return Reference to this iterator after incrementing
   */
  iterator& operator++();

  /**
   * @brief Post-increment operator
   *
   * @return Copy of iterator before incrementing
   */
  iterator operator++(int);

  /**
   * @brief Equality comparison
   *
   * @param other Iterator to compare with
   * @return True if iterators point to the same position
   */
  bool operator==(const iterator& other) const;

  /**
   * @brief Inequality comparison
   *
   * @param other Iterator to compare with
   * @return True if iterators point to different positions
   */
  bool operator!=(const iterator& other) const;

  /**
   * @brief Get i index (west-east)
   *
   * @return Current i index
   */
  size_t i() const { return i_; }

  /**
   * @brief Get j index (south-north)
   *
   * @return Current j index
   */
  size_t j() const { return j_; }

  /**
   * @brief Get k index (bottom-top)
   *
   * @return Current k index
   */
  size_t k() const { return k_; }

 private:
  const Geometry* geometry_;  // Pointer to parent geometry
  size_t i_;                  // West-east index
  size_t j_;                  // South-north index
  size_t k_;                  // Bottom-top index
  size_t index_;              // Linear index into the grid
};

/**
 * @brief Const iterator for the WRF geometry grid
 */
class WRFGeometry::const_iterator : public WRFGeometry::iterator {
 public:
  // Inherit constructors from base class
  using iterator::iterator;
};

}  // namespace metada::backends::wrf