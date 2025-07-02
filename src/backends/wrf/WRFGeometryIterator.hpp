/**
 * @file WRFGeometryIterator.hpp
 * @brief WRF geometry iterator implementation for grid traversal
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details This file contains the iterator implementation for traversing
 * WRF geometry grids. The iterator provides a forward iterator interface
 * for accessing grid points sequentially and supports standard STL iterator
 * operations for integration with algorithms and range-based loops.
 */

#pragma once

#include <cstddef>

namespace metada::backends::wrf {

// Forward declaration
template <typename ConfigBackend>
class WRFGeometry;

/**
 * @brief Forward iterator for traversing WRF geometry grid points
 *
 * @details This class implements a forward iterator that traverses through
 * all grid points in a WRF geometry. The iterator follows STL conventions
 * and can be used with standard algorithms and range-based for loops.
 *
 * The iterator traverses the grid in row-major order:
 * - First dimension: X (west-east)
 * - Second dimension: Y (south-north)
 * - Third dimension: Z (bottom-top)
 *
 * Each iteration returns a Location object containing the geographic
 * coordinates (latitude, longitude, vertical level) of the current grid point.
 *
 * @tparam ConfigBackend Configuration backend type for WRF geometry
 *
 * @note This iterator is designed for read-only traversal of geometry grids
 * @see WRFGeometry
 * @see framework::Location
 */
template <typename ConfigBackend>
class WRFGeometryIterator {
 public:
  ///@{ @name Iterator Type Definitions
  /**
   * @brief Standard iterator type aliases for STL compatibility
   *
   * These type definitions enable the iterator to be used with STL algorithms
   * and provide proper iterator traits for template metaprogramming.
   */
  using iterator_category =
      std::forward_iterator_tag;           ///< Forward iterator category
  using value_type = framework::Location;  ///< Type of objects pointed to
  using difference_type = std::ptrdiff_t;  ///< Type for iterator differences
  using pointer = value_type*;             ///< Pointer to value type
  using reference = value_type&;           ///< Reference to value type
  ///@}

  ///@{ @name Construction and Destruction
  /**
   * @brief Default constructor is deleted
   *
   * @details Default construction is disabled because a valid iterator
   * requires a geometry reference and position information.
   */
  WRFGeometryIterator() = delete;

  /**
   * @brief Construct iterator for specific geometry and position
   *
   * @details Creates an iterator positioned at the specified linear index
   * within the given WRF geometry grid. The constructor automatically
   * calculates 3D grid indices from the linear index.
   *
   * @param[in] geometry Pointer to the WRF geometry backend
   * @param[in] index Linear index into the grid (0 for begin, grid size for
   * end)
   *
   * @note For end iterators, use the total grid size as the index
   * @note The geometry pointer must remain valid for the iterator's lifetime
   */
  WRFGeometryIterator(const WRFGeometry<ConfigBackend>* geometry, size_t index)
      : geometry_(geometry), index_(index) {
    if (geometry_ && index_ < geometry_->totalGridSize()) {
      // Calculate 3D indices from linear index
      const size_t nx = geometry_->x_dim();
      const size_t ny = geometry_->y_dim();

      // Using row-major order: index = k*nx*ny + j*nx + i
      k_ = index_ / (nx * ny);
      const size_t remainder = index_ % (nx * ny);
      j_ = remainder / nx;
      i_ = remainder % nx;
    } else {
      // End iterator or invalid
      i_ = 0;
      j_ = 0;
      k_ = 0;
    }
  }
  ///@}

  ///@{ @name Iterator Operations
  /**
   * @brief Dereference operator to access current grid point
   *
   * @details Returns a Location object containing the geographic coordinates
   * of the current grid point. The location includes latitude, longitude,
   * and vertical level information.
   *
   * @return value_type Location object at the current grid position
   * @throws std::out_of_range If iterator is at end position
   *
   * @note The returned Location object is created on-demand from grid
   * coordinates
   */
  value_type operator*() const { return geometry_->getLocation(index_); }

  /**
   * @brief Pre-increment operator to advance to next grid point
   *
   * @details Advances the iterator to the next grid point in row-major order.
   * Updates both the linear index and the 3D grid indices (i, j, k).
   * The traversal order is: X (west-east) first, then Y (south-north),
   * then Z (bottom-top).
   *
   * @return WRFGeometryIterator& Reference to this iterator after advancement
   *
   * @note Incrementing an end iterator results in undefined behavior
   */
  WRFGeometryIterator& operator++() {
    if (geometry_ && index_ < geometry_->totalGridSize()) {
      // Increment linear index
      ++index_;

      if (index_ < geometry_->totalGridSize()) {
        // Recalculate 3D indices
        const size_t nx = geometry_->x_dim();
        const size_t ny = geometry_->y_dim();

        // Increment i (west-east) first
        ++i_;
        if (i_ >= nx) {
          // Wrap to next j (south-north)
          i_ = 0;
          ++j_;
          if (j_ >= ny) {
            // Wrap to next k (bottom-top)
            j_ = 0;
            ++k_;
          }
        }
      }
    }
    return *this;
  }

  /**
   * @brief Post-increment operator to advance to next grid point
   *
   * @details Advances the iterator to the next grid point but returns
   * a copy of the iterator before advancement. Less efficient than
   * pre-increment due to the copy operation.
   *
   * @return WRFGeometryIterator Copy of iterator before advancement
   *
   * @note Prefer pre-increment (++it) over post-increment (it++) for efficiency
   */
  WRFGeometryIterator operator++(int) {
    WRFGeometryIterator temp = *this;
    ++(*this);
    return temp;
  }
  ///@}

  ///@{ @name Comparison Operations
  /**
   * @brief Equality comparison operator
   *
   * @details Compares two iterators for equality. Iterators are considered
   * equal if they reference the same geometry and are at the same position.
   *
   * @param[in] other Iterator to compare with
   * @return bool True if iterators are equal, false otherwise
   */
  bool operator==(const WRFGeometryIterator& other) const {
    return (geometry_ == other.geometry_ && index_ == other.index_);
  }

  /**
   * @brief Inequality comparison operator
   *
   * @details Compares two iterators for inequality. Convenience function
   * that returns the negation of the equality comparison.
   *
   * @param[in] other Iterator to compare with
   * @return bool True if iterators are different, false otherwise
   */
  bool operator!=(const WRFGeometryIterator& other) const {
    return !(*this == other);
  }
  ///@}

  ///@{ @name Grid Index Access
  /**
   * @brief Get current X-dimension index (west-east)
   *
   * @details Returns the current position in the X dimension (west-east
   * direction). This corresponds to the longitude dimension in geographic
   * coordinates.
   *
   * @return size_t Current X-dimension index (0-based)
   */
  size_t i() const { return i_; }

  /**
   * @brief Get current Y-dimension index (south-north)
   *
   * @details Returns the current position in the Y dimension (south-north
   * direction). This corresponds to the latitude dimension in geographic
   * coordinates.
   *
   * @return size_t Current Y-dimension index (0-based)
   */
  size_t j() const { return j_; }

  /**
   * @brief Get current Z-dimension index (bottom-top)
   *
   * @details Returns the current position in the Z dimension (bottom-top
   * direction). This corresponds to the vertical level in meteorological
   * coordinates.
   *
   * @return size_t Current Z-dimension index (0-based)
   */
  size_t k() const { return k_; }
  ///@}

 private:
  ///@{ @name Member Variables
  const WRFGeometry<ConfigBackend>*
      geometry_;  ///< Pointer to parent geometry backend
  size_t i_;      ///< Current X-dimension index (west-east)
  size_t j_;      ///< Current Y-dimension index (south-north)
  size_t k_;      ///< Current Z-dimension index (bottom-top)
  size_t index_;  ///< Linear index into the grid
  ///@}
};

}  // namespace metada::backends::wrf