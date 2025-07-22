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

// Forward declaration of GridDimensionInfo
struct GridDimensionInfo;

namespace metada::backends::wrf {

// Forward declaration
template <typename ConfigBackend>
class WRFGeometry;

/**
 * @brief Forward iterator for traversing WRF geometry grid points
 *
 * @details This class implements a forward iterator that traverses through
 * all grid points in a WRF geometry across all grid types (unstaggered,
 * U-staggered, V-staggered, W-staggered). The iterator follows STL conventions
 * and can be used with standard algorithms and range-based for loops.
 *
 * The iterator traverses all grids in the following order:
 * 1. Unstaggered grid (mass points)
 * 2. U-staggered grid (U-wind component)
 * 3. V-staggered grid (V-wind component)
 * 4. W-staggered grid (W-wind component)
 *
 * Within each grid, the traversal order is row-major:
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
   * within the given WRF geometry across all grids. The constructor
   * automatically calculates which grid type and local position within that
   * grid.
   *
   * @param[in] geometry Pointer to the WRF geometry backend
   * @param[in] global_index Linear index across all grids (0 for begin, total
   * grid size for end)
   *
   * @note For end iterators, use the total grid size as the index
   * @note The geometry pointer must remain valid for the iterator's lifetime
   */
  WRFGeometryIterator(const WRFGeometry<ConfigBackend>* geometry,
                      size_t global_index)
      : geometry_(geometry), global_index_(global_index) {
    if (geometry_ && global_index_ < geometry_->totalGridSize()) {
      // Calculate which grid and local position
      auto [grid_type, local_index] =
          geometry_->getGridTypeAndIndex(global_index_);
      current_grid_type_ = grid_type;
      local_index_ = local_index;

      // Get the appropriate grid info
      const auto* grid_info = getCurrentGridInfo();
      if (grid_info) {
        const size_t nx = grid_info->nx;
        const size_t ny = grid_info->ny;

        // Calculate 3D indices from local index
        k_ = local_index_ / (nx * ny);
        const size_t remainder = local_index_ % (nx * ny);
        j_ = remainder / nx;
        i_ = remainder % nx;
      } else {
        // Invalid grid info
        i_ = j_ = k_ = 0;
      }
    } else {
      // End iterator or invalid
      current_grid_type_ = 0;
      local_index_ = 0;
      i_ = j_ = k_ = 0;
    }
  }
  ///@}

  ///@{ @name Iterator Operations
  /**
   * @brief Dereference operator to access current grid point
   *
   * @details Returns a Location object containing the geographic coordinates
   * of the current grid point from the appropriate grid type (unstaggered,
   * U-staggered, V-staggered, or W-staggered).
   *
   * @return value_type Location object at the current grid position
   * @throws std::out_of_range If iterator is at end position
   *
   * @note The returned Location object is created on-demand from grid
   * coordinates
   */
  value_type operator*() const { return geometry_->getLocation(global_index_); }

  /**
   * @brief Pre-increment operator to advance to next grid point
   *
   * @details Advances the iterator to the next grid point across all grids.
   * Updates the global index and recalculates the current grid type and
   * local position. When reaching the end of one grid, automatically
   * transitions to the next grid type.
   *
   * @return WRFGeometryIterator& Reference to this iterator after advancement
   *
   * @note Incrementing an end iterator results in undefined behavior
   */
  WRFGeometryIterator& operator++() {
    if (geometry_ && global_index_ < geometry_->totalGridSize()) {
      // Increment global index
      ++global_index_;

      if (global_index_ < geometry_->totalGridSize()) {
        // Calculate new grid type and local position
        auto [grid_type, local_index] =
            geometry_->getGridTypeAndIndex(global_index_);
        current_grid_type_ = grid_type;
        local_index_ = local_index;

        // Get the appropriate grid info
        const auto* grid_info = getCurrentGridInfo();
        if (grid_info) {
          const size_t nx = grid_info->nx;
          const size_t ny = grid_info->ny;

          // Calculate 3D indices from local index
          k_ = local_index_ / (nx * ny);
          const size_t remainder = local_index_ % (nx * ny);
          j_ = remainder / nx;
          i_ = remainder % nx;
        } else {
          // Invalid grid info
          i_ = j_ = k_ = 0;
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
   * equal if they reference the same geometry and are at the same global
   * position.
   *
   * @param[in] other Iterator to compare with
   * @return bool True if iterators are equal, false otherwise
   */
  bool operator==(const WRFGeometryIterator& other) const {
    return (geometry_ == other.geometry_ &&
            global_index_ == other.global_index_);
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

  ///@{ @name Grid Position Access
  /**
   * @brief Get current X-dimension index within current grid
   *
   * @details Returns the current position in the X dimension (west-east
   * direction) within the current grid type. This corresponds to the
   * longitude dimension in geographic coordinates.
   *
   * @return size_t Current X-dimension index (0-based)
   */
  size_t i() const { return i_; }

  /**
   * @brief Get current Y-dimension index within current grid
   *
   * @details Returns the current position in the Y dimension (south-north
   * direction) within the current grid type. This corresponds to the
   * latitude dimension in geographic coordinates.
   *
   * @return size_t Current Y-dimension index (0-based)
   */
  size_t j() const { return j_; }

  /**
   * @brief Get current Z-dimension index within current grid
   *
   * @details Returns the current position in the Z dimension (bottom-top
   * direction) within the current grid type. This corresponds to the
   * vertical level in meteorological coordinates.
   *
   * @return size_t Current Z-dimension index (0-based)
   */
  size_t k() const { return k_; }

  /**
   * @brief Get current grid type
   *
   * @details Returns the grid type that the iterator is currently positioned
   * in. Grid types: 0=unstaggered, 1=U-staggered, 2=V-staggered, 3=W-staggered
   *
   * @return int Current grid type (0-3)
   */
  int currentGridType() const { return current_grid_type_; }

  /**
   * @brief Get current global index across all grids
   *
   * @details Returns the global linear index across all grids.
   *
   * @return size_t Current global index
   */
  size_t globalIndex() const { return global_index_; }

  /**
   * @brief Get current local index within current grid
   *
   * @details Returns the local linear index within the current grid type.
   *
   * @return size_t Current local index
   */
  size_t localIndex() const { return local_index_; }
  ///@}

 private:
  ///@{ @name Private Helper Methods
  /**
   * @brief Get grid info for current grid type
   *
   * @details Returns a pointer to the GridDimensionInfo structure for
   * the current grid type. Used internally to access grid dimensions
   * and coordinate information.
   *
   * @return const GridDimensionInfo* Pointer to current grid info
   */
  const GridDimensionInfo* getCurrentGridInfo() const {
    if (!geometry_) return nullptr;

    switch (current_grid_type_) {
      case 0:
        return &geometry_->unstaggered_info();
      case 1:
        return &geometry_->u_staggered_info();
      case 2:
        return &geometry_->v_staggered_info();
      case 3:
        return &geometry_->w_staggered_info();
      default:
        return nullptr;
    }
  }
  ///@}

  ///@{ @name Member Variables
  const WRFGeometry<ConfigBackend>*
      geometry_;           ///< Pointer to parent geometry backend
  size_t global_index_;    ///< Global linear index across all grids
  int current_grid_type_;  ///< Current grid type (0-3)
  size_t local_index_;     ///< Local index within current grid
  size_t i_;               ///< Current X-dimension index within current grid
  size_t j_;               ///< Current Y-dimension index within current grid
  size_t k_;               ///< Current Z-dimension index within current grid
  ///@}
};

}  // namespace metada::backends::wrf