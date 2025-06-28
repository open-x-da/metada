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
template <typename ConfigBackend>
class WRFGeometry;

// Define standalone iterator classes
template <typename ConfigBackend>
class WRFGeometryIterator {
 public:
  // Iterator traits
  using iterator_category = std::forward_iterator_tag;
  using value_type = framework::Location;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  /**
   * @brief Default constructor - creates an invalid iterator
   */
  WRFGeometryIterator() = delete;

  /**
   * @brief Constructor for creating a valid iterator
   *
   * @param geometry Pointer to the WRF geometry backend
   * @param index Linear index into the grid (0 for begin, grid size for end)
   */
  WRFGeometryIterator(const WRFGeometry<ConfigBackend>* geometry, size_t index)
      : geometry_(geometry), index_(index) {
    if (geometry_ && index_ < geometry_->totalGridSize()) {
      // Calculate 3D indices from linear index
      const size_t nx = geometry_->nx_;
      const size_t ny = geometry_->ny_;

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

  /**
   * @brief Dereference operator
   *
   * @return Location object at the current position
   */
  value_type operator*() const {
    return framework::Location(static_cast<int>(i_), static_cast<int>(j_),
                               static_cast<int>(k_));
  }

  /**
   * @brief Pre-increment operator
   *
   * @return Reference to this iterator after incrementing
   */
  WRFGeometryIterator& operator++() {
    if (geometry_ && index_ < geometry_->totalGridSize()) {
      // Increment linear index
      ++index_;

      if (index_ < geometry_->totalGridSize()) {
        // Recalculate 3D indices
        const size_t nx = geometry_->nx_;
        const size_t ny = geometry_->ny_;

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
   * @brief Post-increment operator
   *
   * @return Copy of iterator before incrementing
   */
  WRFGeometryIterator operator++(int) {
    WRFGeometryIterator temp = *this;
    ++(*this);
    return temp;
  }

  /**
   * @brief Equality comparison
   *
   * @param other Iterator to compare with
   * @return True if iterators point to the same position
   */
  bool operator==(const WRFGeometryIterator& other) const {
    return (geometry_ == other.geometry_ && index_ == other.index_);
  }

  /**
   * @brief Inequality comparison
   *
   * @param other Iterator to compare with
   * @return True if iterators point to different positions
   */
  bool operator!=(const WRFGeometryIterator& other) const {
    return !(*this == other);
  }

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
  const WRFGeometry<ConfigBackend>* geometry_;  // Pointer to parent geometry
  size_t i_;                                    // West-east index
  size_t j_;                                    // South-north index
  size_t k_;                                    // Bottom-top index
  size_t index_;                                // Linear index into the grid
};

template <typename ConfigBackend>
class WRFGeometryConstIterator : public WRFGeometryIterator<ConfigBackend> {
 public:
  // Inherit constructors from base class
  using WRFGeometryIterator<ConfigBackend>::WRFGeometryIterator;
};

}  // namespace metada::backends::wrf