/**
 * @file WRFGeometryIterator.cpp
 * @brief Implementation of WRF geometry iterator
 * @ingroup backends
 * @author Metada Framework Team
 */

#include <tuple>

#include "WRFGeometry.hpp"

namespace metada::backends::wrf {

// Constructor for creating a valid iterator
WRFGeometry::iterator::iterator(const WRFGeometry* geometry, size_t index)
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

// Dereference operator
WRFGeometry::iterator::value_type WRFGeometry::iterator::operator*() const {
  return std::make_tuple(i_, j_, k_);
}

// Pre-increment operator
WRFGeometry::iterator& WRFGeometry::iterator::operator++() {
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

// Post-increment operator
WRFGeometry::iterator WRFGeometry::iterator::operator++(int) {
  iterator temp = *this;
  ++(*this);
  return temp;
}

// Equality comparison
bool WRFGeometry::iterator::operator==(const iterator& other) const {
  // Two iterators are equal if they have the same index and geometry
  return (geometry_ == other.geometry_ && index_ == other.index_);
}

// Inequality comparison
bool WRFGeometry::iterator::operator!=(const iterator& other) const {
  return !(*this == other);
}

}  // namespace metada::backends::wrf