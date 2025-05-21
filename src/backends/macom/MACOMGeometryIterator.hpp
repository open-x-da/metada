/**
 * @file MACOMGeometryIterator.hpp
 * @brief MACOM geometry iterator implementation (simplified)
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <cstddef>
#include <stdexcept>
#include <tuple>

namespace metada::backends::macom {

// Forward declaration
template <typename ConfigBackend>
class MACOMGeometry;

// Define simplified iterator class
template <typename ConfigBackend>
class MACOMGeometryIterator {
 public:
  // Iterator traits
  using iterator_category = std::forward_iterator_tag;
  using value_type = std::tuple<size_t, size_t, size_t>;  // i, j, k indices
  using difference_type = std::ptrdiff_t;
  // Using value_type for reference as operator* will return a temporary tuple
  // by value. Pointer can be value_type* or const value_type* if it were to
  // point to internal state. Since operator* returns by value, a traditional
  // pointer-to-element isn't directly applicable unless we store the tuple
  // inside the iterator and return its address, which adds complexity. For a
  // read-only iterator concept based on returning a temporary, pointer is often
  // value_type*.
  using pointer = value_type*;   // Or `const value_type*`
  using reference = value_type;  // Changed from value_type&

  /**
   * @brief Default constructor is deleted
   */
  MACOMGeometryIterator() = delete;

  /**
   * @brief Constructor for creating a valid iterator
   *
   * @param geometry Pointer to the MACOM geometry backend
   * @param index Linear index into the grid (0 for begin, grid size for end)
   */
  MACOMGeometryIterator(const MACOMGeometry<ConfigBackend>* geometry,
                        size_t index)
      : geometry_(geometry), index_(index) {
    calculate_indices();
  }

  /**
   * @brief Dereference operator
   *
   * @return Tuple containing (i, j, k) indices of current position
   */
  value_type operator*() const {
    if (!geometry_ || index_ >= geometry_->totalGridSize()) {
      throw std::out_of_range(
          "Cannot dereference end or invalid MACOMGeometryIterator");
    }
    return std::make_tuple(i_, j_, k_);
  }

  /**
   * @brief Pre-increment operator
   *
   * @return Reference to this iterator after incrementing
   */
  MACOMGeometryIterator& operator++() {
    if (geometry_ && index_ < geometry_->totalGridSize()) {
      ++index_;
      calculate_indices();
    } else if (geometry_) {                 // If at or past end
      index_ = geometry_->totalGridSize();  // Lock to end state
      set_indices_to_end_state();
    }
    return *this;
  }

  /**
   * @brief Post-increment operator
   *
   * @return Copy of iterator before incrementing
   */
  MACOMGeometryIterator operator++(int) {
    MACOMGeometryIterator temp = *this;
    ++(*this);
    return temp;
  }

  /**
   * @brief Equality comparison
   *
   * @param other Iterator to compare with
   * @return True if iterators point to the same position
   */
  bool operator==(const MACOMGeometryIterator& other) const {
    return (geometry_ == other.geometry_ && index_ == other.index_);
  }

  /**
   * @brief Inequality comparison
   *
   * @param other Iterator to compare with
   * @return True if iterators point to different positions
   */
  bool operator!=(const MACOMGeometryIterator& other) const {
    return !(*this == other);
  }

  /**
   * @brief Get i index
   *
   * @return Current i index
   */
  size_t i() const { return i_; }

  /**
   * @brief Get j index
   *
   * @return Current j index
   */
  size_t j() const { return j_; }

  /**
   * @brief Get k index
   *
   * @return Current k index
   */
  size_t k() const { return k_; }

 protected:
  void calculate_indices() {
    if (geometry_ && index_ < geometry_->totalGridSize()) {
      // These need access to geometry_->nx_, ny_, nz_
      // Ensure MACOMGeometry declares MACOMGeometryIterator as a friend
      // or provides public accessors for nx_, ny_, nz_.
      const size_t nx = geometry_->nx_;
      const size_t ny = geometry_->ny_;
      // const size_t nz = geometry_->nz_; // nz might not be needed if index is
      // 2D slice based

      if (nx == 0 || ny == 0) {      // Protect against division by zero for
                                     // uninitialized geometry
        set_indices_to_end_state();  // Or handle as error
        return;
      }
      // Example for row-major order: index = k * (nx * ny) + j * nx + i
      // Assuming iteration is over a 3D grid. If it's a surface, k might always
      // be 0.
      k_ = index_ / (nx * ny);
      size_t remainder = index_ % (nx * ny);
      j_ = remainder / nx;
      i_ = remainder % nx;
    } else {
      set_indices_to_end_state();
    }
  }

  void set_indices_to_end_state() {
    // Define what an "end" iterator's indices look like.
    // This is mostly for consistency if someone inspects i_,j_,k_ of an end
    // iterator. Often, it doesn't matter as long as operator== works correctly.
    i_ = static_cast<size_t>(-1);
    j_ = static_cast<size_t>(-1);
    k_ = static_cast<size_t>(-1);
  }

 private:
  const MACOMGeometry<ConfigBackend>* geometry_;  // Pointer to parent geometry
  size_t i_;                                      // X index
  size_t j_;                                      // Y index
  size_t k_;                                      // Z index
  size_t index_;                                  // Linear index into the grid
};

template <typename ConfigBackend>
class MACOMGeometryConstIterator : public MACOMGeometryIterator<ConfigBackend> {
 public:
  // Inherit constructors from base class
  using MACOMGeometryIterator<ConfigBackend>::MACOMGeometryIterator;
};

}  // namespace metada::backends::macom