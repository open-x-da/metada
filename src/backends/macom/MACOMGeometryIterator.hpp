/**
 * @file MACOMGeometryIterator.hpp
 * @brief MACOM geometry iterator implementation (simplified)
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <cstddef>

#include "include/MACOMlogging.hpp"
// #include <stdexcept>
// #include <tuple>

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
  using pointer = value_type*;    // Or `const value_type*`
  using reference = value_type&;  // Changed from value_type&

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
    MACOM_LOG_DEBUG(
        "MACOMGeometryIterator",
        "Created iterator for geometry with index " + std::to_string(index));
  }

  /**
   * @brief Dereference operator
   *
   * @return Tuple containing (i, j, k) indices of current position
   */
  value_type operator*() const { return std::make_tuple(i_, j_, k_); }

  /**
   * @brief Pre-increment operator
   *
   * @return Reference to this iterator after incrementing
   */
  MACOMGeometryIterator& operator++() {
    // TODO: Implement
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