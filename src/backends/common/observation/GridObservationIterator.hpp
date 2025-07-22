#pragma once

#include <iterator>
#include <vector>

#include "PointObservation.hpp"

namespace metada::backends::common::observation {

using framework::ObservationPoint;

/**
 * @brief Iterator for grid-based observations
 *
 * @details
 * This iterator provides forward iterator functionality for traversing
 * observation points stored in a vector. It satisfies the standard
 * iterator requirements and can be used with STL algorithms and
 * range-based for loops.
 */
class GridObservationIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = ObservationPoint;
  using difference_type = std::ptrdiff_t;
  using pointer = ObservationPoint*;
  using reference = ObservationPoint&;

  /** @brief Default constructor for end/sentinel iterators */
  GridObservationIterator() = default;

  /**
   * @brief Constructor from data pointer and index
   * @param data Pointer to vector of observation points
   * @param index Current index position
   */
  GridObservationIterator(const std::vector<ObservationPoint>* data,
                          size_t index)
      : data_(data), index_(index) {}

  /**
   * @brief Dereference operator
   * @return Reference to current observation point
   */
  reference operator*() { return const_cast<reference>((*data_)[index_]); }

  /**
   * @brief Arrow operator
   * @return Pointer to current observation point
   */
  pointer operator->() { return &const_cast<reference>((*data_)[index_]); }

  /**
   * @brief Pre-increment operator
   * @return Reference to this iterator after incrementing
   */
  GridObservationIterator& operator++() {
    ++index_;
    return *this;
  }

  /**
   * @brief Post-increment operator
   * @return Copy of iterator before incrementing
   */
  GridObservationIterator operator++(int) {
    GridObservationIterator tmp = *this;
    ++index_;
    return tmp;
  }

  /**
   * @brief Equality comparison operator
   * @param other Iterator to compare with
   * @return True if iterators are equal
   */
  bool operator==(const GridObservationIterator& other) const {
    return data_ == other.data_ && index_ == other.index_;
  }

  /**
   * @brief Inequality comparison operator
   * @param other Iterator to compare with
   * @return True if iterators are not equal
   */
  bool operator!=(const GridObservationIterator& other) const {
    return !(*this == other);
  }

 private:
  const std::vector<ObservationPoint>* data_{
      nullptr};      ///< Pointer to observation data
  size_t index_{0};  ///< Current index position
};

}  // namespace metada::backends::common::observation