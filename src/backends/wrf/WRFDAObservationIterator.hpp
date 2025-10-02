#pragma once

#include <iterator>
#include <memory>

#include "ObsRecord.hpp"

namespace metada::backends::wrf {

/**
 * @brief Iterator for WRFDAObservation class
 *
 * This class provides a standard C++ iterator interface for traversing
 * observations in a WRFDAObservation container. It supports:
 * - Forward iteration
 * - Random access
 * - Comparison operations
 * - Standard iterator traits
 *
 * The iterator provides access to ObsRecord objects that represent
 * individual observations from the WRFDA iv_type/y_type structures.
 *
 * @tparam GeometryBackend The WRF geometry backend type
 */
template <typename GeometryBackend>
class WRFDAObservationIterator {
 public:
  // Standard iterator type definitions
  using iterator_category = std::random_access_iterator_tag;
  using value_type = framework::ObsRecord;
  using difference_type = std::ptrdiff_t;
  using pointer = const value_type*;
  using reference = const value_type&;

  // Forward declare WRFDAObservation to allow friend declaration
  template <typename T>
  friend class WRFDAObservation;

  /**
   * @brief Default constructor creates an invalid iterator
   */
  WRFDAObservationIterator() : obs_(nullptr), index_(0) {}

  /**
   * @brief Copy constructor
   * @param other Iterator to copy from
   */
  WRFDAObservationIterator(const WRFDAObservationIterator& other) = default;

  /**
   * @brief Assignment operator
   * @param other Iterator to assign from
   * @return Reference to this iterator
   */
  WRFDAObservationIterator& operator=(const WRFDAObservationIterator& other) =
      default;

  /**
   * @brief Equality comparison
   * @param other Iterator to compare with
   * @return True if iterators are equal
   */
  bool operator==(const WRFDAObservationIterator& other) const {
    return obs_ == other.obs_ && index_ == other.index_;
  }

  /**
   * @brief Inequality comparison
   * @param other Iterator to compare with
   * @return True if iterators are not equal
   */
  bool operator!=(const WRFDAObservationIterator& other) const {
    return !(*this == other);
  }

  /**
   * @brief Less than comparison
   * @param other Iterator to compare with
   * @return True if this iterator points to earlier element
   */
  bool operator<(const WRFDAObservationIterator& other) const {
    return obs_ == other.obs_ && index_ < other.index_;
  }

  /**
   * @brief Greater than comparison
   * @param other Iterator to compare with
   * @return True if this iterator points to later element
   */
  bool operator>(const WRFDAObservationIterator& other) const {
    return other < *this;
  }

  /**
   * @brief Less than or equal comparison
   * @param other Iterator to compare with
   * @return True if this iterator points to same or earlier element
   */
  bool operator<=(const WRFDAObservationIterator& other) const {
    return !(other < *this);
  }

  /**
   * @brief Greater than or equal comparison
   * @param other Iterator to compare with
   * @return True if this iterator points to same or later element
   */
  bool operator>=(const WRFDAObservationIterator& other) const {
    return !(*this < other);
  }

  /**
   * @brief Pre-increment operator
   * @return Reference to incremented iterator
   */
  WRFDAObservationIterator& operator++() {
    ++index_;
    return *this;
  }

  /**
   * @brief Post-increment operator
   * @return Copy of iterator before increment
   */
  WRFDAObservationIterator operator++(int) {
    WRFDAObservationIterator tmp(*this);
    ++index_;
    return tmp;
  }

  /**
   * @brief Pre-decrement operator
   * @return Reference to decremented iterator
   */
  WRFDAObservationIterator& operator--() {
    --index_;
    return *this;
  }

  /**
   * @brief Post-decrement operator
   * @return Copy of iterator before decrement
   */
  WRFDAObservationIterator operator--(int) {
    WRFDAObservationIterator tmp(*this);
    --index_;
    return tmp;
  }

  /**
   * @brief Addition assignment operator
   * @param n Number of positions to advance
   * @return Reference to advanced iterator
   */
  WRFDAObservationIterator& operator+=(difference_type n) {
    index_ += n;
    return *this;
  }

  /**
   * @brief Subtraction assignment operator
   * @param n Number of positions to retreat
   * @return Reference to retreated iterator
   */
  WRFDAObservationIterator& operator-=(difference_type n) {
    index_ -= n;
    return *this;
  }

  /**
   * @brief Addition operator
   * @param n Number of positions to advance
   * @return New iterator advanced by n
   */
  WRFDAObservationIterator operator+(difference_type n) const {
    WRFDAObservationIterator tmp(*this);
    return tmp += n;
  }

  /**
   * @brief Subtraction operator
   * @param n Number of positions to retreat
   * @return New iterator retreated by n
   */
  WRFDAObservationIterator operator-(difference_type n) const {
    WRFDAObservationIterator tmp(*this);
    return tmp -= n;
  }

  /**
   * @brief Iterator difference operator
   * @param other Iterator to subtract
   * @return Distance between iterators
   */
  difference_type operator-(const WRFDAObservationIterator& other) const {
    return index_ - other.index_;
  }

  /**
   * @brief Dereference operator
   * @return Reference to current observation record
   */
  reference operator*() const {
    if (!obs_) {
      throw std::runtime_error("Cannot dereference invalid iterator");
    }
    return getCurrentRecord();
  }

  /**
   * @brief Arrow operator
   * @return Pointer to current observation record
   */
  pointer operator->() const {
    if (!obs_) {
      throw std::runtime_error("Cannot dereference invalid iterator");
    }
    return &getCurrentRecord();
  }

  /**
   * @brief Array subscript operator
   * @param n Offset from current position
   * @return Reference to observation record at offset
   */
  reference operator[](difference_type n) const { return *(*this + n); }

 private:
  /**
   * @brief Constructor used by WRFDAObservation
   * @param obs Parent observation container
   * @param index Initial position
   */
  template <typename T>
  WRFDAObservationIterator(const WRFDAObservation<T>& obs, size_t index)
      : obs_(&obs), index_(index) {}

  /**
   * @brief Get current observation record
   * @return Reference to current record
   */
  reference getCurrentRecord() const {
    // TODO: Implement access to WRFDA observation at current index
    // This should:
    // 1. Access the appropriate iv_type/y_type data
    // 2. Convert to ObsRecord format
    // 3. Handle different observation types (SYNOP, SOUND, etc.)
    // 4. Include proper error handling
    throw std::runtime_error("Not implemented");
  }

  const void* obs_;  ///< Parent observation container (type-erased)
  size_t index_;     ///< Current position in observation sequence
};

// Non-member operator+ for symmetry
template <typename GeometryBackend>
WRFDAObservationIterator<GeometryBackend> operator+(
    typename WRFDAObservationIterator<GeometryBackend>::difference_type n,
    const WRFDAObservationIterator<GeometryBackend>& it) {
  return it + n;
}

}  // namespace metada::backends::wrf
