/**
 * @file GeometryIterator.hpp
 * @brief Iterator for traversing geometry grid points
 * @ingroup repr
 *
 * @details
 * This header provides an implementation of a forward iterator for traversing
 * geometry grid points in N-dimensional space.
 */

#pragma once

#include <iterator>
#include <vector>

namespace metada::framework {

/**
 * @brief Forward iterator for traversing geometry grid points
 *
 * @tparam T Type of the coordinate value (typically double)
 */
template <typename T>
class GeometryIterator {
 public:
  // Iterator traits
  using iterator_category = std::forward_iterator_tag;
  using value_type = std::vector<T>;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  /**
   * @brief Default constructor creates an end iterator
   */
  GeometryIterator() = default;

  /**
   * @brief Construct iterator with specific position
   */
  GeometryIterator(const std::vector<size_t>& position,
                   const std::vector<size_t>& dimensions,
                   const std::vector<T>& coordinates)
      : position_(position),
        dimensions_(dimensions),
        value_(),
        valid_(true),
        coordinates_(coordinates) {}

  /**
   * @brief Dereference operator returns current point coordinates
   */
  const value_type& operator*() const { return coordinates_; }

  /**
   * @brief Arrow operator for accessing point coordinates
   */
  const value_type* operator->() const { return &coordinates_; }

  /**
   * @brief Pre-increment operator
   */
  GeometryIterator& operator++() {
    // Increment position along the dimensions
    for (size_t i = 0; i < position_.size(); ++i) {
      ++position_[i];
      if (position_[i] < dimensions_[i]) {
        // Update coordinates based on new position
        // This is a placeholder - actual implementation would be provided by
        // derived classes
        break;
      }
      position_[i] = 0;
      if (i == position_.size() - 1) {
        // We've reached the end
        valid_ = false;
      }
    }
    return *this;
  }

  /**
   * @brief Post-increment operator
   */
  GeometryIterator operator++(int) {
    GeometryIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  /**
   * @brief Equality comparison
   */
  bool operator==(const GeometryIterator& other) const {
    if (!valid_ && !other.valid_) return true;
    if (valid_ != other.valid_) return false;
    return position_ == other.position_;
  }

  /**
   * @brief Inequality comparison
   */
  bool operator!=(const GeometryIterator& other) const {
    return !(*this == other);
  }

  /**
   * @brief Check if the iterator has reached the end
   */
  bool isDone() const { return !valid_; }

  /**
   * @brief Get current position indices
   */
  const std::vector<size_t>& getPosition() const { return position_; }

 protected:
  std::vector<size_t> position_;    ///< Current position indices
  std::vector<size_t> dimensions_;  ///< Grid dimensions
  std::vector<T> value_;  ///< Physical coordinates at current position
  bool valid_ = false;    ///< Iterator validity flag

 private:
  std::vector<T> coordinates_;  ///< Copy of coordinates for this position
};

}  // namespace metada::framework