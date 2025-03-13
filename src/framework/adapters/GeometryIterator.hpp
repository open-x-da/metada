/**
 * @file GeometryIterator.hpp
 * @brief Iterator implementation for traversing geometry grid points
 * @ingroup adapters
 *
 * @details
 * This header provides a concrete implementation of the IGeometryIterator
 * interface for traversing geometry grid points in N-dimensional space.
 */

#pragma once

// Include dependencies
#include "IGeometryIterator.hpp"

// Standard library includes
#include <iterator>
#include <memory>
#include <vector>

namespace metada::framework {

/**
 * @brief Concrete implementation of the IGeometryIterator for traversing
 * geometry grid points
 *
 * @tparam T Type of the coordinate value (typically double)
 */
template <typename T>
class GeometryIterator : public IGeometryIterator<T> {
 public:
  // Inherit iterator traits from the interface
  using typename IGeometryIterator<T>::value_type;
  using typename IGeometryIterator<T>::pointer;
  using typename IGeometryIterator<T>::reference;

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
        coordinates_(coordinates),
        valid_(true) {}

  /**
   * @brief Copy constructor
   */
  GeometryIterator(const GeometryIterator& other)
      : position_(other.position_),
        dimensions_(other.dimensions_),
        coordinates_(other.coordinates_),
        valid_(other.valid_) {}

  /**
   * @brief Dereference operator returns current point coordinates
   */
  const value_type& operator*() const override { return coordinates_; }

  /**
   * @brief Arrow operator for accessing point coordinates
   */
  const value_type* operator->() const override { return &coordinates_; }

  /**
   * @brief Pre-increment operator
   */
  IGeometryIterator<T>& operator++() override {
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
   * @brief Post-increment operator (non-virtual convenience method)
   */
  GeometryIterator operator++(int) {
    GeometryIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  /**
   * @brief Equality comparison
   */
  bool operator==(const IGeometryIterator<T>& other) const override {
    // Try to cast to our concrete type
    const GeometryIterator<T>* otherIterator =
        dynamic_cast<const GeometryIterator<T>*>(&other);

    if (!otherIterator) {
      return false;  // Not the same type, can't be equal
    }

    if (!valid_ && !otherIterator->valid_) return true;
    if (valid_ != otherIterator->valid_) return false;
    return position_ == otherIterator->position_;
  }

  /**
   * @brief Inequality comparison
   */
  bool operator!=(const IGeometryIterator<T>& other) const override {
    return !(*this == other);
  }

  /**
   * @brief Clone this iterator
   * @return Unique pointer to a new iterator instance
   */
  std::unique_ptr<IGeometryIterator<T>> clone() const override {
    return std::make_unique<GeometryIterator<T>>(*this);
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
  std::vector<T> coordinates_;      ///< Current coordinates
  bool valid_ = false;              ///< Iterator validity flag
};

}  // namespace metada::framework