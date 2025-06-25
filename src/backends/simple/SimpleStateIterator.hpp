#pragma once

#include <cstddef>
#include <iterator>
#include <utility>

#include "Location.hpp"
#include "SimpleGeometry.hpp"
#include "SimpleGeometryIterator.hpp"

namespace metada::backends::simple {

class SimpleState;  // Forward declaration

class SimpleStateIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = std::pair<framework::Location, double>;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  SimpleStateIterator() = default;
  SimpleStateIterator(
      SimpleState* state,
      SimpleGeometryIterator<std::vector<SimpleGeometry::Coord>::const_iterator>
          it)
      : state_(state),
        geom_it_(it),
        current_value_({framework::Location(0, 0), 0.0}) {}

  // Copy/move constructors/assignments
  SimpleStateIterator(const SimpleStateIterator&) = default;
  SimpleStateIterator& operator=(const SimpleStateIterator&) = default;
  SimpleStateIterator(SimpleStateIterator&&) noexcept = default;
  SimpleStateIterator& operator=(SimpleStateIterator&&) noexcept = default;

  reference operator*();

  SimpleStateIterator& operator++() {
    ++geom_it_;
    return *this;
  }

  SimpleStateIterator operator++(int) {
    SimpleStateIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator==(const SimpleStateIterator& other) const {
    return geom_it_ == other.geom_it_;
  }

  bool operator!=(const SimpleStateIterator& other) const {
    return geom_it_ != other.geom_it_;
  }

 private:
  SimpleState* state_ = nullptr;
  SimpleGeometryIterator<std::vector<SimpleGeometry::Coord>::const_iterator>
      geom_it_;
  value_type current_value_;
};

}  // namespace metada::backends::simple