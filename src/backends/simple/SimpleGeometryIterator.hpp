#pragma once

#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#include "Location.hpp"

namespace metada::backends::simple {

// Iter should be an iterator over std::pair<int, int>
template <typename Iter>
class SimpleGeometryIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = metada::framework::Location;
  using difference_type = typename std::iterator_traits<Iter>::difference_type;
  using pointer = value_type*;   // Not used
  using reference = value_type;  // Return by value

  SimpleGeometryIterator() = default;
  explicit SimpleGeometryIterator(Iter it) : it_(it) {}

  // Copy/move constructors/assignments
  SimpleGeometryIterator(const SimpleGeometryIterator&) = default;
  SimpleGeometryIterator& operator=(const SimpleGeometryIterator&) = default;
  SimpleGeometryIterator(SimpleGeometryIterator&&) noexcept = default;
  SimpleGeometryIterator& operator=(SimpleGeometryIterator&&) noexcept =
      default;

  // Dereference: return Location object
  value_type operator*() const {
    auto [x, y] = *it_;
    return framework::Location(x, y);
  }
  // No operator-> needed (Location is returned by value)

  // Pre-increment
  SimpleGeometryIterator& operator++() {
    ++it_;
    return *this;
  }
  // Post-increment
  SimpleGeometryIterator operator++(int) {
    SimpleGeometryIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  // Equality/inequality
  bool operator==(const SimpleGeometryIterator& other) const {
    return it_ == other.it_;
  }
  bool operator!=(const SimpleGeometryIterator& other) const {
    return it_ != other.it_;
  }

 private:
  Iter it_;
};

}  // namespace metada::backends::simple