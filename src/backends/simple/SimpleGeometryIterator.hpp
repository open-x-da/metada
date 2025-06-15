#pragma once

#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

namespace metada::backends::simple {

template <typename Iter>
class SimpleGeometryIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = typename std::iterator_traits<Iter>::value_type;
  using difference_type = typename std::iterator_traits<Iter>::difference_type;
  using pointer = typename std::iterator_traits<Iter>::pointer;
  using reference = typename std::iterator_traits<Iter>::reference;

  SimpleGeometryIterator() = default;
  explicit SimpleGeometryIterator(Iter it) : it_(it) {}

  // Copy/move constructors/assignments
  SimpleGeometryIterator(const SimpleGeometryIterator&) = default;
  SimpleGeometryIterator& operator=(const SimpleGeometryIterator&) = default;
  SimpleGeometryIterator(SimpleGeometryIterator&&) noexcept = default;
  SimpleGeometryIterator& operator=(SimpleGeometryIterator&&) noexcept =
      default;

  // Dereference
  reference operator*() const { return *it_; }
  pointer operator->() const { return it_.operator->(); }

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