#pragma once

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "SimpleGeometryIterator.hpp"

namespace metada::backends::simple {

class SimpleGeometry {
 public:
  using Coord = std::pair<int, int>;
  using container_type = std::vector<Coord>;
  using value_type = Coord;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using size_type = typename container_type::size_type;
  using difference_type = typename container_type::difference_type;
  using iterator = metada::backends::simple::SimpleGeometryIterator<
      typename container_type::iterator>;
  using const_iterator = metada::backends::simple::SimpleGeometryIterator<
      typename container_type::const_iterator>;

  // Prevent default construction and copy operations
  SimpleGeometry() = delete;
  SimpleGeometry(const SimpleGeometry&) = delete;
  SimpleGeometry& operator=(const SimpleGeometry&) = delete;

  // Move constructor/assignment
  SimpleGeometry(SimpleGeometry&&) noexcept = default;
  SimpleGeometry& operator=(SimpleGeometry&&) noexcept = default;

  ~SimpleGeometry() = default;

  // Construct from config (assume config provides "x_dim" and "y_dim")
  template <typename ConfigBackend>
  explicit SimpleGeometry(const ConfigBackend& config) {
    x_dim_ = config.Get("x_dim").asInt();
    y_dim_ = config.Get("y_dim").asInt();
    if (x_dim_ <= 0 || y_dim_ <= 0) {
      throw std::invalid_argument(
          "SimpleGeometry: x_dim and y_dim must be positive");
    }
    coordinates_.reserve(x_dim_ * y_dim_);
    for (int y = 0; y < y_dim_; ++y) {
      for (int x = 0; x < x_dim_; ++x) {
        coordinates_.emplace_back(x, y);
      }
    }
  }

  // Iterators
  iterator begin() { return iterator(coordinates_.begin()); }
  iterator end() { return iterator(coordinates_.end()); }
  const_iterator begin() const { return const_iterator(coordinates_.cbegin()); }
  const_iterator end() const { return const_iterator(coordinates_.cend()); }
  const_iterator cbegin() const {
    return const_iterator(coordinates_.cbegin());
  }
  const_iterator cend() const { return const_iterator(coordinates_.cend()); }

  // Size information
  size_type size() const { return coordinates_.size(); }
  bool empty() const { return coordinates_.empty(); }
  size_type max_size() const { return coordinates_.max_size(); }

  // Element access
  reference operator[](size_type idx) { return coordinates_[idx]; }
  const_reference operator[](size_type idx) const { return coordinates_[idx]; }
  reference at(size_type idx) { return coordinates_.at(idx); }
  const_reference at(size_type idx) const { return coordinates_.at(idx); }
  reference front() { return coordinates_.front(); }
  const_reference front() const { return coordinates_.front(); }
  reference back() { return coordinates_.back(); }
  const_reference back() const { return coordinates_.back(); }

  // Cloning
  SimpleGeometry clone() const {
    SimpleGeometry copy(x_dim_, y_dim_, coordinates_);
    return copy;
  }

  // Accessors
  int x_dim() const { return x_dim_; }
  int y_dim() const { return y_dim_; }
  const container_type& coordinates() const { return coordinates_; }

 private:
  int x_dim_ = 0;
  int y_dim_ = 0;
  container_type coordinates_;

  SimpleGeometry(int x_dim, int y_dim, const container_type& coords)
      : x_dim_(x_dim), y_dim_(y_dim), coordinates_(coords) {}
};

}  // namespace metada::backends::simple