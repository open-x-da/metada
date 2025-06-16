#pragma once

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "SimpleGeometry.hpp"
#include "SimpleGeometryIterator.hpp"
#include "SimpleStateIterator.hpp"

namespace metada::backends::simple {

class SimpleState {
 public:
  // Prevent default construction and copy operations
  SimpleState() = delete;
  SimpleState(const SimpleState&) = delete;
  SimpleState& operator=(const SimpleState&) = delete;

  // Move constructor/assignment
  SimpleState(SimpleState&& other) noexcept
      : x_dim_(other.x_dim_),
        y_dim_(other.y_dim_),
        data_(std::move(other.data_)),
        dimensions_(std::move(other.dimensions_)),
        variable_names_(std::move(other.variable_names_)),
        geometry_(other.geometry_) {}

  SimpleState& operator=(SimpleState&& other) noexcept {
    if (this != &other) {
      data_ = std::move(other.data_);
      x_dim_ = other.x_dim_;
      y_dim_ = other.y_dim_;
      dimensions_ = std::move(other.dimensions_);
      variable_names_ = std::move(other.variable_names_);
      // Don't move geometry_ as it's a reference
    }
    return *this;
  }

  ~SimpleState() = default;

  // Construct from config and geometry
  template <typename ConfigBackend>
  SimpleState(const ConfigBackend& config, const SimpleGeometry& geometry)
      : geometry_(geometry) {
    // Get dimensions from geometry
    x_dim_ = static_cast<size_t>(geometry.x_dim());
    y_dim_ = static_cast<size_t>(geometry.y_dim());

    // Read data from file specified in config
    std::string filename = config.Get("file").asString();
    readFromFile(filename);
  }

  // Data access
  void* getData() { return data_.data(); }
  const void* getData() const { return data_.data(); }

  // Coordinate-based access
  double& at(const SimpleGeometry::Coord& coord) {
    return data_[coord.second * x_dim_ + coord.first];
  }
  const double& at(const SimpleGeometry::Coord& coord) const {
    return data_[coord.second * x_dim_ + coord.first];
  }

  // Iterators
  SimpleStateIterator begin() {
    return SimpleStateIterator(this, geometry_.begin());
  }
  SimpleStateIterator end() {
    return SimpleStateIterator(this, geometry_.end());
  }

  const SimpleStateIterator begin() const {
    return SimpleStateIterator(const_cast<SimpleState*>(this),
                               geometry_.begin());
  }
  const SimpleStateIterator end() const {
    return SimpleStateIterator(const_cast<SimpleState*>(this), geometry_.end());
  }

  // Variable information
  const std::vector<std::string>& getVariableNames() const {
    return variable_names_;
  }
  const std::vector<size_t>& getDimensions(const std::string& name) const {
    if (name != "state") {
      throw std::invalid_argument("SimpleState only supports 'state' variable");
    }
    return dimensions_;
  }

  // Vector operations
  void zero() { std::fill(data_.begin(), data_.end(), 0.0); }

  void add(const SimpleState& other) {
    if (data_.size() != other.data_.size()) {
      throw std::runtime_error("Cannot add states of different sizes");
    }
    for (size_t i = 0; i < data_.size(); ++i) {
      data_[i] += other.data_[i];
    }
  }

  void subtract(const SimpleState& other) {
    if (data_.size() != other.data_.size()) {
      throw std::runtime_error("Cannot subtract states of different sizes");
    }
    for (size_t i = 0; i < data_.size(); ++i) {
      data_[i] -= other.data_[i];
    }
  }

  void multiply(double scalar) {
    for (size_t i = 0; i < data_.size(); ++i) {
      data_[i] *= scalar;
    }
  }

  double dot(const SimpleState& other) const {
    if (data_.size() != other.data_.size()) {
      throw std::runtime_error(
          "Cannot compute dot product of states of different sizes");
    }
    double result = 0.0;
    for (size_t i = 0; i < data_.size(); ++i) {
      result += data_[i] * other.data_[i];
    }
    return result;
  }

  double norm() const { return std::sqrt(dot(*this)); }

  bool equals(const SimpleState& other) const {
    if (data_.size() != other.data_.size()) {
      return false;
    }
    for (size_t i = 0; i < data_.size(); ++i) {
      if (data_[i] != other.data_[i]) {
        return false;
      }
    }
    return true;
  }

  // Cloning
  std::unique_ptr<SimpleState> clone() const {
    return std::unique_ptr<SimpleState>(new SimpleState(*this, true));
  }

 private:
  // Private constructor for cloning
  SimpleState(const SimpleState& other, bool)
      : x_dim_(other.x_dim_),
        y_dim_(other.y_dim_),
        data_(other.data_),
        dimensions_(other.dimensions_),
        variable_names_(other.variable_names_),
        geometry_(other.geometry_) {}

  void readFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
      throw std::runtime_error("Cannot open state file: " + filename);
    }

    // Initialize dimensions
    dimensions_ = {y_dim_, x_dim_};
    variable_names_ = {"state"};

    // Read data
    data_.resize(x_dim_ * y_dim_);
    for (size_t y = 0; y < y_dim_; ++y) {
      for (size_t x = 0; x < x_dim_; ++x) {
        if (!(file >> data_[y * x_dim_ + x])) {
          throw std::runtime_error("Error reading state data from file");
        }
      }
    }
  }

  size_t x_dim_ = 0;
  size_t y_dim_ = 0;
  std::vector<double> data_;
  std::vector<size_t> dimensions_;
  std::vector<std::string> variable_names_;
  const SimpleGeometry& geometry_;
};

// Output operator
inline std::ostream& operator<<(std::ostream& os, const SimpleState& state) {
  for (const auto& [coord, value] : state) {
    os << "(" << coord.first << "," << coord.second << ")=" << value << " ";
  }
  return os;
}

}  // namespace metada::backends::simple