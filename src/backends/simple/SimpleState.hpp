#pragma once

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Location.hpp"
#include "SimpleGeometry.hpp"

namespace metada::backends::simple {

class SimpleState {
 public:
  // Prevent default construction and copy operations
  SimpleState() = delete;
  SimpleState(const SimpleState&) = delete;
  SimpleState& operator=(const SimpleState&) = delete;

  // Move constructor/assignment
  SimpleState(SimpleState&& other) noexcept
      : data_(std::move(other.data_)),
        variable_names_(std::move(other.variable_names_)),
        geometry_(other.geometry_) {}

  SimpleState& operator=(SimpleState&& other) noexcept {
    if (this != &other) {
      data_ = std::move(other.data_);
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
    // Read variables from config
    try {
      auto variables_config = config.Get("variables");
      if (variables_config.isString()) {
        // Single variable as string
        variable_names_ = {variables_config.asString()};
      } else {
        // Multiple variables as vector
        variable_names_ = variables_config.asVectorString();
      }
    } catch (...) {
      // Fallback to default if variables not specified
      variable_names_ = {"state"};
    }

    // Read data from file specified in config
    std::string filename = config.Get("file").asString();
    readFromFile(filename);
  }

  // Data access - REQUIRED by StateBackendImpl
  void* getData() { return data_.data(); }
  const void* getData() const { return data_.data(); }

  // Coordinate-based access
  double& at(const SimpleGeometry::Coord& coord) {
    size_t x_dim = geometry_.x_dim();
    return data_[coord.second * x_dim + coord.first];
  }
  const double& at(const SimpleGeometry::Coord& coord) const {
    size_t x_dim = geometry_.x_dim();
    return data_[coord.second * x_dim + coord.first];
  }

  // Variable information - REQUIRED by StateBackendImpl
  const std::vector<std::string>& getVariableNames() const {
    return variable_names_;
  }

  // Size - REQUIRED by StateBackendImpl
  size_t size() const { return data_.size(); }

  // Geometry access - REQUIRED by IdentityObsOperator
  const SimpleGeometry& geometry() const { return geometry_; }

  // Vector operations - REQUIRED by StateBackendImpl
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

  // Increment operations
  template <typename IncrementBackend>
  void addIncrement(const IncrementBackend& increment) {
    std::vector<double> temp(data_.size());
    increment.extract(temp.data(), nullptr, nullptr, nullptr, nullptr);
    for (size_t i = 0; i < data_.size(); ++i) {
      data_[i] += temp[i];
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

  // Cloning - REQUIRED by StateBackendImpl
  std::unique_ptr<SimpleState> clone() const {
    return std::unique_ptr<SimpleState>(new SimpleState(*this, true));
  }

  // File I/O - REQUIRED by StateBackendImpl
  void saveToFile(const std::string& filename) const {
    // Create directory if it doesn't exist
    std::filesystem::path file_path(filename);
    std::filesystem::path dir_path = file_path.parent_path();
    if (!dir_path.empty() && !std::filesystem::exists(dir_path)) {
      std::filesystem::create_directories(dir_path);
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open file for writing: " + filename);
    }

    // Write data in grid format (y rows, x columns)
    size_t x_dim = geometry_.x_dim();
    size_t y_dim = geometry_.y_dim();
    for (size_t y = 0; y < y_dim; ++y) {
      for (size_t x = 0; x < x_dim; ++x) {
        file << std::fixed << std::setprecision(6) << std::setw(12)
             << data_[y * x_dim + x];
        if (x < x_dim - 1) {
          file << " ";
        }
      }
      file << "\n";
    }
  }

  // Location-based access - REQUIRED by framework::State adapter
  double& at(const framework::Location& loc) {
    auto grid = loc.getGridCoords2D();
    return at(grid);
  }
  const double& at(const framework::Location& loc) const {
    auto grid = loc.getGridCoords2D();
    return at(grid);
  }

  // Index-based access - REQUIRED by BackgroundErrorCovariance
  double& operator[](size_t index) {
    if (index >= data_.size()) {
      throw std::out_of_range("Index out of range");
    }
    return data_[index];
  }
  const double& operator[](size_t index) const {
    if (index >= data_.size()) {
      throw std::out_of_range("Index out of range");
    }
    return data_[index];
  }

  // Output operator
  friend std::ostream& operator<<(std::ostream& os, const SimpleState& state) {
    os << "SimpleState{";
    os << "variables: [";
    for (size_t i = 0; i < state.variable_names_.size(); ++i) {
      if (i > 0) os << ", ";
      os << "\"" << state.variable_names_[i] << "\"";
    }
    os << "]";
    os << ", size: " << state.data_.size();

    // Add statistics if data is available
    if (!state.data_.empty()) {
      auto min_it = std::min_element(state.data_.begin(), state.data_.end());
      auto max_it = std::max_element(state.data_.begin(), state.data_.end());
      double sum = std::accumulate(state.data_.begin(), state.data_.end(), 0.0);
      double mean = sum / state.data_.size();

      os << ", min: " << *min_it;
      os << ", max: " << *max_it;
      os << ", mean: " << mean;
    }

    os << "}";
    return os;
  }

 private:
  // Private constructor for cloning
  SimpleState(const SimpleState& other, bool)
      : data_(other.data_),
        variable_names_(other.variable_names_),
        geometry_(other.geometry_) {}

  void readFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
      throw std::runtime_error("Cannot open state file: " + filename);
    }

    // Read data line by line without assuming dimensions
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      double value;
      while (iss >> value) {
        data_.push_back(value);
      }
    }
  }

  std::vector<double> data_;
  std::vector<std::string> variable_names_;
  const SimpleGeometry& geometry_;
};

}  // namespace metada::backends::simple