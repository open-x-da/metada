#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace metada::backends::lite {

/**
 * @brief Lite state backend for concrete testing
 *
 * Implements a 3-dimensional state vector with basic operations
 */
class LiteState {
 public:
  LiteState() = default;

  // State data storage
  std::vector<double> state_data_;
  std::vector<std::string> variable_names_;

  // Constructor with config and geometry (for compatibility)
  template <typename ConfigBackend, typename GeometryBackend>
  LiteState(const ConfigBackend& config, const GeometryBackend& geometry) {
    // Initialize with 3 variables for testing
    variable_names_ = {"temperature", "pressure", "humidity"};
    state_data_.resize(3, 0.0);
  }

  // Core operations
  void zero() { std::fill(state_data_.begin(), state_data_.end(), 0.0); }

  double dot(const LiteState& other) const {
    double result = 0.0;
    for (size_t i = 0; i < state_data_.size(); ++i) {
      result += state_data_[i] * other.state_data_[i];
    }
    return result;
  }

  double norm() const {
    double result = 0.0;
    for (double val : state_data_) {
      result += val * val;
    }
    return std::sqrt(result);
  }

  // Arithmetic operations
  void add(const LiteState& other) {
    for (size_t i = 0; i < state_data_.size(); ++i) {
      state_data_[i] += other.state_data_[i];
    }
  }

  void subtract(const LiteState& other) {
    for (size_t i = 0; i < state_data_.size(); ++i) {
      state_data_[i] -= other.state_data_[i];
    }
  }

  void multiply(double scalar) {
    for (double& val : state_data_) {
      val *= scalar;
    }
  }

  // Comparison
  bool equals(const LiteState& other) const {
    return state_data_ == other.state_data_;
  }

  // Data access
  void* getData() { return state_data_.empty() ? nullptr : state_data_.data(); }

  const void* getData() const {
    return state_data_.empty() ? nullptr : state_data_.data();
  }

  // Variable information
  const std::vector<std::string>& getVariableNames() const {
    return variable_names_;
  }

  size_t size() const { return state_data_.size(); }

  // Clone
  std::unique_ptr<LiteState> clone() const {
    auto cloned = std::make_unique<LiteState>();
    cloned->state_data_ = state_data_;
    cloned->variable_names_ = variable_names_;
    return cloned;
  }

  // Test helper methods
  void setData(const std::vector<double>& data) { state_data_ = data; }

  void setVariables(const std::vector<std::string>& vars) {
    variable_names_ = vars;
  }

  // File I/O (stub implementations)
  void saveToFile([[maybe_unused]] const std::string& filename) const {
    // Stub implementation for testing
  }

  // Initialize
  void initialize() {
    // Stub implementation for testing
  }
};

}  // namespace metada::backends::lite