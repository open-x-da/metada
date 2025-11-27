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
 *
 * This class is designed to comply with the StateBackendImpl concept.
 */
class LiteState {
 public:
  // Deleted default constructor (concept compliance)
  LiteState() = delete;

  // Deleted copy constructor and copy assignment operator (concept compliance)
  LiteState(const LiteState&) = delete;
  LiteState& operator=(const LiteState&) = delete;

  LiteState(LiteState&&) noexcept = default;
  LiteState& operator=(LiteState&&) noexcept = default;

  // Explicit templated config/geometry constructor (concept compliance)
  template <typename ConfigBackend, typename GeometryBackend>
  explicit LiteState([[maybe_unused]] const ConfigBackend& config,
                     [[maybe_unused]] const GeometryBackend& geometry) {
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
    auto cloned = std::make_unique<LiteState>((void*)nullptr, (void*)nullptr);
    cloned->state_data_ = state_data_;
    cloned->variable_names_ = variable_names_;
    return cloned;
  }

  // Test helper methods
  void setData(const std::vector<double>& data) { state_data_ = data; }

  void setVariables(const std::vector<std::string>& vars) {
    variable_names_ = vars;
  }

  // Increment operations
  template <typename IncrementBackend>
  void addIncrement(const IncrementBackend& increment) {
    // For LiteState, just add all increment data to state data
    auto inc_data = increment.getData();
    for (size_t i = 0; i < std::min(state_data_.size(), inc_data.size()); ++i) {
      state_data_[i] += inc_data[i];
    }
  }

  // File I/O (stub implementations)
  void saveToFile([[maybe_unused]] const std::string& filename) const {
    // Stub implementation for testing
  }

  // Initialize
  void initialize() {
    // Stub implementation for testing
  }

  // Output operator
  friend std::ostream& operator<<(std::ostream& os, const LiteState& state) {
    os << "LiteState[";
    for (size_t i = 0; i < state.state_data_.size(); ++i) {
      if (i > 0) os << ", ";
      os << state.variable_names_[i] << "=" << state.state_data_[i];
    }
    os << "]";
    return os;
  }

 private:
  std::vector<double> state_data_;
  std::vector<std::string> variable_names_;
};

}  // namespace metada::backends::lite