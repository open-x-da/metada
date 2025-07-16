#pragma once

#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace metada::backends::lite {

/**
 * @brief Lite observation backend for concrete testing
 *
 * Implements a linear observation operator: y = H*x where H is a simple matrix
 * For testing purposes, we use a 2x3 observation operator:
 * H = [1.0  0.5  0.0]
 *     [0.0  1.0  0.5]
 *
 * This class is designed to comply with the ObservationBackendImpl concept.
 */
class LiteObs {
 public:
  // Type aliases for concept compliance
  using value_type = double;
  using iterator_type = std::vector<double>::const_iterator;

  // Deleted default constructor (stricter compliance)
  LiteObs() = delete;

  // Deleted copy constructor and copy assignment operator (stricter compliance)
  LiteObs(const LiteObs&) = delete;
  LiteObs& operator=(const LiteObs&) = delete;

  LiteObs(LiteObs&&) noexcept = default;
  LiteObs& operator=(LiteObs&&) noexcept = default;

  // Explicit config constructor (for concept compliance and flexibility)
  template <typename ConfigBackend>
  explicit LiteObs(const ConfigBackend& /*config*/) {
    // For test, initialize with 2 obs and diagonal covariance
    observations_ = {0.0, 0.0};
    covariance_ = {1.0, 1.0};
  }

  // Iteration capabilities
  iterator_type begin() const { return observations_.begin(); }
  iterator_type end() const { return observations_.end(); }
  size_t size() const { return observations_.size(); }
  const value_type& operator[](size_t idx) const { return observations_[idx]; }

  // Data access
  void* getData() {
    return observations_.empty() ? nullptr : observations_.data();
  }
  const void* getData() const {
    return observations_.empty() ? nullptr : observations_.data();
  }
  template <typename T>
  T getData() const {
    return T(observations_.begin(), observations_.end());
  }

  // Variable names
  std::vector<std::string> getTypeNames() const {
    return {"temperature", "pressure"};
  }
  std::vector<std::string> getVariableNames(const std::string& typeName) const {
    return {typeName};
  }

  // Covariance matrix
  std::vector<double> getCovariance() const { return covariance_; }

  // Cloning
  std::unique_ptr<LiteObs> clone() const {
    auto cloned =
        std::make_unique<LiteObs>((void*)nullptr);  // Use dummy config
    cloned->observations_ = observations_;
    cloned->covariance_ = covariance_;
    return cloned;
  }

  // Vector arithmetic
  void add(const LiteObs& other) {
    for (size_t i = 0; i < observations_.size(); ++i) {
      observations_[i] += other.observations_[i];
    }
  }
  void subtract(const LiteObs& other) {
    for (size_t i = 0; i < observations_.size(); ++i) {
      observations_[i] -= other.observations_[i];
    }
  }
  void multiply(double scalar) {
    for (double& obs : observations_) {
      obs *= scalar;
    }
  }

  // Comparison
  bool equals(const LiteObs& other) const {
    return observations_ == other.observations_ &&
           covariance_ == other.covariance_;
  }

  // Lifecycle management
  void initialize() {
    // Stub: could reset or reinitialize data
  }
  void applyQC() {
    // Stub: could apply quality control
  }

  // Observation operator (for test)
  static constexpr double H[2][3] = {{1.0, 0.5, 0.0}, {0.0, 1.0, 0.5}};
  std::vector<double> apply(const std::vector<double>& state) const {
    std::vector<double> result(2);
    for (int i = 0; i < 2; ++i) {
      result[i] = 0.0;
      for (int j = 0; j < 3; ++j) {
        result[i] += H[i][j] * state[j];
      }
    }
    return result;
  }
  std::vector<double> applyTangentLinear(
      const std::vector<double>& state_increment,
      const std::vector<double>& /*reference_state*/,
      const std::vector<double>& /*reference_obs*/) const {
    return apply(state_increment);
  }
  void applyAdjoint(const std::vector<double>& obs_increment,
                    const std::vector<double>& /*reference_state*/,
                    std::vector<double>& state_increment,
                    const std::vector<double>& /*reference_obs*/) const {
    state_increment.resize(3);
    for (int i = 0; i < 3; ++i) {
      state_increment[i] = 0.0;
      for (int j = 0; j < 2; ++j) {
        state_increment[i] += H[j][i] * obs_increment[j];
      }
    }
  }

  // Quadratic form: dy^T * R^(-1) * dy
  double quadraticForm(const std::vector<double>& innovation) const {
    double result = 0.0;
    for (size_t i = 0; i < innovation.size(); ++i) {
      result += innovation[i] * innovation[i] / covariance_[i];
    }
    return result;
  }
  // Apply inverse covariance: R^(-1) * dy
  std::vector<double> applyInverseCovariance(
      const std::vector<double>& innovation) const {
    std::vector<double> result(innovation.size());
    for (size_t i = 0; i < innovation.size(); ++i) {
      result[i] = innovation[i] / covariance_[i];
    }
    return result;
  }
  // Get inverse covariance diagonal
  std::vector<double> getInverseCovarianceDiagonal() const {
    std::vector<double> result(covariance_.size());
    for (size_t i = 0; i < covariance_.size(); ++i) {
      result[i] = 1.0 / covariance_[i];
    }
    return result;
  }
  // Check if covariance is diagonal
  bool isDiagonalCovariance() const { return true; }

  // Test helper methods
  void setObservations(const std::vector<double>& obs) { observations_ = obs; }
  void setCovariance(const std::vector<double>& cov) { covariance_ = cov; }

 private:
  std::vector<double> observations_;
  std::vector<double> covariance_;
};

}  // namespace metada::backends::lite