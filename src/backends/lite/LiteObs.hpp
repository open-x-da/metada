#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
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

  /**
   * @brief Stream insertion operator for LiteObs summary
   *
   * @details Outputs a summary of the lite observation data including:
   * - Total number of observations
   * - Observation types and variables
   * - Quality control statistics
   * - Covariance information
   *
   * @param os Output stream to write to
   * @param obs LiteObs object to summarize
   * @return Reference to the output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const LiteObs& obs) {
    os << "=== Lite Observation Summary ===\n";
    os << "Total observations: " << obs.observations_.size() << "\n\n";

    // Show observation types and variables
    const auto& type_names = obs.getTypeNames();
    os << "Observation Types (" << type_names.size() << "):\n";
    os << std::string(50, '-') << "\n";

    for (const auto& type_name : type_names) {
      const auto& var_names = obs.getVariableNames(type_name);
      os << "Type: " << std::setw(15) << std::left << type_name;
      os << " | Variables: " << std::setw(3) << std::right << var_names.size();
      os << " | Total obs: " << std::setw(6) << std::right
         << obs.observations_.size() << "\n";

      if (!var_names.empty()) {
        os << "  Variables: ";
        for (size_t i = 0; i < var_names.size(); ++i) {
          if (i > 0) os << ", ";
          os << var_names[i] << "(" << obs.observations_.size() << ")";
        }
        os << "\n";
      }
      os << "\n";
    }

    // Show covariance information
    os << "Covariance Information:\n";
    os << std::string(50, '-') << "\n";
    os << "Diagonal: " << (obs.isDiagonalCovariance() ? "Yes" : "No") << "\n";

    const auto& cov = obs.getCovariance();
    if (!cov.empty()) {
      double min_error = *std::min_element(cov.begin(), cov.end());
      double max_error = *std::max_element(cov.begin(), cov.end());
      os << "Error variance range: [" << std::scientific << std::setprecision(2)
         << min_error << ", " << max_error << "]\n";
    }

    return os;
  }

 private:
  std::vector<double> observations_;
  std::vector<double> covariance_;
};

}  // namespace metada::backends::lite