#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "LiteGeometry.hpp"

namespace metada::backends::lite {

/**
 * @brief Lite increment backend for concrete testing
 *
 * @details Simple vector-based increment implementation for testing.
 * Independent of LiteState, operates on its own storage.
 */
class LiteIncrement {
 public:
  // Default constructor for convenience in tests/operators
  LiteIncrement() : data_(3, 0.0) {}
  explicit LiteIncrement(const LiteGeometry& geometry)
      : data_(3, 0.0) {  // 3 variables: temperature, pressure, humidity
    (void)geometry;      // Unused for now
  }

  void zero() { std::fill(data_.begin(), data_.end(), 0.0); }

  void scale(double alpha) {
    for (auto& val : data_) val *= alpha;
  }

  void axpy(double alpha, const LiteIncrement& other) {
    for (size_t i = 0; i < data_.size(); ++i) {
      data_[i] += alpha * other.data_[i];
    }
  }

  double dot(const LiteIncrement& other) const {
    double result = 0.0;
    for (size_t i = 0; i < data_.size(); ++i) {
      result += data_[i] * other.data_[i];
    }
    return result;
  }

  double norm() const {
    double sum_sq = 0.0;
    for (double val : data_) {
      sum_sq += val * val;
    }
    return std::sqrt(sum_sq);
  }

  LiteIncrement& operator+=(const LiteIncrement& other) {
    axpy(1.0, other);
    return *this;
  }

  LiteIncrement& operator-=(const LiteIncrement& other) {
    axpy(-1.0, other);
    return *this;
  }

  LiteIncrement& operator*=(double scalar) {
    scale(scalar);
    return *this;
  }

  LiteIncrement& operator/=(double scalar) {
    scale(1.0 / scalar);
    return *this;
  }

  size_t size() const { return data_.size(); }

  // Stub methods for compatibility with WRFIncrement interface
  int getNx() const { return static_cast<int>(data_.size()); }
  int getNy() const { return 1; }
  int getNz() const { return 1; }

  void extract(double* u, double* v, double* t, double* q, double* psfc) const {
    // For Lite: just copy data to first array
    std::copy(data_.begin(), data_.end(), u);
    (void)v;
    (void)t;
    (void)q;
    (void)psfc;
  }

  void update(const double* u, const double* v, const double* t,
              const double* q, const double* psfc) {
    // For Lite: just copy from first array
    std::copy(u, u + data_.size(), data_.begin());
    (void)v;
    (void)t;
    (void)q;
    (void)psfc;
  }

  const LiteGeometry& geometry() const {
    throw std::runtime_error("LiteIncrement::geometry() not implemented");
  }

  std::vector<double> getData() const { return data_; }

  void randomize() {
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    for (auto& val : data_) {
      val = dist(gen);
    }
  }

  /**
   * @brief Transfer state difference to increment
   */
  template <typename StateBackend>
  void transferFromState(const StateBackend& state_backend) {
    const double* state_data =
        static_cast<const double*>(state_backend.getData());
    std::copy(state_data, state_data + data_.size(), data_.begin());
  }

 private:
  std::vector<double> data_;
};

}  // namespace metada::backends::lite
