#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "SimpleGeometry.hpp"

namespace metada::backends::simple {

/**
 * @brief Simple increment backend for testing
 *
 * @details Vector-based increment implementation for simple backend.
 * Independent of SimpleState, operates on its own storage.
 */
class SimpleIncrement {
 public:
  explicit SimpleIncrement(const SimpleGeometry& geometry)
      : data_(geometry.size(), 0.0), geometry_(&geometry) {}

  void zero() { std::fill(data_.begin(), data_.end(), 0.0); }

  void scale(double alpha) {
    for (auto& val : data_) val *= alpha;
  }

  void axpy(double alpha, const SimpleIncrement& other) {
    for (size_t i = 0; i < data_.size(); ++i) {
      data_[i] += alpha * other.data_[i];
    }
  }

  double dot(const SimpleIncrement& other) const {
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

  SimpleIncrement& operator+=(const SimpleIncrement& other) {
    axpy(1.0, other);
    return *this;
  }

  SimpleIncrement& operator-=(const SimpleIncrement& other) {
    axpy(-1.0, other);
    return *this;
  }

  SimpleIncrement& operator*=(double scalar) {
    scale(scalar);
    return *this;
  }

  SimpleIncrement& operator/=(double scalar) {
    scale(1.0 / scalar);
    return *this;
  }

  size_t size() const { return data_.size(); }

  const SimpleGeometry& geometry() const { return *geometry_; }

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

  // Stub methods for compatibility with WRFIncrement interface
  int getNx() const { return static_cast<int>(data_.size()); }
  int getNy() const { return 1; }
  int getNz() const { return 1; }

  void extract(double* u, double* v, double* t, double* q, double* psfc) const {
    std::copy(data_.begin(), data_.end(), u);
    (void)v;
    (void)t;
    (void)q;
    (void)psfc;
  }

  void update(const double* u, const double* v, const double* t,
              const double* q, const double* psfc) {
    std::copy(u, u + data_.size(), data_.begin());
    (void)v;
    (void)t;
    (void)q;
    (void)psfc;
  }

 private:
  std::vector<double> data_;
  const SimpleGeometry* geometry_;
};

}  // namespace metada::backends::simple
