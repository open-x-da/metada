#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "MockGeometry.hpp"

namespace metada::backends::gmock {

/**
 * @brief Mock increment backend for testing
 *
 * @details Vector-based increment implementation for mock backend.
 * Independent of MockState, operates on its own storage.
 * Can be extended with Google Mock expectations if needed.
 */
class MockIncrement {
 public:
  explicit MockIncrement(const MockGeometry& geometry)
      : data_(10, 0.0), geometry_(&geometry) {}  // Default size for mocking

  void zero() { std::fill(data_.begin(), data_.end(), 0.0); }

  void scale(double alpha) {
    for (auto& val : data_) val *= alpha;
  }

  void axpy(double alpha, const MockIncrement& other) {
    for (size_t i = 0; i < data_.size(); ++i) {
      data_[i] += alpha * other.data_[i];
    }
  }

  double dot(const MockIncrement& other) const {
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

  MockIncrement& operator+=(const MockIncrement& other) {
    axpy(1.0, other);
    return *this;
  }

  MockIncrement& operator-=(const MockIncrement& other) {
    axpy(-1.0, other);
    return *this;
  }

  MockIncrement& operator*=(double scalar) {
    scale(scalar);
    return *this;
  }

  MockIncrement& operator/=(double scalar) {
    scale(1.0 / scalar);
    return *this;
  }

  size_t size() const { return data_.size(); }

  const MockGeometry& geometry() const { return *geometry_; }

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
    if (state_data == nullptr) {
      // State has no data, zero out the increment
      zero();
      return;
    }

    // Get the actual size of the state data
    size_t state_size = state_backend.size();

    // Resize increment to match state size, or copy minimum if sizes differ
    if (data_.size() != state_size) {
      data_.resize(state_size);
    }

    // Copy only the available elements
    std::copy(state_data, state_data + state_size, data_.begin());
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
  const MockGeometry* geometry_;
};

}  // namespace metada::backends::gmock
