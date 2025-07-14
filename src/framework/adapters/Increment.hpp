#pragma once

#include <memory>
#include <random>
#include <string>
#include <vector>

#include "BackendTraits.hpp"
#include "State.hpp"
#include "StateConcepts.hpp"

namespace metada::framework {

template <typename BackendTag>
  requires StateBackendType<BackendTag>
class Increment {
 public:
  using StateType = State<BackendTag>;

  // Constructors
  Increment() = default;
  explicit Increment(const StateType& state) : state_(state.clone()) {}
  Increment(const StateType& a, const StateType& b) : state_(a.clone()) {
    state_ -= b;
  }

  // Factory methods
  static Increment createFromEntity(const StateType& state) {
    return Increment(state);
  }
  static Increment createFromDifference(const StateType& a,
                                        const StateType& b) {
    return Increment(a, b);
  }

  // Core operations
  void zero() { state_.zero(); }
  void scale(double alpha) { state_ *= alpha; }
  void axpy(double alpha, const Increment& other) {
    state_ += other.state_ * alpha;
  }
  double dot(const Increment& other) const { return state_.dot(other.state_); }
  double norm() const { return state_.norm(); }

  // Randomize (for testing)
  void randomize() {
    auto* data = state_.template getDataPtr<double>();
    size_t n = state_.size();
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    for (size_t i = 0; i < n; ++i) data[i] = dist(gen);
  }

  // Data access (for tests)
  template <typename T>
  T& getData() {
    return state_.template getData<T>();
  }
  template <typename T>
  const T& getData() const {
    return state_.template getData<T>();
  }

  // Underlying state access (if needed)
  StateType& state() { return state_; }
  const StateType& state() const { return state_; }

  // Addition assignment operator
  Increment& operator+=(const Increment& other) {
    state_ += other.state_;
    return *this;
  }

  // Subtraction assignment operator
  Increment& operator-=(const Increment& other) {
    state_ -= other.state_;
    return *this;
  }

  // Scalar multiplication operator
  Increment& operator*=(double scalar) {
    state_ *= scalar;
    return *this;
  }

 private:
  StateType state_;
};

// Non-member scalar multiplication operators
template <typename BackendTag>
Increment<BackendTag> operator*(const Increment<BackendTag>& inc,
                                double scalar) {
  Increment<BackendTag> result(inc.state());
  result *= scalar;
  return result;
}

template <typename BackendTag>
Increment<BackendTag> operator*(double scalar,
                                const Increment<BackendTag>& inc) {
  return inc * scalar;
}

}  // namespace metada::framework
