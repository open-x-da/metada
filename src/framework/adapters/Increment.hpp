#pragma once

#include "BackendTraits.hpp"
#include "GeometryConcepts.hpp"
#include "IncrementConcepts.hpp"

namespace metada::framework {

// Forward declare Geometry to avoid circular dependency
template <typename BackendTag>
  requires GeometryBackendType<BackendTag>
class Geometry;

template <typename BackendTag>
  requires IncrementBackendType<BackendTag>
class Increment {
 public:
  using IncrementBackendType =
      typename traits::BackendTraits<BackendTag>::IncrementBackend;
  using GeometryType = Geometry<BackendTag>;
  using GeometryBackendType =
      typename traits::BackendTraits<BackendTag>::GeometryBackend;

  // Constructors
  Increment() = delete;  // Must construct from geometry

  explicit Increment(const GeometryBackendType& geometry)
      : increment_(geometry) {}

  // Factory method
  static Increment createFromGeometry(const GeometryBackendType& geometry) {
    return Increment(geometry);
  }

  // Core operations
  void zero() { increment_.zero(); }
  void scale(double alpha) { increment_.scale(alpha); }
  void axpy(double alpha, const Increment& other) {
    increment_.axpy(alpha, other.increment_);
  }
  double dot(const Increment& other) const {
    return increment_.dot(other.increment_);
  }
  double norm() const { return increment_.norm(); }

  // Randomize (for testing)
  void randomize() { increment_.randomize(); }

  // Increment backend access
  IncrementBackendType& incrementBackend() { return increment_; }
  const IncrementBackendType& incrementBackend() const { return increment_; }

  // Geometry access
  const GeometryBackendType& geometry() const { return increment_.geometry(); }

  // Data access (for tests and gradient checks)
  template <typename T>
  T getData() const {
    // Delegate to backend implementation
    static_assert(std::is_same_v<T, std::vector<double>>,
                  "Increment::getData only supports std::vector<double>");
    return increment_.getData();
  }

  // Addition assignment operator
  Increment& operator+=(const Increment& other) {
    increment_ += other.increment_;
    return *this;
  }

  // Subtraction assignment operator
  Increment& operator-=(const Increment& other) {
    increment_ -= other.increment_;
    return *this;
  }

  // Scalar multiplication operator
  Increment& operator*=(double scalar) {
    increment_ *= scalar;
    return *this;
  }

  // Scalar division operator
  Increment& operator/=(double scalar) {
    increment_ /= scalar;
    return *this;
  }

 private:
  IncrementBackendType increment_;  ///< Increment backend (operates on grid%xa)
};

// Non-member scalar multiplication operators
template <typename BackendTag>
Increment<BackendTag> operator*(const Increment<BackendTag>& inc,
                                double scalar) {
  auto result = inc;
  result *= scalar;
  return result;
}

template <typename BackendTag>
Increment<BackendTag> operator*(double scalar,
                                const Increment<BackendTag>& inc) {
  return inc * scalar;
}

// Output operator
template <typename BackendTag>
  requires IncrementBackendType<BackendTag>
inline std::ostream& operator<<(std::ostream& os,
                                const Increment<BackendTag>& increment) {
  os << "Increment(norm=" << increment.norm() << ")";
  return os;
}

}  // namespace metada::framework
