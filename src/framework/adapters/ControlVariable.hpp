#pragma once

#include "BackendTraits.hpp"
#include "ControlVariableConcepts.hpp"
#include "GeometryConcepts.hpp"

namespace metada::framework {

// Forward declare Geometry to avoid circular dependency
template <typename BackendTag>
  requires GeometryBackendType<BackendTag>
class Geometry;

/**
 * @brief Adapter class for control variables in incremental variational DA
 *
 * @details In incremental variational data assimilation, the control variable v
 * is the actual optimization variable that is minimized. The control variable
 * is transformed to state-space increments δx through a transformation
 * operator:
 *
 *   δx = U v
 *
 * where U is typically related to B^(1/2) (square root of background error
 * covariance). This formulation provides several advantages:
 * - The cost function in control space is: J(v) = 1/2 ||v||² + 1/2 ||d -
 * H(δx(v))||²_R
 * - The background term is simplified to a simple L2 norm
 * - Better numerical conditioning
 * - Natural preconditioning through the B-matrix transformation
 *
 * The ControlVariable adapter provides a uniform interface for optimization
 * algorithms (L-BFGS, CG, etc.) to work with control variables without
 * exposing backend-specific details.
 *
 * @tparam BackendTag The backend tag type that must satisfy the required
 * concepts
 */
template <typename BackendTag>
class ControlVariable {
 public:
  using ControlVariableBackendType =
      typename traits::BackendTraits<BackendTag>::ControlVariableBackend;
  using GeometryType = Geometry<BackendTag>;
  using GeometryBackendType =
      typename traits::BackendTraits<BackendTag>::GeometryBackend;

  // Constructors
  ControlVariable() = delete;  // Must construct from geometry

  /**
   * @brief Construct control variable from geometry
   *
   * @param geometry Geometry defining the domain
   */
  explicit ControlVariable(const GeometryBackendType& geometry)
      : control_variable_(geometry) {}

  /**
   * @brief Factory method to create control variable from geometry
   *
   * @param geometry Geometry defining the domain
   * @return Newly constructed ControlVariable
   */
  static ControlVariable createFromGeometry(
      const GeometryBackendType& geometry) {
    return ControlVariable(geometry);
  }

  // Core vector space operations (required by optimization algorithms)

  /**
   * @brief Set all elements to zero
   */
  void zero() { control_variable_.zero(); }

  /**
   * @brief Scale by a scalar: v = α * v
   *
   * @param alpha Scaling factor
   */
  void scale(double alpha) { control_variable_.scale(alpha); }

  /**
   * @brief AXPY operation: v = α * other + v
   *
   * @param alpha Scaling factor for other
   * @param other Control variable to add
   */
  void axpy(double alpha, const ControlVariable& other) {
    control_variable_.axpy(alpha, other.control_variable_);
  }

  /**
   * @brief Compute dot product with another control variable
   *
   * @param other Control variable to compute dot product with
   * @return Dot product value
   */
  double dot(const ControlVariable& other) const {
    return control_variable_.dot(other.control_variable_);
  }

  /**
   * @brief Compute L2 norm: ||v||
   *
   * @return L2 norm value
   */
  double norm() const { return control_variable_.norm(); }

  /**
   * @brief Randomize control variable values (for testing)
   */
  void randomize() { control_variable_.randomize(); }

  // Backend access

  /**
   * @brief Access to backend implementation
   *
   * @return Reference to backend control variable
   */
  ControlVariableBackendType& backend() { return control_variable_; }

  /**
   * @brief Const access to backend implementation
   *
   * @return Const reference to backend control variable
   */
  const ControlVariableBackendType& backend() const {
    return control_variable_;
  }

  /**
   * @brief Access to geometry
   *
   * @return Const reference to geometry
   */
  const GeometryBackendType& geometry() const {
    return control_variable_.geometry();
  }

  // Data access (for optimization algorithms and testing)

  /**
   * @brief Get control variable data as a vector
   *
   * @details This method is used by optimization algorithms to access the
   * underlying control variable data as a flat vector. The specific format
   * depends on the backend implementation.
   *
   * @tparam T Return type (typically std::vector<double>)
   * @return Control variable data
   */
  template <typename T>
  T getData() const {
    static_assert(std::is_same_v<T, std::vector<double>>,
                  "ControlVariable::getData only supports std::vector<double>");
    return control_variable_.getData();
  }

  /**
   * @brief Set control variable data from a vector
   *
   * @param data Vector containing control variable data
   */
  void setFromVector(const std::vector<double>& data) {
    control_variable_.setFromVector(data);
  }

  // Arithmetic operators (for convenient expressions)

  /**
   * @brief Addition assignment operator: v += other
   *
   * @param other Control variable to add
   * @return Reference to this control variable
   */
  ControlVariable& operator+=(const ControlVariable& other) {
    control_variable_ += other.control_variable_;
    return *this;
  }

  /**
   * @brief Subtraction assignment operator: v -= other
   *
   * @param other Control variable to subtract
   * @return Reference to this control variable
   */
  ControlVariable& operator-=(const ControlVariable& other) {
    control_variable_ -= other.control_variable_;
    return *this;
  }

  /**
   * @brief Scalar multiplication operator: v *= scalar
   *
   * @param scalar Scalar value to multiply by
   * @return Reference to this control variable
   */
  ControlVariable& operator*=(double scalar) {
    control_variable_ *= scalar;
    return *this;
  }

  /**
   * @brief Scalar division operator: v /= scalar
   *
   * @param scalar Scalar value to divide by
   * @return Reference to this control variable
   */
  ControlVariable& operator/=(double scalar) {
    control_variable_ /= scalar;
    return *this;
  }

 private:
  ControlVariableBackendType
      control_variable_;  ///< Control variable backend implementation
};

// Non-member arithmetic operators

/**
 * @brief Scalar multiplication: result = v * scalar
 *
 * @param cv Control variable
 * @param scalar Scalar value
 * @return Scaled control variable
 */
template <typename BackendTag>
ControlVariable<BackendTag> operator*(const ControlVariable<BackendTag>& cv,
                                      double scalar) {
  auto result = cv;
  result *= scalar;
  return result;
}

/**
 * @brief Scalar multiplication: result = scalar * v
 *
 * @param scalar Scalar value
 * @param cv Control variable
 * @return Scaled control variable
 */
template <typename BackendTag>
ControlVariable<BackendTag> operator*(double scalar,
                                      const ControlVariable<BackendTag>& cv) {
  return cv * scalar;
}

/**
 * @brief Output stream operator for control variables
 *
 * @param os Output stream
 * @param control_variable Control variable to output
 * @return Reference to output stream
 */
template <typename BackendTag>
inline std::ostream& operator<<(
    std::ostream& os, const ControlVariable<BackendTag>& control_variable) {
  os << "ControlVariable(norm=" << control_variable.norm() << ")";
  return os;
}

}  // namespace metada::framework
