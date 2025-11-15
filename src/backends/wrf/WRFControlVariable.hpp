#pragma once

#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

#include "WRFGeometry.hpp"

namespace metada::backends::wrf {

/**
 * @brief WRF control variable backend implementation
 *
 * @details This class represents control variables for WRF. In the identity
 * backend case, the control variable is identical to the state-space increment
 * (i.e., grid%xa directly). In the WRFDA CV5 case, this would represent the
 * control vector with balance relationships and transforms.
 *
 * For now, this implementation supports the identity case where the control
 * variable has the same structure as the increment, but it is designed to be
 * extended to support WRFDA CV5 control vectors in the future.
 *
 * @tparam ConfigBackend Configuration backend type
 * @tparam GeometryBackend Geometry backend type
 */
template <typename ConfigBackend, typename GeometryBackend>
class WRFControlVariable {
 public:
  using GeometryType = GeometryBackend;

  /**
   * @brief Constructor from geometry
   *
   * @param geometry WRF geometry defining the domain
   */
  explicit WRFControlVariable(const GeometryBackend& geometry)
      : geometry_(&geometry), data_(computeSize(geometry), 0.0) {}

  /**
   * @brief Set all elements to zero
   */
  void zero() { std::fill(data_.begin(), data_.end(), 0.0); }

  /**
   * @brief Scale by a scalar: v = α * v
   *
   * @param alpha Scaling factor
   */
  void scale(double alpha) {
    for (auto& val : data_) {
      val *= alpha;
    }
  }

  /**
   * @brief AXPY operation: v = α * other + v
   *
   * @param alpha Scaling factor for other
   * @param other Control variable to add
   */
  void axpy(double alpha, const WRFControlVariable& other) {
    if (data_.size() != other.data_.size()) {
      throw std::runtime_error("WRFControlVariable::axpy: size mismatch");
    }
    for (size_t i = 0; i < data_.size(); ++i) {
      data_[i] += alpha * other.data_[i];
    }
  }

  /**
   * @brief Compute dot product with another control variable
   *
   * @param other Control variable to compute dot product with
   * @return Dot product value
   */
  double dot(const WRFControlVariable& other) const {
    if (data_.size() != other.data_.size()) {
      throw std::runtime_error("WRFControlVariable::dot: size mismatch");
    }
    double result = 0.0;
    for (size_t i = 0; i < data_.size(); ++i) {
      result += data_[i] * other.data_[i];
    }
    return result;
  }

  /**
   * @brief Compute L2 norm
   *
   * @return L2 norm value
   */
  double norm() const { return std::sqrt(dot(*this)); }

  /**
   * @brief Randomize control variable values (for testing)
   */
  void randomize() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, 1.0);
    for (auto& val : data_) {
      val = dist(gen);
    }
  }

  /**
   * @brief Get control variable data as a vector
   *
   * @return Copy of control variable data
   */
  std::vector<double> getData() const { return data_; }

  /**
   * @brief Set control variable data from a vector
   *
   * @param new_data Vector containing control variable data
   */
  void setFromVector(const std::vector<double>& new_data) { data_ = new_data; }

  /**
   * @brief Resize control variable storage and zero the new contents.
   *
   * @param new_size Desired storage size.
   */
  void resize(size_t new_size) { data_.assign(new_size, 0.0); }

  /**
   * @brief Access to geometry
   *
   * @return Const reference to geometry
   */
  const GeometryBackend& geometry() const { return *geometry_; }

  // Arithmetic operators

  /**
   * @brief Addition assignment operator
   */
  WRFControlVariable& operator+=(const WRFControlVariable& other) {
    axpy(1.0, other);
    return *this;
  }

  /**
   * @brief Subtraction assignment operator
   */
  WRFControlVariable& operator-=(const WRFControlVariable& other) {
    axpy(-1.0, other);
    return *this;
  }

  /**
   * @brief Scalar multiplication operator
   */
  WRFControlVariable& operator*=(double scalar) {
    scale(scalar);
    return *this;
  }

  /**
   * @brief Scalar division operator
   */
  WRFControlVariable& operator/=(double scalar) {
    if (std::abs(scalar) < 1e-15) {
      throw std::runtime_error("WRFControlVariable: division by zero");
    }
    scale(1.0 / scalar);
    return *this;
  }

 private:
  /**
   * @brief Compute the size of the control variable vector
   *
   * @details For identity backend, this is the same as the increment size.
   * For WRFDA CV5, this would be the size of the control vector (may be
   * different from increment size).
   *
   * @param geometry WRF geometry
   * @return Size of control variable vector
   */
  static size_t computeSize(const GeometryBackend& geometry) {
    // For now, use same size as increment
    // TODO: For WRFDA CV5, compute appropriate control vector size
    const size_t nx = geometry.x_dim();
    const size_t ny = geometry.y_dim();
    const size_t nz = geometry.z_dim();
    size_t total_size = 0;

    // 3D fields: U, V, T, QVAPOR (4 fields)
    total_size += 4 * nx * ny * nz;

    // 2D field: PSFC (1 field)
    total_size += nx * ny;

    return total_size;
  }

  const GeometryBackend* geometry_;  ///< Pointer to geometry
  std::vector<double> data_;         ///< Control variable data
};

}  // namespace metada::backends::wrf
