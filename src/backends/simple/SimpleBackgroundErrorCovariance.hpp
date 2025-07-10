#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "SimpleState.hpp"

namespace metada::backends::simple {

/**
 * @brief Simple backend implementation for background error covariance matrix
 *
 * @details This class provides a basic implementation of background error
 * covariance operations suitable for testing and development. It supports
 * diagonal covariance matrices with optional localization.
 */
class SimpleBackgroundErrorCovariance {
 public:
  /**
   * @brief Constructor
   *
   * @param size Size of the state vector (matrix dimension)
   * @param variance Default variance value for diagonal elements
   */
  explicit SimpleBackgroundErrorCovariance(size_t size = 100,
                                           double variance = 1.0)
      : size_(size),
        variance_(variance),
        localization_radius_(0.0),
        initialized_(true) {}

  // Initialization and state
  bool isInitialized() const { return initialized_; }

  // Representation support queries
  bool supportsDiagonal() const { return true; }
  bool supportsEnsemble() const { return false; }
  bool supportsParametric() const { return false; }
  bool supportsHybrid() const { return false; }
  bool supportsFull() const { return false; }

  // Localization support
  bool supportsLocalization() const { return true; }
  void setLocalizationRadius(double radius) { localization_radius_ = radius; }
  double getLocalizationRadius() const { return localization_radius_; }

  // Matrix size information
  size_t getMatrixSize() const { return size_; }
  size_t getRank() const { return size_; }

  // Diagonal covariance operations
  double getVariance() const { return variance_; }
  void setVariance(double variance) { variance_ = variance; }

  /**
   * @brief Apply B^{-1} to a vector (for diagonal case: divide by variance)
   *
   * @param input Input vector
   * @return Result vector
   */
  std::vector<double> applyInverse(const std::vector<double>& input) const {
    if (input.size() != size_) {
      throw std::invalid_argument("Vector size mismatch");
    }

    std::vector<double> result(size_);
    for (size_t i = 0; i < size_; ++i) {
      result[i] = input[i] / variance_;
    }
    return result;
  }

  /**
   * @brief Apply B to a vector (for diagonal case: multiply by variance)
   *
   * @param input Input vector
   * @return Result vector
   */
  std::vector<double> apply(const std::vector<double>& input) const {
    if (input.size() != size_) {
      throw std::invalid_argument("Vector size mismatch");
    }

    std::vector<double> result(size_);
    for (size_t i = 0; i < size_; ++i) {
      result[i] = input[i] * variance_;
    }
    return result;
  }

  /**
   * @brief Compute quadratic form x^T B^{-1} x
   *
   * @param vector Input vector x
   * @return Quadratic form value
   */
  double quadraticForm(const std::vector<double>& vector) const {
    if (vector.size() != size_) {
      throw std::invalid_argument("Vector size mismatch");
    }

    double result = 0.0;
    for (size_t i = 0; i < size_; ++i) {
      result += vector[i] * vector[i] / variance_;
    }
    return result;
  }

  /**
   * @brief Apply square root B^{1/2} to a vector
   *
   * @param input Input vector
   * @return Result vector
   */
  std::vector<double> applySqrt(const std::vector<double>& input) const {
    if (input.size() != size_) {
      throw std::invalid_argument("Vector size mismatch");
    }

    std::vector<double> result(size_);
    double sqrt_variance = std::sqrt(variance_);
    for (size_t i = 0; i < size_; ++i) {
      result[i] = input[i] * sqrt_variance;
    }
    return result;
  }

  /**
   * @brief Compute quadratic form for diagonal covariance
   */
  double computeQuadraticFormDiagonal(const SimpleState& increment) const {
    if (increment.size() != size_) {
      throw std::invalid_argument("Increment size mismatch");
    }

    double result = 0.0;

    // For diagonal covariance: x^T B^-1 x = sum(x_i^2 / sigma_i^2)
    for (size_t i = 0; i < size_; ++i) {
      double value = increment[i];  // Use operator[] instead of data()
      result += (value * value) / variance_;
    }

    return result;
  }

  /**
   * @brief Compute quadratic form for ensemble covariance (fallback to
   * diagonal)
   */
  double computeQuadraticFormEnsemble(const SimpleState& increment) const {
    // Fallback to diagonal for simple backend
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for parametric covariance (fallback to
   * diagonal)
   */
  double computeQuadraticFormParametric(const SimpleState& increment) const {
    // Fallback to diagonal for simple backend
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for hybrid covariance (fallback to diagonal)
   */
  double computeQuadraticFormHybrid(const SimpleState& increment) const {
    // Fallback to diagonal for simple backend
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for full covariance (fallback to diagonal)
   */
  double computeQuadraticFormFull(const SimpleState& increment) const {
    // Fallback to diagonal for simple backend
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Apply inverse diagonal covariance
   */
  void applyInverseDiagonal(const SimpleState& increment,
                            SimpleState& result) const {
    if (increment.size() != size_) {
      throw std::invalid_argument("Increment size mismatch");
    }

    // Clone the increment to get the same structure
    result = std::move(*(increment.clone()));

    // For diagonal covariance: B^-1 = diag(1/sigma_i^2)
    for (size_t i = 0; i < size_; ++i) {
      result[i] /= variance_;  // Use operator[] instead of data()
    }
  }

  /**
   * @brief Apply inverse ensemble covariance (fallback to diagonal)
   */
  void applyInverseEnsemble(const SimpleState& increment,
                            SimpleState& result) const {
    // Fallback to diagonal for simple backend
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse parametric covariance (fallback to diagonal)
   */
  void applyInverseParametric(const SimpleState& increment,
                              SimpleState& result) const {
    // Fallback to diagonal for simple backend
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse hybrid covariance (fallback to diagonal)
   */
  void applyInverseHybrid(const SimpleState& increment,
                          SimpleState& result) const {
    // Fallback to diagonal for simple backend
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse full covariance (fallback to diagonal)
   */
  void applyInverseFull(const SimpleState& increment,
                        SimpleState& result) const {
    // Fallback to diagonal for simple backend
    applyInverseDiagonal(increment, result);
  }

 private:
  size_t size_;      ///< Matrix size (state vector dimension)
  double variance_;  ///< Diagonal variance value
  double
      localization_radius_;  ///< Localization radius (not used in simple case)
  bool initialized_;         ///< Initialization status
};

}  // namespace metada::backends::simple