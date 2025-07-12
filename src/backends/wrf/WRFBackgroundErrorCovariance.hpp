#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace metada::backends::wrf {

/**
 * @brief WRF backend implementation for background error covariance matrix
 *
 * @details This class provides a WRF-specific implementation of background
 * error covariance operations. For now, it provides diagonal covariance
 * functionality that can be extended to support WRF-specific covariance
 * structures.
 *
 * @tparam ConfigBackend Configuration backend type
 */
template <typename ConfigBackend>
class WRFBackgroundErrorCovariance {
 public:
  /**
   * @brief Constructor from configuration
   *
   * @param config Configuration object containing covariance parameters
   */
  explicit WRFBackgroundErrorCovariance(const ConfigBackend& config)
      : initialized_(false) {
    // Parse WRF-specific covariance configuration with safe defaults
    try {
      // Read configuration values if they exist, otherwise use defaults
      size_ = config.HasKey("size")
                  ? static_cast<size_t>(config.Get("size").asInt())
                  : 1000;

      variance_ =
          config.HasKey("variance") ? config.Get("variance").asFloat() : 1.0;

      localization_radius_ = config.HasKey("localization_radius")
                                 ? config.Get("localization_radius").asFloat()
                                 : 0.0;

      initialized_ = true;
    } catch (const std::exception& e) {
      // If configuration parsing fails, use safe defaults
      size_ = 1000;
      variance_ = 1.0;
      localization_radius_ = 0.0;
      initialized_ = true;
    }
  }

  /**
   * @brief Constructor with explicit parameters
   *
   * @param size Size of the state vector (matrix dimension)
   * @param variance Default variance value for diagonal elements
   */
  explicit WRFBackgroundErrorCovariance(size_t size = 1000,
                                        double variance = 1.0)
      : size_(size),
        variance_(variance),
        localization_radius_(0.0),
        initialized_(true) {}

  // Initialization and state
  bool isInitialized() const { return initialized_; }

  // Representation support queries
  bool supportsDiagonal() const { return true; }
  bool supportsEnsemble() const { return false; }    // TODO: Implement for WRF
  bool supportsParametric() const { return false; }  // TODO: Implement for WRF
  bool supportsHybrid() const { return false; }      // TODO: Implement for WRF
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
  template <typename StateType>
  double computeQuadraticFormDiagonal(const StateType& increment) const {
    // For WRF, implement a basic diagonal quadratic form
    // This is a simplified implementation that treats WRF state as a vector
    size_t state_size = increment.size();
    const double* data = increment.template getDataPtr<double>();

    double result = 0.0;
    for (size_t i = 0; i < state_size; ++i) {
      result += (data[i] * data[i]) / variance_;
    }

    return result;
  }

  /**
   * @brief Compute quadratic form for ensemble covariance (fallback to
   * diagonal)
   */
  template <typename StateType>
  double computeQuadraticFormEnsemble(const StateType& increment) const {
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for parametric covariance (fallback to
   * diagonal)
   */
  template <typename StateType>
  double computeQuadraticFormParametric(const StateType& increment) const {
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for hybrid covariance (fallback to diagonal)
   */
  template <typename StateType>
  double computeQuadraticFormHybrid(const StateType& increment) const {
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for full covariance (fallback to diagonal)
   */
  template <typename StateType>
  double computeQuadraticFormFull(const StateType& increment) const {
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Apply inverse diagonal covariance
   */
  template <typename StateType>
  void applyInverseDiagonal(const StateType& increment,
                            StateType& result) const {
    // Clone the increment to get the same structure
    result = increment.clone();

    size_t state_size = increment.size();
    double* result_data = result.template getDataPtr<double>();

    // For diagonal covariance: B^-1 = diag(1/sigma_i^2)
    for (size_t i = 0; i < state_size; ++i) {
      result_data[i] /= variance_;
    }
  }

  /**
   * @brief Apply inverse ensemble covariance (fallback to diagonal)
   */
  template <typename StateType>
  void applyInverseEnsemble(const StateType& increment,
                            StateType& result) const {
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse parametric covariance (fallback to diagonal)
   */
  template <typename StateType>
  void applyInverseParametric(const StateType& increment,
                              StateType& result) const {
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse hybrid covariance (fallback to diagonal)
   */
  template <typename StateType>
  void applyInverseHybrid(const StateType& increment, StateType& result) const {
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse full covariance (fallback to diagonal)
   */
  template <typename StateType>
  void applyInverseFull(const StateType& increment, StateType& result) const {
    applyInverseDiagonal(increment, result);
  }

 private:
  size_t size_;                 ///< Matrix size (state vector dimension)
  double variance_;             ///< Diagonal variance value
  double localization_radius_;  ///< Localization radius
  bool initialized_;            ///< Initialization status
};

}  // namespace metada::backends::wrf