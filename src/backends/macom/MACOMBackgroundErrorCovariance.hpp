/**
 * @file MACOMBackgroundErrorCovariance.hpp
 * @brief MACOM background error covariance backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This class provides a MACOM-specific implementation of background error
 * covariance operations. It supports diagonal covariance matrices with optional
 * localization and can be extended for ensemble-based covariance in the future.
 */

#pragma once

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "MACOMState.hpp"

namespace metada::backends::macom {

/**
 * @brief MACOM backend implementation for background error covariance matrix
 *
 * @details This class provides a MACOM-specific implementation of background
 * error covariance operations. It supports diagonal covariance matrices with
 * optional localization and can be extended for ensemble-based covariance in
 * the future.
 */
template <typename ConfigBackend>
class MACOMBackgroundErrorCovariance {
 public:
  // =============================================================================
  // FRAMEWORK CONCEPTS REQUIRED INTERFACES
  // Required by BackgroundErrorCovarianceBackendType concept
  // =============================================================================

  // --- Resource management (required by framework) ---

  /**
   * @brief Default constructor is deleted (required by framework)
   */
  MACOMBackgroundErrorCovariance() = delete;

  /**
   * @brief Copy constructor is deleted (required by framework)
   */
  MACOMBackgroundErrorCovariance(const MACOMBackgroundErrorCovariance&) =
      delete;

  /**
   * @brief Copy assignment operator is deleted (required by framework)
   */
  MACOMBackgroundErrorCovariance& operator=(
      const MACOMBackgroundErrorCovariance&) = delete;

  /**
   * @brief Constructor that initializes from configuration (required by
   * framework)
   * @param config Configuration backend containing B matrix parameters
   */
  explicit MACOMBackgroundErrorCovariance(const ConfigBackend& config)
      : size_(config.Get("size").asInt()),
        variance_(config.Get("variance").asFloat()),
        localization_radius_(config.Get("localization_radius").asFloat()),
        initialized_(true) {}

  /**
   * @brief Move constructor (required by framework)
   */
  MACOMBackgroundErrorCovariance(
      MACOMBackgroundErrorCovariance&& other) noexcept = default;

  /**
   * @brief Move assignment operator (required by framework)
   */
  MACOMBackgroundErrorCovariance& operator=(
      MACOMBackgroundErrorCovariance&& other) noexcept = default;

  /**
   * @brief Destructor (required by framework)
   */
  ~MACOMBackgroundErrorCovariance() = default;

  // --- Initialization and state interface (required by framework) ---

  /**
   * @brief Check if the background error covariance is initialized
   * @return true if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  // --- Representation support queries (required by framework) ---

  /**
   * @brief Check if diagonal covariance is supported
   * @return true if diagonal covariance is supported
   */
  bool supportsDiagonal() const { return true; }

  /**
   * @brief Check if ensemble covariance is supported
   * @return true if ensemble covariance is supported
   */
  bool supportsEnsemble() const { return false; }

  /**
   * @brief Check if parametric covariance is supported
   * @return true if parametric covariance is supported
   */
  bool supportsParametric() const { return false; }

  /**
   * @brief Check if hybrid covariance is supported
   * @return true if hybrid covariance is supported
   */
  bool supportsHybrid() const { return false; }

  /**
   * @brief Check if full covariance is supported
   * @return true if full covariance is supported
   */
  bool supportsFull() const { return false; }

  // --- Localization support (required by framework) ---

  /**
   * @brief Check if localization is supported
   * @return true if localization is supported
   */
  bool supportsLocalization() const { return true; }

  /**
   * @brief Set localization radius
   * @param radius Localization radius
   */
  void setLocalizationRadius(double radius) { localization_radius_ = radius; }

  /**
   * @brief Get localization radius
   * @return Localization radius
   */
  double getLocalizationRadius() const { return localization_radius_; }

  // --- Matrix size information (required by framework) ---

  /**
   * @brief Get matrix size
   * @return Size of the covariance matrix
   */
  size_t getMatrixSize() const { return size_; }

  /**
   * @brief Get matrix rank
   * @return Rank of the covariance matrix
   */
  size_t getRank() const { return size_; }

  // --- Core covariance operations (required by framework) ---

  /**
   * @brief Apply B^{-1} to a vector (for diagonal case: divide by variance)
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

  // --- Framework-specific covariance operations (required by framework) ---

  /**
   * @brief Compute quadratic form for diagonal covariance
   * @param increment State increment
   * @return Quadratic form value
   */
  template <typename StateType>
  double computeQuadraticFormDiagonal(const StateType& increment) const {
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
   * @param increment State increment
   * @return Quadratic form value
   */
  template <typename StateType>
  double computeQuadraticFormEnsemble(const StateType& increment) const {
    // Fallback to diagonal for MACOM backend
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for parametric covariance (fallback to
   * diagonal)
   * @param increment State increment
   * @return Quadratic form value
   */
  template <typename StateType>
  double computeQuadraticFormParametric(const StateType& increment) const {
    // Fallback to diagonal for MACOM backend
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for hybrid covariance (fallback to diagonal)
   * @param increment State increment
   * @return Quadratic form value
   */
  template <typename StateType>
  double computeQuadraticFormHybrid(const StateType& increment) const {
    // Fallback to diagonal for MACOM backend
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for full covariance (fallback to diagonal)
   * @param increment State increment
   * @return Quadratic form value
   */
  template <typename StateType>
  double computeQuadraticFormFull(const StateType& increment) const {
    // Fallback to diagonal for MACOM backend
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Apply inverse diagonal covariance
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applyInverseDiagonal(const StateType& increment,
                            StateType& result) const {
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
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applyInverseEnsemble(const StateType& increment,
                            StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse parametric covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applyInverseParametric(const StateType& increment,
                              StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse hybrid covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applyInverseHybrid(const StateType& increment, StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse full covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applyInverseFull(const StateType& increment, StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply diagonal covariance
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applyDiagonal(const StateType& increment, StateType& result) const {
    if (increment.size() != size_) {
      throw std::invalid_argument("Increment size mismatch");
    }

    // Clone the increment to get the same structure
    result = std::move(*(increment.clone()));

    // For diagonal covariance: B = diag(sigma_i^2)
    for (size_t i = 0; i < size_; ++i) {
      result[i] *= variance_;  // Use operator[] instead of data()
    }
  }

  /**
   * @brief Apply ensemble covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applyEnsemble(const StateType& increment, StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applyDiagonal(increment, result);
  }

  /**
   * @brief Apply parametric covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applyParametric(const StateType& increment, StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applyDiagonal(increment, result);
  }

  /**
   * @brief Apply hybrid covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applyHybrid(const StateType& increment, StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applyDiagonal(increment, result);
  }

  /**
   * @brief Apply full covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applyFull(const StateType& increment, StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applyDiagonal(increment, result);
  }

  /**
   * @brief Apply square root diagonal covariance
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applySquareRootDiagonal(const StateType& increment,
                               StateType& result) const {
    if (increment.size() != size_) {
      throw std::invalid_argument("Increment size mismatch");
    }

    // Clone the increment to get the same structure
    result = std::move(*(increment.clone()));

    // For diagonal covariance: B^{1/2} = diag(sigma_i)
    double sqrt_variance = std::sqrt(variance_);
    for (size_t i = 0; i < size_; ++i) {
      result[i] *= sqrt_variance;  // Use operator[] instead of data()
    }
  }

  /**
   * @brief Apply square root ensemble covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applySquareRootEnsemble(const StateType& increment,
                               StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applySquareRootDiagonal(increment, result);
  }

  /**
   * @brief Apply square root parametric covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applySquareRootParametric(const StateType& increment,
                                 StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applySquareRootDiagonal(increment, result);
  }

  /**
   * @brief Apply square root hybrid covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applySquareRootHybrid(const StateType& increment,
                             StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applySquareRootDiagonal(increment, result);
  }

  /**
   * @brief Apply square root full covariance (fallback to diagonal)
   * @param increment Input state increment
   * @param result Output state increment
   */
  template <typename StateType>
  void applySquareRootFull(const StateType& increment,
                           StateType& result) const {
    // Fallback to diagonal for MACOM backend
    applySquareRootDiagonal(increment, result);
  }

  // =============================================================================
  // MACOM SPECIFIC FUNCTIONALITY
  // These are MACOM-specific methods beyond framework requirements
  // =============================================================================

  /**
   * @brief Get the variance value for diagonal covariance
   * @return Variance value
   */
  double getVariance() const { return variance_; }

  /**
   * @brief Set the variance value for diagonal covariance
   * @param variance New variance value
   */
  void setVariance(double variance) { variance_ = variance; }

  /**
   * @brief Get the covariance matrix size
   * @return Size of the covariance matrix
   */
  size_t getSize() const { return size_; }

  /**
   * @brief Check if the covariance matrix is diagonal
   * @return true if diagonal, false otherwise
   */
  bool isDiagonal() const { return true; }

  /**
   * @brief Get the covariance matrix type
   * @return String describing the covariance type
   */
  std::string getCovarianceType() const { return "diagonal"; }

  /**
   * @brief Get the covariance matrix description
   * @return String describing the covariance matrix
   */
  std::string getDescription() const {
    return "MACOM diagonal background error covariance with variance " +
           std::to_string(variance_) + " and size " + std::to_string(size_);
  }

  /**
   * @brief Validate covariance parameters
   * @return true if parameters are valid, false otherwise
   */
  bool validateParameters() const {
    return size_ > 0 && variance_ > 0.0 && initialized_;
  }

  /**
   * @brief Get covariance statistics
   * @return String with covariance statistics
   */
  std::string getStatistics() const {
    std::stringstream ss;
    ss << "Covariance Statistics:\n";
    ss << "  Type: " << getCovarianceType() << "\n";
    ss << "  Size: " << size_ << "\n";
    ss << "  Variance: " << variance_ << "\n";
    ss << "  Localization radius: " << localization_radius_ << "\n";
    ss << "  Initialized: " << (initialized_ ? "yes" : "no") << "\n";
    return ss.str();
  }

 private:
  size_t size_;      ///< Matrix size (state vector dimension)
  double variance_;  ///< Diagonal variance value
  double
      localization_radius_;  ///< Localization radius (not used in MACOM case)
  bool initialized_;         ///< Initialization status
};

}  // namespace metada::backends::macom
