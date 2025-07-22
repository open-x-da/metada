#pragma once

#include <concepts>
#include <string>

#include "BackendTraits.hpp"
#include "CommonConcepts.hpp"

namespace metada::framework {

/**
 * @brief Concept defining the interface requirements for
 * BackgroundErrorCovariance backend implementations
 *
 * @details This concept ensures that backend types provide the necessary
 * operations for background error covariance matrix operations in variational
 * data assimilation. Backends must support various B matrix representations
 * and operations.
 *
 * Required operations:
 * - Quadratic form computation: x^T B^-1 x
 * - Inverse application: B^-1 x
 * - Forward application: B x
 * - Square root application: B^1/2 x
 * - Various matrix representations (diagonal, full, ensemble-based, etc.)
 *
 * @tparam T The backend type to be checked against this concept
 */
template <typename T>
concept BackgroundErrorCovarianceBackendType = requires(T t) {
  // Initialization and state
  { t.isInitialized() } -> std::convertible_to<bool>;

  // Representation support queries
  { t.supportsDiagonal() } -> std::convertible_to<bool>;
  { t.supportsEnsemble() } -> std::convertible_to<bool>;
  { t.supportsParametric() } -> std::convertible_to<bool>;
  { t.supportsHybrid() } -> std::convertible_to<bool>;
  { t.supportsFull() } -> std::convertible_to<bool>;

  // Localization support
  { t.supportsLocalization() } -> std::convertible_to<bool>;
  { t.setLocalizationRadius(double{}) } -> std::same_as<void>;
  { t.getLocalizationRadius() } -> std::convertible_to<double>;

  // Matrix size information
  { t.getMatrixSize() } -> std::convertible_to<size_t>;
  { t.getRank() } -> std::convertible_to<size_t>;
};

/**
 * @brief Concept for diagonal B matrix backends
 *
 * @details Specialized concept for backends supporting diagonal background
 * error covariance matrices (climatological errors).
 *
 * @tparam T The backend type to be checked
 */
template <typename T>
concept DiagonalCovarianceBackendType =
    BackgroundErrorCovarianceBackendType<T> && requires(T t) {
      // Diagonal-specific operations would be defined by the backend
      // For now, we just require the basic covariance interface
    };

/**
 * @brief Concept for ensemble-based B matrix backends
 *
 * @details Concept for backends supporting ensemble-based background error
 * covariance (NMC method, EnVar, hybrid methods).
 *
 * @tparam T The backend type to be checked
 */
template <typename T>
concept EnsembleCovarianceBackendType =
    BackgroundErrorCovarianceBackendType<T> && requires(T t) {
      // Ensemble size management
      { t.getEnsembleSize() } -> std::convertible_to<size_t>;
      { t.setEnsembleSize(size_t{}) } -> std::same_as<void>;

      // Ensemble data management
      { t.hasEnsembleData() } -> std::convertible_to<bool>;
      { t.computeEnsembleCovariance() } -> std::same_as<void>;

      // Inflation support
      { t.setInflationFactor(double{}) } -> std::same_as<void>;
      { t.getInflationFactor() } -> std::convertible_to<double>;
    };

/**
 * @brief Concept for parametric B matrix backends
 *
 * @details Concept for backends supporting parametric background error
 * covariance (recursive filters, spectral methods).
 *
 * @tparam T The backend type to be checked
 */
template <typename T>
concept ParametricCovarianceBackendType =
    BackgroundErrorCovarianceBackendType<T> && requires(T t) {
      // Correlation length scale management
      { t.setCorrelationLengthScale(double{}) } -> std::same_as<void>;
      { t.getCorrelationLengthScale() } -> std::convertible_to<double>;

      // Variance field management
      { t.setVarianceField(const void*) } -> std::same_as<void>;
      { t.hasVarianceField() } -> std::convertible_to<bool>;

      // Filter support
      { t.supportsRecursiveFilter() } -> std::convertible_to<bool>;
      { t.supportsSpectralFilter() } -> std::convertible_to<bool>;
    };

/**
 * @brief Concept for hybrid B matrix backends
 *
 * @details Concept for backends supporting hybrid background error covariance
 * (combination of static and ensemble-based covariance).
 *
 * @tparam T The backend type to be checked
 */
template <typename T>
concept HybridCovarianceBackendType =
    BackgroundErrorCovarianceBackendType<T> &&
    EnsembleCovarianceBackendType<T> && requires(T t) {
      // Hybrid weight management
      { t.setStaticWeight(double{}) } -> std::same_as<void>;
      { t.setEnsembleWeight(double{}) } -> std::same_as<void>;
      { t.getStaticWeight() } -> std::convertible_to<double>;
      { t.getEnsembleWeight() } -> std::convertible_to<double>;

      // Static covariance component
      { t.hasStaticCovariance() } -> std::convertible_to<bool>;
      { t.setStaticCovariance(const void*) } -> std::same_as<void>;
    };

/**
 * @brief Concept for full B matrix backends
 *
 * @details Concept for backends supporting full background error covariance
 * matrices (dense storage, for small problems).
 *
 * @tparam T The backend type to be checked
 */
template <typename T>
concept FullCovarianceBackendType =
    BackgroundErrorCovarianceBackendType<T> && requires(T t) {
      // Matrix storage queries
      { t.getStorageSize() } -> std::convertible_to<size_t>;
      { t.isSymmetric() } -> std::convertible_to<bool>;
      { t.isPositiveDefinite() } -> std::convertible_to<bool>;

      // Direct matrix access (for small problems)
      { t.getMatrixElement(size_t{}, size_t{}) } -> std::convertible_to<double>;
      {
        t.setMatrixElement(size_t{}, size_t{}, double{})
      } -> std::same_as<void>;

      // Decomposition support
      { t.supportsCholesky() } -> std::convertible_to<bool>;
      { t.supportsEigendecomposition() } -> std::convertible_to<bool>;
    };

}  // namespace metada::framework