#pragma once

#include <memory>
#include <string>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "Increment.hpp"
#include "Logger.hpp"
#include "NonCopyable.hpp"
#include "State.hpp"

namespace metada::framework {

// Forward declarations
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

/**
 * @brief Background error covariance adapter for variational data assimilation
 *
 * @details This class provides an interface for background error covariance
 * matrix operations (B matrix) used in variational data assimilation. It
 * supports various representations and operations including:
 * - Quadratic form evaluation: x^T B^-1 x
 * - Inverse covariance application: B^-1 x
 * - Forward covariance application: B x (for preconditioning)
 * - Square root operations: B^1/2 x
 *
 * The implementation can handle different B matrix representations:
 * - Full matrix (for small problems)
 * - Diagonal matrix (climatological errors)
 * - Ensemble-based covariance (NMC method, ensemble of differences)
 * - Parametric covariance (recursive filters, spectral methods)
 * - Hybrid covariance (combination of static and ensemble-based)
 *
 * @tparam BackendTag The backend tag type
 */
template <typename BackendTag>
class BackgroundErrorCovariance : public NonCopyable {
 public:
  using BackgroundErrorCovarianceBackend = typename traits::BackendTraits<
      BackendTag>::BackgroundErrorCovarianceBackend;

  /** @brief Default constructor is deleted */
  BackgroundErrorCovariance() = delete;

  /**
   * @brief Constructor with configuration
   *
   * @param config Configuration object containing B matrix parameters
   */
  explicit BackgroundErrorCovariance(const Config<BackendTag>& config)
      : backend_(config.backend()),
        representation_(determineRepresentation(config)),
        localization_enabled_(config.Get("localization_enabled").asBool()),
        localization_radius_(config.Get("localization_radius").asFloat()) {
    logger_.Info() << "BackgroundErrorCovariance constructed with "
                   << getRepresentationName() << " representation";

    if (localization_enabled_) {
      logger_.Info() << "Localization enabled with radius: "
                     << localization_radius_;
    }
  }

  /**
   * @brief Move constructor
   */
  BackgroundErrorCovariance(BackgroundErrorCovariance&& other) noexcept =
      default;

  /**
   * @brief Move assignment operator
   */
  BackgroundErrorCovariance& operator=(
      BackgroundErrorCovariance&& other) noexcept = default;

  /**
   * @brief Compute quadratic form x^T B^-1 x
   *
   * @param increment The increment vector x
   * @return Quadratic form value
   */
  double quadraticForm(const Increment<State<BackendTag>>& increment) const {
    logger_.Debug() << "Computing quadratic form x^T B^-1 x";

    switch (representation_) {
      case Representation::Diagonal:
        return computeQuadraticFormDiagonal(increment);
      case Representation::EnsembleBased:
        return computeQuadraticFormEnsemble(increment);
      case Representation::Parametric:
        return computeQuadraticFormParametric(increment);
      case Representation::Hybrid:
        return computeQuadraticFormHybrid(increment);
      case Representation::Full:
      default:
        return computeQuadraticFormFull(increment);
    }
  }

  /**
   * @brief Apply inverse covariance: B^-1 x
   *
   * @param increment Input increment x
   * @return Result of B^-1 x
   */
  Increment<State<BackendTag>> applyInverse(
      const Increment<State<BackendTag>>& increment) const {
    logger_.Debug() << "Applying B^-1";

    switch (representation_) {
      case Representation::Diagonal:
        return applyInverseDiagonal(increment);
      case Representation::EnsembleBased:
        return applyInverseEnsemble(increment);
      case Representation::Parametric:
        return applyInverseParametric(increment);
      case Representation::Hybrid:
        return applyInverseHybrid(increment);
      case Representation::Full:
      default:
        return applyInverseFull(increment);
    }
  }

  /**
   * @brief Apply forward covariance: B x
   *
   * @param increment Input increment x
   * @return Result of B x
   */
  Increment<State<BackendTag>> apply(
      const Increment<State<BackendTag>>& increment) const {
    logger_.Debug() << "Applying B";

    switch (representation_) {
      case Representation::Diagonal:
        return applyDiagonal(increment);
      case Representation::EnsembleBased:
        return applyEnsemble(increment);
      case Representation::Parametric:
        return applyParametric(increment);
      case Representation::Hybrid:
        return applyHybrid(increment);
      case Representation::Full:
      default:
        return applyFull(increment);
    }
  }

  /**
   * @brief Apply square root: B^1/2 x
   *
   * @param increment Input increment x
   * @return Result of B^1/2 x
   */
  Increment<State<BackendTag>> applySquareRoot(
      const Increment<State<BackendTag>>& increment) const {
    logger_.Debug() << "Applying B^1/2";

    switch (representation_) {
      case Representation::Diagonal:
        return applySquareRootDiagonal(increment);
      case Representation::EnsembleBased:
        return applySquareRootEnsemble(increment);
      case Representation::Parametric:
        return applySquareRootParametric(increment);
      case Representation::Hybrid:
        return applySquareRootHybrid(increment);
      case Representation::Full:
      default:
        return applySquareRootFull(increment);
    }
  }

  /**
   * @brief Get the representation type name
   */
  std::string getRepresentationName() const {
    switch (representation_) {
      case Representation::Diagonal:
        return "Diagonal";
      case Representation::EnsembleBased:
        return "Ensemble-based";
      case Representation::Parametric:
        return "Parametric";
      case Representation::Hybrid:
        return "Hybrid";
      case Representation::Full:
        return "Full";
      default:
        return "Unknown";
    }
  }

  /**
   * @brief Check if localization is enabled
   */
  bool isLocalizationEnabled() const { return localization_enabled_; }

  /**
   * @brief Get localization radius
   */
  double getLocalizationRadius() const { return localization_radius_; }

 private:
  enum class Representation {
    Full,           ///< Full covariance matrix
    Diagonal,       ///< Diagonal covariance matrix
    EnsembleBased,  ///< Ensemble-based covariance (NMC, EnVar)
    Parametric,     ///< Parametric covariance (recursive filters)
    Hybrid          ///< Hybrid static-ensemble covariance
  };

  /**
   * @brief Determine representation from configuration
   */
  Representation determineRepresentation(const Config<BackendTag>& config) {
    std::string type = config.Get("background_covariance_type").asString();
    if (type == "diagonal") return Representation::Diagonal;
    if (type == "ensemble") return Representation::EnsembleBased;
    if (type == "parametric") return Representation::Parametric;
    if (type == "hybrid") return Representation::Hybrid;
    if (type == "full") return Representation::Full;

    // Default to diagonal for simplicity
    return Representation::Diagonal;
  }

  // Diagonal representation methods
  double computeQuadraticFormDiagonal(
      const Increment<State<BackendTag>>& increment) const {
    // For diagonal B: x^T B^-1 x = sum_i (x_i^2 / sigma_i^2)
    return backend_.computeQuadraticFormDiagonal(increment.entity());
  }

  Increment<State<BackendTag>> applyInverseDiagonal(
      const Increment<State<BackendTag>>& increment) const {
    // For diagonal B: B^-1 x = x_i / sigma_i^2
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applyInverseDiagonal(increment.entity(), result.entity());
    return result;
  }

  Increment<State<BackendTag>> applyDiagonal(
      const Increment<State<BackendTag>>& increment) const {
    // For diagonal B: B x = x_i * sigma_i^2
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applyDiagonal(increment.entity(), result.entity());
    return result;
  }

  Increment<State<BackendTag>> applySquareRootDiagonal(
      const Increment<State<BackendTag>>& increment) const {
    // For diagonal B: B^1/2 x = x_i * sigma_i
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applySquareRootDiagonal(increment.entity(), result.entity());
    return result;
  }

  // Ensemble-based representation methods
  double computeQuadraticFormEnsemble(
      const Increment<State<BackendTag>>& increment) const {
    // For ensemble B: x^T B^-1 x using ensemble perturbations
    return backend_.computeQuadraticFormEnsemble(increment.entity());
  }

  Increment<State<BackendTag>> applyInverseEnsemble(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applyInverseEnsemble(increment.entity(), result.entity());
    return result;
  }

  Increment<State<BackendTag>> applyEnsemble(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applyEnsemble(increment.entity(), result.entity());
    return result;
  }

  Increment<State<BackendTag>> applySquareRootEnsemble(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applySquareRootEnsemble(increment.entity(), result.entity());
    return result;
  }

  // Parametric representation methods
  double computeQuadraticFormParametric(
      const Increment<State<BackendTag>>& increment) const {
    return backend_.computeQuadraticFormParametric(increment.entity());
  }

  Increment<State<BackendTag>> applyInverseParametric(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applyInverseParametric(increment.entity(), result.entity());
    return result;
  }

  Increment<State<BackendTag>> applyParametric(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applyParametric(increment.entity(), result.entity());
    return result;
  }

  Increment<State<BackendTag>> applySquareRootParametric(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applySquareRootParametric(increment.entity(), result.entity());
    return result;
  }

  // Hybrid representation methods
  double computeQuadraticFormHybrid(
      const Increment<State<BackendTag>>& increment) const {
    return backend_.computeQuadraticFormHybrid(increment.entity());
  }

  Increment<State<BackendTag>> applyInverseHybrid(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applyInverseHybrid(increment.entity(), result.entity());
    return result;
  }

  Increment<State<BackendTag>> applyHybrid(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applyHybrid(increment.entity(), result.entity());
    return result;
  }

  Increment<State<BackendTag>> applySquareRootHybrid(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applySquareRootHybrid(increment.entity(), result.entity());
    return result;
  }

  // Full matrix representation methods
  double computeQuadraticFormFull(
      const Increment<State<BackendTag>>& increment) const {
    return backend_.computeQuadraticFormFull(increment.entity());
  }

  Increment<State<BackendTag>> applyInverseFull(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applyInverseFull(increment.entity(), result.entity());
    return result;
  }

  Increment<State<BackendTag>> applyFull(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applyFull(increment.entity(), result.entity());
    return result;
  }

  Increment<State<BackendTag>> applySquareRootFull(
      const Increment<State<BackendTag>>& increment) const {
    auto result =
        Increment<State<BackendTag>>::createFromEntity(increment.entity());
    backend_.applySquareRootFull(increment.entity(), result.entity());
    return result;
  }

  BackgroundErrorCovarianceBackend backend_;
  Representation representation_;
  bool localization_enabled_;
  double localization_radius_;
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework