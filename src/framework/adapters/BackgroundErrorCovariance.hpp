#pragma once

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
        localization_enabled_(config.HasKey("localization_enabled")
                                  ? config.Get("localization_enabled").asBool()
                                  : false),
        localization_radius_(config.HasKey("localization_radius")
                                 ? config.Get("localization_radius").asFloat()
                                 : 0.0) {
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
  double quadraticForm(const Increment<BackendTag>& increment) const {
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
  Increment<BackendTag> applyInverse(
      const Increment<BackendTag>& increment) const {
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
  Increment<BackendTag> apply(const Increment<BackendTag>& increment) const {
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
  Increment<BackendTag> applySquareRoot(
      const Increment<BackendTag>& increment) const {
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
    if (!config.HasKey("background_covariance_type")) {
      // Default to diagonal for simplicity
      return Representation::Diagonal;
    }

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
      const Increment<BackendTag>& increment) const {
    // For diagonal B: x^T B^-1 x = sum_i (x_i^2 / sigma_i^2)
    return backend_.computeQuadraticFormDiagonal(increment.backend());
  }

  Increment<BackendTag> applyInverseDiagonal(
      const Increment<BackendTag>& increment) const {
    // For diagonal B: B^-1 x = x_i / sigma_i^2
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applyInverseDiagonal(increment.backend(), result.backend());
    return result;
  }

  Increment<BackendTag> applyDiagonal(
      const Increment<BackendTag>& increment) const {
    // For diagonal B: B x = x_i * sigma_i^2
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applyDiagonal(increment.backend(), result.backend());
    return result;
  }

  Increment<BackendTag> applySquareRootDiagonal(
      const Increment<BackendTag>& increment) const {
    // For diagonal B: B^1/2 x = x_i * sigma_i
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applySquareRootDiagonal(increment.backend(), result.backend());
    return result;
  }

  // Ensemble-based representation methods
  double computeQuadraticFormEnsemble(
      const Increment<BackendTag>& increment) const {
    // For ensemble B: x^T B^-1 x using ensemble perturbations
    return backend_.computeQuadraticFormEnsemble(increment.backend());
  }

  Increment<BackendTag> applyInverseEnsemble(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applyInverseEnsemble(increment.backend(), result.backend());
    return result;
  }

  Increment<BackendTag> applyEnsemble(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applyEnsemble(increment.backend(), result.backend());
    return result;
  }

  Increment<BackendTag> applySquareRootEnsemble(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applySquareRootEnsemble(increment.backend(), result.backend());
    return result;
  }

  // Parametric representation methods
  double computeQuadraticFormParametric(
      const Increment<BackendTag>& increment) const {
    return backend_.computeQuadraticFormParametric(increment.backend());
  }

  Increment<BackendTag> applyInverseParametric(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applyInverseParametric(increment.backend(), result.backend());
    return result;
  }

  Increment<BackendTag> applyParametric(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applyParametric(increment.backend(), result.backend());
    return result;
  }

  Increment<BackendTag> applySquareRootParametric(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applySquareRootParametric(increment.backend(), result.backend());
    return result;
  }

  // Hybrid representation methods
  double computeQuadraticFormHybrid(
      const Increment<BackendTag>& increment) const {
    return backend_.computeQuadraticFormHybrid(increment.backend());
  }

  Increment<BackendTag> applyInverseHybrid(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applyInverseHybrid(increment.backend(), result.backend());
    return result;
  }

  Increment<BackendTag> applyHybrid(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applyHybrid(increment.backend(), result.backend());
    return result;
  }

  Increment<BackendTag> applySquareRootHybrid(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applySquareRootHybrid(increment.backend(), result.backend());
    return result;
  }

  // Full matrix representation methods
  double computeQuadraticFormFull(
      const Increment<BackendTag>& increment) const {
    return backend_.computeQuadraticFormFull(increment.backend());
  }

  Increment<BackendTag> applyInverseFull(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applyInverseFull(increment.backend(), result.backend());
    return result;
  }

  Increment<BackendTag> applyFull(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applyFull(increment.backend(), result.backend());
    return result;
  }

  Increment<BackendTag> applySquareRootFull(
      const Increment<BackendTag>& increment) const {
    auto result =
        Increment<BackendTag>::createFromGeometry(increment.geometry());
    backend_.applySquareRootFull(increment.backend(), result.backend());
    return result;
  }

  BackgroundErrorCovarianceBackend backend_;
  Representation representation_;
  bool localization_enabled_;
  double localization_radius_;
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework