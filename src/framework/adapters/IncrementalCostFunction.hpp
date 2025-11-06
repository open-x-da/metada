#pragma once

#include "BackendTraits.hpp"
#include "BackgroundErrorCovariance.hpp"
#include "ConfigConcepts.hpp"
#include "Increment.hpp"
#include "Logger.hpp"
#include "Model.hpp"
#include "NonCopyable.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "State.hpp"

namespace metada::framework {

// Forward declarations
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

/**
 * @brief Incremental cost function for variational data assimilation
 *
 * @details This class implements the incremental variational formulation:
 * J(δx) = 1/2 * δx^T B^-1 δx + 1/2 * Σᵢ (dᵢ - Hᵢ(Mᵢ(xb + δx)))^T Rᵢ^-1 (dᵢ -
 * Hᵢ(Mᵢ(xb + δx)))
 *
 * Where:
 * - δx is the analysis increment (δx = x - xb)
 * - xb is the background state (first guess)
 * - B is the background error covariance
 * - dᵢ = yᵢ - Hᵢ(Mᵢ(xb)) is the innovation vector (O-B)
 * - Hᵢ is the observation operator at time i
 * - Mᵢ is the model propagated to time i
 * - Rᵢ is the observation error covariance at time i
 *
 * This formulation offers better numerical conditioning and is the standard
 * approach used in operational systems like WRFDA.
 *
 * @tparam BackendTag The backend tag type
 */
template <typename BackendTag>
  requires StateBackendType<BackendTag>
class IncrementalCostFunction : public NonCopyable {
 public:
  /** @brief Default constructor is deleted */
  IncrementalCostFunction() = delete;

  /**
   * @brief Constructor for incremental cost function
   *
   * @param config Configuration object
   * @param background Background state (first guess)
   * @param observations Vector of observations at different times
   * @param obs_operators Vector of observation operators for each time
   * @param model Model for forward propagation
   * @param bg_error_cov Background error covariance
   */
  IncrementalCostFunction(
      const Config<BackendTag>& config, const State<BackendTag>& background,
      const std::vector<Observation<BackendTag>>& observations,
      const std::vector<ObsOperator<BackendTag>>& obs_operators,
      Model<BackendTag>& model,
      const BackgroundErrorCovariance<BackendTag>& bg_error_cov)
      : background_(background),
        observations_(observations),
        obs_operators_(obs_operators),
        model_(model),
        bg_error_cov_(bg_error_cov),
        time_windows_(observations.size()),
        var_type_(determineVariationalType(config)) {
    if (observations_.size() != obs_operators_.size()) {
      throw std::runtime_error(
          "Number of observations must match number of observation operators");
    }

    // Pre-compute innovation vectors dᵢ = yᵢ - Hᵢ(Mᵢ(xb))
    precomputeInnovations();

    logger_.Info() << "IncrementalCostFunction constructed for "
                   << getVariationalTypeName() << " with " << time_windows_
                   << " time windows";
  }

  /**
   * @brief Move constructor
   */
  IncrementalCostFunction(IncrementalCostFunction&& other) noexcept = default;

  /**
   * @brief Move assignment operator
   */
  IncrementalCostFunction& operator=(IncrementalCostFunction&& other) noexcept =
      default;

  /**
   * @brief Evaluate the incremental cost function at increment δx
   *
   * @param increment Analysis increment δx
   * @return Cost function value J(δx)
   */
  double evaluate(const Increment<BackendTag>& increment) const {
    double total_cost = 0.0;

    // Background term: 1/2 * δx^T B^-1 δx
    double bg_cost = 0.5 * bg_error_cov_.quadraticForm(increment);
    total_cost += bg_cost;

    // Observation terms: 1/2 * Σᵢ (dᵢ - Hᵢ(Mᵢ(xb + δx)))^T Rᵢ^-1 (dᵢ - Hᵢ(Mᵢ(xb
    // + δx)))
    double obs_cost = 0.0;
    if (var_type_ == VariationalType::ThreeDVAR) {
      obs_cost = evaluateObservationCost3DVAR(increment);
    } else if (var_type_ == VariationalType::FGAT) {
      obs_cost = evaluateObservationCostFGAT(increment);
    } else {  // 4DVAR
      obs_cost = evaluateObservationCost4DVAR(increment);
    }

    total_cost += obs_cost;

    return total_cost;
  }

  /**
   * @brief Compute the gradient of the incremental cost function
   *
   * @param increment Current increment δx
   * @param gradient Output gradient vector ∇J(δx)
   */
  void gradient(const Increment<BackendTag>& increment,
                Increment<BackendTag>& gradient) const {
    // Initialize gradient to zero
    gradient.zero();

    // Background term gradient: B^-1 δx
    auto bg_gradient = bg_error_cov_.applyInverse(increment);
    gradient += bg_gradient;

    // Observation term gradients
    if (var_type_ == VariationalType::ThreeDVAR) {
      computeObservationGradient3DVAR(increment, gradient);
    } else if (var_type_ == VariationalType::FGAT) {
      computeObservationGradientFGAT(increment, gradient);
    } else {  // 4DVAR
      computeObservationGradient4DVAR(increment, gradient);
    }
  }

  /**
   * @brief Get the number of time windows
   */
  size_t getTimeWindows() const { return time_windows_; }

  /**
   * @brief Get the variational type (3DVAR, FGAT, 4DVAR)
   */
  std::string getVariationalTypeName() const {
    switch (var_type_) {
      case VariationalType::ThreeDVAR:
        return "3DVAR";
      case VariationalType::FGAT:
        return "FGAT";
      case VariationalType::FourDVAR:
        return "4DVAR";
      default:
        return "Unknown";
    }
  }

  /**
   * @brief Get the background state
   */
  const State<BackendTag>& getBackground() const { return background_; }

  /**
   * @brief Get the pre-computed innovation vectors
   */
  const std::vector<std::vector<double>>& getInnovations() const {
    return innovations_;
  }

 private:
  enum class VariationalType {
    ThreeDVAR,  ///< 3DVAR: Single time, all obs at analysis time
    FGAT,       ///< FGAT: Obs at proper times, single trajectory
    FourDVAR    ///< 4DVAR: Full time window with adjoint
  };

  /**
   * @brief Add increment to state: state = state + increment
   */
  void addIncrementToState(const Increment<BackendTag>& increment,
                           State<BackendTag>& state) const {
    state += increment;
  }

  /**
   * @brief Determine variational type from configuration
   */
  VariationalType determineVariationalType(const Config<BackendTag>& config) {
    std::string type = config.Get("variational_type").asString();
    if (type == "3DVAR") return VariationalType::ThreeDVAR;
    if (type == "FGAT") return VariationalType::FGAT;
    if (type == "4DVAR") return VariationalType::FourDVAR;

    // Default to 4DVAR if multiple time windows, 3DVAR if single
    return (time_windows_ > 1) ? VariationalType::FourDVAR
                               : VariationalType::ThreeDVAR;
  }

  /**
   * @brief Pre-compute innovation vectors dᵢ = yᵢ - Hᵢ(Mᵢ(xb))
   */
  void precomputeInnovations() {
    innovations_.resize(observations_.size());

    if (var_type_ == VariationalType::ThreeDVAR) {
      // 3DVAR: All observations at analysis time
      const auto& obs = observations_[0];
      const auto& obs_op = obs_operators_[0];

      auto simulated_obs = obs_op.apply(background_, obs);
      auto obs_data = obs.getObservationValues();
      innovations_[0] = computeInnovation(obs_data, simulated_obs);

    } else if (var_type_ == VariationalType::FGAT ||
               var_type_ == VariationalType::FourDVAR) {
      // FGAT/4DVAR: Observations at proper times
      auto current_state = background_.clone();

      for (size_t i = 0; i < observations_.size(); ++i) {
        const auto& obs = observations_[i];
        const auto& obs_op = obs_operators_[i];

        // Propagate state to observation time (if needed)
        if (i > 0) {
          auto next_state = current_state.clone();
          model_.run(current_state, next_state);
          current_state = std::move(next_state);
        }

        auto simulated_obs = obs_op.apply(current_state, obs);
        auto obs_data = obs.getObservationValues();
        innovations_[i] = computeInnovation(obs_data, simulated_obs);
      }
    }

    logger_.Info() << "Pre-computed innovation vectors for " << time_windows_
                   << " time windows";
  }

  /**
   * @brief Evaluate observation cost for 3DVAR (incremental form)
   */
  double evaluateObservationCost3DVAR(
      const Increment<BackendTag>& increment) const {
    // 3DVAR: All observations assumed at analysis time
    const auto& obs = observations_[0];
    const auto& obs_op = obs_operators_[0];

    // For incremental 3D-Var, we only need H'(δx) where δx is the increment
    // The innovation vector d = H(xb) - y is already pre-computed
    auto simulated_increment =
        obs_op.applyTangentLinear(increment, background_, obs);

    // Compute observation cost using WRFDA's proven formula:
    // Jo = 0.5 * re^T * R^{-1} * re where re = (O-B) - H'(δx)
    // This uses WRFDA's da_calculate_residual + da_calculate_grady internally
    return obs_op.backend().computeObservationCost(background_.backend(),
                                                   obs.backend());
  }

  /**
   * @brief Evaluate observation cost for FGAT (incremental form)
   */
  double evaluateObservationCostFGAT(
      const Increment<BackendTag>& increment) const {
    // FGAT: Observations at proper times but single forward trajectory
    double total_cost = 0.0;

    for (size_t i = 0; i < observations_.size(); ++i) {
      const auto& obs = observations_[i];
      const auto& obs_op = obs_operators_[i];

      // For incremental FGAT, we only need H'(δx) at each observation time
      // The innovation vectors d = H(xb) - y are already pre-computed
      // Note: For FGAT, we use the same increment at all times (no model
      // propagation)
      auto simulated_increment =
          obs_op.applyTangentLinear(increment, background_, obs);

      // Compute observation cost using WRFDA's proven formula
      total_cost += obs_op.backend().computeObservationCost(
          background_.backend(), obs.backend());
    }

    return total_cost;
  }

  /**
   * @brief Evaluate observation cost for 4DVAR (incremental form)
   */
  double evaluateObservationCost4DVAR(
      const Increment<BackendTag>& increment) const {
    // 4DVAR: Full time window evaluation
    return evaluateObservationCostFGAT(increment);  // Same forward evaluation
  }

  /**
   * @brief Compute observation gradient for 3DVAR (incremental form)
   */
  void computeObservationGradient3DVAR(const Increment<BackendTag>& increment,
                                       Increment<BackendTag>& gradient) const {
    const auto& obs = observations_[0];
    const auto& obs_op = obs_operators_[0];

    // For incremental 3D-Var, we only need H'(δx) where δx is the increment
    // The innovation vector d = H(xb) - y is already pre-computed
    auto simulated_increment =
        obs_op.applyTangentLinear(increment, background_, obs);

    // Compute residual: d - H'(δx) where d is pre-computed innovation
    if (innovations_[0].size() != simulated_increment.size()) {
      throw std::runtime_error(
          "Size mismatch in gradient computation: innovations.size()=" +
          std::to_string(innovations_[0].size()) +
          ", simulated_increment.size()=" +
          std::to_string(simulated_increment.size()));
    }

    std::vector<double> residual(innovations_[0].size());
    for (size_t i = 0; i < innovations_[0].size(); ++i) {
      residual[i] = innovations_[0][i] - simulated_increment[i];
    }

    // Compute norm manually for std::vector<double>
    double residual_norm = 0.0;
    for (const auto& val : residual) {
      residual_norm += val * val;
    }
    residual_norm = std::sqrt(residual_norm);

    // Apply R^{-1} weighting
    // NOTE: For WRF backend, this computation is redundant.
    // WRFObsOperator::applyAdjoint ignores weighted_residual and uses WRFDA's
    // proven da_calculate_residual + da_calculate_grady workflow instead. Kept
    // here for backend-agnosticism.
    auto weighted_residual = obs.applyInverseCovariance(residual);

    // Compute norm manually for std::vector<double>
    double weighted_norm = 0.0;
    for (const auto& val : weighted_residual) {
      weighted_norm += val * val;
    }
    weighted_norm = std::sqrt(weighted_norm);

    // Compute observation gradient: -H'^T R^{-1} [d - H'(δx)]
    // This is the adjoint of the observation operator applied to the weighted
    // residual. The NEGATIVE sign comes from WRFDA's da_calculate_grady:
    // jo_grad_y = -R^{-1} · re NOTE: WRF backend ignores weighted_residual
    // parameter and uses WRFDA's internal workflow (da_calculate_residual +
    // da_calculate_grady + da_transform_xtoy_adj)
    auto obs_gradient =
        obs_op.applyAdjoint(weighted_residual, background_, obs);
    gradient += obs_gradient;  // ADD the observation gradient (already has
                               // negative sign)
  }

  /**
   * @brief Compute observation gradient for FGAT (incremental form)
   */
  void computeObservationGradientFGAT(const Increment<BackendTag>& increment,
                                      Increment<BackendTag>& gradient) const {
    // For incremental FGAT, we only need H'(δx) at each observation time
    // The innovation vectors d = H(xb) - y are already pre-computed

    for (size_t i = 0; i < observations_.size(); ++i) {
      const auto& obs = observations_[i];
      const auto& obs_op = obs_operators_[i];

      // For incremental FGAT, we only need H'(δx) at each observation time
      // Note: For FGAT, we use the same increment at all times (no model
      // propagation)
      auto simulated_increment =
          obs_op.applyTangentLinear(increment, background_, obs);
      auto residual = computeInnovation(innovations_[i], simulated_increment);
      auto weighted_residual = obs.applyInverseCovariance(residual);
      auto obs_gradient =
          obs_op.applyAdjoint(weighted_residual, background_, obs);

      // For FGAT, we approximate by not using model adjoint
      // Note: NEGATIVE sign from chain rule!
      gradient -= obs_gradient;
    }
  }

  /**
   * @brief Compute observation gradient for 4DVAR (incremental form)
   */
  void computeObservationGradient4DVAR(const Increment<BackendTag>& increment,
                                       Increment<BackendTag>& gradient) const {
    // 4DVAR requires full adjoint model integration
    // Forward pass: compute trajectory for xb + δx
    std::vector<State<BackendTag>> trajectory;
    trajectory.reserve(observations_.size());

    auto current_state = background_.clone();
    addIncrementToState(increment, current_state);
    trajectory.push_back(current_state.clone());

    for (size_t i = 1; i < observations_.size(); ++i) {
      auto next_state = current_state.clone();
      model_.run(current_state, next_state);
      current_state = std::move(next_state);
      trajectory.push_back(current_state.clone());
    }

    // Backward pass with adjoint
    auto adjoint_forcing = Increment<BackendTag>::createFromGeometry(
        background_.geometry()->backend());
    adjoint_forcing.zero();

    // Process observations in reverse order
    for (int i = observations_.size() - 1; i >= 0; --i) {
      const auto& obs = observations_[i];
      const auto& obs_op = obs_operators_[i];
      const auto& state_at_time = trajectory[i];

      auto simulated_obs = obs_op.apply(state_at_time, obs);
      auto residual = computeInnovation(innovations_[i], simulated_obs);
      auto weighted_residual = obs.applyInverseCovariance(residual);
      auto obs_adjoint =
          obs_op.applyAdjoint(weighted_residual, state_at_time, obs);

      // Note: NEGATIVE sign from chain rule!
      adjoint_forcing -= obs_adjoint;

      // Adjoint model integration (if not at initial time)
      if (i > 0) {
        auto next_adjoint_forcing = Increment<BackendTag>::createFromGeometry(
            trajectory[i - 1].geometry()->backend());
        model_.runAdjoint(trajectory[i - 1], trajectory[i], adjoint_forcing,
                          next_adjoint_forcing);
        adjoint_forcing = std::move(next_adjoint_forcing);
      }
    }

    gradient += adjoint_forcing;
  }

  /**
   * @brief Compute residual vector (d - H(xb + δx))
   */
  std::vector<double> computeInnovation(
      const std::vector<double>& innovation,
      const std::vector<double>& simulated_obs) const {
    // Debug: Check vector sizes before accessing
    if (innovation.size() != simulated_obs.size()) {
      throw std::runtime_error(
          "computeInnovation: Size mismatch - innovation.size() = " +
          std::to_string(innovation.size()) +
          ", simulated_obs.size() = " + std::to_string(simulated_obs.size()));
    }

    std::vector<double> residual(innovation.size());

    for (size_t i = 0; i < residual.size(); ++i) {
      residual[i] = innovation[i] - simulated_obs[i];
    }

    return residual;
  }

  const State<BackendTag>& background_;
  const std::vector<Observation<BackendTag>>& observations_;
  const std::vector<ObsOperator<BackendTag>>& obs_operators_;
  Model<BackendTag>& model_;
  const BackgroundErrorCovariance<BackendTag>& bg_error_cov_;
  size_t time_windows_;
  VariationalType var_type_;
  std::vector<std::vector<double>>
      innovations_;  ///< Pre-computed innovation vectors dᵢ
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework
