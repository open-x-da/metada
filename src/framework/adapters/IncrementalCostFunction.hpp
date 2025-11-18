#pragma once

#include "BackendTraits.hpp"
#include "BackgroundErrorCovariance.hpp"
#include "ConfigConcepts.hpp"
#include "ControlVariable.hpp"
#include "ControlVariableBackend.hpp"
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
      const BackgroundErrorCovariance<BackendTag>& bg_error_cov,
      const ControlVariableBackend<BackendTag>& control_backend)
      : background_(background),
        observations_(observations),
        obs_operators_(obs_operators),
        model_(model),
        bg_error_cov_(bg_error_cov),
        control_backend_(control_backend),
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
   * @brief Move assignment operator (deleted due to reference members)
   */
  IncrementalCostFunction& operator=(IncrementalCostFunction&&) noexcept =
      delete;

  /**
   * @brief Evaluate the incremental cost function at control variable v
   *
   * @details The cost function in control space is:
   *   J(v) = 1/2 ||v||² + 1/2 Σᵢ ||dᵢ - Hᵢ(δx(v))||²_Rᵢ
   *
   * where δx = U v is the transformation from control to state space.
   *
   * @param control Control variable v
   * @return Cost function value J(v)
   */
  double evaluate(const ControlVariable<BackendTag>& control) const {
    double total_cost = 0.0;

    // Background term in control space: 1/2 ||v||²
    // (The B^(-1) is absorbed into the transformation U)
    double bg_cost = 0.5 * control.dot(control);
    total_cost += bg_cost;

    // Observation term: 1/2 Σᵢ ||dᵢ - Hᵢ(U v)||²_Rᵢ
    // The observation cost methods work in control space
    double obs_cost = 0.0;
    if (var_type_ == VariationalType::ThreeDVAR) {
      obs_cost = evaluateObservationCost3DVAR(control);
    } else if (var_type_ == VariationalType::FGAT) {
      obs_cost = evaluateObservationCostFGAT(control);
    } else {
      obs_cost = evaluateObservationCost4DVAR(control);
    }

    total_cost += obs_cost;

    logger_.Info() << "Cost diagnostics: background=" << bg_cost
                   << ", observation=" << obs_cost << ", total=" << total_cost;

    return total_cost;
  }

  /**
   * @brief Compute the gradient of the incremental cost function in control
   * space
   *
   * @details The gradient in control space is:
   *   ∇_v J = v + ∇_v J_obs
   *
   * where:
   * - v is the control variable (background term gradient)
   * - ∇_v J_obs is the observation term gradient in control space
   *
   * All transformations between control and state space are handled internally
   * by the observation gradient computation methods.
   *
   * @param control Current control variable v
   * @param gradient Output gradient vector ∇_v J (in control space)
   */
  void gradient(const ControlVariable<BackendTag>& control,
                ControlVariable<BackendTag>& gradient) const {
    gradient.zero();

    // Background term gradient in control space: ∇_v J_b = v
    // (The B^(-1) is absorbed into the control variable transformation U)
    gradient.axpy(1.0, control);

    // Compute observation term gradient directly in control space
    auto control_gradient_obs =
        control_backend_.createControlVariable(control.geometry());
    control_gradient_obs.zero();

    if (var_type_ == VariationalType::ThreeDVAR) {
      computeObservationGradient3DVAR(control, control_gradient_obs);
    } else if (var_type_ == VariationalType::FGAT) {
      computeObservationGradientFGAT(control, control_gradient_obs);
    } else {
      computeObservationGradient4DVAR(control, control_gradient_obs);
    }

    // Add observation term gradient to total gradient
    gradient.axpy(1.0, control_gradient_obs);

    logGradientDiagnostics(control, control_gradient_obs, gradient);
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
    control_backend_.addIncrementToState(increment, state);
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
   * @brief Evaluate observation cost for 3DVAR in control space
   */
  double evaluateObservationCost3DVAR(
      const ControlVariable<BackendTag>& control) const {
    const auto& obs = observations_[0];
    const auto& obs_op = obs_operators_[0];

    // Apply tangent linear in control space: H' U v
    auto simulated_increment =
        obs_op.applyTangentLinear(control, background_, obs);

    if (obs_op.hasNativeIncrementalCost(background_, obs)) {
      (void)simulated_increment;
      auto jo = obs_op.computeObservationCost(background_, obs);
      if (jo.has_value()) {
        return *jo;
      }
    }

    auto residual = computeInnovation(innovations_[0], simulated_increment);
    auto weighted_residual = obs.applyInverseCovariance(residual);

    return computeObservationCostFromResidual(residual, weighted_residual);
  }

  /**
   * @brief Evaluate observation cost for FGAT in control space
   */
  double evaluateObservationCostFGAT(
      const ControlVariable<BackendTag>& control) const {
    double total_cost = 0.0;

    for (size_t i = 0; i < observations_.size(); ++i) {
      const auto& obs = observations_[i];
      const auto& obs_op = obs_operators_[i];

      // Apply tangent linear in control space: H' U v
      auto simulated_increment =
          obs_op.applyTangentLinear(control, background_, obs);

      auto residual = computeInnovation(innovations_[i], simulated_increment);
      auto weighted_residual = obs.applyInverseCovariance(residual);

      total_cost +=
          computeObservationCostFromResidual(residual, weighted_residual);
    }

    return total_cost;
  }

  /**
   * @brief Evaluate observation cost for 4DVAR in control space
   */
  double evaluateObservationCost4DVAR(
      const ControlVariable<BackendTag>& control) const {
    return evaluateObservationCostFGAT(control);
  }

  /**
   * @brief Compute observation gradient for 3DVAR in control space
   *
   * @details Computes ∇_v J_obs directly in control space. The ObsOperator
   * handles all internal transformations between control and state space.
   *
   * @param control Control variable v
   * @param gradient Output gradient in control space
   */
  void computeObservationGradient3DVAR(
      const ControlVariable<BackendTag>& control,
      ControlVariable<BackendTag>& gradient) const {
    const auto& obs = observations_[0];
    const auto& obs_op = obs_operators_[0];

    // Apply tangent linear in control space: H' U v
    // Use direct da_transform_vtoy when available (matches WRFDA's
    // da_calculate_gradj workflow during CG iterations)
    auto simulated_increment =
        obs_op.applyTangentLinearDirect(control, background_, obs);

    if (obs_op.hasNativeIncrementalCost(background_, obs)) {
      (void)simulated_increment;
      std::vector<double> empty;
      // Apply adjoint in control space: U^T H^T
      auto control_gradient = obs_op.applyAdjoint(empty, background_, obs);
      gradient = control_gradient;
    } else {
      auto residual = computeInnovation(innovations_[0], simulated_increment);
      auto weighted_residual = obs.applyInverseCovariance(residual);
      // Apply pure adjoint in control space: U^T H'^T
      auto control_gradient =
          obs_op.applyAdjointPure(weighted_residual, background_, obs);
      gradient -= control_gradient;
    }
  }

  /**
   * @brief Compute observation gradient for FGAT in control space
   *
   * @details Computes ∇_v J_obs for FGAT by accumulating gradients from all
   * time windows. The ObsOperator handles all internal transformations.
   *
   * @param control Control variable v
   * @param gradient Output gradient in control space
   */
  void computeObservationGradientFGAT(
      const ControlVariable<BackendTag>& control,
      ControlVariable<BackendTag>& gradient) const {
    for (size_t i = 0; i < observations_.size(); ++i) {
      const auto& obs = observations_[i];
      const auto& obs_op = obs_operators_[i];

      // Apply tangent linear in control space: H' U v
      // Use direct da_transform_vtoy when available (matches WRFDA's
      // da_calculate_gradj workflow during CG iterations)
      auto simulated_increment =
          obs_op.applyTangentLinearDirect(control, background_, obs);

      if (obs_op.hasNativeIncrementalCost(background_, obs)) {
        (void)simulated_increment;
        std::vector<double> empty;
        // Apply adjoint in control space: U^T H^T
        auto control_gradient = obs_op.applyAdjoint(empty, background_, obs);
        gradient += control_gradient;
        continue;
      }

      auto residual = computeInnovation(innovations_[i], simulated_increment);
      auto weighted_residual = obs.applyInverseCovariance(residual);

      // Apply pure adjoint in control space: U^T H'^T
      auto control_gradient =
          obs_op.applyAdjointPure(weighted_residual, background_, obs);
      gradient -= control_gradient;
    }
  }

  /**
   * @brief Compute observation gradient for 4DVAR in control space
   *
   * @details For now, 4DVAR uses the same implementation as FGAT.
   * Full 4DVAR would require model adjoint integration.
   *
   * @param control Control variable v
   * @param gradient Output gradient in control space
   */
  void computeObservationGradient4DVAR(
      const ControlVariable<BackendTag>& control,
      ControlVariable<BackendTag>& gradient) const {
    computeObservationGradientFGAT(control, gradient);
  }

  void logGradientDiagnostics(
      const ControlVariable<BackendTag>& control,
      const ControlVariable<BackendTag>& control_obs_gradient,
      const ControlVariable<BackendTag>& total_gradient) const {
    const double control_norm = control.norm();
    const double control_obs_norm = control_obs_gradient.norm();
    const double total_norm = total_gradient.norm();

    logger_.Info() << "Gradient diagnostics: |v|=" << control_norm
                   << ", |grad_obs_control|=" << control_obs_norm
                   << ", |grad_total|=" << total_norm;
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

  double computeObservationCostFromResidual(
      const std::vector<double>& residual,
      const std::vector<double>& weighted_residual) const {
    double cost = 0.0;
    for (size_t i = 0; i < residual.size(); ++i) {
      cost += residual[i] * weighted_residual[i];
    }
    return 0.5 * cost;
  }

  const State<BackendTag>& background_;
  const std::vector<Observation<BackendTag>>& observations_;
  const std::vector<ObsOperator<BackendTag>>& obs_operators_;
  Model<BackendTag>& model_;
  const BackgroundErrorCovariance<BackendTag>& bg_error_cov_;
  const ControlVariableBackend<BackendTag>& control_backend_;
  size_t time_windows_;
  VariationalType var_type_;
  std::vector<std::vector<double>>
      innovations_;  ///< Pre-computed innovation vectors dᵢ
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework
