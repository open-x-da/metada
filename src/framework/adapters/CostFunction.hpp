#pragma once

#include <memory>
#include <vector>

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

template <typename BackendTag>
  requires StateBackendType<BackendTag>
class Increment;

/**
 * @brief Cost function adapter for variational data assimilation
 *
 * @details This class implements the 4DVAR cost function:
 * J(x) = 1/2 * (x-xb)^T B^-1 (x-xb) + 1/2 * sum_i (y_i - H_i(M_i(x)))^T R_i^-1
 * (y_i - H_i(M_i(x)))
 *
 * Where:
 * - x is the analysis state
 * - xb is the background state
 * - B is the background error covariance
 * - y_i are observations at time i
 * - H_i is the observation operator at time i
 * - M_i is the model propagated to time i
 * - R_i is the observation error covariance at time i
 *
 * For 3DVAR: Only one time window (i=0)
 * For FGAT: Observations at proper times but single forward trajectory
 * For 4DVAR: Full time window with multiple observation times
 *
 * @tparam BackendTag The backend tag type
 */
template <typename BackendTag>
  requires StateBackendType<BackendTag>
class CostFunction : public NonCopyable {
 public:
  /** @brief Default constructor is deleted */
  CostFunction() = delete;

  /**
   * @brief Constructor for 4DVAR cost function
   *
   * @param config Configuration object
   * @param background Background state (first guess)
   * @param observations Vector of observations at different times
   * @param obs_operators Vector of observation operators for each time
   * @param model Model for forward propagation
   * @param bg_error_cov Background error covariance
   */
  CostFunction(const Config<BackendTag>& config,
               const State<BackendTag>& background,
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

    logger_.Info() << "CostFunction constructed for "
                   << getVariationalTypeName() << " with " << time_windows_
                   << " time windows";
  }

  /**
   * @brief Move constructor
   */
  CostFunction(CostFunction&& other) noexcept = default;

  /**
   * @brief Move assignment operator
   */
  CostFunction& operator=(CostFunction&& other) noexcept = default;

  /**
   * @brief Evaluate the cost function at state x
   *
   * @param state Current state vector
   * @return Cost function value J(x)
   */
  double evaluate(const State<BackendTag>& state) const {
    logger_.Debug() << "Evaluating cost function";

    double total_cost = 0.0;

    // Background term: 1/2 * (x - xb)^T B^-1 (x - xb)
    auto background_increment =
        Increment<BackendTag>::createFromDifference(state, background_);
    double bg_cost = 0.5 * bg_error_cov_.quadraticForm(background_increment);
    total_cost += bg_cost;

    logger_.Debug() << "Background cost: " << bg_cost;

    // Observation terms
    double obs_cost = 0.0;
    if (var_type_ == VariationalType::ThreeDVAR) {
      obs_cost = evaluateObservationCost3DVAR(state);
    } else if (var_type_ == VariationalType::FGAT) {
      obs_cost = evaluateObservationCostFGAT(state);
    } else {  // 4DVAR
      obs_cost = evaluateObservationCost4DVAR(state);
    }

    total_cost += obs_cost;
    logger_.Debug() << "Observation cost: " << obs_cost;
    logger_.Debug() << "Total cost: " << total_cost;

    return total_cost;
  }

  /**
   * @brief Compute the gradient of the cost function
   *
   * @param state Current state vector
   * @param gradient Output gradient vector
   */
  void gradient(const State<BackendTag>& state,
                Increment<BackendTag>& gradient) const {
    logger_.Debug() << "Computing cost function gradient";

    // Initialize gradient to zero
    gradient.zero();

    // Background term gradient: B^-1 (x - xb)
    auto background_increment =
        Increment<BackendTag>::createFromDifference(state, background_);
    auto bg_gradient = bg_error_cov_.applyInverse(background_increment);
    gradient += bg_gradient;

    // Observation term gradients
    if (var_type_ == VariationalType::ThreeDVAR) {
      computeObservationGradient3DVAR(state, gradient);
    } else if (var_type_ == VariationalType::FGAT) {
      computeObservationGradientFGAT(state, gradient);
    } else {  // 4DVAR
      computeObservationGradient4DVAR(state, gradient);
    }

    logger_.Debug() << "Gradient computed, norm: " << gradient.norm();
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

 private:
  enum class VariationalType {
    ThreeDVAR,  ///< 3DVAR: Single time, all obs at analysis time
    FGAT,       ///< FGAT: Obs at proper times, single trajectory
    FourDVAR    ///< 4DVAR: Full time window with adjoint
  };

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
   * @brief Evaluate observation cost for 3DVAR
   */
  double evaluateObservationCost3DVAR(const State<BackendTag>& state) const {
    // 3DVAR: All observations assumed at analysis time
    const auto& obs = observations_[0];
    const auto& obs_op = obs_operators_[0];

    auto simulated_obs = obs_op.apply(state, obs);
    auto obs_data = obs.template getData<std::vector<double>>();
    auto innovation = computeInnovation(obs, simulated_obs);

    // Debug print: show first 5 values of H(x), y, and innovation
    size_t n_debug = std::min<size_t>(5, simulated_obs.size());
    logger_.Info() << "Obs debug: idx | H(x) | y | y-H(x)";
    for (size_t i = 0; i < n_debug; ++i) {
      logger_.Info() << "Obs debug: " << i << " | " << simulated_obs[i] << " | "
                     << obs_data[i] << " | "
                     << (obs_data[i] - simulated_obs[i]);
    }

    return 0.5 * obs.quadraticForm(innovation);
  }

  /**
   * @brief Evaluate observation cost for FGAT
   */
  double evaluateObservationCostFGAT(const State<BackendTag>& state) const {
    // FGAT: Observations at proper times but single forward trajectory
    double total_cost = 0.0;
    auto current_state = state.clone();

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
      auto innovation = computeInnovation(obs, simulated_obs);
      total_cost += 0.5 * obs.quadraticForm(innovation);
    }

    return total_cost;
  }

  /**
   * @brief Evaluate observation cost for 4DVAR
   */
  double evaluateObservationCost4DVAR(const State<BackendTag>& state) const {
    // 4DVAR: Full time window evaluation
    return evaluateObservationCostFGAT(state);  // Same forward evaluation
  }

  /**
   * @brief Compute observation gradient for 3DVAR
   */
  void computeObservationGradient3DVAR(const State<BackendTag>& state,
                                       Increment<BackendTag>& gradient) const {
    const auto& obs = observations_[0];
    const auto& obs_op = obs_operators_[0];

    auto simulated_obs = obs_op.apply(state, obs);
    auto innovation = computeInnovation(obs, simulated_obs);
    auto weighted_innovation = obs.applyInverseCovariance(innovation);

    // H^T R^-1 (H(x) - y)
    auto obs_gradient = obs_op.applyAdjoint(weighted_innovation, state);
    gradient += obs_gradient;
  }

  /**
   * @brief Compute observation gradient for FGAT
   */
  void computeObservationGradientFGAT(const State<BackendTag>& state,
                                      Increment<BackendTag>& gradient) const {
    // For FGAT, we need to accumulate gradients but don't use full adjoint
    // model This is an approximation where we ignore model trajectory
    // sensitivity
    auto current_state = state.clone();

    for (size_t i = 0; i < observations_.size(); ++i) {
      const auto& obs = observations_[i];
      const auto& obs_op = obs_operators_[i];

      if (i > 0) {
        auto next_state = current_state.clone();
        model_.run(current_state, next_state);
        current_state = std::move(next_state);
      }

      auto simulated_obs = obs_op.apply(current_state, obs);
      auto innovation = computeInnovation(obs, simulated_obs);
      auto weighted_innovation = obs.applyInverseCovariance(innovation);
      auto obs_gradient =
          obs_op.applyAdjoint(weighted_innovation, current_state);

      // For FGAT, we approximate by not using model adjoint
      gradient += obs_gradient;
    }
  }

  /**
   * @brief Compute observation gradient for 4DVAR
   */
  void computeObservationGradient4DVAR(const State<BackendTag>& state,
                                       Increment<BackendTag>& gradient) const {
    // 4DVAR requires full adjoint model integration
    // Forward pass
    std::vector<State<BackendTag>> trajectory;
    trajectory.reserve(observations_.size());

    auto current_state = state.clone();
    trajectory.push_back(current_state.clone());

    for (size_t i = 1; i < observations_.size(); ++i) {
      auto next_state = current_state.clone();
      model_.run(current_state, next_state);
      current_state = std::move(next_state);
      trajectory.push_back(current_state.clone());
    }

    // Backward pass with adjoint
    auto adjoint_forcing = Increment<BackendTag>::createFromEntity(state);
    adjoint_forcing.zero();

    // Process observations in reverse order
    for (int i = observations_.size() - 1; i >= 0; --i) {
      const auto& obs = observations_[i];
      const auto& obs_op = obs_operators_[i];
      const auto& state_at_time = trajectory[i];

      auto simulated_obs = obs_op.apply(state_at_time, obs);
      auto innovation = computeInnovation(obs, simulated_obs);
      auto weighted_innovation = obs.applyInverseCovariance(innovation);
      auto obs_adjoint =
          obs_op.applyAdjoint(weighted_innovation, state_at_time);

      adjoint_forcing += obs_adjoint;

      // Adjoint model integration (if not at initial time)
      if (i > 0) {
        auto next_adjoint_forcing =
            Increment<BackendTag>::createFromEntity(trajectory[i - 1]);
        model_.runAdjoint(trajectory[i - 1], trajectory[i], adjoint_forcing,
                          next_adjoint_forcing);
        adjoint_forcing = std::move(next_adjoint_forcing);
      }
    }

    gradient += adjoint_forcing;
  }

  /**
   * @brief Compute innovation vector (y - H(x))
   */
  std::vector<double> computeInnovation(
      const Observation<BackendTag>& obs,
      const std::vector<double>& simulated_obs) const {
    const auto obs_data = obs.template getData<std::vector<double>>();
    std::vector<double> innovation(obs_data.size());

    for (size_t i = 0; i < innovation.size(); ++i) {
      innovation[i] = obs_data[i] - simulated_obs[i];
    }

    return innovation;
  }

  const State<BackendTag>& background_;
  const std::vector<Observation<BackendTag>>& observations_;
  const std::vector<ObsOperator<BackendTag>>& obs_operators_;
  Model<BackendTag>& model_;
  const BackgroundErrorCovariance<BackendTag>& bg_error_cov_;
  size_t time_windows_;
  VariationalType var_type_;
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework