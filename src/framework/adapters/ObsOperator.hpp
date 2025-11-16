#pragma once

#include <optional>
#include <string>
#include <vector>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "ControlVariable.hpp"
#include "ControlVariableBackend.hpp"
#include "Logger.hpp"
#include "NonCopyable.hpp"
#include "ObsOperatorConcepts.hpp"
#include "ObservationConcepts.hpp"

namespace metada::framework {

// Forward declaration
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

template <typename BackendTag>
  requires ObservationBackendType<BackendTag>
class Observation;

template <typename BackendTag>
class State;

template <typename BackendTag>
class ControlVariableBackend;

/**
 * @brief Adapter class for observation operators in data assimilation systems
 *
 * @details This class provides a type-safe interface for observation operators
 * that map between model state space and observation space. It wraps a backend
 * implementation that satisfies the ObsOperatorBackendType concept, providing
 * consistent access patterns across different backend implementations.
 *
 * The ObsOperator class inherits from NonCopyable to prevent unintended
 * copying, but supports move semantics for efficient resource management.
 *
 * Features:
 * - Forward observation operator (H): Maps state to observation space
 * - Type-safe data access through templates
 * - Delegation to backend implementation
 * - Comprehensive error handling
 * - Move semantics for efficient resource management
 *
 * @tparam BackendTag The backend tag type that must satisfy
 * ObsOperatorBackendType concept
 *
 * @see ObsOperatorBackendType
 * @see NonCopyable
 */
template <typename BackendTag>
  requires ObsOperatorBackendType<BackendTag>
class ObsOperator : public NonCopyable {
 public:
  /** @brief Type alias for the backend implementation type */
  using ObsOperatorBackend =
      typename traits::BackendTraits<BackendTag>::ObsOperatorBackend;

  /** @brief Default constructor is deleted since we need a backend */
  ObsOperator() = delete;

 public:
  /**
   * @brief Constructor that initializes observation operator with configuration
   *
   * @details Creates and initializes the observation operator backend using the
   * provided configuration object. The backend handles its own initialization
   * in its constructor.
   *
   * @tparam ConfigBackend The configuration backend type
   * @param config Configuration object containing initialization parameters
   * @param control_backend Control variable backend for transformations
   */
  explicit ObsOperator(
      const Config<BackendTag>& config,
      const ControlVariableBackend<BackendTag>& control_backend)
      : backend_([&]() {
          using CB = typename traits::BackendTraits<BackendTag>::ConfigBackend;
          using B = ObsOperatorBackend;
          const CB& cfg = config.backend();
          if constexpr (requires { B{cfg, control_backend}; }) {
            return B(cfg, control_backend);
          } else {
            return B(cfg);
          }
        }()),
        control_backend_(&control_backend) {
    logger_.Info() << "ObsOperator constructed with control variable backend";
  }

  /**
   * @brief Constructor that initializes observation operator with configuration
   * and observation
   *
   * @details Creates and initializes the observation operator backend using the
   * provided configuration object and observation. Derives required variables
   * from the observation and augments the configuration before passing to the
   * backend.
   *
   * @param config Configuration object containing initialization parameters
   * @param obs Observation object to derive required variables from
   * @param control_backend Control variable backend for transformations
   */
  ObsOperator(const Config<BackendTag>& config,
              const Observation<BackendTag>& obs,
              const ControlVariableBackend<BackendTag>& control_backend)
      : backend_([&]() {
          // Extract observation types and variables
          const auto type_names = obs.getTypeNames();

          // Collect unique observation variable names across all types
          std::vector<std::string> required_obs_vars;
          auto add_unique = [&](const std::string& v) {
            if (std::find(required_obs_vars.begin(), required_obs_vars.end(),
                          v) == required_obs_vars.end()) {
              required_obs_vars.push_back(v);
            }
          };
          for (const auto& tname : type_names) {
            for (const auto& v : obs.getVariableNames(tname)) add_unique(v);
          }

          // Derive required state vars from observation variables (same names)
          std::vector<std::string> required_state_vars = required_obs_vars;

          // Determine operator family: if single type enabled, use it; else
          // empty
          std::string operator_family =
              (type_names.size() == 1) ? type_names[0] : "";

          // Build augmented config with derived metadata
          framework::ConfigMap config_map;

          // Copy specific configuration keys that are commonly needed across
          // backends This preserves backend-specific configuration (WRFDA, GSI,
          // DART, etc.) while maintaining framework genericity
          try {
            if (config.HasKey("external_root")) {
              config_map["external_root"] = config.Get("external_root");
            }
          } catch (...) {
            // Ignore if key doesn't exist
          }

          try {
            if (config.HasKey("external_system")) {
              config_map["external_system"] = config.Get("external_system");
            }
          } catch (...) {
            // Ignore if key doesn't exist
          }

          try {
            if (config.HasKey("operator_family")) {
              config_map["operator_family"] = config.Get("operator_family");
            }
          } catch (...) {
            // Ignore if key doesn't exist
          }

          // Set derived values - handle both single and multiple operator
          // families
          if (!operator_family.empty()) {
            // If we have a single type, use it as the primary family
            config_map["operator_family"] =
                framework::ConfigValue(operator_family);
          }
          config_map["required_state_vars"] =
              framework::ConfigValue(required_state_vars);
          config_map["required_obs_vars"] =
              framework::ConfigValue(required_obs_vars);

          // Create and return backend with the augmented config map, passing
          // control backend if supported
          using CB = typename traits::BackendTraits<BackendTag>::ConfigBackend;
          CB cfg(config_map);
          using B = ObsOperatorBackend;
          if constexpr (requires { B{cfg, control_backend}; }) {
            return B(cfg, control_backend);
          } else {
            return B(cfg);
          }
        }()),
        control_backend_(&control_backend) {
    logger_.Info()
        << "ObsOperator constructed with observation-derived configuration";
  }

  /**
   * @brief Move constructor
   *
   * @details Transfers ownership of the backend from another observation
   * operator.
   *
   * @param other The observation operator to move from
   */
  ObsOperator(ObsOperator&& other) noexcept
      : backend_(std::move(other.backend_)),
        control_backend_(other.control_backend_) {}

  /**
   * @brief Move assignment operator
   *
   * @details Transfers ownership of the backend from another observation
   * operator.
   *
   * @param other The observation operator to move from
   * @return Reference to this observation operator after assignment
   */
  ObsOperator& operator=(ObsOperator&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
    }
    return *this;
  }

  /**
   * @brief Get const reference to the backend implementation
   *
   * @return Const reference to the backend implementation
   */
  const ObsOperatorBackend& backend() const { return backend_; }

  /**
   * @brief Check if the observation operator is initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return backend_.isInitialized(); }

  /**
   * @brief Apply forward observation operator: H(x)
   *
   * @details Maps state to observation space (y = H(x)) by transforming model
   * state variables into simulated observations that can be compared with
   * actual observations.
   *
   * @tparam StateBackend The state backend type
   * @tparam ObsBackend The observation backend type
   * @param state Model state to transform
   * @param obs Output observation to store the result
   * @throws std::runtime_error If the observation operator is not initialized
   */
  template <typename StateBackend, typename ObsBackend>
  std::vector<double> apply(const State<StateBackend>& state,
                            const Observation<ObsBackend>& obs) const {
    return backend_.apply(state.backend(), obs.backend());
  }

  /**
   * @brief Get the state variables required by this observation operator
   *
   * @return Const reference to vector of required state variable names
   */
  const std::vector<std::string>& getRequiredStateVars() const {
    return backend_.getRequiredStateVars();
  }

  /**
   * @brief Get the observation variables required by this observation operator
   *
   * @return Const reference to vector of required observation variable names
   */
  const std::vector<std::string>& getRequiredObsVars() const {
    return backend_.getRequiredObsVars();
  }

  /**
   * @brief Apply tangent linear observation operator in control space: H' U v
   *
   * @details Applies the tangent linear of the observation operator to a
   * control variable. This method delegates to the backend's control-space
   * method, which handles the transformation from control space to state space
   * internally.
   *
   * @tparam StateBackend The state backend type
   * @tparam ObsBackend The observation backend type
   * @param control Control variable to transform
   * @param reference_state Reference state around which to linearize
   * @param obs Reference observations for context
   * @return Vector containing the transformed result in observation space
   * @throws std::runtime_error If the observation operator is not initialized
   */
  template <typename StateBackend, typename ObsBackend>
  std::vector<double> applyTangentLinear(
      const ControlVariable<StateBackend>& control,
      const State<StateBackend>& reference_state,
      const Observation<ObsBackend>& obs) const {
    // Transform control to state: δx = U v
    auto increment = control_backend_->createIncrement(control.geometry());
    control_backend_->controlToIncrement(control, increment);

    // Apply tangent linear: H'(δx)
    return backend_.applyTangentLinear(
        increment.backend(), reference_state.backend(), obs.backend());
  }

  /**
   * @brief Apply adjoint observation operator in control space: U^T H^T delta_y
   *
   * @details Applies the adjoint of the observation operator to an observation
   * increment and returns the result in control space. This method handles the
   * transformation from state space to control space internally.
   *
   * Steps:
   * 1. Apply adjoint: ∇_δx = H^T(delta_y)
   * 2. Transform to control: ∇_v = U^T ∇_δx
   * 3. Return result in control space
   *
   * @param obs_increment Observation space increment
   * @param reference_state Reference state around which the adjoint is computed
   * @param obs Observations to determine grid coordinates for adjoint mapping
   * @return Control variable containing the adjoint transformation result
   * @throws std::runtime_error If the observation operator is not initialized
   */
  ControlVariable<BackendTag> applyAdjoint(
      const std::vector<double>& obs_increment,
      const State<BackendTag>& reference_state,
      const Observation<BackendTag>& obs) const;

  template <typename StateBackend, typename ObsBackend>
  bool hasNativeIncrementalCost(const State<StateBackend>& reference_state,
                                const Observation<ObsBackend>& obs) const {
    if constexpr (requires {
                    backend_.computeObservationCost(reference_state.backend(),
                                                    obs.backend());
                  }) {
      (void)reference_state;
      (void)obs;
      return true;
    } else {
      (void)reference_state;
      (void)obs;
      return false;
    }
  }

  template <typename StateBackend, typename ObsBackend>
  std::optional<double> computeObservationCost(
      const State<StateBackend>& reference_state,
      const Observation<ObsBackend>& obs) const {
    if constexpr (requires {
                    backend_.computeObservationCost(reference_state.backend(),
                                                    obs.backend());
                  }) {
      return backend_.computeObservationCost(reference_state.backend(),
                                             obs.backend());
    } else {
      (void)reference_state;
      (void)obs;
      return std::nullopt;
    }
  }

  /**
   * @brief Apply pure adjoint observation operator in control space: U^T H'^T
   * delta_y
   *
   * @details This method delegates to the backend's control-space adjoint
   * method, which handles the transformation between control and state space
   * internally. For now, this is equivalent to applyAdjoint since backends
   * handle both weighted and pure adjoints through the same interface.
   */
  template <typename StateBackend, typename ObsBackend>
  ControlVariable<BackendTag> applyAdjointPure(
      const std::vector<double>& obs_increment,
      const State<StateBackend>& reference_state,
      const Observation<ObsBackend>& obs) const {
    // Apply pure adjoint in state space if available; otherwise fallback
    auto state_gradient = control_backend_->createIncrement(
        reference_state.geometry()->backend());
    state_gradient.zero();

    if constexpr (requires {
                    backend_.applyAdjointPure(obs_increment,
                                              reference_state.backend(),
                                              obs.backend());
                  }) {
      auto backend_increment = backend_.applyAdjointPure(
          obs_increment, reference_state.backend(), obs.backend());
      state_gradient.backend() = std::move(backend_increment);
    } else {
      backend_.applyAdjoint(obs_increment, reference_state.backend(),
                            state_gradient.backend(), obs.backend());
    }

    // Transform to control space: ∇_v = U^T ∇_δx
    auto control_gradient = control_backend_->createControlVariable(
        reference_state.geometry()->backend());
    control_gradient.zero();
    control_backend_->incrementAdjointToControl(state_gradient,
                                                control_gradient);

    return control_gradient;
  }

  /**
   * @brief Check if tangent linear and adjoint operators are available
   *
   * @return True if tangent linear/adjoint operators are supported
   */
  bool supportsLinearization() const {
    return backend_.supportsLinearization();
  }

  /**
   * @brief Check if the observation operator is linear
   *
   * @details Linear observation operators have the property that H(x+dx) = H(x)
   * + H(dx). For linear operators, the tangent linear and forward operators are
   * identical.
   *
   * @return True if the observation operator is linear
   */
  bool isLinear() const { return backend_.isLinear(); }

 private:
  ObsOperatorBackend backend_; /**< Backend implementation */
  const ControlVariableBackend<BackendTag>*
      control_backend_; /**< Control variable backend for transformations */
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

// Implementation of applyAdjoint for ObsOperator
template <typename BackendTag>
  requires ObsOperatorBackendType<BackendTag>
ControlVariable<BackendTag> ObsOperator<BackendTag>::applyAdjoint(
    const std::vector<double>& obs_adjoint,
    const State<BackendTag>& reference_state,
    const Observation<BackendTag>& obs) const {
  // If backend provides a control-space adjoint that returns a control vector,
  // delegate to it and wrap the backend result.
  if constexpr (requires {
                  backend_.applyAdjoint(obs_adjoint, reference_state.backend(),
                                        obs.backend());
                }) {
    auto control_gradient_backend = backend_.applyAdjoint(
        obs_adjoint, reference_state.backend(), obs.backend());
    // Wrap backend control result into ControlVariable adapter
    auto control_gradient = control_backend_->createControlVariable(
        reference_state.geometry()->backend());
    control_gradient.zero();
    control_gradient.backend() = std::move(control_gradient_backend);
    return control_gradient;
  } else {
    // Fallback: Apply adjoint in state space and transform to control space
    auto state_gradient = control_backend_->createIncrement(
        reference_state.geometry()->backend());
    state_gradient.zero();
    backend_.applyAdjoint(obs_adjoint, reference_state.backend(),
                          state_gradient.backend(), obs.backend());

    auto control_gradient = control_backend_->createControlVariable(
        reference_state.geometry()->backend());
    control_gradient.zero();
    control_backend_->incrementAdjointToControl(state_gradient,
                                                control_gradient);
    return control_gradient;
  }
}

}  // namespace metada::framework