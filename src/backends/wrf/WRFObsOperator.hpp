#pragma once

#include <string>
#include <vector>

#include "WRFDAObsOperator_c_api.h"

// Forward declarations for WRFDA structure accessors
extern "C" {
void* wrfda_get_head_grid_ptr_();
void* wrfda_get_iv_ptr(void);
void* wrfda_get_y_ptr(void);
void* wrfda_get_y_tl_ptr(void);
}

namespace metada::backends::wrf {

/**
 * @brief WRF observation operator backend with direct WRFDA integration
 *
 * This class serves as the WRF observation operator backend, providing a clean
 * interface between the framework adapter and WRFDA C API calls. It follows
 * the same pattern as other WRF backend classes (WRFObservation, WRFState,
 * WRFGeometry).
 *
 * Key responsibilities:
 * 1. Implement framework ObsOperatorBackendImpl concept requirements
 * 2. Call WRFDA C API functions directly for observation operators
 * 3. Provide WRF-specific configuration and initialization
 * 4. Handle linearization and adjoint operations
 * 5. Bridge between framework types and WRFDA C API
 *
 * @tparam StateBackend The WRF state backend type
 * @tparam ObsBackend The WRF observation backend type
 * @tparam IncrementBackend The WRF increment backend type
 */
template <typename StateBackend, typename ObsBackend, typename IncrementBackend>
class WRFObsOperator {
 public:
  // Delete default constructor and copying (required by framework concept)
  WRFObsOperator() = delete;
  WRFObsOperator(const WRFObsOperator&) = delete;
  WRFObsOperator& operator=(const WRFObsOperator&) = delete;

  /**
   * @brief Constructor with configuration
   * @param config Configuration object containing operator settings
   *
   * @details Initializes the WRF observation operator with the provided
   * configuration and calls WRFDA initialization directly.
   */
  template <typename ConfigBackend>
  explicit WRFObsOperator(const ConfigBackend& config) : initialized_(false) {
    initialize(config);
  }

  /**
   * @brief Move constructor
   * @param other WRFObsOperator to move from
   */
  WRFObsOperator(WRFObsOperator&& other) noexcept
      : initialized_(other.initialized_) {
    other.initialized_ = false;
  }

  /**
   * @brief Move assignment operator
   * @param other WRFObsOperator to move from
   * @return Reference to this object
   */
  WRFObsOperator& operator=(WRFObsOperator&& other) noexcept {
    if (this != &other) {
      initialized_ = other.initialized_;
      other.initialized_ = false;
    }
    return *this;
  }

  /**
   * @brief Initialize the observation operator
   * @param config Configuration object
   *
   * @details Initializes WRFDA directly via C API calls.
   */
  template <typename ConfigBackend>
  void initialize(const ConfigBackend& /* config */) {
    if (initialized_) {
      throw std::runtime_error("WRFObsOperator already initialized");
    }

    // Note: operator_family configuration is no longer needed
    // The refactored extract functions automatically handle all families

    // Initialize WRFDA 3D-Var system directly
    initialize_wrfda_3dvar();
    initialized_ = true;
  }

  /**
   * @brief Check if the operator is initialized
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Apply the forward observation operator: H(x)
   * @param state Model state to transform
   * @param obs Observations defining observation locations
   * @return Vector of simulated observations H(x)
   *
   * @details Maps model state to observation space by applying the WRFDA
   * observation operator. This computes H(x) where H is the observation
   * operator and x is the model state.
   *
   * Implementation: WRFDA computes innovations d = y_obs - H(x), so we extract:
   * - innovations d from iv structure
   * - observations y_obs from y structure
   * - Then compute H(x) = y_obs - d
   */
  std::vector<double> apply(const StateBackend& state,
                            const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }

    // Transfer WRF fields to background state (grid%xb) using WRFDA's standard
    // workflow This is called in da_setup_firstguess BEFORE da_get_innov_vector
    int rc = wrfda_transfer_wrftoxb();
    if (rc != 0) {
      throw std::runtime_error("WRFDA da_transfer_wrftoxb failed with code " +
                               std::to_string(rc));
    }

    // Get observation structures from observation backend
    // iv_type and y_type are already allocated and populated by WRFObservation
    void* ob_ptr = obs.getYTypeData();
    void* iv_ptr = obs.getIVTypeData();

    // Get grid pointer from geometry (allocated during WRFGeometry
    // construction)
    void* grid_ptr = state.geometry().getGridPtr();

    // Call WRFDA to compute innovations: d = y_obs - H(x)
    const int it = 1;
    rc = wrfda_get_innov_vector(&it, ob_ptr, iv_ptr, grid_ptr);
    if (rc != 0) {
      throw std::runtime_error(
          "WRFDA innovation computation failed with code " +
          std::to_string(rc));
    }

    (void)state;  // Silence unused parameter warning

    // Count innovations across ALL observation types
    int num_innovations = 0;
    rc = wrfda_count_innovations(iv_ptr, &num_innovations);
    if (rc != 0) {
      throw std::runtime_error("Failed to count innovations with code " +
                               std::to_string(rc));
    }
    if (num_innovations <= 0) {
      throw std::runtime_error(
          "No innovations computed by WRFDA. Check that observations are "
          "within "
          "the domain and that use_* flags are enabled in namelist.input");
    }

    // Extract innovations for all families at once
    std::vector<double> innovations(num_innovations);
    int actual_innovations = 0;
    rc = wrfda_extract_innovations(iv_ptr, innovations.data(),
                                   &actual_innovations);
    if (rc != 0) {
      throw std::runtime_error("Failed to extract innovations with code " +
                               std::to_string(rc));
    }

    if (actual_innovations != num_innovations) {
      throw std::runtime_error(
          "Mismatch between counted and extracted innovations: " +
          std::to_string(num_innovations) + " vs " +
          std::to_string(actual_innovations));
    }

    // Extract observations for all families at once
    std::vector<double> observations(num_innovations);
    int actual_observations = 0;
    rc = wrfda_extract_observations(iv_ptr, ob_ptr, observations.data(),
                                    &actual_observations);
    if (rc != 0) {
      throw std::runtime_error("Failed to extract observations with code " +
                               std::to_string(rc));
    }

    if (actual_observations != num_innovations) {
      throw std::runtime_error(
          "Mismatch between number of observations and innovations: " +
          std::to_string(actual_observations) + " vs " +
          std::to_string(num_innovations));
    }

    // Compute H(x) = y_obs - d for all observations
    std::vector<double> all_simulated_obs;
    all_simulated_obs.reserve(num_innovations);
    for (int i = 0; i < num_innovations; ++i) {
      all_simulated_obs.emplace_back(observations[i] - innovations[i]);
    }

    return all_simulated_obs;
  }

  /**
   * @brief Apply the tangent linear observation operator: H'(δx)
   * @param state_increment State increment to transform
   * @param reference_state Reference state around which to linearize
   * @param obs Observations defining observation locations
   * @return Vector of observation increments H'(xb)·δx
   *
   * @details Applies the tangent linear of the observation operator to a
   * state increment. This is used in variational data assimilation for
   * computing the linearized observation operator.
   *
   * Implementation: Since apply() now returns H(x), the tangent linear
   * returns H'(xb)·δx directly without sign changes.
   */
  std::vector<double> applyTangentLinear(const IncrementBackend& increment,
                                         const StateBackend& reference_state,
                                         const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }

    // Sync increment's internal CV storage to grid%xa before WRFDA operation
    // This writes the increment data to WRFDA's working memory
    increment.syncToGrid();

    // Note: reference_state parameter is not used because WRF grid already
    // contains the state. Background state (grid%xb) is already set up during
    // the nonlinear apply() call in cost function initialization.
    (void)reference_state;
    (void)obs;

    // Get WRFDA structures directly (grid, iv, y are allocated and managed by
    // WRFDA)
    void* iv_ptr = obs.getIVTypeData();

    // Apply tangent linear operator: H'(xb)·xa
    // Use WRFDA-side structures directly (iv and y are already allocated)
    int rc = wrfda_xtoy_apply_grid();
    if (rc != 0) {
      throw std::runtime_error("WRFDA tangent linear failed with code " +
                               std::to_string(rc));
    }

    // Extract H'·δx from temporary TL y structure (wrfda_y_tl)
    // This preserves wrfda_ob (observation values) from corruption
    int num_observations = 0;
    rc = wrfda_count_innovations(iv_ptr, &num_observations);
    if (rc != 0 || num_observations <= 0) {
      return std::vector<double>();
    }

    // Get TL y pointer (separate from observation y to avoid corruption)
    void* y_tl_ptr = wrfda_get_y_tl_ptr();
    if (!y_tl_ptr) {
      throw std::runtime_error("Tangent linear y_type not allocated");
    }

    std::vector<double> out_y(num_observations);
    int actual_count = 0;
    rc = wrfda_extract_observations(iv_ptr, y_tl_ptr, out_y.data(),
                                    &actual_count);
    if (rc != 0) {
      throw std::runtime_error(
          "Failed to extract tangent linear output with code " +
          std::to_string(rc));
    }

    // WRFDA's da_transform_xtoy computes +H'(xb)·δx and stores in wrfda_y_tl
    // wrfda_ob remains unchanged, preserving observation values
    return out_y;
  }

  /**
   * @brief Apply PURE adjoint for TL/AD testing (without R^{-1} weighting)
   *
   * @param obs_adjoint Raw observation space vector dy
   * @param reference_state Reference state (not used)
   * @param obs Observations
   * @return Increment containing H'^T · dy
   *
   * @details This method applies the pure adjoint operator H'^T directly
   * without R^{-1} weighting or residual computation. It's used specifically
   * for TL/AD checks where we need: <H'dx, dy> = <dx, H'^T dy>
   *
   * This is separate from the full applyAdjoint workflow used in cost function
   * evaluation.
   */
  template <typename StateBackendArg, typename ObsBackendArg>
  IncrementBackend applyAdjointPure(const std::vector<double>& obs_adjoint,
                                    const StateBackendArg& reference_state,
                                    const ObsBackendArg& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }

    (void)reference_state;

    // Get WRFDA structures
    void* grid_ptr = wrfda_get_head_grid_ptr_();
    void* iv_ptr = obs.getIVTypeData();

    if (!grid_ptr || !iv_ptr) {
      throw std::runtime_error("WRFDA structures not available");
    }

    // Create result increment from geometry
    auto x_adjoint = IncrementBackend(reference_state.geometry());
    x_adjoint.zero();

    // Apply pure adjoint: H'^T · dy → grid%xa
    int rc = wrfda_xtoy_adjoint_pure(grid_ptr, iv_ptr, obs_adjoint.data(),
                                     static_cast<int>(obs_adjoint.size()));
    if (rc != 0) {
      throw std::runtime_error("WRFDA pure adjoint failed with code " +
                               std::to_string(rc));
    }

    // Sync grid%xa back to increment's internal CV storage
    x_adjoint.syncFromGrid();

    return x_adjoint;
  }

  /**
   * @brief Apply adjoint using WRFDA's complete proven workflow
   *
   * @param obs_adjoint IGNORED for WRF backend (see details)
   * @param reference_state Reference state
   * @param x_adjoint Increment to receive adjoint gradient
   * @param obs Observations
   *
   * @details For WRF backend, this uses WRFDA's complete proven workflow:
   *
   * **WRFDA Workflow (replaces framework's manual residual + R^{-1}):**
   * 1. `da_calculate_residual`: re = (O-B) - H'(δx)
   * 2. `da_calculate_grady`: jo_grad_y = -R^{-1} · re
   * 3. `da_transform_xtoy_adj`: H^T · jo_grad_y → grid%xa
   *
   * @note obs_adjoint parameter is **IGNORED**. WRFDA computes weighted
   * residuals internally. Framework's manual computation in line 340-363 of
   * IncrementalCostFunction.hpp is redundant for WRF but kept for other
   * backends.
   */
  void applyAdjoint(const std::vector<double>& obs_adjoint,
                    const StateBackend& reference_state,
                    IncrementBackend& x_adjoint, const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }

    // Note: reference_state and obs parameters are not used because WRF grid
    // already contains the state, and WRFDA manages observation structures
    // internally. Background state (grid%xb) is already set up during the
    // nonlinear apply() call.
    (void)reference_state;
    (void)obs_adjoint;  // Using WRFDA's internal workflow instead
    (void)obs;

    // Get WRFDA structures directly (grid, iv, y are managed by WRFDA)
    void* grid_ptr = wrfda_get_head_grid_ptr_();
    void* iv_ptr = wrfda_get_iv_ptr();

    // CRITICAL: Use wrfda_y_tl (TL output), not wrfda_ob (observation values)
    // The adjoint needs re = (O-B) - H'(δx), where H'(δx) is in wrfda_y_tl
    void* y_tl_ptr = wrfda_get_y_tl_ptr();

    if (!grid_ptr || !iv_ptr || !y_tl_ptr) {
      throw std::runtime_error(
          "WRFDA structures not available. Ensure observations are read and "
          "TL operator has been called.");
    }

    // Use WRFDA's proven workflow instead of manual residual computation:
    // 1. da_calculate_residual: re = (O-B) - H'(δx)
    //    where H'(δx) is in wrfda_y_tl from applyTangentLinear
    // 2. da_calculate_grady: jo_grad_y = -R^{-1} · re
    // 3. Jo = -0.5 * Σ(re · jo_grad_y) = 0.5 * re^T · R^{-1} · re
    // This replaces MetaDA's manual residual and R^{-1} application
    double jo_cost = 0.0;
    int rc = wrfda_compute_weighted_residual(iv_ptr, y_tl_ptr, &jo_cost);
    if (rc != 0) {
      throw std::runtime_error(
          "WRFDA compute weighted residual failed with code " +
          std::to_string(rc));
    }

    // Store the observation cost for retrieval by cost function
    last_observation_cost_ = jo_cost;

    // Zero the increment first (using WRFDA's da_zero_x for ALL fields)
    x_adjoint.zero();

    // Call WRFDA adjoint which uses the persistent_jo_grad_y computed above
    // This computes H^T · jo_grad_y → grid%xa for ALL observation types
    rc = wrfda_xtoy_adjoint_grid(grid_ptr, iv_ptr);
    if (rc != 0) {
      throw std::runtime_error("WRFDA adjoint failed with code " +
                               std::to_string(rc));
    }

    // Sync grid%xa back to increment's internal CV storage
    // WRFDA's adjoint has populated grid%xa with ALL fields (35+ fields)
    // using the complete workflow: residual → R^{-1} → adjoint H^T
    x_adjoint.syncFromGrid();

    // WRFDA computes the mathematically correct observation gradient:
    // H^T · jo_grad_y where jo_grad_y = -R^{-1} · [d - H'(δx)]
    // This gives: -H^T R^{-1} [d - H'(δx)]
    // The framework will ADD this to the background gradient (see
    // IncrementalCostFunction) No negation needed - WRFDA's sign convention
    // matches the mathematics!
  }

  /**
   * @brief Compute observation cost function without full adjoint
   * @return Observation cost Jo = 0.5 * re^T * R^{-1} * re
   *
   * @details This method computes only the observation cost using WRFDA's
   * proven workflow (da_calculate_residual + da_calculate_grady). It's used
   * by the cost function evaluation before the adjoint is computed. The
   * computed jo_grad_y is stored for later use by applyAdjoint.
   *
   * @note This requires that applyTangentLinear was called first to populate
   * wrfda_y_tl with H'(δx).
   */
  template <typename StateBackendArg, typename ObsBackendArg>
  double computeObservationCost(const StateBackendArg& /*reference_state*/,
                                const ObsBackendArg& /*obs*/) const {
    // Get WRFDA structures (y_tl contains H'(δx) from applyTangentLinear)
    void* iv_ptr = wrfda_get_iv_ptr();
    void* y_tl_ptr = wrfda_get_y_tl_ptr();

    if (!iv_ptr || !y_tl_ptr) {
      throw std::runtime_error(
          "WRFDA structures not available. Ensure applyTangentLinear was "
          "called first.");
    }

    // Compute ONLY the observation cost (without overwriting
    // persistent_jo_grad_y) This is critical: gradient checks evaluate costs
    // AFTER computing the adjoint, so we must not overwrite
    // persistent_jo_grad_y which is needed by the adjoint
    double jo_cost = 0.0;
    int rc = wrfda_compute_jo_only(iv_ptr, y_tl_ptr, &jo_cost);
    if (rc != 0) {
      throw std::runtime_error("WRFDA compute cost failed with code " +
                               std::to_string(rc));
    }

    // Store for later retrieval
    last_observation_cost_ = jo_cost;

    return jo_cost;
  }

  /**
   * @brief Get required state variables
   * @return Reference to vector of required state variable names
   */
  const std::vector<std::string>& getRequiredStateVars() const {
    static const std::vector<std::string> required_vars = {"U", "V", "T",
                                                           "QVAPOR", "PSFC"};
    return required_vars;
  }

  /**
   * @brief Get required observation variables
   * @return Reference to vector of required observation variable names
   */
  const std::vector<std::string>& getRequiredObsVars() const {
    static const std::vector<std::string> required_obs_vars = {"u", "v", "t",
                                                               "p", "q"};
    return required_obs_vars;
  }

  /**
   * @brief Check if linearization is supported
   * @return True if tangent linear and adjoint operators are available
   */
  bool supportsLinearization() const { return true; }

  /**
   * @brief Check if the operator is linear
   * @return True if the observation operator is linear
   *
   * @details Linear operators have the property that H(x+dx) = H(x) + H(dx).
   * For linear operators, the tangent linear and forward operators are
   * identical.
   */
  bool isLinear() const { return true; }

 private:
  /// Initialization status flag
  bool initialized_;

  /// Last computed observation cost (from wrfda_compute_weighted_residual)
  mutable double last_observation_cost_ = 0.0;

 public:
  /**
   * @brief Get the last computed observation cost
   * @return Observation cost from last applyAdjoint call
   *
   * @details Returns the observation cost Jo = 0.5 * re^T * R^{-1} * re
   * computed by wrfda_compute_weighted_residual during the last applyAdjoint
   * call. This follows WRFDA's proven formula: Jo = -0.5 * Σ(re · jo_grad_y)
   */
  double getLastObservationCost() const { return last_observation_cost_; }
};

}  // namespace metada::backends::wrf
