#pragma once

#include <string>
#include <vector>

#include "WRFDAObsOperator_c_api.h"

// Forward declaration for WRFDA grid accessor
extern "C" {
void* wrfda_get_head_grid_ptr_();
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
 */
template <typename StateBackend, typename ObsBackend>
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
      all_simulated_obs.push_back(observations[i] - innovations[i]);
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
  std::vector<double> applyTangentLinear(const StateBackend& state_increment,
                                         const StateBackend& reference_state,
                                         const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }

    // Note: reference_state parameter is not used because WRF grid already
    // contains the state
    (void)reference_state;

    // Transfer WRF fields to background state (grid%xb) using WRFDA's standard
    // workflow
    int rc = wrfda_transfer_wrftoxb();
    if (rc != 0) {
      throw std::runtime_error(
          "WRFDA da_transfer_wrftoxb failed in tangent linear with code " +
          std::to_string(rc));
    }

    // Get grid pointer from geometry (allocated during WRFGeometry
    // construction)
    void* grid_ptr = reference_state.geometry().getGridPtr();

    // Extract state increment data and update analysis increments (grid%xa)
    const auto state_data = state_increment.getObsOperatorData();
    wrfda_update_analysis_increments(state_data.u, state_data.v, state_data.t,
                                     state_data.q, state_data.psfc, grid_ptr);

    // Apply tangent linear operator: H'(xb)·xa
    // Get observation structures from observation backend (already allocated)
    void* ob_ptr = obs.getYTypeData();
    void* iv_ptr = obs.getIVTypeData();
    rc = wrfda_xtoy_apply_grid(grid_ptr, ob_ptr, iv_ptr);
    if (rc != 0) {
      throw std::runtime_error("WRFDA tangent linear failed with code " +
                               std::to_string(rc));
    }

    // Get output count and extract values
    int num_innovations = 0;
    rc = wrfda_get_tangent_linear_count(&num_innovations);
    if (rc != 0 || num_innovations <= 0) {
      return std::vector<double>();
    }

    std::vector<double> out_y(num_innovations);
    rc = wrfda_extract_tangent_linear_output(out_y.data(), &num_innovations);
    if (rc != 0) {
      throw std::runtime_error(
          "Failed to extract tangent linear output with code " +
          std::to_string(rc));
    }

    // WRFDA's da_transform_xtoy computes +H'(xb)·xa, which is exactly what we
    // need since apply() now returns H(x), so the derivative is +H'(x)
    return out_y;
  }

  /**
   * @brief Apply the adjoint observation operator: H^T(δy)
   * @param obs_increment Observation space increment
   * @param reference_state Reference state around which adjoint is computed
   * @param result_state State to accumulate adjoint result in
   * @param obs Observations defining observation locations
   *
   * @details Applies the adjoint of the observation operator to map from
   * observation space back to state space. This is essential for computing
   * gradients in variational data assimilation.
   *
   * Implementation: Since apply() now returns H(x), the adjoint returns
   * H'^T(xb) directly without sign changes.
   */
  void applyAdjoint(const std::vector<double>& obs_increment,
                    const StateBackend& reference_state,
                    StateBackend& result_state, const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }

    // Note: reference_state parameter is not used because WRF grid already
    // contains the state
    (void)reference_state;

    // Transfer WRF fields to background state (grid%xb) using WRFDA's standard
    // workflow
    int rc = wrfda_transfer_wrftoxb();
    if (rc != 0) {
      throw std::runtime_error(
          "WRFDA da_transfer_wrftoxb failed in adjoint with code " +
          std::to_string(rc));
    }

    // Set delta_y input for adjoint operator
    rc = wrfda_set_delta_y(obs_increment.data(),
                           static_cast<int>(obs_increment.size()));
    if (rc != 0) {
      throw std::runtime_error("WRFDA set delta_y failed with code " +
                               std::to_string(rc));
    }

    // Get observation structures from observation backend (already allocated)
    void* grid_ptr = wrfda_get_head_grid_ptr_();
    void* iv_ptr = obs.getIVTypeData();

    // Call WRFDA adjoint directly
    rc = wrfda_xtoy_adjoint_grid(grid_ptr, iv_ptr);
    if (rc != 0) {
      throw std::runtime_error("WRFDA adjoint failed with code " +
                               std::to_string(rc));
    }

    // Get state variable pointers
    double* u = static_cast<double*>(result_state.getData("U"));
    double* v = static_cast<double*>(result_state.getData("V"));
    double* t = static_cast<double*>(result_state.getData("T"));
    double* q = static_cast<double*>(result_state.getData("QVAPOR"));
    double* psfc = static_cast<double*>(result_state.getData("PSFC"));

    // Retrieve adjoint gradients from persistent arrays
    rc = wrfda_get_adjoint_gradients(u, v, t, q, psfc);
    if (rc != 0) {
      throw std::runtime_error("WRFDA get adjoint gradients failed with code " +
                               std::to_string(rc));
    }

    // WRFDA's adjoint computes +H'^T(xb), which is exactly what we need
    // since apply() now returns H(x), so the adjoint is +H'^T(x)
    // No sign changes needed!
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
};

}  // namespace metada::backends::wrf
