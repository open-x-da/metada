#pragma once

#include <string>
#include <vector>

#include "WRFDAObsOperator_c_api.h"

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
  explicit WRFObsOperator(const ConfigBackend& config)
      : initialized_(false), iv_type_data_(nullptr) {
    initialize(config);
  }

  /**
   * @brief Move constructor
   * @param other WRFObsOperator to move from
   */
  WRFObsOperator(WRFObsOperator&& other) noexcept
      : initialized_(other.initialized_),
        iv_type_data_(other.iv_type_data_),
        operator_families_(std::move(other.operator_families_)) {
    other.initialized_ = false;
    other.iv_type_data_ = nullptr;
  }

  /**
   * @brief Move assignment operator
   * @param other WRFObsOperator to move from
   * @return Reference to this object
   */
  WRFObsOperator& operator=(WRFObsOperator&& other) noexcept {
    if (this != &other) {
      initialized_ = other.initialized_;
      iv_type_data_ = other.iv_type_data_;
      operator_families_ = std::move(other.operator_families_);
      other.initialized_ = false;
      other.iv_type_data_ = nullptr;
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
  void initialize(const ConfigBackend& config) {
    if (initialized_) {
      throw std::runtime_error("WRFObsOperator already initialized");
    }

    // Extract operator family configuration
    try {
      auto family_value = config.Get("operator_family");
      if (family_value.isVectorString()) {
        operator_families_ = family_value.asVectorString();
      } else if (family_value.isString()) {
        operator_families_ = {family_value.asString()};
      } else {
        operator_families_ = {"synop"};  // Default to synop
      }
    } catch (...) {
      operator_families_ = {"synop"};  // Default fallback
    }

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
   * @return Vector of simulated observations
   *
   * @details Maps model state to observation space by applying the WRFDA
   * observation operator. This computes H(x) where H is the observation
   * operator and x is the model state.
   */
  std::vector<double> apply([[maybe_unused]] const StateBackend& state,
                            const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }

    // Ensure IV type is constructed for this observation set
    ensureIVTypeConstructed(obs);

    // Get observation structure from observation backend
    void* ob_ptr = obs.getYTypeData();
    void* iv_ptr = iv_type_data_;

    // Call WRFDA to compute innovations directly
    const int it = 1;
    int rc = wrfda_get_innov_vector(&it, ob_ptr, iv_ptr);
    if (rc != 0) {
      throw std::runtime_error(
          "WRFDA innovation computation failed with code " +
          std::to_string(rc));
    }

    // Extract innovations from iv structure
    std::string family =
        operator_families_.empty() ? "synop" : operator_families_[0];

    int num_innovations = 0;
    rc = wrfda_count_innovations(family.c_str(), &num_innovations);
    if (rc != 0 || num_innovations <= 0) {
      return std::vector<double>(obs.size(), 0.0);
    }

    std::vector<double> innovations(num_innovations);
    rc = wrfda_extract_innovations(family.c_str(), innovations.data(),
                                   &num_innovations);
    if (rc != 0) {
      throw std::runtime_error("Failed to extract innovations with code " +
                               std::to_string(rc));
    }

    return innovations;
  }

  /**
   * @brief Apply the tangent linear observation operator: H'(δx)
   * @param state_increment State increment to transform
   * @param reference_state Reference state around which to linearize
   * @param obs Observations defining observation locations
   * @return Vector of observation increments
   *
   * @details Applies the tangent linear of the observation operator to a
   * state increment. This is used in variational data assimilation for
   * computing the linearized observation operator.
   */
  std::vector<double> applyTangentLinear(
      const StateBackend& state_increment,
      [[maybe_unused]] const StateBackend& reference_state,
      const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }

    // Ensure IV type is constructed for TL/AD checks
    ensureIVTypeConstructed(obs);

    // Extract state increment data
    const auto state_data = state_increment.getObsOperatorData();

    // Update WRFDA analysis increments directly
    wrfda_update_analysis_increments(state_data.u, state_data.v, state_data.t,
                                     state_data.q, state_data.psfc);

    // Apply tangent linear operator directly
    int rc = wrfda_xtoy_apply_grid();
    if (rc != 0) {
      throw std::runtime_error("WRFDA tangent linear failed with code " +
                               std::to_string(rc));
    }

    // Get output count and extract values directly
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
   */
  void applyAdjoint(const std::vector<double>& obs_increment,
                    [[maybe_unused]] const StateBackend& reference_state,
                    StateBackend& result_state, const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }

    // Ensure IV type is constructed for TL/AD checks
    ensureIVTypeConstructed(obs);

    // Extract grid dimensions and state data
    const int nx = static_cast<int>(result_state.geometry().x_dim());
    const int ny = static_cast<int>(result_state.geometry().y_dim());
    const int nz = static_cast<int>(result_state.geometry().z_dim());

    // Get state variable pointers (or create zero arrays)
    double* u = static_cast<double*>(result_state.getData("U"));
    double* v = static_cast<double*>(result_state.getData("V"));
    double* t = static_cast<double*>(result_state.getData("T"));
    double* q = static_cast<double*>(result_state.getData("QVAPOR"));
    double* psfc = static_cast<double*>(result_state.getData("PSFC"));

    // Create zero arrays if variables don't exist
    std::vector<double> u_zero, v_zero, t_zero, q_zero, psfc_zero;
    if (!u) {
      u_zero.assign(nx * ny * nz, 0.0);
      u = u_zero.data();
    }
    if (!v) {
      v_zero.assign(nx * ny * nz, 0.0);
      v = v_zero.data();
    }
    if (!t) {
      t_zero.assign(nx * ny * nz, 0.0);
      t = t_zero.data();
    }
    if (!q) {
      q_zero.assign(nx * ny * nz, 0.0);
      q = q_zero.data();
    }
    if (!psfc) {
      psfc_zero.assign(nx * ny, 0.0);
      psfc = psfc_zero.data();
    }

    // Call WRFDA adjoint directly
    int rc = wrfda_xtoy_adjoint_grid(nx, ny, nz, obs_increment.data(),
                                     static_cast<int>(obs_increment.size()), u,
                                     v, t, q, psfc);
    if (rc != 0) {
      throw std::runtime_error("WRFDA adjoint failed with code " +
                               std::to_string(rc));
    }
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

  /**
   * @brief Get configured operator families
   * @return Reference to vector of operator family names
   */
  const std::vector<std::string>& getOperatorFamilies() const {
    return operator_families_;
  }

 private:
  /**
   * @brief Ensure IV type is constructed for the given observation set
   * @param obs Observation backend containing observation data
   *
   * @details Constructs WRFDA iv_type structure once and caches it for reuse
   * in TL/AD checks. This ensures the IV structure is ready before any
   * operator methods are called.
   */
  void ensureIVTypeConstructed(const ObsBackend& obs) const {
    if (iv_type_data_ == nullptr) {
      constructIVType(obs);
    }
  }

  /**
   * @brief Construct IV type (innovation vector) from observation data
   * @param obs Observation backend containing observation data
   *
   * @details Constructs WRFDA iv_type structure by calling
   * wrfda_construct_iv_type from the WRFDA bridge. The IV type contains
   * innovation vectors (O-B) and is used by WRFDA for data assimilation
   * computations.
   */
  void constructIVType(const ObsBackend& obs) const {
    const auto obs_data = obs.getObsOperatorData();

    if (obs_data.num_obs == 0) {
      return;
    }

    // Determine the operator family (default to "synop" for surface
    // observations)
    std::string family = "synop";

    // Extract structured data for WRFDA iv_type construction
    std::vector<double> u_values, v_values, t_values, p_values, q_values;
    std::vector<double> u_errors, v_errors, t_errors, p_errors, q_errors;
    std::vector<double> lats, lons, levels;
    std::vector<int> u_qc, v_qc, t_qc, p_qc, q_qc;
    std::vector<int> u_available, v_available, t_available, p_available,
        q_available;

    // Prepare observation types as a flat string
    std::string obs_types_flat;
    for (size_t i = 0; i < obs_data.num_obs; ++i) {
      std::string obs_type =
          (i < obs_data.types.size()) ? obs_data.types[i] : "ADPSFC";
      obs_types_flat += obs_type + std::string(20 - obs_type.length(), '\0');
    }

    // Extract data from the hierarchical structure
    for (size_t station = 0; station < obs_data.num_obs; ++station) {
      lats.push_back(obs_data.lats[station]);
      lons.push_back(obs_data.lons[station]);

      // Process each level for this station
      for (const auto& level : obs_data.levels[station]) {
        levels.push_back(level.level);

        // Use WRFDA standard missing value (-888888.0) for unavailable
        // observations
        const double missing_r = -888888.0;
        u_values.push_back(level.u.available ? level.u.value : missing_r);
        v_values.push_back(level.v.available ? level.v.value : missing_r);
        t_values.push_back(level.t.available ? level.t.value : missing_r);
        p_values.push_back(level.p.available ? level.p.value : missing_r);
        q_values.push_back(level.q.available ? level.q.value : missing_r);

        u_errors.push_back(level.u.available ? level.u.error : 1.0);
        v_errors.push_back(level.v.available ? level.v.error : 1.0);
        t_errors.push_back(level.t.available ? level.t.error : 1.0);
        p_errors.push_back(level.p.available ? level.p.error : 1.0);
        q_errors.push_back(level.q.available ? level.q.error : 1.0);

        const int missing_data = -88;
        u_qc.push_back(level.u.available ? level.u.qc : missing_data);
        v_qc.push_back(level.v.available ? level.v.qc : missing_data);
        t_qc.push_back(level.t.available ? level.t.qc : missing_data);
        p_qc.push_back(level.p.available ? level.p.qc : missing_data);
        q_qc.push_back(level.q.available ? level.q.qc : missing_data);

        // Extract availability flags
        u_available.push_back(level.u.available ? 1 : 0);
        v_available.push_back(level.v.available ? 1 : 0);
        t_available.push_back(level.t.available ? 1 : 0);
        p_available.push_back(level.p.available ? 1 : 0);
        q_available.push_back(level.q.available ? 1 : 0);
      }
    }

    // Call the WRFDA bridge function to construct iv_type
    int num_obs_int = static_cast<int>(obs_data.num_obs);
    int num_levels_int =
        static_cast<int>(u_values.size());  // Total number of level records

    void* iv_ptr = wrfda_construct_iv_type(
        &num_obs_int, &num_levels_int, const_cast<double*>(u_values.data()),
        const_cast<double*>(v_values.data()),
        const_cast<double*>(t_values.data()),
        const_cast<double*>(p_values.data()),
        const_cast<double*>(q_values.data()),
        const_cast<double*>(u_errors.data()),
        const_cast<double*>(v_errors.data()),
        const_cast<double*>(t_errors.data()),
        const_cast<double*>(p_errors.data()),
        const_cast<double*>(q_errors.data()), const_cast<int*>(u_qc.data()),
        const_cast<int*>(v_qc.data()), const_cast<int*>(t_qc.data()),
        const_cast<int*>(p_qc.data()), const_cast<int*>(q_qc.data()),
        const_cast<int*>(u_available.data()),
        const_cast<int*>(v_available.data()),
        const_cast<int*>(t_available.data()),
        const_cast<int*>(p_available.data()),
        const_cast<int*>(q_available.data()), const_cast<double*>(lats.data()),
        const_cast<double*>(lons.data()), const_cast<double*>(levels.data()),
        const_cast<double*>(obs_data.elevations.data()),
        const_cast<char*>(obs_types_flat.c_str()),
        const_cast<char*>(family.c_str()));

    // Store the pointer (memory is managed by Fortran persistent structures)
    if (iv_ptr) {
      iv_type_data_ = iv_ptr;
    }
  }

  /// Initialization status flag
  bool initialized_;

  /// IV type data (non-owning pointer to Fortran-managed memory)
  mutable void* iv_type_data_;

  /// Configured operator families
  std::vector<std::string> operator_families_;
};

}  // namespace metada::backends::wrf
