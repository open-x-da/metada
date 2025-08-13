/**
 * @file WRFDAObsOperator.hpp
 * @brief Generic WRFDA observation operator backend (header-only)
 * @ingroup backends
 *
 * @details
 * This operator provides a WRFDA-style observation operator that can be used by
 * any real model backend. It currently delegates to the identity/grid
 * interpolation operator for functionality. It reserves configuration keys for
 * integrating native WRFDA Fortran routines via a C/Fortran bridge in the
 * future, without coupling to any specific model backend.
 *
 * Config keys (optional):
 * - wrfda_root: Path to WRFDA sources (e.g., D:/linux/WRF/var/da)
 * - wrfda_operator_family: Observation operator family (e.g., "metar", "gpspw")
 * - required_state_vars: [array of strings]
 * - required_obs_vars: [array of strings]
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "IdentityObsOperator.hpp"
#include "Location.hpp"
#include "WRFDAObsOperator_c_api.h"

namespace metada::backends::common::obsoperator {

template <typename StateBackend, typename ObsBackend>
class WRFDAObsOperator {
 public:
  WRFDAObsOperator() = delete;
  WRFDAObsOperator(const WRFDAObsOperator&) = delete;
  WRFDAObsOperator& operator=(const WRFDAObsOperator&) = delete;

  template <typename ConfigBackend>
  explicit WRFDAObsOperator(const ConfigBackend& config) {
    initialize(config);
  }

  WRFDAObsOperator(WRFDAObsOperator&& other) noexcept
      : initialized_(other.initialized_),
        required_state_vars_(std::move(other.required_state_vars_)),
        required_obs_vars_(std::move(other.required_obs_vars_)),
        wrfda_root_(std::move(other.wrfda_root_)),
        operator_family_(std::move(other.operator_family_)),
        delegate_(std::move(other.delegate_)) {
    other.initialized_ = false;
  }

  WRFDAObsOperator& operator=(WRFDAObsOperator&& other) noexcept {
    if (this != &other) {
      initialized_ = other.initialized_;
      required_state_vars_ = std::move(other.required_state_vars_);
      required_obs_vars_ = std::move(other.required_obs_vars_);
      wrfda_root_ = std::move(other.wrfda_root_);
      operator_family_ = std::move(other.operator_family_);
      delegate_ = std::move(other.delegate_);
      other.initialized_ = false;
    }
    return *this;
  }

  template <typename ConfigBackend>
  void initialize(const ConfigBackend& config) {
    if (isInitialized()) {
      throw std::runtime_error("WRFDAObsOperator already initialized");
    }

    try {
      wrfda_root_ = config.Get("wrfda_root").asString();
    } catch (...) {
      wrfda_root_.clear();
    }

    try {
      operator_family_ = config.Get("wrfda_operator_family").asString();
    } catch (...) {
      operator_family_.clear();
    }

    try {
      required_state_vars_ = config.Get("required_state_vars").asVectorString();
    } catch (...) {
      required_state_vars_.clear();
    }

    try {
      required_obs_vars_ = config.Get("required_obs_vars").asVectorString();
    } catch (...) {
      required_obs_vars_.clear();
    }

    // Current implementation: delegate to identity/grid interpolation
    delegate_ =
        std::make_unique<IdentityObsOperator<StateBackend, ObsBackend>>(config);

    initialized_ = true;
  }

  bool isInitialized() const { return initialized_; }

  std::vector<double> apply(const StateBackend& state,
                            const ObsBackend& obs) const {
    ensureInitialized();
    // Collect observation coordinates (geographic)
    const size_t num_observations = obs.size();
    std::vector<double> obs_lats(num_observations);
    std::vector<double> obs_lons(num_observations);
    std::vector<double> obs_levels(num_observations, 0.0);
    for (size_t i = 0; i < num_observations; ++i) {
      const auto& loc = obs[i].location;
      if (loc.getCoordinateSystem() !=
          framework::CoordinateSystem::GEOGRAPHIC) {
        throw std::runtime_error(
            "WRFDAObsOperator requires geographic observation locations");
      }
      auto [lat, lon, level] = loc.getGeographicCoords();
      obs_lats[i] = lat;
      obs_lons[i] = lon;
      obs_levels[i] = level;
    }

    // Helper to obtain variable buffers or zero-filled fallbacks
    auto get_or_zero = [&](const std::string& var, std::vector<double>& storage,
                           const double*& ptr, int& nx, int& ny, int& nz) {
      ptr = static_cast<const double*>(state.getData(var));
      if (ptr == nullptr) {
        nx = static_cast<int>(state.geometry().x_dim());
        ny = static_cast<int>(state.geometry().y_dim());
        nz = static_cast<int>(state.geometry().z_dim());
        const size_t n = static_cast<size_t>(std::max(1, nx) * std::max(1, ny) *
                                             std::max(1, nz));
        storage.assign(n, 0.0);
        ptr = storage.data();
      } else {
        try {
          const auto& dims = state.getVariableDimensions(var);
          if (dims.size() == 3) {
            nz = static_cast<int>(dims[0]);
            ny = static_cast<int>(dims[1]);
            nx = static_cast<int>(dims[2]);
          } else if (dims.size() == 2) {
            ny = static_cast<int>(dims[0]);
            nx = static_cast<int>(dims[1]);
            nz = 1;
          }
        } catch (...) {
        }
      }
    };

    const double *u = nullptr, *v = nullptr, *t = nullptr, *q = nullptr,
                 *psfc = nullptr;
    std::vector<double> u_zero, v_zero, t_zero, q_zero, psfc_zero;
    int nx_u = 1, ny_u = 1, nz_u = 1;
    int nx_v = 1, ny_v = 1, nz_v = 1;
    int nx_t = 1, ny_t = 1, nz_t = 1;
    int nx_q = 1, ny_q = 1, nz_q = 1;
    int nx_ps = 1, ny_ps = 1, nz_ps = 1;
    get_or_zero("U", u_zero, u, nx_u, ny_u, nz_u);
    get_or_zero("V", v_zero, v, nx_v, ny_v, nz_v);
    get_or_zero("T", t_zero, t, nx_t, ny_t, nz_t);
    get_or_zero("Q", q_zero, q, nx_q, ny_q, nz_q);
    get_or_zero("PSFC", psfc_zero, psfc, nx_ps, ny_ps, nz_ps);

    // Canonical unstaggered grid dims
    const int nx = static_cast<int>(state.geometry().x_dim());
    const int ny = static_cast<int>(state.geometry().y_dim());
    const int nz = static_cast<int>(state.geometry().z_dim());

    // Geometry arrays
    const auto& gi = state.geometry().unstaggered_info();
    const double* lats2d = gi.latitude_2d.data();
    const double* lons2d = gi.longitude_2d.data();
    const double* levels_ptr =
        gi.vertical_coords.empty() ? nullptr : gi.vertical_coords.data();
    std::vector<double> dummy_levels(1, 0.0);
    if (levels_ptr == nullptr) levels_ptr = dummy_levels.data();

    std::vector<double> out_y(num_observations, 0.0);

    const char* family_cstr =
        (operator_family_.empty() ? "" : operator_family_.c_str());

    const int rc = wrfda_xtoy_apply_grid(
        family_cstr, nx, ny, nz, u, v, t, q, psfc, lats2d, lons2d, levels_ptr,
        static_cast<int>(num_observations), obs_lats.data(), obs_lons.data(),
        obs_levels.data(), out_y.data());
    if (rc != 0) {
      throw std::runtime_error("wrfda_xtoy_apply failed with code " +
                               std::to_string(rc));
    }

    return out_y;
  }

  const std::vector<std::string>& getRequiredStateVars() const {
    return required_state_vars_;
  }

  const std::vector<std::string>& getRequiredObsVars() const {
    return required_obs_vars_;
  }

  std::vector<double> applyTangentLinear(
      const StateBackend& state_increment,
      [[maybe_unused]] const StateBackend& reference_state,
      const ObsBackend& obs) const {
    // For linearized operator, use same array-based call on increment
    ensureInitialized();
    return apply(state_increment, obs);
  }

  void applyAdjoint(const std::vector<double>& obs_increment,
                    [[maybe_unused]] const StateBackend& reference_state,
                    StateBackend& result_state, const ObsBackend& obs) const {
    ensureInitialized();

    const size_t num_observations = obs.size();
    if (obs_increment.size() != num_observations) {
      throw std::runtime_error(
          "Adjoint increment size does not match number of observations");
    }

    std::vector<double> obs_lats(num_observations);
    std::vector<double> obs_lons(num_observations);
    std::vector<double> obs_levels(num_observations, 0.0);
    for (size_t i = 0; i < num_observations; ++i) {
      const auto& loc = obs[i].location;
      if (loc.getCoordinateSystem() !=
          framework::CoordinateSystem::GEOGRAPHIC) {
        throw std::runtime_error(
            "WRFDAObsOperator requires geographic observation locations");
      }
      auto [lat, lon, level] = loc.getGeographicCoords();
      obs_lats[i] = lat;
      obs_lons[i] = lon;
      obs_levels[i] = level;
    }

    // Prepare inout arrays for u,v,t,q,psfc accumulation
    auto get_or_zero_inout = [&](const std::string& var,
                                 std::vector<double>& storage, double*& ptr,
                                 int& nx, int& ny, int& nz) {
      ptr = static_cast<double*>(result_state.getData(var));
      if (ptr == nullptr) {
        nx = static_cast<int>(result_state.geometry().x_dim());
        ny = static_cast<int>(result_state.geometry().y_dim());
        nz = static_cast<int>(result_state.geometry().z_dim());
        const size_t n = static_cast<size_t>(std::max(1, nx) * std::max(1, ny) *
                                             std::max(1, nz));
        storage.assign(n, 0.0);
        ptr = storage.data();
      } else {
        try {
          const auto& dims = result_state.getVariableDimensions(var);
          if (dims.size() == 3) {
            nz = static_cast<int>(dims[0]);
            ny = static_cast<int>(dims[1]);
            nx = static_cast<int>(dims[2]);
          } else if (dims.size() == 2) {
            ny = static_cast<int>(dims[0]);
            nx = static_cast<int>(dims[1]);
            nz = 1;
          }
        } catch (...) {
        }
      }
    };

    double *u = nullptr, *v = nullptr, *t = nullptr, *q = nullptr,
           *psfc = nullptr;
    std::vector<double> u_zero, v_zero, t_zero, q_zero, psfc_zero;
    int nx_u = 1, ny_u = 1, nz_u = 1;
    int nx_v = 1, ny_v = 1, nz_v = 1;
    int nx_t = 1, ny_t = 1, nz_t = 1;
    int nx_q = 1, ny_q = 1, nz_q = 1;
    int nx_ps = 1, ny_ps = 1, nz_ps = 1;
    get_or_zero_inout("U", u_zero, u, nx_u, ny_u, nz_u);
    get_or_zero_inout("V", v_zero, v, nx_v, ny_v, nz_v);
    get_or_zero_inout("T", t_zero, t, nx_t, ny_t, nz_t);
    get_or_zero_inout("Q", q_zero, q, nx_q, ny_q, nz_q);
    get_or_zero_inout("PSFC", psfc_zero, psfc, nx_ps, ny_ps, nz_ps);

    const int nx = static_cast<int>(result_state.geometry().x_dim());
    const int ny = static_cast<int>(result_state.geometry().y_dim());
    const int nz = static_cast<int>(result_state.geometry().z_dim());

    const auto& gi = result_state.geometry().unstaggered_info();
    const double* lats2d = gi.latitude_2d.data();
    const double* lons2d = gi.longitude_2d.data();
    const double* levels_ptr =
        gi.vertical_coords.empty() ? nullptr : gi.vertical_coords.data();
    std::vector<double> dummy_levels(1, 0.0);
    if (levels_ptr == nullptr) levels_ptr = dummy_levels.data();

    const char* family_cstr =
        (operator_family_.empty() ? "" : operator_family_.c_str());

    const int rc = wrfda_xtoy_adjoint_grid(
        family_cstr, nx, ny, nz, obs_increment.data(), lats2d, lons2d,
        levels_ptr, static_cast<int>(num_observations), obs_lats.data(),
        obs_lons.data(), obs_levels.data(), u, v, t, q, psfc);
    if (rc != 0) {
      throw std::runtime_error("wrfda_xtoy_adjoint failed with code " +
                               std::to_string(rc));
    }
  }

  bool supportsLinearization() const { return true; }
  bool isLinear() const { return true; }

 private:
  void ensureInitialized() const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFDAObsOperator not initialized");
    }
  }

  bool initialized_ = false;
  std::vector<std::string> required_state_vars_;
  std::vector<std::string> required_obs_vars_;
  std::string wrfda_root_;
  std::string operator_family_;
  std::unique_ptr<IdentityObsOperator<StateBackend, ObsBackend>> delegate_;
};

}  // namespace metada::backends::common::obsoperator
