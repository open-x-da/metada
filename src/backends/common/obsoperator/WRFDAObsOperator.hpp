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

#include "Location.hpp"
#include "IdentityObsOperator.hpp"

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
    delegate_ = std::make_unique<IdentityObsOperator<StateBackend, ObsBackend>>(config);

    initialized_ = true;
  }

  bool isInitialized() const { return initialized_; }

  std::vector<double> apply(const StateBackend& state, const ObsBackend& obs) const {
    ensureInitialized();
    return delegate_->apply(state, obs);
  }

  const std::vector<std::string>& getRequiredStateVars() const {
    return required_state_vars_;
  }

  const std::vector<std::string>& getRequiredObsVars() const { return required_obs_vars_; }

  std::vector<double> applyTangentLinear(const StateBackend& state_increment,
                                         const StateBackend& reference_state,
                                         const ObsBackend& obs) const {
    ensureInitialized();
    return delegate_->applyTangentLinear(state_increment, reference_state, obs);
  }

  void applyAdjoint(const std::vector<double>& obs_increment,
                    const StateBackend& reference_state,
                    StateBackend& result_state,
                    const ObsBackend& obs) const {
    ensureInitialized();
    delegate_->applyAdjoint(obs_increment, reference_state, result_state, obs);
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


