#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "ControlVariable.hpp"
#include "ControlVariableBackend.hpp"
#include "Increment.hpp"
#include "State.hpp"
#include "WRFControlVariable.hpp"
#include "WRFIncrement.hpp"
#include "wrfda/WRFDAControlBackendBridge.hpp"

namespace metada::backends::wrf {

/**
 * @brief WRFDA CV5 control variable backend (stub for future implementation)
 *
 * @details This backend will implement the WRFDA CV5 control variable
 * transformation, which includes:
 * - Balance relationships (geostrophic, hydrostatic, etc.)
 * - Vertical EOF decomposition
 * - Horizontal recursive filter transforms
 *
 * The transformation U relates control variables to increments:
 *   δx = U v
 *
 * where v is the CV5 control vector and δx is the state-space increment.
 *
 * This corresponds to WRFDA's da_transform_vtox (forward) and
 * da_transform_vtox_adj (adjoint) routines.
 *
 * @tparam BackendTag Backend tag satisfying framework concepts
 */
template <typename BackendTag>
class WRFDACV5ControlVariableBackend
    : public metada::framework::ControlVariableBackend<BackendTag> {
 public:
  using typename metada::framework::ControlVariableBackend<
      BackendTag>::IncrementType;
  using typename metada::framework::ControlVariableBackend<
      BackendTag>::ControlVariableType;
  using
      typename metada::framework::ControlVariableBackend<BackendTag>::StateType;
  using ConfigType = typename metada::framework::Config<BackendTag>;

  /**
   * @brief Constructor from configuration
   *
   * @param config Configuration containing WRFDA CV5 parameters
   */
  explicit WRFDACV5ControlVariableBackend(const ConfigType& config) {
    (void)config;  // Suppress unused warning for now
    // TODO: Initialize WRFDA CV5 parameters from config
    //   - Balance transformation parameters
    //   - Vertical EOF parameters
    //   - Horizontal recursive filter parameters
  }

  ~WRFDACV5ControlVariableBackend() override {
    try {
      finalizeBackend();
    } catch (...) {
      // Destructors must not throw; swallow errors during finalization.
    }
  }

  /**
   * @brief Get backend name
   *
   * @return "wrfda_cv5"
   */
  std::string name() const override { return "wrfda_cv5"; }

  /**
   * @brief Create a control variable from geometry
   *
   * @param geometry Geometry backend defining the domain
   * @return Newly created control variable
   */
  ControlVariableType createControlVariable(
      const typename ControlVariableType::GeometryBackendType& geometry)
      const override {
    auto control = ControlVariableType::createFromGeometry(geometry);
    auto* grid_ptr = geometry.getGridPtr();
    if (grid_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDACV5ControlVariableBackend::createControlVariable: grid pointer "
          "is null");
    }
    ensureBackendInitialized(grid_ptr);
    control.backend().resize(static_cast<size_t>(cv_size_));
    control.zero();
    return control;
  }

  /**
   * @brief Create an increment from geometry
   *
   * @param geometry Geometry backend defining the domain
   * @return Newly created increment
   */
  IncrementType createIncrement(
      const typename IncrementType::GeometryBackendType& geometry)
      const override {
    return IncrementType::createFromGeometry(geometry);
  }

  /**
   * @brief Forward transformation: δx = U v (da_transform_vtox)
   *
   * @details Transforms WRFDA CV5 control variable to state-space increment.
   * This involves:
   * 1. Horizontal recursive filter inverse
   * 2. Vertical EOF reconstruction
   * 3. Balance relationship application
   *
   * @param control CV5 control variable v
   * @param increment Output state-space increment δx (will be overwritten)
   */
  void controlToIncrement(const ControlVariableType& control,
                          IncrementType& increment) const override {
    auto* grid_ptr = control.geometry().getGridPtr();
    if (grid_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDACV5ControlVariableBackend::controlToIncrement: grid pointer is "
          "null");
    }

    ensureBackendInitialized(grid_ptr);

    auto control_data = control.template getData<std::vector<double>>();
    if (static_cast<int>(control_data.size()) != cv_size_) {
      throw std::runtime_error(
          "WRFDACV5ControlVariableBackend::controlToIncrement: control vector "
          "size mismatch");
    }

    wrfda::WRFDAControlBackendBridge::controlToState(
        grid_ptr, be_ptr_, control_data.data(), cv_size_);

    increment.backend().syncFromGrid();
  }

  /**
   * @brief Adjoint transformation: ∇_v J = U^T ∇_δx J (da_transform_vtox_adj)
   *
   * @details Transforms state-space gradient to CV5 control-space gradient.
   * This is the adjoint of the forward transformation.
   *
   * @param state_gradient Gradient in state space (∇_δx J)
   * @param control_gradient Output gradient in CV5 control space (∇_v J)
   */
  void incrementAdjointToControl(
      const IncrementType& state_gradient,
      ControlVariableType& control_gradient) const override {
    auto* grid_ptr = state_gradient.geometry().getGridPtr();
    if (grid_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDACV5ControlVariableBackend::incrementAdjointToControl: grid "
          "pointer is null");
    }

    ensureBackendInitialized(grid_ptr);

    state_gradient.backend().syncToGrid();

    std::vector<double> control_buffer(cv_size_, 0.0);
    wrfda::WRFDAControlBackendBridge::stateGradientToControl(
        grid_ptr, be_ptr_, control_buffer.data(), cv_size_);

    control_gradient.backend().setFromVector(control_buffer);
  }

  /**
   * @brief Apply increment to state: xa = xb + δx
   *
   * @param increment State-space increment δx
   * @param state State to be updated (input: xb, output: xa = xb + δx)
   */
  void addIncrementToState(const IncrementType& increment,
                           StateType& state) const override {
    state += increment;
  }

 private:
  using GeometryBackendType = typename ControlVariableType::GeometryBackendType;

  void ensureBackendInitialized(void* grid_ptr) const {
    if (grid_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDACV5ControlVariableBackend::ensureBackendInitialized: grid "
          "pointer is null");
    }

    if (be_ptr_ != nullptr && grid_ptr == active_grid_ptr_) {
      return;
    }

    finalizeBackend();
    cv_size_ = 0;
    be_ptr_ = wrfda::WRFDAControlBackendBridge::setup(grid_ptr, cv_size_);
    active_grid_ptr_ = grid_ptr;
  }

  void finalizeBackend() const {
    if (be_ptr_ != nullptr) {
      wrfda::WRFDAControlBackendBridge::finalize(be_ptr_);
      be_ptr_ = nullptr;
      active_grid_ptr_ = nullptr;
      cv_size_ = 0;
    }
  }

  mutable void* be_ptr_ = nullptr;
  mutable void* active_grid_ptr_ = nullptr;
  mutable int cv_size_ = 0;
};

}  // namespace metada::backends::wrf
