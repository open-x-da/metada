#pragma once

#include <stdexcept>
#include <string>

#include "Config.hpp"
#include "ControlVariableBackend.hpp"
#include "Logger.hpp"
#include "backends/wrf/WRFGeometry.hpp"

extern "C" {
void wrfda_control_backend_setup(void* grid_ptr, void** be_ptr,
                                 int* cv_size_out, int* error_code);
void wrfda_control_backend_finalize(void* be_ptr, int* error_code);
void wrfda_control_backend_control_to_state(void* grid_ptr, void* be_ptr,
                                            const double* control, int cv_size,
                                            int* error_code);
void wrfda_control_backend_state_to_control(void* grid_ptr, void* be_ptr,
                                            double* control, int cv_size,
                                            int* error_code);
void wrfda_control_backend_state_gradient_to_control(void* grid_ptr,
                                                     void* be_ptr,
                                                     double* control_gradient,
                                                     int cv_size,
                                                     int* error_code);
void wrfda_control_backend_control_gradient_to_state(
    void* grid_ptr, void* be_ptr, const double* control_gradient, int cv_size,
    int* error_code);
}

namespace metada::backends::wrf {

/**
 * @brief Control-variable backend that wraps WRFDA CV5 control space.
 *
 * @details Initializes WRFDA background-error structures via the Fortran
 *          bridge and exposes control increments through the framework
 *          interface. Placeholder transforms remain until later steps wire in
 *          da_transform_vtox / da_transform_vptox.
 *
 * @tparam BackendTag Backend tag satisfying framework concepts.
 */
template <typename BackendTag>
class WRFDACV5ControlVariableBackend
    : public framework::ControlVariableBackend<BackendTag> {
  using Base = framework::ControlVariableBackend<BackendTag>;

 public:
  using typename Base::ControlVector;
  using typename Base::GeometryBackendType;
  using typename Base::IncrementType;
  using typename Base::StateType;

  explicit WRFDACV5ControlVariableBackend(
      const framework::Config<BackendTag>& config)
      : wrfda_config_(extractControlConfig(config)),
        be_statistics_path_(parseStringOption("be_statistics")),
        cv_options_(parseIntOption("cv_options", 5)) {
    auto& logger = framework::Logger<BackendTag>::Instance();
    logger.Info() << "WRFDA CV5 control backend configured";
    logger.Info() << "  be_statistics: "
                  << (be_statistics_path_.empty() ? "<not provided>"
                                                  : be_statistics_path_);
    logger.Info() << "  cv_options: " << cv_options_;
  }

  ~WRFDACV5ControlVariableBackend() override { finalizeBackend(); }

  IncrementType createIncrement(
      const GeometryBackendType& geometry) const override {
    ensureInitialized(geometry);
    return IncrementType::createFromGeometry(geometry);
  }

  void addIncrementToState(const IncrementType& increment,
                           StateType& state) const override {
    if (auto* geom = state.geometry()) {
      ensureInitialized(geom->backend());
    }

    ControlVector control = createControlVector(increment.geometry());
    convertStateIncrementToControl(increment, control);

    IncrementType transformed(increment.geometry());
    convertControlToStateIncrement(control, transformed);

    state += transformed;
  }

  std::string name() const override { return "wrfda_cv5"; }

  const std::string& beStatisticsPath() const { return be_statistics_path_; }
  int cvOptions() const { return cv_options_; }

  ControlVector createControlVector(
      const GeometryBackendType& geometry) const override {
    ensureInitialized(geometry);
    return ControlVector(static_cast<size_t>(cv_size_), 0.0);
  }

  void convertStateIncrementToControl(const IncrementType& increment,
                                      ControlVector& control) const override {
    ensureInitialized(increment.geometry());
    ensureControlVectorSize(control);
    increment.backend().syncToGrid();
    int error_code = 0;
    wrfda_control_backend_state_to_control(
        const_cast<void*>(increment.geometry().getGridPtr()), be_ptr_,
        control.data(), cv_size_, &error_code);
    if (error_code != 0) {
      throw std::runtime_error(
          "wrfda_control_backend_state_to_control failed with error code " +
          std::to_string(error_code));
    }
  }

  void convertControlToStateIncrement(const ControlVector& control,
                                      IncrementType& increment) const override {
    ensureInitialized(increment.geometry());
    if (static_cast<int>(control.size()) != cv_size_) {
      throw std::runtime_error(
          "Control vector size mismatch for WRFDA CV5 backend");
    }
    int error_code = 0;
    wrfda_control_backend_control_to_state(
        const_cast<void*>(increment.geometry().getGridPtr()), be_ptr_,
        control.data(), cv_size_, &error_code);
    if (error_code != 0) {
      throw std::runtime_error(
          "wrfda_control_backend_control_to_state failed with error code " +
          std::to_string(error_code));
    }
    increment.backend().syncFromGrid();
  }

  void convertStateGradientToControl(
      const IncrementType& state_gradient,
      ControlVector& control_gradient) const override {
    ensureInitialized(state_gradient.geometry());
    ensureControlVectorSize(control_gradient);
    state_gradient.backend().syncToGrid();
    int error_code = 0;
    wrfda_control_backend_state_gradient_to_control(
        const_cast<void*>(state_gradient.geometry().getGridPtr()), be_ptr_,
        control_gradient.data(), cv_size_, &error_code);
    if (error_code != 0) {
      throw std::runtime_error(
          "wrfda_control_backend_state_gradient_to_control failed with error "
          "code " +
          std::to_string(error_code));
    }
  }

  void convertControlGradientToState(
      const ControlVector& control_gradient,
      IncrementType& state_gradient) const override {
    ensureInitialized(state_gradient.geometry());
    if (static_cast<int>(control_gradient.size()) != cv_size_) {
      throw std::runtime_error(
          "Control gradient size mismatch for WRFDA CV5 backend");
    }
    int error_code = 0;
    wrfda_control_backend_control_gradient_to_state(
        const_cast<void*>(state_gradient.geometry().getGridPtr()), be_ptr_,
        control_gradient.data(), cv_size_, &error_code);
    if (error_code != 0) {
      throw std::runtime_error(
          "wrfda_control_backend_control_gradient_to_state failed with error "
          "code " +
          std::to_string(error_code));
    }
    state_gradient.backend().syncFromGrid();
  }

 private:
  void ensureInitialized(const GeometryBackendType& geometry) const {
    if (backend_initialized_) {
      return;
    }

    void* grid_ptr = geometry.getGridPtr();
    if (grid_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDA CV5 backend initialization failed: null grid pointer");
    }

    void* local_be_ptr = nullptr;
    int local_cv_size = 0;
    int error_code = 0;
    wrfda_control_backend_setup(grid_ptr, &local_be_ptr, &local_cv_size,
                                &error_code);

    if (error_code != 0) {
      throw std::runtime_error(
          "wrfda_control_backend_setup failed with error code " +
          std::to_string(error_code));
    }

    if (local_be_ptr == nullptr) {
      throw std::runtime_error(
          "wrfda_control_backend_setup returned null be pointer");
    }

    be_ptr_ = local_be_ptr;
    cv_size_ = local_cv_size;
    backend_initialized_ = true;

    framework::Logger<BackendTag>::Instance().Info()
        << "WRFDA CV5 backend initialized (cv_size=" << cv_size_ << ")";
  }

  void finalizeBackend() noexcept {
    if (!backend_initialized_ || be_ptr_ == nullptr) {
      return;
    }

    int error_code = 0;
    wrfda_control_backend_finalize(be_ptr_, &error_code);
    if (error_code != 0) {
      framework::Logger<BackendTag>::Instance().Warning()
          << "wrfda_control_backend_finalize returned error code "
          << error_code;
    }

    be_ptr_ = nullptr;
    backend_initialized_ = false;
  }

  static framework::Config<BackendTag> extractControlConfig(
      const framework::Config<BackendTag>& config) {
    if (!config.HasKey("wrfda_control")) {
      throw std::runtime_error(
          "wrfda_cv5 backend requires a 'wrfda_control' configuration section");
    }
    return config.GetSubsection("wrfda_control");
  }

  std::string parseStringOption(const std::string& key) const {
    if (wrfda_config_.HasKey(key)) {
      return wrfda_config_.Get(key).asString();
    }
    return {};
  }

  int parseIntOption(const std::string& key, int default_value) const {
    if (wrfda_config_.HasKey(key)) {
      return wrfda_config_.Get(key).asInt();
    }
    return default_value;
  }

  framework::Config<BackendTag> wrfda_config_;
  std::string be_statistics_path_;
  int cv_options_;
  mutable bool backend_initialized_ = false;
  mutable void* be_ptr_ = nullptr;
  mutable int cv_size_ = 0;

  void ensureControlVectorSize(ControlVector& control) const {
    if (static_cast<int>(control.size()) != cv_size_) {
      control.assign(static_cast<size_t>(cv_size_), 0.0);
    }
  }
};

}  // namespace metada::backends::wrf
