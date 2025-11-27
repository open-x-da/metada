#pragma once

#include <stdexcept>
#include <string>

namespace metada::backends::wrf::wrfda {

extern "C" {
void wrfda_control_backend_setup(void* grid_ptr, void** be_ptr,
                                 int* cv_size_out, int* error_code);
void wrfda_control_backend_finalize(void* be_ptr, int* error_code);
void wrfda_control_backend_control_to_state(void* grid_ptr, void* be_ptr,
                                            const double* control,
                                            int cv_size_in, int* error_code);
void wrfda_control_backend_state_to_control(void* grid_ptr, void* be_ptr,
                                            double* control, int cv_size_out,
                                            int* error_code);
void wrfda_control_backend_state_gradient_to_control(void* grid_ptr,
                                                     void* be_ptr,
                                                     double* control_gradient,
                                                     int cv_size_grad,
                                                     int* error_code);
void wrfda_control_backend_control_gradient_to_state(void* grid_ptr,
                                                     void* be_ptr,
                                                     const double* control_grad,
                                                     int cv_size_grad,
                                                     int* error_code);
void wrfda_control_backend_control_dot(void* grid_ptr, void* be_ptr,
                                       const double* control_a,
                                       const double* control_b, int cv_size,
                                       double* dot_value, int* error_code);
}  // extern "C"

class WRFDAControlBackendBridge {
 public:
  static void* setup(void* grid_ptr, int& cv_size) {
    if (grid_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::setup: grid_ptr is null");
    }
    void* be_ptr = nullptr;
    int error_code = 0;
    int cv_size_out = 0;
    wrfda_control_backend_setup(grid_ptr, &be_ptr, &cv_size_out, &error_code);
    if (error_code != 0 || be_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::setup failed with error code " +
          std::to_string(error_code));
    }
    cv_size = cv_size_out;
    return be_ptr;
  }

  static void finalize(void* be_ptr) {
    if (be_ptr == nullptr) {
      return;
    }
    int error_code = 0;
    wrfda_control_backend_finalize(be_ptr, &error_code);
    if (error_code != 0) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::finalize failed with error code " +
          std::to_string(error_code));
    }
  }

  static void controlToState(void* grid_ptr, void* be_ptr,
                             const double* control, int cv_size) {
    if (grid_ptr == nullptr || be_ptr == nullptr || control == nullptr) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::controlToState received null pointer");
    }
    int error_code = 0;
    wrfda_control_backend_control_to_state(grid_ptr, be_ptr, control, cv_size,
                                           &error_code);
    if (error_code != 0) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::controlToState failed with error code " +
          std::to_string(error_code));
    }
  }

  static void stateGradientToControl(void* grid_ptr, void* be_ptr,
                                     double* control_gradient, int cv_size) {
    if (grid_ptr == nullptr || be_ptr == nullptr ||
        control_gradient == nullptr) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::stateGradientToControl received null "
          "pointer");
    }
    int error_code = 0;
    wrfda_control_backend_state_gradient_to_control(
        grid_ptr, be_ptr, control_gradient, cv_size, &error_code);
    if (error_code != 0) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::stateGradientToControl failed with error "
          "code " +
          std::to_string(error_code));
    }
  }

  static void controlGradientToState(void* grid_ptr, void* be_ptr,
                                     const double* control_gradient,
                                     int cv_size) {
    if (grid_ptr == nullptr || be_ptr == nullptr ||
        control_gradient == nullptr) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::controlGradientToState received null "
          "pointer");
    }
    int error_code = 0;
    wrfda_control_backend_control_gradient_to_state(
        grid_ptr, be_ptr, control_gradient, cv_size, &error_code);
    if (error_code != 0) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::controlGradientToState failed with error "
          "code " +
          std::to_string(error_code));
    }
  }

  static double controlDot(void* grid_ptr, void* be_ptr,
                           const double* control_a, const double* control_b,
                           int cv_size) {
    if (grid_ptr == nullptr || be_ptr == nullptr || control_a == nullptr ||
        control_b == nullptr) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::controlDot received null pointer");
    }
    double dot_value = 0.0;
    int error_code = 0;
    wrfda_control_backend_control_dot(grid_ptr, be_ptr, control_a, control_b,
                                      cv_size, &dot_value, &error_code);
    if (error_code != 0) {
      throw std::runtime_error(
          "WRFDAControlBackendBridge::controlDot failed with error code " +
          std::to_string(error_code));
    }
    return dot_value;
  }
};

}  // namespace metada::backends::wrf::wrfda
