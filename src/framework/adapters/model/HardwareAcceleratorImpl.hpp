/**
 * @file HardwareAcceleratorImpl.hpp
 * @brief Implementation of the hardware acceleration capability for models
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains the implementation of the IHardwareAccelerator interface
 * that delegates to a backend model implementation. It provides hardware
 * acceleration functionality for both physical and AI models that support
 * specialized hardware.
 */

#pragma once

#include <string>
#include <vector>

#include "../../interfaces/IHardwareAccelerator.hpp"

namespace metada::framework {

/**
 * @brief Implementation of the hardware acceleration capability
 *
 * @details
 * This class implements the IHardwareAccelerator interface by delegating
 * to the backend model implementation. It provides hardware acceleration
 * functionality for both physical and AI models that support specialized
 * hardware.
 *
 * Note: This implementation assumes the backend supports all methods.
 * The Model class is responsible for checking capability support before
 * creating this implementation.
 *
 * @tparam Backend The backend type implementing the model
 */
template <typename Backend>
class HardwareAcceleratorImpl : public IHardwareAccelerator {
 private:
  Backend& backend_;  ///< Reference to the backend model

 public:
  /**
   * @brief Construct a new HardwareAcceleratorImpl
   *
   * @param backend Reference to the backend model
   */
  explicit HardwareAcceleratorImpl(Backend& backend) : backend_(backend) {}

  /**
   * @brief Set the computational device for the model
   *
   * @param device The device specification string
   */
  void setDevice(const std::string& device) override {
    try {
      backend_.setDevice(device);
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Failed to set device: ") +
                               e.what());
    }
  }

  /**
   * @brief Get the current computational device
   *
   * @return The current device specification string
   */
  std::string getDevice() const override { return backend_.getDevice(); }

  /**
   * @brief Check if the model supports a specific device type
   *
   * @param deviceType The device type to check
   * @return true if the model supports the device type, false otherwise
   */
  bool supportsDeviceType(const std::string& deviceType) const override {
    return backend_.supportsDeviceType(deviceType);
  }

  /**
   * @brief Get a list of available devices for this model
   *
   * @return Vector of device specification strings
   */
  std::vector<std::string> getAvailableDevices() const override {
    return backend_.getAvailableDevices();
  }

  /**
   * @brief Get the amount of memory available on the current device
   *
   * @return Available memory in bytes
   */
  size_t getAvailableMemory() const override {
    return backend_.getAvailableMemory();
  }

  /**
   * @brief Set the memory limit for the model on the current device
   *
   * @param memoryLimitBytes Maximum memory to use, in bytes
   */
  void setMemoryLimit(size_t memoryLimitBytes) override {
    backend_.setMemoryLimit(memoryLimitBytes);
  }

  /**
   * @brief Get the current memory limit
   *
   * @return Current memory limit in bytes
   */
  size_t getMemoryLimit() const override { return backend_.getMemoryLimit(); }

  /**
   * @brief Enable or disable asynchronous execution
   *
   * @param enable true to enable asynchronous execution, false to disable
   */
  void setAsynchronousExecution(bool enable) override {
    backend_.setAsynchronousExecution(enable);
  }

  /**
   * @brief Check if asynchronous execution is enabled
   *
   * @return true if asynchronous execution is enabled, false otherwise
   */
  bool isAsynchronousExecutionEnabled() const override {
    return backend_.isAsynchronousExecutionEnabled();
  }

  /**
   * @brief Synchronize the device
   */
  void synchronize() override { backend_.synchronize(); }
};

}  // namespace metada::framework