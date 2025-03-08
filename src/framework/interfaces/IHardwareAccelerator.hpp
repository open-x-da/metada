/**
 * @file IHardwareAccelerator.hpp
 * @brief Interface defining the capability for hardware-accelerated models
 * @ingroup repr
 *
 * @details
 * This header provides the capability interface for models that support
 * hardware acceleration, including GPU acceleration, multi-node computation,
 * and other specialized hardware. This interface can be used by both
 * physical numerical models and AI-based models.
 */

#pragma once

#include <string>
#include <vector>

namespace metada::framework {

/**
 * @brief Capability interface for hardware-accelerated models
 *
 * @details
 * This interface defines the capabilities required for models that
 * can run on specialized hardware like GPUs, TPUs, or distributed
 * computing resources. It provides methods for device management,
 * resource allocation, and hardware capability discovery.
 *
 * This capability is typically used by:
 * - GPU-accelerated numerical models (physics-based)
 * - AI models running on GPUs, TPUs, etc.
 * - Models using MPI for distributed computing
 * - Any model that can benefit from specialized hardware
 *
 * Implementations should ensure:
 * - Proper device selection and management
 * - Resource allocation and deallocation
 * - Error handling for device-specific issues
 */
class IHardwareAccelerator {
 public:
  /**
   * @brief Virtual destructor
   */
  virtual ~IHardwareAccelerator() = default;

  /**
   * @brief Set the computational device for the model
   *
   * @param device The device specification string (format depends on backend)
   *               Examples: "cpu", "cuda:0", "cuda:1", "mpi:4", "tpu:0", etc.
   * @throws std::runtime_error if the device is not available or not supported
   */
  virtual void setDevice(const std::string& device) = 0;

  /**
   * @brief Get the current computational device
   *
   * @return The current device specification string
   */
  virtual std::string getDevice() const = 0;

  /**
   * @brief Check if the model supports a specific device type
   *
   * @param deviceType The device type to check ("cpu", "cuda", "mpi", "tpu",
   * etc.)
   * @return true if the model supports the device type, false otherwise
   */
  virtual bool supportsDeviceType(const std::string& deviceType) const = 0;

  /**
   * @brief Get a list of available devices for this model
   *
   * @return Vector of device specification strings that can be used with
   * setDevice()
   */
  virtual std::vector<std::string> getAvailableDevices() const = 0;

  /**
   * @brief Get the amount of memory available on the current device
   *
   * @return Available memory in bytes
   */
  virtual size_t getAvailableMemory() const = 0;

  /**
   * @brief Set the memory limit for the model on the current device
   *
   * @param memoryLimitBytes Maximum memory to use, in bytes
   * @throws std::invalid_argument if the memory limit is invalid
   */
  virtual void setMemoryLimit(size_t memoryLimitBytes) = 0;

  /**
   * @brief Get the current memory limit
   *
   * @return Current memory limit in bytes
   */
  virtual size_t getMemoryLimit() const = 0;

  /**
   * @brief Enable or disable asynchronous execution
   *
   * Asynchronous execution can improve performance but may affect determinism.
   *
   * @param enable true to enable asynchronous execution, false to disable
   */
  virtual void setAsynchronousExecution(bool enable) = 0;

  /**
   * @brief Check if asynchronous execution is enabled
   *
   * @return true if asynchronous execution is enabled, false otherwise
   */
  virtual bool isAsynchronousExecutionEnabled() const = 0;

  /**
   * @brief Synchronize the device (wait for all pending operations to complete)
   *
   * This is particularly important when using asynchronous execution.
   */
  virtual void synchronize() = 0;
};

}  // namespace metada::framework