/**
 * @file IModel.hpp
 * @brief Interface defining the contract for model implementations
 * @ingroup repr
 *
 * @details
 * This header provides the abstract interface that all model implementations
 * must follow. It defines a unified API for handling numerical models in
 * scientific computations and data assimilation, as well as AI-based models
 * like GraphCast and Pangu-Weather.
 *
 * The interface follows a component-based design where:
 * - Core initialization and parameter management are required for all models
 * - Specific capabilities (time stepping, AI prediction, etc.) are provided
 *   through optional capability interfaces
 *
 * This design allows support for both traditional numerical models and
 * modern AI-based forecasting systems within the same framework.
 */

#pragma once

#include <string>
#include <vector>

#include "../common/utils/NonCopyable.hpp"

namespace metada::framework {

// Forward declarations
class IState;
class IConfig;
class ITimeStepper;
class IAIPredictor;
class IBatchProcessor;
class IHardwareAccelerator;

/**
 * @brief Abstract interface for model implementations
 *
 * @details
 * This interface defines the core contract that all model implementations must
 * follow, supporting both physical numerical models and AI-based models.
 *
 * The interface uses a component-based design where models can advertise
 * specific capabilities through capability interface queries.
 *
 * Core functionality required by all models includes:
 * - Model initialization from configuration
 * - Parameter management
 * - Basic state handling
 *
 * Optional capabilities include:
 * - Time stepping (for numerical models)
 * - AI prediction (for pre-trained models)
 * - Batch processing
 *
 * Implementations should ensure:
 * - Thread safety for all operations
 * - Proper exception handling and error reporting
 * - Efficient memory management
 */
class IModel : private NonCopyable {
 protected:
  /**
   * @brief Protected constructor requiring configuration
   *
   * All models must be initialized with a configuration object.
   * Derived classes should provide their own constructors that take
   * appropriate configuration parameters.
   *
   * @param config Configuration object with model parameters
   */
  explicit IModel([[maybe_unused]] const IConfig& config) {}

 public:
  /**
   * @brief Virtual destructor
   */
  virtual ~IModel() = default;

  /**
   * @brief Initialize the model with the given configuration
   *
   * @param config Configuration object with model parameters
   * @throws std::runtime_error if initialization fails
   */
  virtual void initialize(const IConfig& config) = 0;

  /**
   * @brief Reset the model to its initial state
   */
  virtual void reset() = 0;

  /**
   * @brief Check if the model is initialized
   *
   * @return true if the model is initialized, false otherwise
   */
  virtual bool isInitialized() const = 0;

  /**
   * @brief Get a parameter from the model
   *
   * @param name The name of the parameter
   * @return The value of the parameter as a string
   * @throws std::out_of_range if the parameter doesn't exist
   */
  virtual std::string getParameter(const std::string& name) const = 0;

  /**
   * @brief Set a parameter in the model
   *
   * @param name The name of the parameter
   * @param value The value of the parameter
   * @throws std::invalid_argument if the parameter is invalid
   */
  virtual void setParameter(const std::string& name,
                            const std::string& value) = 0;

  /**
   * @brief Run the model from initial state to final state
   *
   * This is the primary model execution method that all models must implement.
   * The implementation may utilize any of the capabilities the model provides
   * (time stepping, AI prediction, etc.) to produce the final state.
   *
   * @param initialState The initial state of the model
   * @param finalState The final state after model execution (output parameter)
   * @param startTime The start time for the model run
   * @param endTime The end time for the model run
   * @throws std::runtime_error if the model run fails
   */
  virtual void run(const IState& initialState, IState& finalState,
                   double startTime, double endTime) = 0;

  /**
   * @brief Get the model's time stepping capability
   *
   * This method returns a pointer to the model's time stepping capability
   * interface if the model supports time stepping. If the model does not
   * support time stepping, this method returns nullptr.
   *
   * @return Pointer to the model's time stepping capability, or nullptr if not
   * supported
   */
  virtual ITimeStepper* getTimeStepper() { return nullptr; }

  /**
   * @brief Get the model's AI prediction capability
   *
   * This method returns a pointer to the model's AI prediction capability
   * interface if the model is an AI-based model. If the model does not
   * support AI prediction, this method returns nullptr.
   *
   * @return Pointer to the model's AI prediction capability, or nullptr if not
   * supported
   */
  virtual IAIPredictor* getAIPredictor() { return nullptr; }

  /**
   * @brief Get the model's batch processing capability
   *
   * This method returns a pointer to the model's batch processing capability
   * interface if the model supports batch processing. If the model does not
   * support batch processing, this method returns nullptr.
   *
   * @return Pointer to the model's batch processing capability, or nullptr if
   * not supported
   */
  virtual IBatchProcessor* getBatchProcessor() { return nullptr; }

  /**
   * @brief Get the model's hardware acceleration capability
   *
   * This method returns a pointer to the model's hardware acceleration
   * capability if the model supports hardware acceleration. If the model does
   * not support hardware acceleration, this method returns nullptr.
   *
   * @return Pointer to the model's hardware acceleration capability, or nullptr
   */
  virtual IHardwareAccelerator* getHardwareAccelerator() { return nullptr; }
};

}  // namespace metada::framework