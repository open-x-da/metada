#ifndef METADA_FRAMEWORK_REPR_IMODEL_HPP_
#define METADA_FRAMEWORK_REPR_IMODEL_HPP_

#include <memory>
#include <string>
#include "IState.hpp"

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Abstract interface for model implementations
 *
 * This interface defines the contract that all model implementations must follow.
 * It provides a unified API for model operations including initialization,
 * time stepping, state management, and parameter configuration.
 *
 * Key features:
 * - Model lifecycle management (initialize, step, finalize)
 * - State handling
 * - Parameter configuration
 * - Model metadata access
 *
 * Example usage:
 * @code
 * template<typename StateType>
 * class MyModel : public IModel<StateType> {
 *   void initialize() override;
 *   void step(double dt) override;
 *   void finalize() override;
 *   // ... implement other interface methods
 * };
 * @endcode
 */
class IModel {
public:
    /**
     * @brief Virtual destructor
     *
     * Ensures proper cleanup of derived classes through base pointer
     */
    virtual ~IModel() = default;

    // Core model operations
    /**
     * @brief Initialize the model
     *
     * Sets up internal data structures and prepares the model for time stepping.
     * Must be called before any step() operations.
     */
    virtual void initialize() = 0;

    /**
     * @brief Advance the model state by one time step
     *
     * @param dt Time step size
     */
    virtual void step(double dt) = 0;

    /**
     * @brief Cleanup model resources
     *
     * Performs any necessary cleanup operations when the model is no longer needed.
     */
    virtual void finalize() = 0;

    // State management
    /**
     * @brief Set the model state
     *
     * @param state New state to be used by the model
     */
    virtual void setState(const IState& state) = 0;

    /**
     * @brief Get the current model state
     *
     * @return Current state of the model
     */
    virtual const IState& getState() const = 0;

    // Model configuration
    /**
     * @brief Set a model parameter value
     *
     * @param name Parameter name
     * @param value Parameter value
     */
    virtual void setParameter(const std::string& name, double value) = 0;

    /**
     * @brief Get a model parameter value
     *
     * @param name Parameter name
     * @return Current value of the parameter
     */
    virtual double getParameter(const std::string& name) const = 0;

    // Model metadata
    /**
     * @brief Get the model name
     *
     * @return Name of the model
     */
    virtual std::string getName() const = 0;

    /**
     * @brief Get the model version
     *
     * @return Version string of the model
     */
    virtual std::string getVersion() const = 0;

    // Model state
    /**
     * @brief Get the current model time
     *
     * @return Current time of the model
     */
    virtual double getCurrentTime() const = 0;

    /**
     * @brief Check if the model is initialized
     *
     * @return True if the model is initialized, false otherwise
     */
    virtual bool isInitialized() const = 0;
};

}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_REPR_IMODEL_HPP_ 