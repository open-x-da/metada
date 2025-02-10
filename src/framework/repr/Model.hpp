#ifndef METADA_FRAMEWORK_REPR_MODEL_HPP_
#define METADA_FRAMEWORK_REPR_MODEL_HPP_

#include "State.hpp"
#include <memory>
#include <string>

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Base class for scientific model representations
 * 
 * Provides a generic interface for scientific models, supporting:
 * - Model state initialization and management
 * - Time integration/stepping
 * - Model parameters and configuration
 * - Input/Output operations
 * 
 * @tparam StateType Type of state used by the model
 */
template<typename StateType>
class Model {
public:
    virtual ~Model() = default;

    // Core model operations
    virtual void initialize() = 0;
    virtual void step(double dt) = 0;
    virtual void finalize() = 0;

    // State management
    virtual void setState(std::shared_ptr<StateType> state) = 0;
    virtual std::shared_ptr<StateType> getState() const = 0;

    // Model configuration
    virtual void setParameter(const std::string& name, double value) = 0;
    virtual double getParameter(const std::string& name) const = 0;

    // Model metadata
    virtual std::string getName() const = 0;
    virtual std::string getVersion() const = 0;

protected:
    std::shared_ptr<StateType> state_;
    double current_time_{0.0};
    bool initialized_{false};
};

}}} // namespace metada::framework::repr

#endif // METADA_FRAMEWORK_REPR_MODEL_HPP_ 