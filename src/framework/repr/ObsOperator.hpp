#ifndef METADA_FRAMEWORK_REPR_OBSOPERATOR_HPP_
#define METADA_FRAMEWORK_REPR_OBSOPERATOR_HPP_

#include "State.hpp"
#include "Observation.hpp"
#include "Increment.hpp"
#include <memory>
#include <vector>

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Base class for observation operators
 * 
 * Provides a generic interface for mapping between model space and observation space.
 * Supports:
 * - Forward operator (model -> observation)
 * - Tangent linear operator
 * - Adjoint operator
 * - Observation error handling
 * 
 * @tparam StateType Type of model state
 * @tparam ObsType Type of observation
 */
template<typename StateType, typename ObsType>
class ObsOperator {
public:
    virtual ~ObsOperator() = default;

    // Core operations
    virtual void initialize() = 0;
    virtual void finalize() = 0;

    // Forward operator: model state -> observation space
    virtual void apply(const State<StateType>& state, 
                      Observation<ObsType>& observation) const = 0;

    // Tangent linear operator: increment -> observation space
    virtual void applyTangentLinear(const Increment<StateType>& increment,
                                  Observation<ObsType>& observation) const = 0;

    // Adjoint operator: observation -> increment space
    virtual void applyAdjoint(const Observation<ObsType>& observation,
                            Increment<StateType>& increment) const = 0;

    // Error handling
    virtual void setObservationError(const Observation<ObsType>& obs) = 0;
    virtual double getObservationError(const Observation<ObsType>& obs) const = 0;

    // Configuration
    virtual void setParameter(const std::string& name, double value) = 0;
    virtual double getParameter(const std::string& name) const = 0;

protected:
    bool initialized_{false};
    std::vector<std::string> required_state_variables_;
    std::vector<std::string> required_obs_variables_;
};

}}} // namespace metada::framework::repr

#endif // METADA_FRAMEWORK_REPR_OBSOPERATOR_HPP_ 