#ifndef METADA_FRAMEWORK_REPR_IOBSOPERATOR_HPP_
#define METADA_FRAMEWORK_REPR_IOBSOPERATOR_HPP_

#include <memory>
#include <vector>

#include "IIncrement.hpp"
#include "IObservation.hpp"
#include "IState.hpp"

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Abstract interface for observation operator implementations
 *
 * This interface defines the contract that all observation operator
 * implementations must follow. It provides a unified API for mapping between
 * model space and observation space.
 *
 * Key features:
 * - Forward operator (model -> observation)
 * - Tangent linear operator
 * - Adjoint operator
 * - Observation error handling
 */
class IObsOperator {
 public:
  virtual ~IObsOperator() = default;

  // Core operations
  virtual void initialize() = 0;
  virtual void finalize() = 0;

  // Forward operator: model state -> observation space
  virtual void apply(const IState& state, IObservation& observation) const = 0;

  // Tangent linear operator: increment -> observation space
  virtual void applyTangentLinear(const IIncrement& increment,
                                  IObservation& observation) const = 0;

  // Adjoint operator: observation -> increment space
  virtual void applyAdjoint(const IObservation& observation,
                            IIncrement& increment) const = 0;

  // Error handling
  virtual void setObservationError(const IObservation& obs) = 0;
  virtual double getObservationError(const IObservation& obs) const = 0;

  // Configuration
  virtual void setParameter(const std::string& name, double value) = 0;
  virtual double getParameter(const std::string& name) const = 0;

  // Required variables
  virtual const std::vector<std::string>& getRequiredStateVariables() const = 0;
  virtual const std::vector<std::string>& getRequiredObsVariables() const = 0;
  virtual bool isInitialized() const = 0;
};

}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_REPR_IOBSOPERATOR_HPP_