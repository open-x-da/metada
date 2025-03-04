#pragma once

#include <memory>
#include <vector>

#include "IIncrement.hpp"
#include "IObservation.hpp"
#include "IState.hpp"
#include "utils/NonCopyable.hpp"
#include "utils/config/IConfig.hpp"

namespace metada::framework {

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
class IObsOperator : public NonCopyable {
 public:
  virtual ~IObsOperator() = default;

  // Core operations
  virtual void initialize(const IConfig& config) = 0;
  virtual bool isInitialized() const = 0;

  // Forward operator
  virtual void apply(const IState& state, IObservation& obs) const = 0;

  // Tangent linear and adjoint
  virtual void applyTangentLinear(const IIncrement& dx,
                                  IObservation& dy) const = 0;
  virtual void applyAdjoint(const IObservation& dy, IIncrement& dx) const = 0;

  // Metadata
  virtual const std::vector<std::string>& getRequiredStateVars() const = 0;
  virtual const std::vector<std::string>& getRequiredObsVars() const = 0;
};

}  // namespace metada::framework