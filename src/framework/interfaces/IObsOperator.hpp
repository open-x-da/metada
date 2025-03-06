#pragma once

#include <vector>

#include "IObservation.hpp"
#include "utils/NonCopyable.hpp"

namespace metada::framework {

// Forward declaration
class IState;
class IConfig;

/**
 * @brief Interface for observation operator implementations
 *
 * Defines the contract for mapping between model space and observation space
 * in data assimilation algorithms. Implementations must provide:
 *
 * 1. Forward operator (H): Map model state to observation space
 * 2. Tangent linear (H'): Linearized approximation of H
 * 3. Adjoint operator (H'*): Transpose of the tangent linear
 */
class IObsOperator : public NonCopyable {
 public:
  virtual ~IObsOperator() = default;

  // Core operations
  virtual void initialize(const IConfig& config) = 0;
  virtual bool isInitialized() const = 0;

  // Forward operator: H(x)
  virtual void apply(const IState& state, IObservation& obs) const = 0;

  // Tangent linear: H'(x)δx - Working with raw data pointers
  virtual void applyTangentLinear(const void* dx, void* dy) const = 0;

  // Adjoint: H'*(x)δy - Working with raw data pointers
  virtual void applyAdjoint(const void* dy, void* dx) const = 0;

  // Metadata
  virtual const std::vector<std::string>& getRequiredStateVars() const = 0;
  virtual const std::vector<std::string>& getRequiredObsVars() const = 0;
};

}  // namespace metada::framework