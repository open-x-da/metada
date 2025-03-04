#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "IObsOperator.hpp"
#include "Increment.hpp"
#include "Observation.hpp"
#include "State.hpp"
#include "utils/NonCopyable.hpp"
#include "utils/config/IConfig.hpp"

namespace metada::framework {

/**
 * @brief Main observation operator class template providing a generic interface
 *
 * This class template provides a static interface for observation operators
 * using a backend specified by the ObsOperatorBackend template parameter. The
 * backend must implement the IObsOperator interface.
 *
 * The observation operator is responsible for:
 * - Mapping model state variables to observation space (H operator)
 * - Computing tangent linear approximations (H' operator)
 * - Computing adjoint operations (H'* operator)
 * - Managing observation errors and uncertainties
 *
 * @tparam ObsOperatorBackend The observation operator backend type
 * @tparam StateType The state type used by this operator
 * @tparam ObsType The observation type used by this operator
 */
template <typename Backend, typename StateType, typename IncrementType,
          typename ObsType>
class ObsOperator : public NonCopyable {
 private:
  Backend& backend_;
  bool initialized_{false};

 public:
  /** @brief Default constructor - deleted to prevent usage without backend */
  ObsOperator() = delete;

  /** @brief Constructor with backend reference */
  explicit ObsOperator(Backend& backend) : backend_(backend) {}

  /** @brief Get direct access to the backend instance */
  Backend& backend() { return backend_; }

  /** @brief Get const access to the backend instance */
  const Backend& backend() const { return backend_; }

  void initialize(const IConfig& config) {
    backend_.initialize(config);
    initialized_ = true;
  }

  bool isInitialized() const { return initialized_; }

  // Forward operator: model state -> observation space
  void apply(const State<StateType>& state,
             const Observation<ObsType>& obs) const {
    if (!initialized_) throw std::runtime_error("ObsOperator not initialized");
    backend_.apply(state.backend(), obs.backend());
  }

  // Tangent linear operator: increment -> observation space
  void applyTangentLinear(const Increment<IncrementType>& dx,
                          const Observation<ObsType>& dy) const {
    if (!initialized_) throw std::runtime_error("ObsOperator not initialized");
    backend_.applyTangentLinear(dx.backend(), dy.backend());
  }

  // Adjoint operator: observation -> increment space
  void applyAdjoint(const Observation<ObsType>& dy,
                    Increment<IncrementType>& dx) const {
    if (!initialized_) throw std::runtime_error("ObsOperator not initialized");
    backend_.applyAdjoint(dy.backend(), dx.backend());
  }

  // Required variables
  const std::vector<std::string>& getRequiredStateVars() const {
    return backend_.getRequiredStateVars();
  }

  const std::vector<std::string>& getRequiredObsVars() const {
    return backend_.getRequiredObsVars();
  }
};

}  // namespace metada::framework