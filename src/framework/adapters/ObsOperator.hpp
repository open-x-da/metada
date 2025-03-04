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
template <typename Backend>
class ObsOperator : public NonCopyable {
 private:
  std::unique_ptr<Backend> backend_;
  bool initialized_{false};

 public:
  /** @brief Default constructor - deleted to prevent usage without config */
  ObsOperator() = delete;

  /** @brief Constructor with configuration */
  template <typename Config>
  explicit ObsOperator(const Config& config)
      : backend_(std::make_unique<Backend>()) {
    initialize(config);
  }

  /** @brief Move constructor */
  ObsOperator(ObsOperator&& other) noexcept
      : backend_(std::move(other.backend_)), initialized_(other.initialized_) {
    other.initialized_ = false;
  }

  /** @brief Move assignment operator */
  ObsOperator& operator=(ObsOperator&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
      initialized_ = other.initialized_;
      other.initialized_ = false;
    }
    return *this;
  }

  /** @brief Get const access to the backend instance */
  const Backend& backend() const { return *backend_; }

  template <typename Config>
  void initialize(const Config& config) {
    if (initialized_) {
      throw std::runtime_error("ObsOperator already initialized");
    }
    backend_->initialize(config);
    initialized_ = true;
  }

  bool isInitialized() const { return initialized_; }

  // Forward operator: model state -> observation space
  template <typename StateType, typename ObsType>
  void apply(const State<StateType>& state,
             const Observation<ObsType>& obs) const {
    if (!initialized_) throw std::runtime_error("ObsOperator not initialized");
    backend_->apply(state.backend(), obs.backend());
  }

  // Tangent linear operator: increment -> observation space
  template <typename IncrementType, typename ObsType>
  void applyTangentLinear(const Increment<IncrementType>& dx,
                          const Observation<ObsType>& dy) const {
    if (!initialized_) throw std::runtime_error("ObsOperator not initialized");
    backend_->applyTangentLinear(dx.backend(), dy.backend());
  }

  // Adjoint operator: observation -> increment space
  template <typename IncrementType, typename ObsType>
  void applyAdjoint(const Observation<ObsType>& dy,
                    Increment<IncrementType>& dx) const {
    if (!initialized_) throw std::runtime_error("ObsOperator not initialized");
    backend_->applyAdjoint(dy.backend(), dx.backend());
  }

  // Required variables
  const std::vector<std::string>& getRequiredStateVars() const {
    return backend_->getRequiredStateVars();
  }

  const std::vector<std::string>& getRequiredObsVars() const {
    return backend_->getRequiredObsVars();
  }
};

}  // namespace metada::framework