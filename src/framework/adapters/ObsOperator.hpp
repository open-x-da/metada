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

namespace metada::framework {

/**
 * @brief Observation operator for data assimilation algorithms
 *
 * Maps between model space and observation space, providing essential
 * operations for variational and ensemble data assimilation:
 * - H: Forward observation operator (state → observation)
 * - H': Tangent linear operator (increment → observation)
 * - H'*: Adjoint operator (observation → increment)
 *
 * @tparam Backend Backend implementation type
 */
template <typename Backend>
class ObsOperator : public NonCopyable {
 private:
  std::unique_ptr<Backend> backend_;
  bool initialized_{false};

  /** @brief Check initialization status and throw if not initialized */
  void checkInitialized() const {
    if (!isInitialized()) {
      throw std::runtime_error("ObsOperator not initialized");
    }
  }

 public:
  /** @brief Default constructor deleted - must initialize with config */
  ObsOperator() = delete;

  /**
   * @brief Constructor with configuration
   *
   * @tparam Config Configuration type
   * @param config Configuration object
   */
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

  /** @brief Get const access to backend */
  const Backend& backend() const { return *backend_; }

  /**
   * @brief Initialize with configuration
   *
   * @tparam Config Configuration type
   * @param config Configuration object
   * @throws std::runtime_error If already initialized
   */
  template <typename Config>
  void initialize(const Config& config) {
    if (initialized_) {
      throw std::runtime_error("ObsOperator already initialized");
    }
    backend_->initialize(config);
    initialized_ = true;
  }

  /** @brief Check if initialized */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Apply forward operator: H(x)
   *
   * Maps state to observation space (y = H(x))
   *
   * @tparam StateType State type
   * @tparam ObsType Observation type
   * @param state Model state
   * @param obs Output observation
   * @throws std::runtime_error If not initialized
   */
  template <typename StateType, typename ObsType>
  void apply(const State<StateType>& state, Observation<ObsType>& obs) const {
    checkInitialized();
    backend_->apply(state.backend(), obs.backend());
  }

  /**
   * @brief Apply tangent linear operator: H'(x)δx
   *
   * Maps increment to observation space (δy = H'(x)δx)
   *
   * @tparam IncrementType Increment type
   * @tparam ObsType Observation type
   * @param dx State increment
   * @param dy Output observation increment
   * @throws std::runtime_error If not initialized
   */
  template <typename IncrementType, typename ObsType>
  void applyTangentLinear(const Increment<IncrementType>& dx,
                          Observation<ObsType>& dy) const {
    checkInitialized();
    backend_->applyTangentLinear(dx.backend().getData(), dy.backend());
  }

  /**
   * @brief Apply adjoint operator: H'*(x)δy
   *
   * Maps observation to increment space (δx = H'*(x)δy)
   *
   * @tparam IncrementType Increment type
   * @tparam ObsType Observation type
   * @param dy Observation increment
   * @param dx Output state increment
   * @throws std::runtime_error If not initialized
   */
  template <typename IncrementType, typename ObsType>
  void applyAdjoint(const Observation<ObsType>& dy,
                    Increment<IncrementType>& dx) const {
    checkInitialized();
    backend_->applyAdjoint(dy.backend(), dx.backend().getData());
  }

  /** @brief Get required state variables */
  const std::vector<std::string>& getRequiredStateVars() const {
    return backend_->getRequiredStateVars();
  }

  /** @brief Get required observation variables */
  const std::vector<std::string>& getRequiredObsVars() const {
    return backend_->getRequiredObsVars();
  }
};

}  // namespace metada::framework