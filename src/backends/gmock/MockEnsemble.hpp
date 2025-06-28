/**
 * @file MockEnsemble.hpp
 * @brief Mock implementation of ensemble backend for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the ensemble backend using Google
 * Mock. It allows testing code that depends on ensemble operations by providing
 * mock implementations of all interface methods that can be configured with
 * expectations and behaviors.
 *
 * The mock implementation supports:
 * - Setting expectations on method calls
 * - Configuring return values and behaviors
 * - Verifying interaction patterns
 * - Testing error conditions
 * - Member access and manipulation
 * - Mean and perturbation computations
 * - Move semantics
 *
 * @see Ensemble
 * @see MockState
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <memory>
#include <vector>

#include "MockState.hpp"

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of ensemble backend for testing
 *
 * @details
 * Provides mock methods for all ensemble operations, organized into
 * the following categories:
 *
 * @par Core Operations
 * - Initialize() - Initialize ensemble from configuration and geometry
 * - GetMember() - Get access to an ensemble member
 * - Size() - Get the number of ensemble members
 *
 * @par Ensemble Statistics
 * - ComputeMean() - Compute the ensemble mean
 * - GetMean() - Get access to the ensemble mean
 * - ComputePerturbations() - Compute member perturbations
 * - GetPerturbation() - Get access to a perturbation
 *
 * @par Memory Management
 * - Move constructor and assignment operator
 * - Deleted copy constructor and assignment operator (non-copyable)
 *
 * @note All mock methods use Google Mock's MOCK_METHOD macro to enable
 * setting expectations and verifying calls.
 */
template <typename ConfigBackend, typename GeometryBackend>
class MockEnsemble {
 public:
  using StateType = MockState<ConfigBackend, GeometryBackend>;

  // Disable default constructor
  MockEnsemble() = delete;

  // Destructor
  ~MockEnsemble() = default;

  // Copy constructor
  MockEnsemble(const MockEnsemble& other) = delete;

  // Copy assignment operator
  MockEnsemble& operator=(const MockEnsemble& other) = delete;

  // Move constructor
  MockEnsemble(MockEnsemble&& other) noexcept
      : config_(other.config_),
        geometry_(other.geometry_),
        members_(std::move(other.members_)),
        mean_(std::move(other.mean_)),
        perturbations_(std::move(other.perturbations_)),
        size_(other.size_) {
    other.members_.clear();
    other.perturbations_.clear();
    other.size_ = 0;
  }

  // Move assignment operator
  MockEnsemble& operator=(MockEnsemble&& other) noexcept {
    if (this != &other) {
      members_ = std::move(other.members_);
      mean_ = std::move(other.mean_);
      perturbations_ = std::move(other.perturbations_);
      size_ = other.size_;
      other.members_.clear();
      other.perturbations_.clear();
      other.size_ = 0;
    }
    return *this;
  }

  // Constructor that initializes ensemble from config and geometry
  MockEnsemble(const ConfigBackend& config, const GeometryBackend& geometry)
      : config_(config), geometry_(geometry) {
    Initialize();
  }

  // Core ensemble operations
  MOCK_METHOD(void, Initialize, ());
  MOCK_METHOD(size_t, Size, (), (const));
  MOCK_METHOD(StateType&, GetMember, (size_t index));
  MOCK_METHOD(const StateType&, GetMember, (size_t index), (const));

  // Ensemble statistics
  MOCK_METHOD(void, ComputeMean, ());
  MOCK_METHOD(StateType&, GetMean, ());
  MOCK_METHOD(const StateType&, GetMean, (), (const));
  MOCK_METHOD(void, ComputePerturbations, ());
  MOCK_METHOD(StateType&, GetPerturbation, (size_t index));
  MOCK_METHOD(const StateType&, GetPerturbation, (size_t index), (const));

  // Test helper methods
  void SetSize(size_t size) { size_ = size; }
  void AddMember(std::unique_ptr<StateType> member) {
    members_.push_back(std::move(member));
  }
  void SetMean(std::unique_ptr<StateType> mean) { mean_ = std::move(mean); }
  void AddPerturbation(std::unique_ptr<StateType> perturbation) {
    perturbations_.push_back(std::move(perturbation));
  }

  // Get the config
  const ConfigBackend& Config() const { return config_; }
  // Get the geometry
  const GeometryBackend& Geometry() const { return geometry_; }

 private:
  const ConfigBackend& config_;
  const GeometryBackend& geometry_;
  std::vector<std::unique_ptr<StateType>> members_;
  std::unique_ptr<StateType> mean_;
  std::vector<std::unique_ptr<StateType>> perturbations_;
  size_t size_{0};
};

}  // namespace metada::backends::gmock