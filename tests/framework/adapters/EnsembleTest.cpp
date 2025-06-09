#include <gmock/gmock.h>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "MockBackendTraits.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "State.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::MockBackendTag;
using StateType = fwk::State<BackendTag>;
using Ensemble = fwk::Ensemble<BackendTag>;
using ConfigType = fwk::Config<BackendTag>;

namespace metada::tests {

using ::testing::InSequence;
using ::testing::Invoke;
using ::testing::NiceMock;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;

class EnsembleTest : public ::testing::Test {
 protected:
  void SetUp() override {
    size_ = 3;
    dimensions_ = {10, 20};
    test_data_ = std::vector<double>(dimensions_[0] * dimensions_[1]);
    for (size_t i = 0; i < dimensions_[0] * dimensions_[1]; i++)
      test_data_[i] = i;
  }

  /**
   * @brief Clean up test data after each test
   *
   * Frees allocated test data array
   */
  void TearDown() override { test_data_.clear(); }

  ConfigType config_;
  std::vector<size_t> dimensions_;
  std::vector<double> test_data_;
  size_t size_;
};

// Test ensemble initialization
TEST_F(EnsembleTest, InitializationCreatesCorrectNumberOfMembers) {
  Ensemble ensemble(config_, size_);
  EXPECT_EQ(ensemble.Size(), size_);
}

// Test member access
TEST_F(EnsembleTest, MemberAccessIsValid) {
  Ensemble ensemble(config_, size_);

  // Test access to all members
  for (size_t i = 0; i < size_; ++i) {
    // Here you would mock or check member access as needed
    EXPECT_NO_THROW(ensemble.GetMember(i));
  }

  // Test const access
  const auto& const_ensemble = ensemble;
  EXPECT_NO_THROW(const_ensemble.GetMember(0));
}

// Test mean computation
TEST_F(EnsembleTest, ComputeMean) {
  Ensemble ensemble(config_, size_);
  EXPECT_NO_THROW(ensemble.ComputeMean(config_));
}

// Test out of bounds access
TEST_F(EnsembleTest, OutOfBoundsAccessThrows) {
  Ensemble ensemble(config_, size_);
  EXPECT_THROW(ensemble.GetMember(size_), std::out_of_range);
  EXPECT_THROW(ensemble.GetPerturbation(size_), std::out_of_range);
}

// Test perturbation computation
TEST_F(EnsembleTest, ComputePerturbationsCreatesValidPerturbations) {
  Ensemble ensemble(config_, size_);
  EXPECT_NO_THROW(ensemble.ComputePerturbations());
  EXPECT_NO_THROW(ensemble.GetPerturbation(0));
}

}  // namespace metada::tests