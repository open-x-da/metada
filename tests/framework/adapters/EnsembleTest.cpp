#include <gmock/gmock.h>

#include "Ensemble.hpp"
#include "MockConfig.hpp"
#include "MockState.hpp"
#include "State.hpp"
#include "utils/config/Config.hpp"

namespace metada::tests {

using ::testing::InSequence;
using ::testing::Invoke;
using ::testing::NiceMock;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;

using StateType = framework::State<MockState>;
using Ensemble = framework::Ensemble<StateType>;

using ConfigType = framework::Config<MockConfig>;

class EnsembleTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Setup test data
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
/*TEST_F(EnsembleTest, InitializationCreatesCorrectNumberOfMembers) {
  Ensemble ensemble(config_, size_);
  EXPECT_EQ(ensemble.getSize(), size_);
}*/

// Test member access
/*TEST_F(EnsembleTest, MemberAccessIsValid) {
  Ensemble ensemble(config_, size_);

  // Test access to all members
  for (size_t i = 0; i < size_; ++i) {
    EXPECT_CALL(ensemble.getMember(i).backend(), getDimensions())
        .WillOnce(ReturnRef(dimensions_));  // Return reference to vector
    EXPECT_EQ(ensemble.getMember(i).getDimensions(), dimensions_);
  }

  // Test const access
  const auto& const_ensemble = ensemble;
  EXPECT_NO_THROW(const_ensemble.getMember(0));
}*/

// Test mean computation
/*TEST_F(EnsembleTest, ComputeMean) {
  Ensemble ensemble(config_, size_);

  {
    InSequence seq;

    // Set up expectations for mean state
    EXPECT_CALL(ensemble.getMean().backend(), reset()).Times(1);

    for (size_t i = 0; i < size_; ++i) {
      EXPECT_CALL(ensemble.getMean().backend(),
                  add(testing::Ref(ensemble.getMember(i).backend())))
          .Times(1);
    }

    EXPECT_CALL(ensemble.getMean().backend(),
                multiply(1.0 / static_cast<double>(size_)))
        .Times(1);
  }

  EXPECT_NO_THROW(ensemble.computeMean());
}*/

// Test out of bounds access
/*TEST_F(EnsembleTest, OutOfBoundsAccessThrows) {
  Ensemble ensemble(config_, size_);
  EXPECT_THROW(ensemble.getMember(size_), std::out_of_range);
  EXPECT_THROW(ensemble.getPerturbation(size_), std::out_of_range);
}*/

// Test perturbation computation
/*TEST_F(EnsembleTest, ComputePerturbationsCreatesValidPerturbations) {
  Ensemble ensemble(config_, size_);

  // Test perturbation computation
  EXPECT_NO_THROW(ensemble.computePerturbations());
  EXPECT_NO_THROW(ensemble.getPerturbation(0));
}*/

}  // namespace metada::tests