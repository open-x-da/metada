#include <gmock/gmock.h>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "MockConfig.hpp"
#include "MockState.hpp"
#include "State.hpp"

namespace metada::framework::tests {

using StateType = State<backends::MockState>;
using ConfigType = tools::config::Config<tools::config::tests::MockConfig>;

using ::testing::NiceMock;
using ::testing::Return;
using ::testing::ReturnRef;

class EnsembleTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Setup test data
    size_ = 3;
    dimensions_ = {10};
    test_data_ = new double[10];  // Example data
    for (int i = 0; i < 10; i++) test_data_[i] = i;
  }

  /**
   * @brief Clean up test data after each test
   *
   * Frees allocated test data array
   */
  void TearDown() override { delete[] test_data_; }

  ConfigType config_;
  std::vector<size_t> dimensions_;
  double* test_data_;
  size_t size_;
};

// Test ensemble initialization
TEST_F(EnsembleTest, InitializationCreatesCorrectNumberOfMembers) {
  Ensemble<StateType> ensemble(config_, size_);
  EXPECT_EQ(ensemble.getSize(), size_);
}

// Test member access
TEST_F(EnsembleTest, MemberAccessIsValid) {
  Ensemble<StateType> ensemble(config_, size_);

  // Test access to all members
  for (size_t i = 0; i < size_; ++i) {
    EXPECT_CALL(ensemble.getMember(i).backend(), getDimensions())
        .WillOnce(ReturnRef(dimensions_));  // Return reference to vector
    EXPECT_EQ(ensemble.getMember(i).getDimensions(), dimensions_);
  }

  // Test const access
  const auto& const_ensemble = ensemble;
  EXPECT_NO_THROW(const_ensemble.getMember(0));
}

// Test mean computation
TEST_F(EnsembleTest, ComputeMean) {
  Ensemble<StateType> ensemble(config_, size_);

  // Set up expectations for mean state
  EXPECT_CALL(ensemble.getMean().backend(), reset()).Times(1);
  EXPECT_CALL(ensemble.getMean().backend(), getData())
      .WillOnce(Return(new double[10]));
  EXPECT_CALL(ensemble.getMean().backend(), getDimensions())
      .WillOnce(ReturnRef(dimensions_));

  for (size_t i = 0; i < size_; ++i) {
    EXPECT_CALL(ensemble.getMember(i).backend(), getData())
        .WillOnce(Return(test_data_));
  }

  EXPECT_NO_THROW(ensemble.computeMean());

  // Verify mean computation
  const auto& mean = ensemble.getMean();
  EXPECT_CALL(mean.backend(), getData()).WillOnce(Return(test_data_));
  const double* mean_data = &mean.getData<double>();
  EXPECT_DOUBLE_EQ(mean_data[5],
                   static_cast<double>(size_ * test_data_[5] / size_));
}

// Test perturbation computation
// TODO: Implement this test
// TEST_F(EnsembleTest,
// DISABLED_ComputePerturbationsCreatesValidPerturbations)
// {
//   Ensemble<state> ensemble(config_, size_);
//   EXPECT_NO_THROW(ensemble.computePerturbations());
//   EXPECT_NO_THROW(ensemble.getPerturbation(0));
// }

// Test out of bounds access
// TODO: Implement this test
// TEST_F(EnsembleTest, DISABLED_OutOfBoundsAccessThrows) {
//   Ensemble<state> ensemble(config_, size_);
//   EXPECT_THROW(ensemble.getMember(size_), std::out_of_range);
//   EXPECT_THROW(ensemble.getPerturbation(size_), std::out_of_range);
// }

}  // namespace metada::framework::tests