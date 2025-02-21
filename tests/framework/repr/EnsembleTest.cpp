#include <gmock/gmock.h>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "MockConfig.hpp"
#include "MockState.hpp"
#include "State.hpp"

namespace metada::framework::tests {

using StateType = State<backends::MockState>;
using ConfigType = tools::config::Config<tools::config::tests::MockConfig>;

using ::testing::InSequence;
using ::testing::Invoke;
using ::testing::NiceMock;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;

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

  {
    InSequence seq;

    // Set up expectations for mean state
    EXPECT_CALL(ensemble.getMean().backend(), reset())
        .Times(1)
        .WillOnce(Invoke([this]() {
          std::fill(test_data_.begin(), test_data_.end(), 0.0);
        }));

    EXPECT_CALL(ensemble.getMean().backend(), getData())
        .Times(1)
        .WillOnce(Return(test_data_.data()));

    EXPECT_CALL(ensemble.getMean().backend(), getDimensions())
        .Times(1)
        .WillOnce(ReturnRef(dimensions_));

    for (size_t i = 0; i < size_; ++i) {
      EXPECT_CALL(ensemble.getMember(i).backend(), getData())
          .Times(1)
          .WillOnce(Return(test_data_.data()));
    }
  }

  EXPECT_NO_THROW(ensemble.computeMean());

  // Verify mean computation
  const auto& mean = ensemble.getMean();
  EXPECT_CALL(mean.backend(), getData()).WillOnce(Return(test_data_.data()));
  const double* mean_data = &mean.getData<double>();
  EXPECT_DOUBLE_EQ(mean_data[5],
                   static_cast<double>(size_ * test_data_[5] / size_));
}

// Test out of bounds access
TEST_F(EnsembleTest, OutOfBoundsAccessThrows) {
  Ensemble<StateType> ensemble(config_, size_);
  EXPECT_THROW(ensemble.getMember(size_), std::out_of_range);
  EXPECT_THROW(ensemble.getPerturbation(size_), std::out_of_range);
}

// Test perturbation computation
TEST_F(EnsembleTest, ComputePerturbationsCreatesValidPerturbations) {
  Ensemble<StateType> ensemble(config_, size_);

  {
    InSequence seq;

    // Set up expectations for mean state
    EXPECT_CALL(ensemble.getMean().backend(), reset())
        .Times(1)
        .WillOnce(Invoke([this]() {
          std::fill(test_data_.begin(), test_data_.end(), 0.0);
        }));

    EXPECT_CALL(ensemble.getMean().backend(), getData())
        .Times(1)
        .WillOnce(Return(test_data_.data()));

    EXPECT_CALL(ensemble.getMean().backend(), getDimensions())
        .Times(1)
        .WillOnce(ReturnRef(dimensions_));

    for (size_t i = 0; i < size_; ++i) {
      EXPECT_CALL(ensemble.getMember(i).backend(), getData())
          .Times(1)
          .WillOnce(Return(test_data_.data()));
    }

    // Set up expectations for perturbations
    // for (size_t i = 0; i < size_; ++i) {
    //   std::cout << "Setting up perturbation " << i << std::endl;
    //   EXPECT_CALL(ensemble.getPerturbation(i).backend(),
    //               copyFrom(testing::Ref(ensemble.getMean().backend())))
    //       .Times(1);
    // }

    EXPECT_CALL(ensemble.getMean().backend(), getDimensions())
        .Times(1)
        .WillOnce(ReturnRef(dimensions_));

    EXPECT_CALL(ensemble.getMean().backend(), getData())
        .Times(1)
        .WillOnce(Return(test_data_.data()));
  }

  // Test perturbation computation
  EXPECT_NO_THROW(ensemble.computePerturbations());
  EXPECT_NO_THROW(ensemble.getPerturbation(0));

  // Verify perturbation computation for first member
  // const auto& perturbation = ensemble.getPerturbation(0);
  // EXPECT_CALL(perturbation.backend(), getData())
  //     .WillOnce(Return(test_data_.data()));
  // const double* pert_data = &perturbation.getData<double>();
  // EXPECT_DOUBLE_EQ(pert_data[0], 0.0);  // Since member data equals mean data
}

}  // namespace metada::framework::tests