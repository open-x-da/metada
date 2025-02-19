#include <gmock/gmock.h>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "MockConfig.hpp"
#include "MockState.hpp"
#include "State.hpp"

namespace metada::framework::tests {

using StateType = State<backends::MockState>;
using ConfigType = tools::config::Config<tools::config::tests::MockConfig>;

using ::testing::ReturnRef;

class EnsembleTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Setup test data
    size_ = 3;
  }

  ConfigType config_;
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
  std::vector<size_t> dims{20};  // Create vector with dimension

  // Test access to all members
  for (size_t i = 0; i < size_; ++i) {
    EXPECT_CALL(ensemble.getMember(i).backend(), getDimensions())
        .WillOnce(ReturnRef(dims));  // Return reference to vector
    EXPECT_EQ(ensemble.getMember(i).getDimensions(), dims);
  }

  // Test const access
  const auto& const_ensemble = ensemble;
  EXPECT_NO_THROW(const_ensemble.getMember(0));
}

// Test mean computation
// TEST_F(EnsembleTest, ComputeMeanUpdatesEnsembleMean) {
//   Ensemble<state> ensemble(config_, size_);
//   ensemble.computeMean();
//   //  EXPECT_NO_THROW(ensemble.computeMean());
//   //  EXPECT_NO_THROW(ensemble.getMean());
// }

// Test perturbation computation
// TODO: Implement this test
// TEST_F(EnsembleTest, DISABLED_ComputePerturbationsCreatesValidPerturbations)
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