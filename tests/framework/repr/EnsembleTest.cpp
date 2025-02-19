#include "Config.hpp"
#include "Ensemble.hpp"
#include "MockConfig.hpp"
#include "MockState.hpp"

namespace metada::framework::tests {

class EnsembleTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Setup test data
    ensemble_size_ = 3;
    transform_matrix_ = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    localization_weights_ = {1.0, 0.8, 0.5};
  }

  tools::config::Config<tools::config::tests::MockConfig> config_;
  size_t ensemble_size_;
  std::vector<std::vector<double>> transform_matrix_;
  std::vector<double> localization_weights_;
};

// Test ensemble initialization
TEST_F(EnsembleTest, InitializationCreatesCorrectNumberOfMembers) {
  Ensemble<metada::backends::MockState> ensemble(config_, ensemble_size_);
  EXPECT_EQ(ensemble.getSize(), ensemble_size_);
}

// Test member access
TEST_F(EnsembleTest, MemberAccessIsValid) {
  Ensemble<metada::backends::MockState> ensemble(config_, ensemble_size_);

  // Test access to all members
  for (size_t i = 0; i < ensemble_size_; ++i) {
    EXPECT_NO_THROW(ensemble.getMember(i));
  }

  // Test const access
  const auto& const_ensemble = ensemble;
  EXPECT_NO_THROW(const_ensemble.getMember(0));
}

// Test mean computation
TEST_F(EnsembleTest, ComputeMeanUpdatesEnsembleMean) {
  Ensemble<metada::backends::MockState> ensemble(config_, ensemble_size_);
  ensemble.computeMean();
  //  EXPECT_NO_THROW(ensemble.computeMean());
  //  EXPECT_NO_THROW(ensemble.getMean());
}

// Test perturbation computation
// TODO: Implement this test
// TEST_F(EnsembleTest, DISABLED_ComputePerturbationsCreatesValidPerturbations)
// {
//   Ensemble<metada::backends::repr::MockState> ensemble(config_,
//   ensemble_size_); EXPECT_NO_THROW(ensemble.computePerturbations());
//   EXPECT_NO_THROW(ensemble.getPerturbation(0));
// }

// Test inflation
// TODO: Implement this test
// TEST_F(EnsembleTest, DISABLED_InflationModifiesEnsembleMembers) {
//   Ensemble<metada::backends::repr::MockState> ensemble(config_,
//   ensemble_size_); double inflation_factor = 1.1;
//   EXPECT_NO_THROW(ensemble.inflate(inflation_factor));
// }

// Test transformation
// TODO: Implement this test
// TEST_F(EnsembleTest, DISABLED_TransformationAppliesMatrix) {
//   Ensemble<metada::backends::repr::MockState> ensemble(config_,
//   ensemble_size_); EXPECT_NO_THROW(ensemble.transform(transform_matrix_));
// }

// Test localization
// TODO: Implement this test
// TEST_F(EnsembleTest, DISABLED_LocalizationAppliesWeights) {
//   Ensemble<metada::backends::repr::MockState> ensemble(config_,
//   ensemble_size_);
//   EXPECT_NO_THROW(ensemble.localizeCovariance(localization_weights_));
// }

// Test out of bounds access
// TODO: Implement this test
// TEST_F(EnsembleTest, DISABLED_OutOfBoundsAccessThrows) {
//   Ensemble<metada::backends::repr::MockState> ensemble(config_,
//   ensemble_size_); EXPECT_THROW(ensemble.getMember(ensemble_size_),
//   std::out_of_range); EXPECT_THROW(ensemble.getPerturbation(ensemble_size_),
//   std::out_of_range);
// }

}  // namespace metada::framework::tests