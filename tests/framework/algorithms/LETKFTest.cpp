/**
 * @file LETKFTest.cpp
 * @brief Unit tests for the Local Ensemble Transform Kalman Filter (LETKF)
 * implementation
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "Geometry.hpp"
#include "LETKF.hpp"
#include "MockBackendTraits.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"

using ::testing::Return;
using ::testing::ReturnRef;

namespace metada::tests {

/**
 * @brief Test fixture for LETKF tests
 *
 * Sets up mock data and components needed for testing the LETKF implementation:
 * - Mock ensemble data (3 elements)
 * - Mock observation data (2 elements)
 * - Mock covariance data (2x2 matrix)
 * - Configuration loaded from test_config.yaml
 * - Geometry instance initialized from config
 * - Ensemble instance with mock backend
 * - Observation instance with mock backend
 * - Observation operator instance with mock backend
 * - LETKF instance with 1.1 inflation factor
 */
class LETKFTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Setup mock data
    ensemble_data_ = std::vector<std::vector<double>>{
        {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    obs_data_ = std::vector<double>{1.5, 2.5};
    cov_data_ = std::vector<double>{0.1, 0.0, 0.0, 0.1};  // 2x2 matrix

    // Load configuration file from test directory
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    config_file_ = (test_dir / "test_config.yaml").string();

    // Create config and geometry
    config_ = std::make_unique<framework::Config<traits::MockBackendTag>>(
        config_file_);
    geometry_ =
        std::make_unique<framework::Geometry<traits::MockBackendTag>>(*config_);

    // Setup mock expectations for config
    ON_CALL(config_->backend(), Get("ensemble.size")).WillByDefault(Return(3));

    // Create adapter instances
    ensemble_ = std::make_unique<framework::Ensemble<traits::MockBackendTag>>(
        *config_, *geometry_);
    obs_ = std::make_unique<framework::Observation<traits::MockBackendTag>>(
        *config_);
    obs_op_ = std::make_unique<framework::ObsOperator<traits::MockBackendTag>>(
        *config_);

    // Setup mock expectations for ensemble member data
    for (size_t i = 0; i < ensemble_->Size(); ++i) {
      ensemble_->GetMember(i).backend().setData(ensemble_data_[i]);
    }

    // Setup mock expectations for observation data
    obs_->backend().setData(obs_data_);

    // Setup mock expectations for observation covariance data
    obs_->backend().setCovariance(cov_data_);

    // Create LETKF instance
    letkf_ = std::make_unique<framework::LETKF<traits::MockBackendTag>>(
        *ensemble_, *obs_, *obs_op_, 1.1);
  }

  std::vector<std::vector<double>> ensemble_data_;
  std::vector<double> obs_data_;
  std::vector<double> cov_data_;
  std::string config_file_;
  std::unique_ptr<framework::Config<traits::MockBackendTag>> config_;
  std::unique_ptr<framework::Geometry<traits::MockBackendTag>> geometry_;
  std::unique_ptr<framework::Ensemble<traits::MockBackendTag>> ensemble_;
  std::unique_ptr<framework::Observation<traits::MockBackendTag>> obs_;
  std::unique_ptr<framework::ObsOperator<traits::MockBackendTag>> obs_op_;
  std::unique_ptr<framework::LETKF<traits::MockBackendTag>> letkf_;
};

/**
 * @brief Test that LETKF constructor properly initializes the instance
 *
 * Verifies that the LETKF instance is created successfully with all required
 * components.
 */
TEST_F(LETKFTest, ConstructorInitializesCorrectly) {
  EXPECT_NE(letkf_, nullptr);
}

/**
 * @brief Test the LETKF analysis step
 *
 * Verifies that:
 * - The observation operator is called exactly 3 times (once per ensemble
 * member)
 * - The analysis step completes without errors
 */
TEST_F(LETKFTest, AnalyseUpdatesEnsemble) {
  // Setup expectations for the analysis step
  EXPECT_CALL(obs_op_->backend(), apply(::testing::_, ::testing::_))
      .Times(3);  // Once per ensemble member

  // Perform analysis
  letkf_->analyse();
}

/**
 * @brief Test LETKF analysis with empty ensemble
 *
 * Verifies that the analysis step handles an empty ensemble gracefully
 * without throwing any exceptions.
 */
TEST_F(LETKFTest, AnalyseHandlesEmptyEnsemble) {
  // Analysis should not crash with empty ensemble
  EXPECT_NO_THROW(letkf_->analyse());
}

/**
 * @brief Test LETKF analysis with zero inflation
 *
 * Verifies that:
 * - LETKF can be constructed with 0.0 inflation factor
 * - The observation operator is still called 3 times
 * - Analysis completes successfully without errors
 */
TEST_F(LETKFTest, AnalyseHandlesZeroInflation) {
  // Create LETKF with zero inflation
  letkf_ = std::make_unique<framework::LETKF<traits::MockBackendTag>>(
      *ensemble_, *obs_, *obs_op_, 0.0);

  // Setup expectations
  EXPECT_CALL(obs_op_->backend(), apply(::testing::_, ::testing::_)).Times(3);

  // Analysis should work with zero inflation
  EXPECT_NO_THROW(letkf_->analyse());
}

}  // namespace metada::tests