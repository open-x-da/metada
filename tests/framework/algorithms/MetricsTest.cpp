/**
 * @file MetricsTest.cpp
 * @brief Unit tests for the Metrics class
 */

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "Metrics.hpp"

namespace metada::tests {

using ::testing::Test;

/**
 * @brief Test fixture for Metrics class tests
 *
 * Sets up test data including:
 * - True state vector
 * - Ensemble data with known properties
 * - Expected metric values
 */
class MetricsTest : public Test {
 protected:
  void SetUp() override {
    // Setup test dimensions
    state_dim_ = 4;
    ens_size_ = 3;

    // Create true state: [1.0, 2.0, 3.0, 4.0]
    true_state_ = {1.0, 2.0, 3.0, 4.0};

    // Create ensemble data with known properties:
    // Member 1: true state + [0.1, 0.1, 0.1, 0.1]
    // Member 2: true state + [-0.1, -0.1, -0.1, -0.1]
    // Member 3: true state + [0.0, 0.0, 0.0, 0.0]
    ensemble_data_ = {
        {1.1, 2.1, 3.1, 4.1},  // Member 1
        {0.9, 1.9, 2.9, 3.9},  // Member 2
        {1.0, 2.0, 3.0, 4.0}   // Member 3
    };

    // Expected values
    expected_mean_ = {1.0, 2.0, 3.0, 4.0};  // Should match true state
    expected_spread_ = {0.1, 0.1, 0.1,
                        0.1};     // Standard deviation of perturbations
    expected_rmse_ = 0.0;         // Mean matches true state
    expected_bias_ = 0.0;         // No bias
    expected_correlation_ = 1.0;  // Perfect correlation
    expected_crps_ = 0.022222222222222227;  // CRPS for this ensemble
  }

  size_t state_dim_;
  size_t ens_size_;
  std::vector<double> true_state_;
  std::vector<std::vector<double>> ensemble_data_;
  std::vector<double> expected_mean_;
  std::vector<double> expected_spread_;
  double expected_rmse_;
  double expected_bias_;
  double expected_correlation_;
  double expected_crps_;
};

/**
 * @brief Test ensemble mean calculation
 *
 * Verifies that:
 * - Mean is calculated correctly
 * - Mean matches expected values
 * - Mean has correct dimensions
 */
TEST_F(MetricsTest, CalculateEnsembleMean) {
  auto mean = framework::Metrics<>::calculateEnsembleMean(
      ensemble_data_, state_dim_, ens_size_);

  ASSERT_EQ(mean.size(), state_dim_);
  for (size_t i = 0; i < state_dim_; ++i) {
    EXPECT_NEAR(mean[i], expected_mean_[i], 1e-6);
  }
}

/**
 * @brief Test ensemble spread calculation
 *
 * Verifies that:
 * - Spread is calculated correctly
 * - Spread matches expected values
 * - Spread has correct dimensions
 */
TEST_F(MetricsTest, CalculateEnsembleSpread) {
  auto mean = framework::Metrics<>::calculateEnsembleMean(
      ensemble_data_, state_dim_, ens_size_);
  auto spread = framework::Metrics<>::calculateEnsembleSpread(
      ensemble_data_, mean, state_dim_, ens_size_);

  ASSERT_EQ(spread.size(), state_dim_);
  for (size_t i = 0; i < state_dim_; ++i) {
    EXPECT_NEAR(spread[i], expected_spread_[i], 1e-6);
  }
}

/**
 * @brief Test bias calculation
 *
 * Verifies that:
 * - Bias is calculated correctly
 * - Bias matches expected value
 * - Bias is zero for unbiased ensemble
 */
TEST_F(MetricsTest, CalculateBias) {
  auto mean = framework::Metrics<>::calculateEnsembleMean(
      ensemble_data_, state_dim_, ens_size_);
  double bias =
      framework::Metrics<>::calculateBias(mean, true_state_, state_dim_);

  EXPECT_NEAR(bias, expected_bias_, 1e-6);
}

/**
 * @brief Test correlation calculation
 *
 * Verifies that:
 * - Correlation is calculated correctly
 * - Correlation matches expected value
 * - Correlation is 1.0 for perfect correlation
 */
TEST_F(MetricsTest, CalculateCorrelation) {
  auto mean = framework::Metrics<>::calculateEnsembleMean(
      ensemble_data_, state_dim_, ens_size_);
  double correlation =
      framework::Metrics<>::calculateCorrelation(mean, true_state_, state_dim_);

  EXPECT_NEAR(correlation, expected_correlation_, 1e-6);
}

/**
 * @brief Test CRPS calculation
 *
 * Verifies that:
 * - CRPS is calculated correctly
 * - CRPS matches expected value
 * - CRPS is positive
 */
TEST_F(MetricsTest, CalculateCRPS) {
  double crps = framework::Metrics<>::calculateCRPS(ensemble_data_, true_state_,
                                                    state_dim_, ens_size_);

  EXPECT_NEAR(crps, expected_crps_, 1e-6);
  EXPECT_GT(crps, 0.0);
}

/**
 * @brief Test RMSE calculation
 *
 * Verifies that:
 * - RMSE is calculated correctly
 * - RMSE matches expected value
 * - RMSE is zero for perfect match
 */
TEST_F(MetricsTest, CalculateRMSE) {
  auto mean = framework::Metrics<>::calculateEnsembleMean(
      ensemble_data_, state_dim_, ens_size_);
  double rmse =
      framework::Metrics<>::calculateRMSE(mean, true_state_, state_dim_);

  EXPECT_NEAR(rmse, expected_rmse_, 1e-6);
}

/**
 * @brief Test average spread calculation
 *
 * Verifies that:
 * - Average spread is calculated correctly
 * - Average spread matches expected value
 * - Average spread is positive
 */
TEST_F(MetricsTest, CalculateAverageSpread) {
  auto mean = framework::Metrics<>::calculateEnsembleMean(
      ensemble_data_, state_dim_, ens_size_);
  auto spread = framework::Metrics<>::calculateEnsembleSpread(
      ensemble_data_, mean, state_dim_, ens_size_);
  double avg_spread =
      framework::Metrics<>::calculateAverageSpread(spread, state_dim_);

  EXPECT_DOUBLE_EQ(avg_spread, 0.1);  // Average of expected_spread_
  EXPECT_GT(avg_spread, 0.0);
}

/**
 * @brief Test metrics with biased ensemble
 *
 * Verifies that:
 * - Bias is detected correctly
 * - RMSE increases with bias
 * - Correlation decreases with bias
 */
TEST_F(MetricsTest, BiasedEnsemble) {
  // Create biased ensemble by adding 0.5 to all members
  std::vector<std::vector<double>> biased_ensemble = ensemble_data_;
  for (auto& member : biased_ensemble) {
    for (auto& val : member) {
      val += 0.5;
    }
  }

  auto mean = framework::Metrics<>::calculateEnsembleMean(
      biased_ensemble, state_dim_, ens_size_);
  double bias =
      framework::Metrics<>::calculateBias(mean, true_state_, state_dim_);
  double rmse =
      framework::Metrics<>::calculateRMSE(mean, true_state_, state_dim_);
  double correlation =
      framework::Metrics<>::calculateCorrelation(mean, true_state_, state_dim_);

  EXPECT_NEAR(bias, 0.5, 1e-6);
  EXPECT_NEAR(rmse, 0.5, 1e-6);
  EXPECT_NEAR(correlation, 1.0, 1e-6);  // Still perfect correlation
}

/**
 * @brief Test metrics with uncorrelated ensemble
 *
 * Verifies that:
 * - Correlation is zero for uncorrelated data
 * - RMSE is large for uncorrelated data
 * - CRPS is large for uncorrelated data
 */
TEST_F(MetricsTest, UncorrelatedEnsemble) {
  // Create uncorrelated ensemble
  std::vector<std::vector<double>> uncorrelated_ensemble = {
      {10.0, 2.0, 300.0, 4.0},  // Reversed true state
      {4.0, 3.0, 20.0, 1.0},    // Shuffled true state
      {2.0, 30.0, 1.0, 43.0}    // Another shuffle
  };

  auto mean = framework::Metrics<>::calculateEnsembleMean(
      uncorrelated_ensemble, state_dim_, ens_size_);
  double correlation =
      framework::Metrics<>::calculateCorrelation(mean, true_state_, state_dim_);
  double rmse =
      framework::Metrics<>::calculateRMSE(mean, true_state_, state_dim_);
  double crps = framework::Metrics<>::calculateCRPS(
      uncorrelated_ensemble, true_state_, state_dim_, ens_size_);

  EXPECT_NEAR(correlation, 0.34105511573854358, 1e-6);
  EXPECT_GT(rmse, 1.0);
  EXPECT_GT(crps, 1.0);
}

}  // namespace metada::tests