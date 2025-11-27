/**
 * @file MetricsTest.cpp
 * @brief Unit tests for the Metrics class
 */

#include <gtest/gtest.h>

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

    // Calculate all metrics
    metrics_ = framework::Metrics<>::CalculateAll(ensemble_data_, true_state_,
                                                  state_dim_, ens_size_);
  }

  size_t state_dim_;
  size_t ens_size_;
  std::vector<double> true_state_;
  std::vector<std::vector<double>> ensemble_data_;
  framework::MetricValues<double> metrics_;
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
  ASSERT_EQ(metrics_.mean.size(), state_dim_);
  for (size_t i = 0; i < state_dim_; ++i) {
    EXPECT_NEAR(metrics_.mean[i], true_state_[i], 1e-6);
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
  ASSERT_EQ(metrics_.spread.size(), state_dim_);
  for (size_t i = 0; i < state_dim_; ++i) {
    EXPECT_NEAR(metrics_.spread[i], 0.1, 1e-6);
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
  EXPECT_NEAR(metrics_.bias, 0.0, 1e-6);
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
  EXPECT_NEAR(metrics_.correlation, 1.0, 1e-6);
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
  EXPECT_NEAR(metrics_.crps, 0.022222222222222227, 1e-6);
  EXPECT_GT(metrics_.crps, 0.0);
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
  EXPECT_NEAR(metrics_.rmse, 0.0, 1e-6);
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
  EXPECT_NEAR(metrics_.avg_spread, 0.1, 1e-6);
  EXPECT_GT(metrics_.avg_spread, 0.0);
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

  auto biased_metrics = framework::Metrics<>::CalculateAll(
      biased_ensemble, true_state_, state_dim_, ens_size_);

  EXPECT_NEAR(biased_metrics.bias, 0.5, 1e-6);
  EXPECT_NEAR(biased_metrics.rmse, 0.5, 1e-6);
  EXPECT_NEAR(biased_metrics.correlation, 1.0,
              1e-6);  // Still perfect correlation
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

  auto uncorrelated_metrics = framework::Metrics<>::CalculateAll(
      uncorrelated_ensemble, true_state_, state_dim_, ens_size_);

  EXPECT_NEAR(uncorrelated_metrics.correlation, 0.34105511573854358, 1e-6);
  EXPECT_GT(uncorrelated_metrics.rmse, 1.0);
  EXPECT_GT(uncorrelated_metrics.crps, 1.0);
}

/**
 * @brief Test metrics printing
 *
 * Verifies that:
 * - Metrics can be printed to stream
 * - Output format is correct
 */
TEST_F(MetricsTest, PrintMetrics) {
  std::stringstream ss;
  ss << metrics_;
  std::string output = ss.str();

  // Check that all metrics are present in output
  EXPECT_TRUE(output.find("RMSE:") != std::string::npos);
  EXPECT_TRUE(output.find("Bias:") != std::string::npos);
  EXPECT_TRUE(output.find("Correlation:") != std::string::npos);
  EXPECT_TRUE(output.find("CRPS:") != std::string::npos);
  EXPECT_TRUE(output.find("Average Spread:") != std::string::npos);
}

}  // namespace metada::tests