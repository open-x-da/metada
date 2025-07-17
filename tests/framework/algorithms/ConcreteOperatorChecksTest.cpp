/**
 * @file ConcreteOperatorChecksTest.cpp
 * @brief Concrete unit tests for OperatorChecks algorithms using actual data
 *
 * Tests the mathematical correctness of:
 * - Tangent linear and adjoint consistency checks
 * - Tangent linear implementation verification using Taylor expansion
 * - Cost function gradient checks
 *
 * Uses lite backends with known mathematical properties for verification.
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "BackgroundErrorCovariance.hpp"
#include "Config.hpp"
#include "CostFunction.hpp"
#include "Ensemble.hpp"
#include "Geometry.hpp"
#include "Increment.hpp"
#include "Logger.hpp"
#include "Model.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "OperatorChecks.hpp"
#include "State.hpp"

// Include lite backend implementations
#include "LiteBackendTraits.hpp"

using ::testing::DoubleNear;
using ::testing::Each;
using ::testing::Ge;
using ::testing::Le;

namespace metada::tests {

/**
 * @brief Test fixture for concrete OperatorChecks tests
 */
class ConcreteOperatorChecksTest : public ::testing::Test {
 protected:
  void SetUp() override {
    try {
      // Load configuration
      auto test_dir = std::filesystem::path(__FILE__).parent_path();
      config_file_ = (test_dir / "test_config.yaml").string();
      config_ = std::make_unique<framework::Config<traits::LiteBackendTag>>(
          config_file_);

      // Initialize logger
      framework::Logger<traits::LiteBackendTag>::Init(
          config_->GetSubsection("logger"));

      // Create basic components
      geometry_ = std::make_unique<framework::Geometry<traits::LiteBackendTag>>(
          *config_);
      state_ = std::make_unique<framework::State<traits::LiteBackendTag>>(
          *config_, *geometry_);
      obs_ = std::make_unique<framework::Observation<traits::LiteBackendTag>>(
          *config_);

      // Set up test data
      std::vector<double> test_data = {1.0, 2.0, 3.0};
      state_->backend().setData(test_data);

      // Set up observation data
      std::vector<double> obs_data = {2.0, 3.5};  // H * [1, 2, 3] = [2.0, 3.5]
      obs_->backend().setObservations(obs_data);
      obs_->backend().setCovariance({1.0, 1.0});

    } catch (const std::exception& e) {
      std::cerr << "Setup failed: " << e.what() << std::endl;
    }
  }

  void TearDown() override {
    state_.reset();
    obs_.reset();
    geometry_.reset();
    config_.reset();
    framework::Logger<traits::LiteBackendTag>::Reset();
  }

  std::string config_file_;
  std::unique_ptr<framework::Config<traits::LiteBackendTag>> config_;
  std::unique_ptr<framework::Geometry<traits::LiteBackendTag>> geometry_;
  std::unique_ptr<framework::State<traits::LiteBackendTag>> state_;
  std::unique_ptr<framework::Observation<traits::LiteBackendTag>> obs_;
};

/**
 * @brief Test cost function gradient using finite differences
 */
TEST_F(ConcreteOperatorChecksTest, CostFunctionGradientCheck) {
  if (!state_ || !obs_ || !config_) {
    GTEST_SKIP() << "Setup failed - skipping test";
  }

  // Create background state
  auto background =
      framework::State<traits::LiteBackendTag>(*config_, *geometry_);
  background.backend().setData({0.0, 0.0, 0.0});

  // Create background error covariance
  auto bg_error_cov =
      framework::BackgroundErrorCovariance<traits::LiteBackendTag>(*config_);

  // Create model
  auto model = framework::Model<traits::LiteBackendTag>(*config_);

  // Create observations and obs operators
  std::vector<framework::Observation<traits::LiteBackendTag>> observations;
  std::vector<framework::ObsOperator<traits::LiteBackendTag>> obs_operators;

  auto obs1 = framework::Observation<traits::LiteBackendTag>(*config_);
  auto obs_op1 = framework::ObsOperator<traits::LiteBackendTag>(*config_);

  obs1.backend().setObservations({2.0, 3.5});
  obs1.backend().setCovariance({1.0, 1.0});

  observations.push_back(std::move(obs1));
  obs_operators.push_back(std::move(obs_op1));

  // Create cost function
  auto cost_func = framework::CostFunction<traits::LiteBackendTag>(
      *config_, background, observations, obs_operators, model, bg_error_cov);

  // Test state
  auto test_state =
      framework::State<traits::LiteBackendTag>(*config_, *geometry_);
  test_state.backend().setData({1.0, 2.0, 3.0});

  // Compute gradient
  auto gradient =
      framework::Increment<traits::LiteBackendTag>::createFromEntity(
          test_state);
  cost_func.gradient(test_state, gradient);

  // Compute cost at test state
  double cost_at_x = cost_func.evaluate(test_state);

  // Test gradient using finite differences
  double epsilon = 1e-6;
  for (size_t i = 0; i < 3; ++i) {
    // Create perturbed state
    auto perturbed_state = test_state.clone();
    auto* data = perturbed_state.template getDataPtr<double>();
    data[i] += epsilon;

    // Compute cost at perturbed state
    double cost_at_x_plus_dx = cost_func.evaluate(perturbed_state);

    // Finite difference approximation of gradient
    double finite_diff_grad = (cost_at_x_plus_dx - cost_at_x) / epsilon;

    // Compare with computed gradient
    auto* grad_data = gradient.state().template getDataPtr<double>();
    EXPECT_NEAR(grad_data[i], finite_diff_grad, 1e-5);
  }
}

/**
 * @brief Test observation operator TL/AD consistency check
 */
TEST_F(ConcreteOperatorChecksTest, ObsOperatorTLADCheck) {
  if (!state_ || !obs_) {
    GTEST_SKIP() << "Setup failed - skipping test";
  }

  // Create observation operators
  std::vector<framework::ObsOperator<traits::LiteBackendTag>> obs_operators;
  auto obs_op = framework::ObsOperator<traits::LiteBackendTag>(*config_);
  obs_operators.push_back(std::move(obs_op));

  // Create observations
  std::vector<framework::Observation<traits::LiteBackendTag>> observations;
  auto obs1 = framework::Observation<traits::LiteBackendTag>(*config_);
  obs1.backend().setObservations({2.0, 3.5});
  obs1.backend().setCovariance({1.0, 1.0});
  observations.push_back(std::move(obs1));

  // Test TL/AD consistency check
  bool result = framework::checkObsOperatorTLAD<traits::LiteBackendTag>(
      obs_operators, *state_, observations, 1e-10);

  // Should pass for linear observation operator
  EXPECT_TRUE(result);
}

/**
 * @brief Test cost function gradient consistency using framework check
 */
TEST_F(ConcreteOperatorChecksTest, CostFunctionGradientConsistency) {
  if (!state_ || !obs_ || !config_) {
    GTEST_SKIP() << "Setup failed - skipping test";
  }

  // Create background state
  auto background =
      framework::State<traits::LiteBackendTag>(*config_, *geometry_);
  background.backend().setData({0.0, 0.0, 0.0});

  // Create background error covariance
  auto bg_error_cov =
      framework::BackgroundErrorCovariance<traits::LiteBackendTag>(*config_);

  // Create model
  auto model = framework::Model<traits::LiteBackendTag>(*config_);

  // Create observations and obs operators
  std::vector<framework::Observation<traits::LiteBackendTag>> observations;
  std::vector<framework::ObsOperator<traits::LiteBackendTag>> obs_operators;

  auto obs1 = framework::Observation<traits::LiteBackendTag>(*config_);
  auto obs_op1 = framework::ObsOperator<traits::LiteBackendTag>(*config_);

  obs1.backend().setObservations({2.0, 3.5});
  obs1.backend().setCovariance({1.0, 1.0});

  observations.push_back(std::move(obs1));
  obs_operators.push_back(std::move(obs_op1));

  // Create cost function
  auto cost_func = framework::CostFunction<traits::LiteBackendTag>(
      *config_, background, observations, obs_operators, model, bg_error_cov);

  // Test state
  auto test_state =
      framework::State<traits::LiteBackendTag>(*config_, *geometry_);
  test_state.backend().setData({1.0, 2.0, 3.0});

  // Use the framework's gradient check function
  bool result = framework::checkCostFunctionGradient<traits::LiteBackendTag>(
      cost_func, test_state, 1e-6, 1e-6);
  EXPECT_TRUE(result);
}

/**
 * @brief Test observation operator linearity
 */
TEST_F(ConcreteOperatorChecksTest, ObsOperatorLinearity) {
  if (!state_ || !obs_) {
    GTEST_SKIP() << "Setup failed - skipping test";
  }

  // Create observation operator
  auto obs_op = framework::ObsOperator<traits::LiteBackendTag>(*config_);

  // Create two states
  auto state1 = framework::State<traits::LiteBackendTag>(*config_, *geometry_);
  auto state2 = framework::State<traits::LiteBackendTag>(*config_, *geometry_);

  state1.backend().setData({1.0, 2.0, 3.0});
  state2.backend().setData({4.0, 5.0, 6.0});

  // Test linearity: H(x1 + x2) = H(x1) + H(x2)
  auto state_sum = state1 + state2;
  auto obs_sum = obs_op.apply(state_sum, *obs_);

  auto obs1 = obs_op.apply(state1, *obs_);
  auto obs2 = obs_op.apply(state2, *obs_);

  // Add observation results
  std::vector<double> obs1_plus_obs2(obs1.size());
  for (size_t i = 0; i < obs1.size(); ++i) {
    obs1_plus_obs2[i] = obs1[i] + obs2[i];
  }

  // Check linearity
  for (size_t i = 0; i < obs_sum.size(); ++i) {
    EXPECT_NEAR(obs_sum[i], obs1_plus_obs2[i], 1e-10);
  }
}

/**
 * @brief Test cost function gradient using multiple random directions
 */
TEST_F(ConcreteOperatorChecksTest, CostFunctionGradientMultipleDirections) {
  if (!state_ || !obs_ || !config_) {
    GTEST_SKIP() << "Setup failed - skipping test";
  }

  // Create background state
  auto background =
      framework::State<traits::LiteBackendTag>(*config_, *geometry_);
  background.backend().setData({0.0, 0.0, 0.0});

  // Create background error covariance
  auto bg_error_cov =
      framework::BackgroundErrorCovariance<traits::LiteBackendTag>(*config_);

  // Create model
  auto model = framework::Model<traits::LiteBackendTag>(*config_);

  // Create observations and obs operators
  std::vector<framework::Observation<traits::LiteBackendTag>> observations;
  std::vector<framework::ObsOperator<traits::LiteBackendTag>> obs_operators;

  auto obs1 = framework::Observation<traits::LiteBackendTag>(*config_);
  auto obs_op1 = framework::ObsOperator<traits::LiteBackendTag>(*config_);

  obs1.backend().setObservations({2.0, 3.5});
  obs1.backend().setCovariance({1.0, 1.0});

  observations.push_back(std::move(obs1));
  obs_operators.push_back(std::move(obs_op1));

  // Create cost function
  auto cost_func = framework::CostFunction<traits::LiteBackendTag>(
      *config_, background, observations, obs_operators, model, bg_error_cov);

  // Test state
  auto test_state =
      framework::State<traits::LiteBackendTag>(*config_, *geometry_);
  test_state.backend().setData({1.0, 2.0, 3.0});

  // Run the multiple directions gradient check
  bool result = framework::checkCostFunctionGradientMultipleDirections<
      traits::LiteBackendTag>(cost_func, test_state, 10, 1e-4);
  EXPECT_TRUE(result);
}

/**
 * @brief Test cost function gradient using unit vector directions
 */
TEST_F(ConcreteOperatorChecksTest, CostFunctionGradientUnitDirections) {
  if (!state_ || !obs_ || !config_) {
    GTEST_SKIP() << "Setup failed - skipping test";
  }

  // Create background state
  auto background =
      framework::State<traits::LiteBackendTag>(*config_, *geometry_);
  background.backend().setData({0.0, 0.0, 0.0});

  // Create background error covariance
  auto bg_error_cov =
      framework::BackgroundErrorCovariance<traits::LiteBackendTag>(*config_);

  // Create model
  auto model = framework::Model<traits::LiteBackendTag>(*config_);

  // Create observations and obs operators
  std::vector<framework::Observation<traits::LiteBackendTag>> observations;
  std::vector<framework::ObsOperator<traits::LiteBackendTag>> obs_operators;

  auto obs1 = framework::Observation<traits::LiteBackendTag>(*config_);
  auto obs_op1 = framework::ObsOperator<traits::LiteBackendTag>(*config_);

  obs1.backend().setObservations({2.0, 3.5});
  obs1.backend().setCovariance({1.0, 1.0});

  observations.push_back(std::move(obs1));
  obs_operators.push_back(std::move(obs_op1));

  // Create cost function
  auto cost_func = framework::CostFunction<traits::LiteBackendTag>(
      *config_, background, observations, obs_operators, model, bg_error_cov);

  // Test state
  auto test_state =
      framework::State<traits::LiteBackendTag>(*config_, *geometry_);
  test_state.backend().setData({1.0, 2.0, 3.0});

  // Run the unit directions gradient check
  bool result = framework::checkCostFunctionGradientUnitDirections<
      traits::LiteBackendTag>(cost_func, test_state, 1e-4);
  EXPECT_TRUE(result);
}

/**
 * @brief Test observation operator tangent linear implementation using direct
 * call to checkObsOperatorTangentLinear (single-epsilon mode via unified
 * interface)
 */
TEST_F(ConcreteOperatorChecksTest, ObsOperatorTangentLinearCheck) {
  if (!state_ || !obs_) {
    GTEST_SKIP() << "Setup failed - skipping test";
  }

  // Create observation operators
  std::vector<framework::ObsOperator<traits::LiteBackendTag>> obs_operators;
  auto obs_op = framework::ObsOperator<traits::LiteBackendTag>(*config_);
  obs_operators.push_back(std::move(obs_op));

  // Create observations
  std::vector<framework::Observation<traits::LiteBackendTag>> observations;
  auto obs1 = framework::Observation<traits::LiteBackendTag>(*config_);
  obs1.backend().setObservations({2.0, 3.5});
  obs1.backend().setCovariance({1.0, 1.0});
  observations.push_back(std::move(obs1));

  // Run the tangent linear check in single-epsilon mode using the unified
  // interface
  bool result =
      framework::checkObsOperatorTangentLinear<traits::LiteBackendTag>(
          obs_operators, *state_, observations, 1e-6, {0.1});
  EXPECT_TRUE(result);

  // Run the tangent linear check in multiple-epsilon mode using the unified
  // interface
  bool result2 =
      framework::checkObsOperatorTangentLinear<traits::LiteBackendTag>(
          obs_operators, *state_, observations, 1e-6, {1.0, 0.1, 0.01, 0.001});
  EXPECT_TRUE(result2);
}

}  // namespace metada::tests