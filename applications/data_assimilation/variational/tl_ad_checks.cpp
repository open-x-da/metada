/**
 * @file tl_ad_checks.cpp
 * @brief Tangent Linear and Adjoint Checks Application for ObsOperator
 * @details This application performs tangent linear and adjoint checks for
 *          observation operators to verify the correctness of their TL/AD
 * implementations. It reads configuration from a file and performs systematic
 * checks on the mathematical consistency between tangent linear and adjoint
 * operators.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *            Expected format: tl_ad_checks <config_file>
 * @return 0 on success, 1 on failure
 */

#include <iomanip>
#include <iostream>
#include <vector>

#include "ApplicationContext.hpp"
#include "BackgroundErrorCovariance.hpp"
#include "Config.hpp"
#include "CostFunction.hpp"
#include "Geometry.hpp"
#include "Model.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "OperatorChecks.hpp"
#include "State.hpp"
#include "WRFBackendTraits.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::WRFBackendTag;

// Structure to hold check results
struct CheckResult {
  std::string test_name;
  bool passed;
  double relative_error;
  double tolerance;
  std::string details;
};

// Function to print check results in a formatted way
void printCheckResults(const std::vector<CheckResult>& results) {
  std::cout << "\n" << std::string(80, '=') << "\n";
  std::cout << "TL/AD AND GRADIENT CHECKS SUMMARY\n";
  std::cout << std::string(80, '=') << "\n";

  std::cout << std::left << std::setw(40) << "Test Name" << std::setw(10)
            << "Status" << std::setw(15) << "Rel Error" << std::setw(15)
            << "Tolerance" << "\n";
  std::cout << std::string(80, '-') << "\n";

  for (const auto& result : results) {
    std::cout << std::left << std::setw(40) << result.test_name;
    std::cout << std::setw(10) << (result.passed ? "PASS" : "FAIL");
    std::cout << std::setw(15) << std::scientific << std::setprecision(6)
              << result.relative_error;
    std::cout << std::setw(15) << std::scientific << std::setprecision(6)
              << result.tolerance << "\n";
  }

  std::cout << std::string(80, '=') << "\n";

  // Print detailed results for failed tests
  bool any_failed = false;
  for (const auto& result : results) {
    if (!result.passed) {
      if (!any_failed) {
        std::cout << "\nDETAILED FAILURE INFORMATION:\n";
        std::cout << std::string(50, '-') << "\n";
        any_failed = true;
      }
      std::cout << result.test_name << ": " << result.details << "\n";
    }
  }
}

int main(int argc, char** argv) {
  // Initialize application context
  auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
  auto& logger = context.getLogger();
  auto& config = context.getConfig();

  logger.Info() << "ObsOperator TL/AD Checks application starting...";

  try {
    // Validate command line arguments
    if (argc != 2) {
      logger.Error() << "Usage: tl_ad_checks <config_file>";
      return 1;
    }

    std::vector<CheckResult> results;

    // Get configuration parameters
    double tl_ad_tolerance = config.Get("tl_ad_tolerance").asFloat();
    std::vector<double> epsilons;
    auto eps_floats = config.Get("finite_difference_epsilons").asVectorFloat();
    for (auto f : eps_floats) {
      epsilons.push_back(static_cast<double>(f));
    }

    logger.Info() << "TL/AD tolerance: " << tl_ad_tolerance;
    if (!epsilons.empty()) {
      logger.Info() << "Finite difference epsilons: ";
      for (auto e : epsilons) logger.Info() << e;
    }

    // Initialize geometry
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));
    logger.Info() << "Geometry initialized";

    // Initialize state
    fwk::State<BackendTag> state(config.GetSubsection("state"), geometry);
    logger.Info() << "State initialized";

    // Initialize observations
    fwk::Observation<BackendTag> observation(
        config.GetSubsection("observations"));
    std::vector<fwk::Observation<BackendTag>> observations;
    observations.emplace_back(std::move(observation));
    logger.Info() << "Observations initialized";

    // Initialize observation operator
    fwk::ObsOperator<BackendTag> obs_operator(
        config.GetSubsection("obs_operator"));
    std::vector<fwk::ObsOperator<BackendTag>> obs_operators;
    obs_operators.emplace_back(std::move(obs_operator));
    logger.Info() << "Observation operator initialized";

    // Initialize background error covariance (needed for cost function)
    fwk::BackgroundErrorCovariance<BackendTag> bg_error_cov(
        config.GetSubsection("background_covariance"));
    logger.Info() << "Background error covariance initialized";

    // Initialize model (needed for cost function)
    fwk::Model<BackendTag> model(config.GetSubsection("model"));
    logger.Info() << "Model initialized";

    // Initialize cost function
    fwk::CostFunction<BackendTag> cost_function(
        config, state, observations, obs_operators, model, bg_error_cov);
    logger.Info() << "Cost function initialized";

    // Perform ObsOperator Tangent Linear checks
    logger.Info() << "\n=== Performing ObsOperator Tangent Linear Checks ===";

    try {
      // Check: ObsOperator Tangent Linear correctness
      double final_rel_error = 0.0;
      bool obs_op_tl_passed = false;
      // Use the new signature with configurable epsilons
      obs_op_tl_passed = fwk::checkObsOperatorTangentLinear(
          obs_operators, state, observations, tl_ad_tolerance, epsilons);
      // The function logs the final error, but if you want to extract it, you
      // could modify the function to return it. For now, we set it to 0.0 as
      // before, or you can parse from logs if needed.
      results.push_back(
          {"ObsOperator Tangent Linear", obs_op_tl_passed, final_rel_error,
           tl_ad_tolerance,
           obs_op_tl_passed
               ? "Observation operator tangent linear is correct"
               : "Observation operator tangent linear is incorrect"});

      logger.Info() << "ObsOperator Tangent Linear check: "
                    << (obs_op_tl_passed ? "PASSED" : "FAILED");

    } catch (const std::exception& e) {
      results.push_back({"ObsOperator Tangent Linear", false, 1.0,
                         tl_ad_tolerance,
                         std::string("Exception during check: ") + e.what()});
      logger.Error()
          << "ObsOperator Tangent Linear check failed with exception: "
          << e.what();
    }

    // Perform ObsOperator TL/AD checks
    logger.Info() << "\n=== Performing ObsOperator TL/AD Checks ===";

    try {
      // Check: ObsOperator TL/AD consistency
      bool obs_op_tl_ad_passed = fwk::checkObsOperatorTLAD(
          obs_operators, state, observations, tl_ad_tolerance);

      results.push_back(
          {"ObsOperator TL/AD Consistency", obs_op_tl_ad_passed,
           0.0,  // Will be updated if we have more detailed error info
           tl_ad_tolerance,
           obs_op_tl_ad_passed
               ? "Observation operator TL/AD are consistent"
               : "Observation operator TL/AD are inconsistent"});

      logger.Info() << "ObsOperator TL/AD check: "
                    << (obs_op_tl_ad_passed ? "PASSED" : "FAILED");

    } catch (const std::exception& e) {
      results.push_back({"ObsOperator TL/AD Consistency", false, 1.0,
                         tl_ad_tolerance,
                         std::string("Exception during check: ") + e.what()});
      logger.Error() << "ObsOperator TL/AD check failed with exception: "
                     << e.what();
    }

    // Perform Cost Function Gradient checks
    logger.Info() << "\n=== Performing Cost Function Gradient Checks ===";

    try {
      // Check: Cost Function Gradient correctness (single direction)
      bool cost_func_gradient_passed =
          fwk::checkCostFunctionGradient(cost_function, state, tl_ad_tolerance);

      results.push_back({"Cost Function Gradient (Single)",
                         cost_func_gradient_passed, 0.0, tl_ad_tolerance,
                         cost_func_gradient_passed
                             ? "Cost function gradient is correct"
                             : "Cost function gradient is incorrect"});

      logger.Info() << "Cost Function Gradient check (single): "
                    << (cost_func_gradient_passed ? "PASSED" : "FAILED");

    } catch (const std::exception& e) {
      results.push_back({"Cost Function Gradient (Single)", false, 1.0,
                         tl_ad_tolerance,
                         std::string("Exception during check: ") + e.what()});
      logger.Error() << "Cost Function Gradient check failed with exception: "
                     << e.what();
    }

    // Check: Cost Function Gradient with multiple random directions
    try {
      size_t num_directions = config.Get("gradient_check_directions").asInt();
      bool cost_func_gradient_multi_passed =
          fwk::checkCostFunctionGradientMultipleDirections(
              cost_function, state, num_directions, tl_ad_tolerance);

      results.push_back(
          {"Cost Function Gradient (Multi)", cost_func_gradient_multi_passed,
           0.0, tl_ad_tolerance,
           cost_func_gradient_multi_passed
               ? "Cost function gradient is correct across multiple directions"
               : "Cost function gradient is incorrect across multiple "
                 "directions"});

      logger.Info() << "Cost Function Gradient check (multiple directions): "
                    << (cost_func_gradient_multi_passed ? "PASSED" : "FAILED");

    } catch (const std::exception& e) {
      results.push_back({"Cost Function Gradient (Multi)", false, 1.0,
                         tl_ad_tolerance,
                         std::string("Exception during check: ") + e.what()});
      logger.Error()
          << "Cost Function Gradient check (multiple) failed with exception: "
          << e.what();
    }

    // Check: Cost Function Gradient with unit vector directions
    try {
      bool cost_func_gradient_unit_passed =
          fwk::checkCostFunctionGradientUnitDirections(cost_function, state,
                                                       tl_ad_tolerance);

      results.push_back(
          {"Cost Function Gradient (Unit)", cost_func_gradient_unit_passed, 0.0,
           tl_ad_tolerance,
           cost_func_gradient_unit_passed
               ? "Cost function gradient is correct along unit vectors"
               : "Cost function gradient is incorrect along unit vectors"});

      logger.Info() << "Cost Function Gradient check (unit vectors): "
                    << (cost_func_gradient_unit_passed ? "PASSED" : "FAILED");

    } catch (const std::exception& e) {
      results.push_back({"Cost Function Gradient (Unit)", false, 1.0,
                         tl_ad_tolerance,
                         std::string("Exception during check: ") + e.what()});
      logger.Error() << "Cost Function Gradient check (unit vectors) failed "
                        "with exception: "
                     << e.what();
    }

    // Print summary
    printCheckResults(results);

    // Determine overall success
    bool all_passed = true;
    for (const auto& result : results) {
      if (!result.passed) {
        all_passed = false;
        break;
      }
    }

    if (all_passed) {
      logger.Info() << "\nAll TL/AD and gradient checks PASSED!";
      return 0;
    } else {
      logger.Error() << "\nSome TL/AD and gradient checks FAILED!";
      return 1;
    }
  } catch (const std::exception& e) {
    logger.Error() << "ObsOperator TL/AD checks application failed: "
                   << e.what();
    return 1;
  }
}