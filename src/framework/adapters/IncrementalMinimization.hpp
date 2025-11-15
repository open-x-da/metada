#pragma once

#include <algorithm>
#include <cmath>
#include <format>
#include <iomanip>
#include <numeric>
#include <sstream>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "ControlVariable.hpp"
#include "ControlVariableBackend.hpp"
#include "Logger.hpp"

// Framework-independent optimizers
#include "ConjugateGradientOptimizer.hpp"
#include "LBFGSOptimizer.hpp"
#include "OptimizerBase.hpp"
#include "SteepestDescentOptimizer.hpp"

namespace metada::framework {

// Forward declarations
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

/**
 * @brief Framework adapter for incremental minimization algorithms
 *
 * @details This class bridges between the framework types
 * (ControlVariable<BackendTag>) and the pure mathematical optimization
 * algorithms in base::optimization. It converts framework types to
 * std::vector<double>, calls the appropriate optimizer, and converts results
 * back to framework types.
 *
 * The minimization is performed over control variables v, which are the actual
 * optimization variables in incremental variational data assimilation. The
 * control variables are transformed to state-space increments δx inside the
 * cost function evaluation.
 *
 * @tparam BackendTag The backend tag type
 */
template <typename BackendTag>
  requires StateBackendType<BackendTag>
class IncrementalMinimization : public NonCopyable {
 public:
  /**
   * @brief Cost function type definition for control-space formulation
   *
   * Function that takes a control variable and returns the cost function value
   */
  using CostFunction =
      std::function<double(const ControlVariable<BackendTag>&)>;

  /**
   * @brief Gradient function type definition for control-space formulation
   *
   * Function that takes a control variable and computes the gradient into
   * another control variable
   */
  using GradientFunction = std::function<void(
      const ControlVariable<BackendTag>&, ControlVariable<BackendTag>&)>;

  /**
   * @brief Convergence information structure
   */
  struct ConvergenceInfo {
    bool converged = false;          ///< Whether convergence was achieved
    int iterations = 0;              ///< Number of iterations performed
    double initial_cost = 0.0;       ///< Initial cost function value
    double final_cost = 0.0;         ///< Final cost function value
    double gradient_norm = 0.0;      ///< Final gradient norm
    double cost_reduction = 0.0;     ///< Total cost reduction achieved
    std::string convergence_reason;  ///< Reason for convergence/termination

    // Performance metrics
    int cost_evaluations = 0;  ///< Total number of cost function evaluations
    int gradient_evaluations = 0;  ///< Total number of gradient evaluations
  };

  /** @brief Default constructor is deleted */
  IncrementalMinimization() = delete;

  /**
   * @brief Constructor with configuration
   *
   * @param config Configuration object containing minimization parameters
   */
  IncrementalMinimization(
      const Config<BackendTag>& config,
      const ControlVariableBackend<BackendTag>& control_backend)
      : control_backend_(control_backend), optimizer_(createOptimizer(config)) {
    logger_.Info() << "IncrementalMinimization adapter constructed with "
                   << optimizer_->getName() << " algorithm";
    logger_.Info() << "Max iterations: " << optimizer_->getMaxIterations();
    logger_.Info() << "Tolerance: " << optimizer_->getTolerance();
    logger_.Info() << "Gradient tolerance: "
                   << optimizer_->getGradientTolerance();

    // Set up generic iteration logger for any optimizer that supports it
    optimizer_->setIterationLogger(
        [this](int iter, double previous_cost, double current_cost,
               double gradient_norm, double step_size) {
          logger_.Info() << "Iteration " << iter
                         << ": previous_cost=" << previous_cost
                         << ", current_cost=" << current_cost
                         << ", cost_change=" << (current_cost - previous_cost)
                         << ", gradient_norm=" << gradient_norm
                         << ", step=" << step_size;
        });
  }

  /**
   * @brief Move constructor
   */
  IncrementalMinimization(IncrementalMinimization&& other) noexcept = default;

  /**
   * @brief Move assignment operator
   */
  IncrementalMinimization& operator=(IncrementalMinimization&&) noexcept =
      delete;

  /**
   * @brief Destructor
   */
  ~IncrementalMinimization() = default;

  /**
   * @brief Minimize the cost function over control variables v
   *
   * @param initial_control Initial guess for the control variable v₀
   * @param cost_func Cost function J(v)
   * @param gradient_func Gradient function ∇_v J
   * @param final_control Output: final optimized control variable v*
   * @return Convergence information
   */
  ConvergenceInfo minimize(const ControlVariable<BackendTag>& initial_control,
                           const CostFunction& cost_func,
                           const GradientFunction& gradient_func,
                           ControlVariable<BackendTag>& final_control) const {
    logger_.Info() << "Starting incremental minimization with "
                   << optimizer_->getName();

    // Convert initial control variable to std::vector<double>
    auto initial_vector =
        initial_control.template getData<std::vector<double>>();

    std::vector<double> final_data(initial_vector.size());

    // Cost function wrapper: converts std::vector<double> to ControlVariable
    auto cost_wrapper = [this, &cost_func, &initial_control](
                            const std::vector<double>& x) -> double {
      auto temp_control =
          control_backend_.createControlVariable(initial_control.geometry());
      temp_control.setFromVector(x);
      return cost_func(temp_control);
    };

    // Gradient wrapper: converts std::vector<double> to ControlVariable for
    // input and output
    auto gradient_wrapper =
        [this, &gradient_func, &initial_control](
            const std::vector<double>& x) -> std::vector<double> {
      auto temp_control =
          control_backend_.createControlVariable(initial_control.geometry());
      temp_control.setFromVector(x);

      auto temp_gradient =
          control_backend_.createControlVariable(initial_control.geometry());
      temp_gradient.zero();

      gradient_func(temp_control, temp_gradient);

      return temp_gradient.template getData<std::vector<double>>();
    };

    const double initial_cost = cost_wrapper(initial_vector);
    const auto initial_gradient = gradient_wrapper(initial_vector);
    const double initial_gradient_norm = std::sqrt(
        std::inner_product(initial_gradient.begin(), initial_gradient.end(),
                           initial_gradient.begin(), 0.0));

    const std::string method_prefix =
        "minimize_" + toLowerCopy(optimizer_->getName());
    constexpr int prefix_width = 24;
    constexpr int loop_width = 6;
    constexpr int iter_width = 6;
    constexpr int value_width = 27;

    logSeparator(prefix_width, loop_width, iter_width, value_width);
    logHeader(prefix_width, loop_width, iter_width, value_width);
    logSeparator(prefix_width, loop_width, iter_width, value_width);
    logIteration(method_prefix, prefix_width, loop_width, iter_width,
                 value_width, 1, 0, initial_cost, initial_gradient_norm, 0.0);

    optimizer_->setIterationLogger(
        [method_prefix](int iter, double /*previous_cost*/, double current_cost,
                        double gradient_norm, double step_size) {
          logIteration(method_prefix, prefix_width, loop_width, iter_width,
                       value_width, 1, iter, current_cost, gradient_norm,
                       step_size);
        });

    auto convergence_info = optimizer_->minimize(initial_vector, cost_wrapper,
                                                 gradient_wrapper, final_data);
    logSeparator(prefix_width, loop_width, iter_width, value_width);

    // Convert final vector back to control variable
    final_control =
        control_backend_.createControlVariable(initial_control.geometry());
    final_control.setFromVector(final_data);

    ConvergenceInfo result;
    result.converged = convergence_info.converged;
    result.iterations = convergence_info.iterations;
    result.initial_cost = convergence_info.initial_cost;
    result.final_cost = convergence_info.final_cost;
    result.gradient_norm = convergence_info.gradient_norm;
    result.cost_reduction =
        convergence_info.initial_cost - convergence_info.final_cost;
    result.convergence_reason = convergence_info.convergence_reason;
    result.cost_evaluations = convergence_info.cost_evaluations;
    result.gradient_evaluations = convergence_info.gradient_evaluations;

    return result;
  }

  /**
   * @brief Get the optimizer name
   */
  std::string getAlgorithmName() const { return optimizer_->getName(); }

  /**
   * @brief Get the maximum number of iterations
   */
  int getMaxIterations() const { return optimizer_->getMaxIterations(); }

  /**
   * @brief Get the cost function tolerance
   */
  double getTolerance() const { return optimizer_->getTolerance(); }

  /**
   * @brief Get the gradient tolerance
   */
  double getGradientTolerance() const {
    return optimizer_->getGradientTolerance();
  }

 private:
  /**
   * @brief Create the appropriate optimizer based on configuration
   */
  std::unique_ptr<base::optimization::OptimizerBase> createOptimizer(
      const Config<BackendTag>& config) {
    std::string algorithm = config.Get("minimization_algorithm").asString();
    int max_iterations = config.Get("max_iterations").asInt();
    double tolerance = config.Get("cost_tolerance").asFloat();
    double gradient_tolerance = config.Get("gradient_tolerance").asFloat();

    if (algorithm == "L-BFGS") {
      return std::make_unique<base::optimization::LBFGSOptimizer>(
          max_iterations, tolerance, gradient_tolerance);
    } else if (algorithm == "CG") {
      return std::make_unique<base::optimization::ConjugateGradientOptimizer>(
          max_iterations, tolerance, gradient_tolerance);
    } else if (algorithm == "SteepestDescent") {
      return std::make_unique<base::optimization::SteepestDescentOptimizer>(
          max_iterations, tolerance, gradient_tolerance);
    } else {
      logger_.Warning() << "Unknown minimization algorithm: " << algorithm
                        << ". Defaulting to L-BFGS.";
      return std::make_unique<base::optimization::LBFGSOptimizer>(
          max_iterations, tolerance, gradient_tolerance);
    }
  }

  const ControlVariableBackend<BackendTag>& control_backend_;
  std::unique_ptr<base::optimization::OptimizerBase> optimizer_;
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();

  static std::string toLowerCopy(std::string value) {
    std::transform(
        value.begin(), value.end(), value.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
  }

  static void logSeparator(int prefix_width, int loop_width, int iter_width,
                           int value_width) {
    const int total_width =
        prefix_width + loop_width + iter_width + value_width * 3 + 1;
    Logger<BackendTag>::Instance().Info() << std::string(total_width, '-');
  }

  static void logHeader(int prefix_width, int loop_width, int iter_width,
                        int value_width) {
    auto line =
        std::format(" {0:<{1}}{2:>{3}}{4:>{5}}{6:>{7}}{8:>{9}}{10:>{11}}", "",
                    prefix_width - 1, "Loop", loop_width, "Iter", iter_width,
                    "Cost Function", value_width, "Gradient", value_width,
                    "Step", value_width);
    Logger<BackendTag>::Instance().Info() << line;
  }

  static void logIteration(const std::string& prefix, int prefix_width,
                           int loop_width, int iter_width, int value_width,
                           int loop, int iter, double cost, double gradient,
                           double step) {
    auto line = std::format(
        " {0:<{1}}{2:>{3}}{4:>{5}}{6:>{7}.15e}{8:>{9}.15e}{10:>{11}.15e}",
        prefix, prefix_width - 1, loop, loop_width, iter, iter_width, cost,
        value_width, gradient, value_width, step, value_width);
    Logger<BackendTag>::Instance().Info() << line;
  }
};

}  // namespace metada::framework
