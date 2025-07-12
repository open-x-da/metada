#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "Increment.hpp"
#include "Logger.hpp"
#include "NonCopyable.hpp"
#include "State.hpp"

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
 * @brief Framework adapter for minimization algorithms
 *
 * @details This class bridges between the framework types (State<BackendTag>,
 * Increment<State<BackendTag>>) and the pure mathematical optimization
 * algorithms in base::optimization. It converts framework types to
 * std::vector<double>, calls the appropriate optimizer, and converts results
 * back to framework types.
 *
 * @tparam BackendTag The backend tag type
 */
template <typename BackendTag>
class Minimization : public NonCopyable {
 public:
  /**
   * @brief Cost function type definition for framework types
   *
   * Function that takes a state and returns the cost function value
   */
  using CostFunction = std::function<double(const State<BackendTag>&)>;

  /**
   * @brief Gradient function type definition for framework types
   *
   * Function that takes a state and computes the gradient into an increment
   */
  using GradientFunction = std::function<void(const State<BackendTag>&,
                                              Increment<State<BackendTag>>&)>;

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
  Minimization() = delete;

  /**
   * @brief Constructor with configuration
   *
   * @param config Configuration object containing minimization parameters
   */
  explicit Minimization(const Config<BackendTag>& config)
      : optimizer_(createOptimizer(config)) {
    logger_.Info() << "Minimization adapter constructed with "
                   << optimizer_->getName() << " algorithm";
    logger_.Info() << "Max iterations: " << optimizer_->getMaxIterations();
    logger_.Info() << "Tolerance: " << optimizer_->getTolerance();
    logger_.Info() << "Gradient tolerance: "
                   << optimizer_->getGradientTolerance();
  }

  /**
   * @brief Move constructor
   */
  Minimization(Minimization&& other) noexcept = default;

  /**
   * @brief Move assignment operator
   */
  Minimization& operator=(Minimization&& other) noexcept = default;

  /**
   * @brief Static factory method to create optimizer instances
   */
  static std::unique_ptr<base::optimization::OptimizerBase> createOptimizer(
      const Config<BackendTag>& config) {
    std::string algorithm_name =
        config.Get("minimization_algorithm").asString();

    int max_iterations = config.Get("max_iterations").asInt();
    double tolerance = config.Get("tolerance").asFloat();
    double gradient_tolerance = config.Get("gradient_tolerance").asFloat();
    bool line_search_enabled = config.Get("line_search_enabled").asBool();

    if (algorithm_name == "lbfgs" || algorithm_name == "L-BFGS") {
      size_t memory_size = config.Get("lbfgs_memory").asInt();
      return std::make_unique<base::optimization::LBFGSOptimizer>(
          max_iterations, tolerance, gradient_tolerance, memory_size,
          line_search_enabled);
    } else if (algorithm_name == "cg" ||
               algorithm_name == "conjugate_gradient") {
      int restart_interval = 0;
      if (config.HasKey("cg_restart_interval")) {
        restart_interval = config.Get("cg_restart_interval").asInt();
      }
      return std::make_unique<base::optimization::ConjugateGradientOptimizer>(
          max_iterations, tolerance, gradient_tolerance, line_search_enabled,
          restart_interval);
    } else if (algorithm_name == "steepest_descent" ||
               algorithm_name == "gradient_descent") {
      return std::make_unique<base::optimization::SteepestDescentOptimizer>(
          max_iterations, tolerance, gradient_tolerance, line_search_enabled);
    } else {
      // Default to L-BFGS
      Logger<BackendTag>::Instance().Warning()
          << "Unknown minimization algorithm: " << algorithm_name
          << ". Defaulting to L-BFGS.";
      size_t memory_size = 10;
      if (config.HasKey("lbfgs_memory")) {
        memory_size = config.Get("lbfgs_memory").asInt();
      }
      return std::make_unique<base::optimization::LBFGSOptimizer>(
          max_iterations, tolerance, gradient_tolerance, memory_size,
          line_search_enabled);
    }
  }

  /**
   * @brief Minimize the cost function starting from initial state
   *
   * @param initial_state Initial guess for the optimization
   * @param cost_func Cost function to minimize
   * @param gradient_func Gradient function
   * @param result Output state containing the optimized solution
   * @return Convergence information
   */
  ConvergenceInfo minimize(const State<BackendTag>& initial_state,
                           const CostFunction& cost_func,
                           const GradientFunction& gradient_func,
                           State<BackendTag>& result) {
    logger_.Info() << "Starting minimization with " << optimizer_->getName();

    // Convert initial state to vector
    std::vector<double> x = stateToVector(initial_state);

    // Wrap cost function
    auto vector_cost_func = [&](const std::vector<double>& vec) -> double {
      State<BackendTag> state = vectorToState(vec, initial_state);
      return cost_func(state);
    };

    // Wrap gradient function
    auto vector_gradient_func =
        [&](const std::vector<double>& vec) -> std::vector<double> {
      State<BackendTag> state = vectorToState(vec, initial_state);
      auto increment = Increment<State<BackendTag>>::createFromEntity(state);
      gradient_func(state, increment);
      return incrementToVector(increment);
    };

    // Call the optimizer
    std::vector<double> solution;
    auto optimizer_result = optimizer_->minimize(
        x, vector_cost_func, vector_gradient_func, solution);

    // Convert result back to state
    result = vectorToState(solution, initial_state);

    // Convert result information
    ConvergenceInfo info;
    info.converged = optimizer_result.converged;
    info.iterations = optimizer_result.iterations;
    info.initial_cost = optimizer_result.initial_cost;
    info.final_cost = optimizer_result.final_cost;
    info.gradient_norm = optimizer_result.gradient_norm;
    info.cost_reduction = optimizer_result.cost_reduction;
    info.convergence_reason = optimizer_result.convergence_reason;
    info.cost_evaluations = optimizer_result.cost_evaluations;
    info.gradient_evaluations = optimizer_result.gradient_evaluations;

    logger_.Info() << "Minimization completed: " << info.convergence_reason;
    logger_.Info() << "Iterations: " << info.iterations;
    logger_.Info() << "Final cost: " << info.final_cost;
    logger_.Info() << "Gradient norm: " << info.gradient_norm;
    logger_.Info() << "Cost evaluations: " << info.cost_evaluations;
    logger_.Info() << "Gradient evaluations: " << info.gradient_evaluations;

    return info;
  }

  /**
   * @brief Get the algorithm name
   */
  std::string getAlgorithmName() const { return optimizer_->getName(); }

  /**
   * @brief Get maximum iterations
   */
  int getMaxIterations() const { return optimizer_->getMaxIterations(); }

  /**
   * @brief Get tolerance
   */
  double getTolerance() const { return optimizer_->getTolerance(); }

  /**
   * @brief Get gradient tolerance
   */
  double getGradientTolerance() const {
    return optimizer_->getGradientTolerance();
  }

  /**
   * @brief Get the underlying optimizer instance
   */
  const base::optimization::OptimizerBase& getOptimizer() const {
    return *optimizer_;
  }

 private:
  std::unique_ptr<base::optimization::OptimizerBase> optimizer_;
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();

  /**
   * @brief Convert State to std::vector<double>
   *
   * This is a simplified conversion that assumes the state can be serialized
   * to a vector. In practice, each backend would implement this differently.
   */
  std::vector<double> stateToVector(const State<BackendTag>& state) const {
    // Get raw data pointer and convert to vector
    // This assumes the state data is contiguous double values
    const double* data_ptr = state.template getDataPtr<double>();
    size_t size = state.size();
    return std::vector<double>(data_ptr, data_ptr + size);
  }

  /**
   * @brief Convert std::vector<double> to State
   *
   * This reconstructs a state from a vector, using the reference state
   * for structure/metadata.
   */
  State<BackendTag> vectorToState(
      const std::vector<double>& vec,
      const State<BackendTag>& reference_state) const {
    // Clone the reference state to get the same structure
    auto new_state = reference_state.clone();

    // Update the data with the vector values
    double* data_ptr = new_state.template getDataPtr<double>();
    size_t size = std::min(vec.size(), new_state.size());

    for (size_t i = 0; i < size; ++i) {
      data_ptr[i] = vec[i];
    }

    return new_state;
  }

  /**
   * @brief Convert Increment to std::vector<double>
   */
  std::vector<double> incrementToVector(
      const Increment<State<BackendTag>>& increment) const {
    // Access the underlying state entity and convert to vector
    return stateToVector(increment.entity());
  }
};

}  // namespace metada::framework