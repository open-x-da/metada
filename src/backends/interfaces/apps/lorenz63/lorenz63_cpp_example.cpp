#include <chrono>
#include <stdexcept>
#include <string>

#include "AppTraits.hpp"
#include "ApplicationContext.hpp"
#include "GoogleLogger.hpp"
#include "JsonConfig.hpp"
#include "L63BackendTraits.hpp"
#include "State.hpp"
#include "geometry/lorenz63/Geometry.hpp"
#include "model/lorenz63/Model.hpp"
#include "state/lorenz63/State.hpp"

using state_backend = metada::backends::lorenz63::State;
using geometry_backend = metada::backends::lorenz63::Geometry;
using model_backend = metada::backends::lorenz63::Model;

using metada::backends::lorenz63::Integrator;
using namespace metada::framework;
using namespace metada::backends::logger;
using namespace metada::backends::config;

using Traits = AppTraits<GoogleLogger, JsonConfig, state_backend>;

int main(int argc, char* argv[]) {
  // Check if the config file is provided
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <config_file>" << std::endl;
    return 1;
  }

  // Create the application context
  auto app = ApplicationContext<traits::L63BackendTag>("Lorenz63", argv[1]);

  // Get the logger and config from the application context
  auto& logger = app.getLogger();
  auto& config = app.getConfig();

  try {
    // Create a Lorenz63 state from the config
    State<Traits::StateType> state(config);
    State<Traits::StateType> initial_state(config);

    // Create the Lorenz63 model with standard parameters
    // sigma = 10.0, rho = 28.0, beta = 8/3, dt = 0.01
    // Model model(config);

    // Define the phase space bounds (typical values for the Lorenz attractor)
    // Geometry geometry(config.getPhaseSpaceBounds());
    // logger.Info() << "Phase space geometry: " << geometry;

    // Get the model parameters to verify
    // float sigma, rho, beta, dt;
    // model.getParameters(sigma, rho, beta, dt);
    // logger.Info() << "Model parameters:";
    // logger.Info() << "  sigma = " << sigma;
    // logger.Info() << "  rho = " << rho;
    // logger.Info() << "  beta = " << beta;
    // logger.Info() << "  dt = " << dt;

    // Create an integrator with the model
    // Integrator integrator(model);

    // Check if initial state is within bounds
    // if (geometry.containsPoint(state)) {
    //  logger.Info() << "Initial state is within phase space bounds.";
    //} else {
    //  logger.Warning() << "Initial state is outside phase space bounds!";
    //}

    // Measure performance
    auto start_time = std::chrono::high_resolution_clock::now();

    // Run the integration for 1000 steps
    // const int num_steps = 1000;
    // for (int i = 0; i < num_steps; ++i) {
    //  integrator.step(state);

    // Print every 100 steps
    // if (i % 100 == 0) {
    //  logger.Info() << "Step " << i << ": " << state;

    // Check if current state is within bounds
    // if (!geometry.containsPoint(state)) {
    //    logger.Warning() << "State has left the phase space bounds!";
    //  }
    //}
    //}

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time);

    // Calculate distance from initial to final state
    // float distance = state.distance(initial_state);
    // logger.Info() << "Distance from initial to final state: " << distance;

    // Print performance metrics
    // logger.Info() << "Performance:";
    // logger.Info() << "  Executed " << num_steps << " steps in "
    //              << duration.count() << " ms";
    // logger.Info() << "  Average time per step: "
    //              << static_cast<float>(duration.count()) / num_steps << "
    //              ms";

    // Demonstrate vector access to components
    // auto components = state.getComponents();
    // logger.Info() << "Final state as vector: [" << components[0] << ", "
    //              << components[1] << ", " << components[2] << "]";

    // Final state using stream operator
    // logger.Info() << "Final state using stream operator: " << state;

    // Demonstrate multiple states in one output
    // logger.Info() << "Initial and final states: " << initial_state << " -> "
    //              << state;
  } catch (const std::exception& e) {
    logger.Error() << "Error: " << e.what();
    return 1;
  }

  return 0;
}