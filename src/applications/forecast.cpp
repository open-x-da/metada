/**
 * @file forecast.cpp
 * @brief Forecast application using the Metada framework
 * @author Metada Framework Team
 */

#include <iostream>
#include <string>

#include "ApplicationContext.hpp"
#include "Config.hpp"
#include "Geometry.hpp"
#include "Logger.hpp"
#include "Model.hpp"
#include "State.hpp"
#include "WRFBackendTraits.hpp"

// Default backends, can be changed with template parameters
using BackendTag = metada::traits::WRFBackendTag;
namespace Framework = metada::framework;

int main(int argc, char** argv) {
  try {
    // Initialize application context
    Framework::ApplicationContext<BackendTag> context(argv[0], argv[1]);

    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.setLevel(Framework::Logger::Level::INFO);
    logger.info("Starting forecast application");
    logger.info("Using configuration file: {}", config_file);

    // Load configuration
    Framework::Config<BackendTag> config(config_file);

    // Read forecast parameters
    const int max_steps = config.getInt("forecast.max_steps", 100);
    const double dt = config.getDouble("forecast.time_step", 60.0);  // seconds
    const std::string output_file =
        config.getString("forecast.output_file", "forecast_output.nc");
    const bool write_intermediates =
        config.getBool("forecast.write_intermediates", false);
    const int output_frequency = config.getInt("forecast.output_frequency", 10);

    logger.info("Forecast configuration:");
    logger.info("  - Maximum steps: {}", max_steps);
    logger.info("  - Time step: {} seconds", dt);
    logger.info("  - Output file: {}", output_file);

    // Initialize geometry
    logger.info("Initializing geometry");
    Framework::Geometry<BackendTag> geometry(config);

    // Initialize initial state
    logger.info("Initializing model state");
    Framework::State<BackendTag> state(config);

    // Initialize model
    logger.info("Initializing forecast model");
    Framework::Model<BackendTag> model(config);

    // Validate configuration
    if (!model.isCompatible(state)) {
      logger.error("Model is not compatible with the provided state");
      return 1;
    }
    if (!model.isCompatible(geometry)) {
      logger.error("Model is not compatible with the provided geometry");
      return 1;
    }

    // Run forecast
    logger.info("Starting forecast integration for {} steps", max_steps);
    double current_time = 0.0;

    // Output initial state if needed
    if (write_intermediates) {
      logger.info("Writing initial state to {}", output_file);
      // Code to write initial state
    }

    // Time integration loop
    for (int step = 0; step < max_steps; ++step) {
      logger.debug("Step {}/{}: t = {} seconds", step + 1, max_steps,
                   current_time);

      // Run one model step
      model.runStep(state, dt);
      current_time += dt;

      // Output intermediate results if requested
      if (write_intermediates && (step + 1) % output_frequency == 0) {
        logger.info("Writing intermediate state at step {}", step + 1);
        // Code to write intermediate state
      }
    }

    // Write final output
    logger.info("Forecast completed successfully");
    logger.info("Writing final state to {}", output_file);
    // Code to write final state

    // Advanced usage: Iterate over the geometry to compute diagnostics
    logger.info("Computing diagnostic values");
    double domain_avg = 0.0;
    int count = 0;

    // Example of using geometry iterator
    for (auto it = geometry.begin(); it != geometry.end(); ++it) {
      const auto& point = *it;
      // Example diagnostic calculation
      // domain_avg += state.getValue("temperature", point);
      count++;
    }

    if (count > 0) {
      domain_avg /= count;
      logger.info("Domain average: {}", domain_avg);
    }

    logger.info("Forecast application completed");
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown error occurred" << std::endl;
    return 1;
  }
}