/**
 * @file forecast.cpp
 * @brief Forecast application using the Metada framework
 * @author Metada Framework Team
 */

#include <chrono>
#include <iostream>
#include <string>

#include "ApplicationContext.hpp"
#include "Config.hpp"
#include "DateTime.hpp"
#include "Geometry.hpp"
#include "Logger.hpp"
#include "Model.hpp"
#include "State.hpp"
#include "WRFBackendTraits.hpp"

// Default backends, can be changed with template parameters
using BackendTag = metada::traits::WRFBackendTag;
using namespace metada::framework;
using namespace metada::core;
using namespace metada::framework::runs;

int main([[maybe_unused]] int argc, char** argv) {
  try {
    // Initialize application context
    ApplicationContext<BackendTag> context(argv[0], argv[1]);

    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info() << "Starting forecast application";
    logger.Info() << "Using configuration file: " << argv[1];

    // Read forecast parameters
    const auto start_datetime =
        DateTime(config.Get("forecast.start_time").asString());
    const auto end_datetime =
        DateTime(config.Get("forecast.end_time").asString());
    const size_t max_steps = config.Get("forecast.max_steps").asInt();
    const auto dt =
        std::chrono::seconds(config.Get("forecast.time_step").asInt());
    const std::string output_file =
        config.Get("forecast.output_file").asString();
    const bool write_history = config.Get("forecast.write_history").asBool();
    const int history_frequency =
        config.Get("forecast.history_frequency").asInt();

    logger.Info() << "Forecast configuration:";
    logger.Info() << "  - Maximum steps: " << max_steps;
    logger.Info() << "  - Time step: " << dt << " seconds";
    logger.Info() << "  - Output file: " << output_file;

    // Initialize geometry
    logger.Info() << "Initializing geometry";
    Geometry<BackendTag> geometry(config);

    // Initialize initial state
    logger.Info() << "Initializing model state";
    // State<BackendTag> initialState(config);
    // State<BackendTag> currentState(config);

    // Initialize model
    logger.Info() << "Initializing forecast model";
    // Model<BackendTag> model(config);

    // Run forecast
    // logger.Info() << "Starting forecast integration from " << start_datetime
    //               << " to " << end_datetime;
    auto current_time = start_datetime;
    auto end_time = end_datetime;

    // Output initial state if needed
    if (write_history) {
      logger.Info() << "Writing initial state to " << output_file;
      // Code to write initial state
    }

    // Copy initial state to current state
    // currentState = initialState;  // Assuming State has assignment operator

    // Time integration loop using model.run
    if (write_history && history_frequency > 0) {
      // Use multiple shorter runs if intermediate outputs are needed
      for (int step = 0; step < max_steps; step += history_frequency) {
        // Calculate end time for this segment
        auto segment_end_time =
            std::min(current_time + history_frequency * dt, end_time);

        // logger.Debug() << "Integration from t=" << current_time
        //                << " to t=" << segment_end_time << " seconds";

        // Run model for this time segment
        // State<BackendTag> nextState(config);
        // model.run(currentState, nextState, current_time, segment_end_time);

        // Update current state and time
        // currentState = nextState;
        current_time = segment_end_time;

        // Write intermediate output
        if (write_history) {
          logger.Info() << "Writing history at t=" << current_time;
          // Code to write history
        }
      }
    } else {
      // Single run for the entire time period
      // State<BackendTag> finalState(config);
      logger.Debug() << "Integration from t=0 to t=" << end_time << " seconds";
      // model.run(initialState, finalState, 0.0, end_time);
      // currentState = finalState;
      current_time = end_time;
    }

    // Write final output
    logger.Info() << "Forecast completed successfully";
    logger.Info() << "Writing final state to " << output_file;
    // Code to write final state

    // Advanced usage: Iterate over the geometry to compute diagnostics
    logger.Info() << "Computing diagnostic values";
    double domain_avg = 0.0;
    int count = 0;

    // Example of using geometry iterator
    for (auto it = geometry.begin(); it != geometry.end(); ++it) {
      // const auto& point = *it;
      // Example diagnostic calculation
      // domain_avg += state.getValue("temperature", point);
      count++;
    }

    if (count > 0) {
      domain_avg /= count;
      logger.Info() << "Domain average: " << domain_avg;
    }

    logger.Info() << "Forecast application completed";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown error occurred" << std::endl;
    return 1;
  }
}