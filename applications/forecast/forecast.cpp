/**
 * @file forecast.cpp
 * @brief Forecast application using the Metada framework
 * @author Metada Framework Team
 */

#include <chrono>
#include <string>

#include "ApplicationContext.hpp"
#include "DateTime.hpp"
#include "Geometry.hpp"
#include "Model.hpp"
#include "State.hpp"
#include "WRFBackendTraits.hpp"

// Default backends, can be changed with template parameters
using BackendTag = metada::traits::WRFBackendTag;
using namespace metada::framework;
using namespace metada::core;
using namespace metada::framework::runs;

int main(int argc, char** argv) {
  try {
    // Initialize application context
    ApplicationContext<BackendTag> context(argv[0], argv[1]);

    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info() << "Starting forecast application";
    logger.Info() << "Using configuration file: " << argv[1];

    const auto time_control = config.GetSubsection("time_control");
    // Read forecast parameters
    const auto start_datetime =
        DateTime(time_control.Get("start_datetime").asString());

    // Exit if neither forecast_length nor end_datetime is available
    if (!time_control.HasKey("forecast_length") &&
        !time_control.HasKey("end_datetime")) {
      logger.Error() << "Either forecast_length or end_datetime must be "
                        "specified in the config";
      return 1;
    }

    DateTime end_datetime;

    // If forecast_length exists, use it to calculate end_datetime
    if (time_control.HasKey("forecast_length")) {
      const int forecast_length_hours =
          time_control.Get("forecast_length").asInt();
      end_datetime = start_datetime + std::chrono::hours(forecast_length_hours);
      logger.Info() << "Using forecast_length of " << forecast_length_hours
                    << " hours to set end_datetime to " << end_datetime;
    } else {
      // Otherwise use the provided end_datetime
      end_datetime = DateTime(time_control.Get("end_datetime").asString());
      logger.Info() << "Using provided end_datetime: " << end_datetime;
    }

    const auto time_step =
        std::chrono::seconds(time_control.Get("time_step").asInt());

    const auto output_history =
        time_control.Get("output_history", "false").asBool();
    const std::string output_file = time_control.Get("history_file").asString();
    logger.Info() << "output_file: " << output_file;
    const bool write_history = true;
    logger.Info() << "write_history: " << write_history;
    const int history_frequency = time_control.Get("history_frequency").asInt();
    logger.Info() << "history_frequency: " << history_frequency;
    logger.Info() << "Forecast configuration:";
    logger.Info() << "  - Output file: " << output_file;

    // Initialize geometry
    logger.Info() << "Initializing geometry";
    Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    // Initialize initial state
    logger.Info() << "Initializing model state";
    // State<BackendTag> initialState(config);
    // State<BackendTag> currentState(config);

    // Initialize model
    logger.Info() << "Initializing forecast model";
    // Model<BackendTag> model(config);

    // Copy initial state to current state
    // currentState = initialState;  // Assuming State has assignment operator

    // Time integration loop using model.run
    auto current_datetime = start_datetime;
    while (current_datetime != end_datetime) {
      logger.Debug() << "Time Step on " << current_datetime;

      // Run model for this time segment
      // State<BackendTag> nextState(config);
      // model.run(currentState, nextState, current_time, segment_end_time);

      // Update current state and time
      // currentState = nextState;
      current_datetime += time_step;

      // Write intermediate output
      if (output_history) {
        logger.Info() << "Writing history at " << current_datetime;
        // Code to write history
      }
    }

    // Write final output
    logger.Info() << "Forecast completed successfully";
    logger.Info() << "Writing final state to " << output_file;
    // Code to write final state

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