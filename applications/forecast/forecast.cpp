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
using namespace metada;
using namespace metada::framework;

int main(int argc, char** argv) {
  try {
    // Initialize application context
    ApplicationContext<BackendTag> context(argc, argv);

    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info() << "Starting forecast application";
    logger.Info() << "Using configuration file: " << argv[1];

    // Initialize geometry
    logger.Info() << "Initializing geometry";
    Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    // Initialize initial state
    logger.Info() << "Initializing model state";
    State<BackendTag> initialState(config.GetSubsection("state"));
    auto currentState = initialState.clone();

    // Initialize model
    logger.Info() << "Initializing forecast model";
    Model<BackendTag> model(config.GetSubsection("model"));

    // Create final state
    State<BackendTag> finalState(config.GetSubsection("state"));

    // Run the model
    logger.Info() << "Running forecast model...";
    model.run(currentState, finalState);

    // Update current state
    currentState = std::move(finalState);

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