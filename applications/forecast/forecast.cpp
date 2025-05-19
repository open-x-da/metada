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
#include "MACOMBackendTraits.hpp"
#include "Model.hpp"
#include "State.hpp"
// #include "WRFBackendTraits.hpp"

// Default backends, can be changed with template parameters
// using BackendTag = metada::traits::WRFBackendTag;
using BackendTag = metada::traits::MACOMBackendTag;
using namespace metada;
using namespace metada::framework;

int main(int argc, char** argv) {
  try {
    // Initialize application context
    ApplicationContext<BackendTag> context(argc, argv);

    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    // 1. 先获取 "geometry" 子配置块
    auto geometry_config_subsection = config.GetSubsection("geometry");

    // 2. 然后从 "geometry" 子配置块中获取 "input_file"
    //    进行错误检查是一个好习惯
    if (geometry_config_subsection.HasKey("input_file")) {
      logger.Info() << "Geometry input file: "
                    << geometry_config_subsection.Get("input_file").asString();
    } else {
      logger.Error() << "Key 'input_file' not found in 'geometry' subsection.";
    }

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