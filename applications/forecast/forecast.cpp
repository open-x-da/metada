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
#ifdef USE_MPI
    // Initialize parallel environment
    auto& parallel = backends::macom::MACOMParallel::getInstance();
    parallel.setFortranMode(true);  // Explicitly use Fortran MPI initialization
    if (!parallel.initialize(argc, argv)) {
      std::cerr << "Failed to initialize parallel environment" << std::endl;
      return 1;
    }
#endif

    // Initialize application context
    ApplicationContext<BackendTag> context(argc, argv);
    auto& logger = context.getLogger();
    auto& config = context.getConfig();

#ifdef USE_MPI
    // Log process information
    if (parallel.isParallel()) {
      logger.Info() << "Running as MPI process " << parallel.getRank();
    }
#else
    logger.Info() << "Running in serial mode";
#endif

#ifdef USE_MPI
    parallel.barrier();

    if (parallel.getRank() == 0) {
      logger.Info() << "Starting forecast application";
      logger.Info() << "Using configuration file: " << argv[1];
    } else {
      parallel.finalize();
      return 0;
    }
#endif

    // // Initialize geometry
    // logger.Info() << "Initializing geometry";
    // Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    // // Initialize initial state
    // logger.Info() << "Initializing model state";
    // State<BackendTag> initialState(config.GetSubsection("state"), geometry);
    // auto currentState = initialState.clone();

    // // Initialize model
    // logger.Info() << "Initializing forecast model";
    // Model<BackendTag> model(config.GetSubsection("model"));

    // // Create final state
    // auto finalState = initialState.clone();
    // // State<BackendTag> finalState(config.GetSubsection("state"), geometry);

    // // Run the model
    // logger.Info() << "Running forecast model...";
    // model.run(currentState, finalState);

    // // Update current state
    // currentState = std::move(finalState);

    logger.Info() << "Forecast application completed";

#ifdef USE_MPI
    // Finalize MPI if MPI support is enabled
    parallel.finalize();
#endif
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    // Ensure parallel environment is finalized even in case of error
#ifdef USE_MPI
    backends::macom::MACOMParallel::getInstance().finalize();
#endif
    return 1;
  } catch (...) {
    std::cerr << "Unknown error occurred" << std::endl;
    // Ensure parallel environment is finalized even in case of error
#ifdef USE_MPI
    backends::macom::MACOMParallel::getInstance().finalize();
#endif
    return 1;
  }
}