#include <iostream>
#include <string>
#include <vector>

#include "ApplicationContext.hpp"
#include "BufrObsIO.hpp"
#include "ObsIO.hpp"
#include "ObsRecord.hpp"
#include "WRFBackendTraits.hpp"

using BackendTag = metada::traits::WRFBackendTag;
using namespace metada;
using namespace metada::framework;
using namespace metada::backends::io;

int main(int argc, char* argv[]) {
  try {
    // Initialize application context
    ApplicationContext<BackendTag> context(argc, argv);

    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info() << "Starting obs_processing application";
    logger.Info() << "Using configuration file: " << argv[1];

    // Create an ObsIO object using BufrObsIO as the backend
    // The filename is already provided in the observation config section
    ObsIO<BackendTag> obsIO(config.GetSubsection("observation"));

    // Read observation records from the configured data source
    std::vector<ObsRecord> records = obsIO.read();

    // Iterate over the records and print the required information
    for (const auto& record : records) {
      std::cout << "Type: " << record.type << ", Value: " << record.value
                << ", Location: " << record.location
                << ", DateTime: " << record.datetime.iso8601() << std::endl;
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}