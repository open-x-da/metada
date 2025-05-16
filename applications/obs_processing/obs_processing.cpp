#include <iomanip>  // For formatting output
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
    logger.Info() << "ObsIO created";

    // Read observation records from the configured data source
    std::vector<ObsRecord> records = obsIO.read();
    logger.Info() << "Records read";

    // Iterate over the records and print the information
    // int count = 0;
    // for (const auto& record : records) {
    //   if (count % 100 == 0) {
    //     std::cout << record << std::endl;
    //   }
    //   count++;
    // }

    logger.Info() << "\nTotal records: " << records.size();
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}