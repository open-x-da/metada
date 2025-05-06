#include <iostream>
#include <string>
#include <vector>

#include "ApplicationContext.hpp"
#include "WRFBackendTraits.hpp"
#include "backends/common/io/BufrObsIO.hpp"
#include "framework/adapters/common/io/ObsIO.hpp"
#include "framework/adapters/common/io/ObsRecord.hpp"

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
    ObsIO<BackendTag> obsIO(std::move(config.GetSubsection("observation")));

    auto filename =
        config.GetSubsection("observation").Get("filename").asString();

    // Check if the file can be read
    if (!obsIO.canRead(filename)) {
      std::cerr << "Cannot read the specified file: " << filename << std::endl;
      return 1;
    }

    // Read observation records from the file
    std::vector<ObsRecord> records = obsIO.read(filename);

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