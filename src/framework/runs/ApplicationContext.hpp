#pragma once
#include <iostream>
#include <stdexcept>

#include "BackendTraits.hpp"
#include "common/utils/config/Config.hpp"
#include "common/utils/logger/Logger.hpp"

namespace metada::framework {

/**
 * @brief RAII wrapper managing application-wide runtime context and services
 *
 * @details
 * The ApplicationContext class provides centralized initialization, access, and
 * cleanup of core application services. It follows RAII principles to ensure
 * proper initialization order and cleanup of all managed services.
 *
 * The context maintains singleton instances of critical application services
 * and handles their lifecycle in a thread-safe manner. Services are initialized
 * in a specific order during construction and cleaned up in reverse order
 * during destruction.
 *
 * @par Managed Services
 * - Logging System: Application-wide logging and diagnostics
 * - Configuration: Reading and accessing application settings from YAML/JSON
 * files
 * - Performance Monitoring: Tracking application metrics (planned)
 * - Resource Management: Managing system resources (planned)
 * - Runtime Environment: MPI and threading configuration (planned)
 *
 * @par Key Features
 * - RAII design ensures proper resource initialization and cleanup
 * - Thread-safe access to all managed services
 * - Move semantics support for RAII usage
 * - Non-copyable to maintain service uniqueness
 * - Centralized error handling and logging
 * - Support for multiple backend implementations through traits
 *
 * @par Example Usage
 * @code{.cpp}
 * int main(int argc, char* argv[]) {
 *   // Create context with app name and config file
 *   auto context = ApplicationContext<L63BackendTag>(argc, argv);
 *
 *   // Access managed services
 *   auto& logger = context.getLogger();
 *   auto& config = context.getConfig();
 *
 *   logger.Info() << "Application started";
 *
 *   // Use services...
 *
 *   return 0;
 *   // Context and services cleaned up automatically
 * }
 * @endcode
 *
 * @tparam BackendTag Tag type defining backend service implementations through
 * traits
 *
 * @note The context should be instantiated exactly once at application startup
 * @note Services are initialized in order: Config -> Logger -> Future Services
 * @note Cleanup occurs automatically in reverse order on destruction
 * @see metada::traits::BackendTraits for backend implementation details
 */
template <typename BackendTag>
class ApplicationContext {
 public:
  /**
   * @brief Constructs and initializes the application context
   *
   * @param argc Number of command line arguments
   * @param argv Array of command line arguments
   * @throws std::runtime_error If initialization of any service fails or config
   * file is missing
   *
   * @details Initializes services in the following order:
   * 1. Check if config file path is provided in command line
   * 2. Configuration loading from specified file
   * 3. Logger initialization with configuration
   * 4. Logs successful initialization with application name
   */
  ApplicationContext(const int argc, char** argv)
      : config_(validateAndGetConfigPath(argc, argv)) {
    Logger<BackendTag>::Init(config_.GetSubsection("logger"));
    Logger<BackendTag>::Instance().Info()
        << "Application context initialized: " << argv[0];
  }

  /**
   * @brief Destructor - cleans up application services
   *
   * @details Services are cleaned up in reverse initialization order:
   * 1. Future services (when implemented)
   * 2. Logger
   * 3. Configuration
   */
  ~ApplicationContext() {
    Logger<BackendTag>::Instance().Info() << "Application context destroyed";
  }

  /**
   * @brief Disabled copy constructor
   *
   * @details ApplicationContext instances are not intended to be copied
   * to ensure single instance of services.
   */
  ApplicationContext(const ApplicationContext&) = delete;

  /**
   * @brief Disabled copy assignment operator
   *
   * @details ApplicationContext instances are not intended to be copied
   * to ensure single instance of services.
   */
  ApplicationContext& operator=(const ApplicationContext&) = delete;

  /**
   * @brief Move constructor
   *
   * @details Allows for moving ApplicationContext instances when needed
   * to support RAII pattern.
   */
  ApplicationContext(ApplicationContext&&) = default;

  /**
   * @brief Move assignment operator
   *
   * @details Allows for moving ApplicationContext instances when needed
   * to support RAII pattern.
   *
   * @return Reference to this ApplicationContext instance
   */
  ApplicationContext& operator=(ApplicationContext&&) = default;

  /**
   * @brief Get reference to the logging service
   *
   * @details Provides access to the application's logging system for
   * recording information, warnings, errors, and debug messages.
   *
   * @return Reference to the logger instance
   */
  Logger<BackendTag>& getLogger() { return Logger<BackendTag>::Instance(); }

  /**
   * @brief Get reference to the configuration service
   *
   * @details Provides access to the application's configuration system
   * for retrieving settings and parameters.
   *
   * @return Reference to the config instance
   */
  Config<BackendTag>& getConfig() { return config_; }

  // Timer access will be added later
  // Timer& getTimer() { return timer_; }

 private:
  /**
   * @brief Validates command line arguments and returns config file path
   *
   * @param argc Number of command line arguments
   * @param argv Array of command line arguments
   * @return const char* Path to configuration file
   * @throws std::runtime_error If config file path is missing
   */
  static const char* validateAndGetConfigPath(const int argc, char** argv) {
    if (argc < 2) {
      std::cerr << "Error: No configuration file specified." << std::endl;
      std::cerr << "Usage: " << (argc > 0 ? argv[0] : "application")
                << " <path_to_config_file>" << std::endl;
      throw std::runtime_error(
          "Missing configuration file path in command line arguments");
    }
    return argv[1];
  }

  Config<BackendTag> config_;  ///< Configuration service instance
  // Timer timer_;  // To be implemented
};

}  // namespace metada::framework