#pragma once
#include <stdexcept>
#include <string>

#include "BackendTraits.hpp"
#include "utils/config/Config.hpp"
#include "utils/logger/Logger.hpp"

namespace metada::framework::runs {

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
 * - Configuration: Reading and accessing application settings
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
 *
 * @par Example Usage
 * @code{.cpp}
 * int main(int argc, char* argv[]) {
 *   // Create context with app name and config file
 *   auto context = ApplicationContext<Traits>("my_app", argv[1]);
 *
 *   // Access managed services
 *   auto& logger = context.getLogger();
 *   auto& config = context.getConfig();
 *
 *   logger.Info("Application started");
 *
 *   // Use services...
 *
 *   return 0;
 *   // Context and services cleaned up automatically
 * }
 * @endcode
 *
 * @tparam Traits Configuration traits class defining service types:
 *         - LoggerType: The logging service implementation
 *         - ConfigType: The configuration service implementation
 *
 * @note The context should be instantiated exactly once at application startup
 * @note Services are initialized in order: Logger -> Config -> Future Services
 * @note Cleanup occurs automatically in reverse order on destruction
 */
template <typename BackendTag>
class ApplicationContext {
 public:
  using MyTraits = metada::traits::BackendTraits<BackendTag>;
  using ConfigBackend = typename MyTraits::ConfigBackend;
  using LoggerBackend = typename MyTraits::LoggerBackend;

  /**
   * @brief Constructs and initializes the application context
   *
   * @param app_name Name of the application used for logging and identification
   * @param config_file Optional path to configuration file (empty for defaults)
   * @throws std::runtime_error If initialization of any service fails
   *
   * @details Initializes services in the following order:
   * 1. Logger initialization with application name
   * 2. Configuration loading if config file specified
   * 3. Logs successful initialization
   */
  ApplicationContext(const std::string& app_name,
                     const std::string& config_file = "") {
    initLogger(app_name);
    if (!config_file.empty()) {
      loadConfig(config_file);
    }

    logger_.Info() << "Application context initialized: " << app_name;
  }

  /**
   * @brief Destructor - cleans up application services
   *
   * @details Services are cleaned up in reverse initialization order:
   * 1. Future services (when implemented)
   * 2. Configuration
   * 3. Logger
   */
  ~ApplicationContext() {
    logger_.Info() << "Shutting down application context";
    shutdownLogger();
  }

  // Prevent copying to ensure single instance of services
  ApplicationContext(const ApplicationContext&) = delete;
  ApplicationContext& operator=(const ApplicationContext&) = delete;

  // Allow moving to support RAII
  ApplicationContext(ApplicationContext&&) = default;
  ApplicationContext& operator=(ApplicationContext&&) = default;

  /**
   * @brief Get reference to the logging service
   * @return Reference to the logger instance
   */
  Logger<LoggerBackend>& getLogger() { return logger_; }

  /**
   * @brief Get reference to the configuration service
   * @return Reference to the config instance
   */
  Config<BackendTag>& getConfig() { return config_; }

  // Timer access will be added later
  // Timer& getTimer() { return timer_; }

 private:
  Logger<LoggerBackend> logger_;
  Config<BackendTag> config_;
  // Timer timer_;  // To be implemented

  /**
   * @brief Initialize the logging service
   * @param app_name Application name for logger identification
   */
  void initLogger(const std::string& app_name) {
    logger_.backend().Init(app_name);
  }

  /**
   * @brief Shutdown the logging service cleanly
   */
  void shutdownLogger() { logger_.backend().Shutdown(); }

  /**
   * @brief Load configuration from specified file
   * @param config_file Path to configuration file
   * @throws std::runtime_error If configuration loading fails
   */
  void loadConfig(const std::string& config_file) {
    try {
      if (!config_.LoadFromFile(config_file)) {
        logger_.Error() << "Failed to load configuration from: " << config_file;
        throw std::runtime_error("Failed to load configuration from: " +
                                 config_file);
      }
      logger_.Info() << "Loaded configuration from: " << config_file;
    } catch (const std::exception& e) {
      logger_.Error() << "Error loading configuration: " << e.what();
      throw std::runtime_error("Error loading configuration: " +
                               std::string(e.what()));
    }
  }
};

}  // namespace metada::framework::runs