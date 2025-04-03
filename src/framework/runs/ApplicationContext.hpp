#pragma once
#include <stdexcept>
#include <string>

#include "BackendTraits.hpp"
#include "common/utils/config/Config.hpp"
#include "common/utils/logger/Logger.hpp"

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
 *   auto context = ApplicationContext<L63BackendTag>("letkf_app", argv[1]);
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
   * @param app_name Name of the application used for logging and identification
   * @param config_file Path to configuration file (JSON/YAML)
   * @throws std::runtime_error If initialization of any service fails
   *
   * @details Initializes services in the following order:
   * 1. Configuration loading from specified file
   * 2. Logger initialization with configuration
   * 3. Logs successful initialization with application name
   */
  ApplicationContext(const std::string& app_name,
                     const std::string& config_file)
      : config_(config_file), logger_(config_) {
    logger_.Info() << "Application context initialized: " << app_name;
  }

  /**
   * @brief Destructor - cleans up application services
   *
   * @details Services are cleaned up in reverse initialization order:
   * 1. Future services (when implemented)
   * 2. Logger
   * 3. Configuration
   */
  ~ApplicationContext() { logger_.Info() << "Application context destroyed"; }

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
  Logger<BackendTag>& getLogger() { return logger_; }

  /**
   * @brief Get reference to the configuration service
   * @return Reference to the config instance
   */
  Config<BackendTag>& getConfig() { return config_; }

  // Timer access will be added later
  // Timer& getTimer() { return timer_; }

 private:
  Config<BackendTag> config_;
  Logger<BackendTag> logger_;
  // Timer timer_;  // To be implemented
};

}  // namespace metada::framework::runs