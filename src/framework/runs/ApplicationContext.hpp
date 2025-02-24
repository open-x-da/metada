#pragma once
#include <stdexcept>

#include "utils/config/Config.hpp"
#include "utils/logger/Logger.hpp"
// Add other service headers as needed

namespace metada::framework::runs {

/**
 * @file ApplicationContext.hpp
 * @brief RAII wrapper managing application-wide runtime context and services
 *
 * @details
 * The ApplicationContext class provides centralized initialization, access, and
 * cleanup of core application services. It follows RAII principles to ensure
 * proper initialization order and cleanup of all managed services.
 *
 * Key features:
 * - Single point of initialization and cleanup for application services
 * - RAII design ensures proper resource management
 * - Thread-safe access to core services
 * - Prevents copying to maintain single service instances
 * - Supports moving for RAII usage
 *
 * Managed services:
 * - Logging system - Application-wide logging and diagnostics
 * - Configuration - Reading and accessing application settings
 * - Performance monitoring - Tracking application metrics (planned)
 * - Resource management - Managing system resources (planned)
 * - Runtime environment - MPI and threading configuration (planned)
 *
 * Usage:
 * @code
 * int main() {
 *   ApplicationContext ctx("MyApp", "config.yaml");
 *   auto& logger = ctx.getLogger();
 *   auto& config = ctx.getConfig();
 *
 *   logger.Info("Application started");
 *   // Use services...
 *   return 0;
 * } // Context cleaned up automatically
 * @endcode
 *
 * @note The context should be instantiated once at application startup and
 * destroyed when the application exits.
 */
template <typename Traits>
class ApplicationContext {
  using LoggerType = typename Traits::LoggerType;
  using ConfigType = typename Traits::ConfigType;

 public:
  /**
   * @brief Initialize application context with required services
   * @param app_name Name of the application for logging
   * @param config_file Optional path to configuration file
   * @throws std::runtime_error If initialization of any service fails
   */
  ApplicationContext(const std::string& app_name,
                     const std::string& config_file = "") {
    initLogger(app_name);
    if (!config_file.empty()) {
      loadConfig(config_file);
    }

    logger_.Info("Application context initialized: " + app_name);
  }

  /**
   * @brief Clean up application services in reverse initialization order
   */
  ~ApplicationContext() {
    logger_.Info("Shutting down application context");
    shutdownLogger();
  }

  // Prevent copying to ensure single instance of services
  ApplicationContext(const ApplicationContext&) = delete;
  ApplicationContext& operator=(const ApplicationContext&) = delete;

  // Allow moving to support RAII
  ApplicationContext(ApplicationContext&&) = default;
  ApplicationContext& operator=(ApplicationContext&&) = default;

  // Access to services
  LoggerType& getLogger() { return logger_; }

  ConfigType& getConfig() { return config_; }

  // Timer access will be added later
  // Timer& getTimer() { return timer_; }

 private:
  LoggerType logger_;
  ConfigType config_;
  // Timer timer_;  // To be implemented

  void initLogger(const std::string& app_name) {
    logger_.backend().Init(app_name);
  }

  void shutdownLogger() { logger_.backend().Shutdown(); }

  void loadConfig(const std::string& config_file) {
    if (!config_.LoadFromFile(config_file)) {
      throw std::runtime_error("Failed to load configuration from: " +
                               config_file);
    }
    logger_.Info("Loaded configuration from: " + config_file);
  }
};

}  // namespace metada::framework::runs