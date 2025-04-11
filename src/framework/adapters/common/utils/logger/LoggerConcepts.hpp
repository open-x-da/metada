#pragma once

#include <concepts>
#include <string>

#include "../../../CommonConcepts.hpp"
#include "LogLevel.hpp"

namespace metada::framework {

//-----------------------------------------------------------------------------
// Core logging capability concepts
//-----------------------------------------------------------------------------

/**
 * @brief Concept that checks if a type implements the core logging
 * functionality
 *
 * @details Verifies that a type T provides a LogMessage method accepting a
 * LogLevel and a string message, which is the minimum requirement for any
 * logger backend implementation. This method is responsible for actually
 * writing log messages to the underlying logging system.
 *
 * @tparam T The type to check for logging capability
 */
template <typename T>
concept HasLogMessage =
    requires(T t, LogLevel level, const std::string& message) {
      { t.LogMessage(level, message) } -> std::same_as<void>;
    };

/**
 * @brief Concept that checks if a type provides lifecycle management methods
 *
 * @details Verifies that a type T has static Init and Shutdown methods
 * for initialization and cleanup of logger backends. These methods are
 * typically used to set up and tear down global logging resources.
 *
 * @tparam T The type to check for lifecycle management capability
 */
template <typename T>
concept HasStaticLifecycle = requires(const std::string& app_name) {
  { T::Init(app_name) } -> std::same_as<void>;
  { T::Shutdown() } -> std::same_as<void>;
};

//-----------------------------------------------------------------------------
// Component concepts
//-----------------------------------------------------------------------------

/**
 * @brief Concept that defines requirements for a logger backend implementation
 *
 * @details A valid logger backend must:
 * - Implement LogMessage method for handling log entries
 * - Be constructible from a ConfigBackend instance
 * - Provide static Init/Shutdown methods for lifecycle management
 * - Have deleted default constructor, copy constructor, and copy assignment
 *
 * This concept is used to ensure that backend implementations provide
 * all the necessary functionality required by the Logger class.
 *
 * Example implementations include GoogleLogger and MockLogger.
 *
 * @tparam T The logger backend implementation type
 * @tparam ConfigBackend The configuration backend type
 *
 * @see HasLogMessage
 * @see HasConfigConstructor
 * @see HasStaticLifecycle
 */
template <typename T, typename ConfigBackend>
concept LoggerBackendImpl =
    HasLogMessage<T> && HasConfigConstructor<T, ConfigBackend> &&
    HasStaticLifecycle<T> && HasDeletedDefaultConstructor<T> &&
    HasDeletedCopyConstructor<T> && HasDeletedCopyAssignment<T>;

/**
 * @brief Concept that defines requirements for a logger backend tag type
 *
 * @details A valid backend tag must:
 * - Provide both LoggerBackend and ConfigBackend types via BackendTraits
 * - Ensure the LoggerBackend type satisfies the LoggerBackendImpl concept
 *   when paired with the ConfigBackend type
 *
 * This concept constrains the template parameter of the Logger class,
 * ensuring that only valid backend configurations can be used. It provides
 * compile-time validation of backend compatibility.
 *
 * @tparam T The backend tag type to check
 *
 * @see HasLoggerBackend
 * @see HasConfigBackend
 * @see LoggerBackendImpl
 */
template <typename T>
concept LoggerBackendType =
    HasLoggerBackend<T> && HasConfigBackend<T> &&
    LoggerBackendImpl<typename traits::BackendTraits<T>::LoggerBackend,
                      typename traits::BackendTraits<T>::ConfigBackend>;

}  // namespace metada::framework