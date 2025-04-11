#pragma once

#include <concepts>
#include <string>

#include "../../../CommonConcepts.hpp"
#include "LogLevel.hpp"

namespace metada::framework {

//-----------------------------------------------------------------------------
// Basic trait detection concepts
//-----------------------------------------------------------------------------

/**
 * @brief Checks if a type provides a LoggerBackend type through BackendTraits
 */
template <typename T>
concept HasLoggerBackend =
    requires { typename traits::BackendTraits<T>::LoggerBackend; };

/**
 * @brief Checks if a type provides a ConfigBackend type through BackendTraits
 */
template <typename T>
concept HasConfigBackend =
    requires { typename traits::BackendTraits<T>::ConfigBackend; };

/**
 * @brief Checks if a type implements the core logging functionality
 *
 * @details Verifies that a type T provides a LogMessage method accepting a
 * LogLevel and a string message, which is the minimum requirement for any
 * logger backend implementation.
 */
template <typename T>
concept HasLogMessage =
    requires(T t, LogLevel level, const std::string& message) {
      { t.LogMessage(level, message) } -> std::same_as<void>;
    };

/**
 * @brief Checks if a type provides lifecycle management methods
 *
 * @details Verifies that a type T has static Init and Shutdown methods
 * for initialization and cleanup of logger backends.
 */
template <typename T>
concept HasStaticLifecycle = requires() {
  { T::Init() } -> std::same_as<void>;
  { T::Shutdown() } -> std::same_as<void>;
};

//-----------------------------------------------------------------------------
// Component concepts
//-----------------------------------------------------------------------------

/**
 * @brief Defines requirements for a logger backend implementation
 *
 * @details A valid logger backend must:
 * - Implement LogMessage method for handling log entries
 * - Be constructible from a ConfigBackend instance
 * - Optionally provide static Init/Shutdown methods
 *
 * @tparam T The logger backend implementation type
 * @tparam ConfigBackend The configuration backend type
 */
template <typename T, typename ConfigBackend>
concept LoggerBackend =
    HasLogMessage<T> && HasConfigConstructor<T, ConfigBackend> &&
    (HasStaticLifecycle<T> || true);  // Static lifecycle methods are optional

/**
 * @brief Defines requirements for a logger backend tag type
 *
 * @details A valid backend tag must:
 * - Provide both LoggerBackend and ConfigBackend types via BackendTraits
 * - Ensure the LoggerBackend type satisfies the LoggerBackend concept
 *   when paired with the ConfigBackend type
 *
 * This concept constrains the template parameter of the Logger class,
 * ensuring that only valid backend configurations can be used.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept LoggerBackendType =
    HasLoggerBackend<T> && HasConfigBackend<T> &&
    LoggerBackend<typename traits::BackendTraits<T>::LoggerBackend,
                  typename traits::BackendTraits<T>::ConfigBackend>;

}  // namespace metada::framework