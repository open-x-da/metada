#pragma once

#include <concepts>
#include <string>

#include "LogLevel.hpp"

namespace metada::framework {

/**
 * @brief Concept that checks if a type has a LoggerBackend type member
 *
 * @details This concept verifies that a type T provides a LoggerBackend type
 * through the BackendTraits specialization.
 */
template <typename T>
concept HasLoggerBackend =
    requires { typename traits::BackendTraits<T>::LoggerBackend; };

/**
 * @brief Concept that checks if a type has a ConfigBackend type member
 *
 * @details This concept verifies that a type T provides a ConfigBackend type
 * through the BackendTraits specialization.
 */
template <typename T>
concept HasConfigBackend =
    requires { typename traits::BackendTraits<T>::ConfigBackend; };

/**
 * @brief Concept that checks if a type has a constructor taking a ConfigBackend
 *
 * @details This concept ensures that a type T can be constructed with a
 * ConfigBackend parameter, which is necessary for backend implementations
 * to be initialized with configuration data.
 */
template <typename T, typename ConfigBackend>
concept HasConfigConstructor = requires(const ConfigBackend& config) {
  { T(config) } -> std::same_as<T>;
};

/**
 * @brief Concept that checks if a type has a LogMessage method
 *
 * @details This concept verifies that a type T provides a LogMessage method
 * that takes a LogLevel and a string message, which is the core functionality
 * required for any logger backend implementation.
 */
template <typename T>
concept HasLogMessage =
    requires(T t, LogLevel level, const std::string& message) {
      { t.LogMessage(level, message) } -> std::same_as<void>;
    };

/**
 * @brief Concept that checks if a type has static Init and Shutdown methods
 *
 * @details This concept verifies that a type T provides static Init and
 * Shutdown methods, which are used for lifecycle management of logger backends
 * that require initialization and cleanup.
 */
template <typename T>
concept HasStaticLifecycle = requires() {
  { T::Init() } -> std::same_as<void>;
  { T::Shutdown() } -> std::same_as<void>;
};

/**
 * @brief Concept that defines the requirements for a logger backend
 *
 * @details A logger backend must have a LogMessage method that takes a LogLevel
 * and string message. It must also be constructible from a ConfigBackend
 * instance. It may optionally have static Init and Shutdown methods for
 * lifecycle management.
 *
 * This concept is used to constrain the types that can be used as logger
 * backends in the Logger class template.
 *
 * @tparam T The type to check
 * @tparam ConfigBackend The configuration backend type
 */
template <typename T, typename ConfigBackend>
concept LoggerBackend =
    HasLogMessage<T> && HasConfigConstructor<T, ConfigBackend> &&
    (HasStaticLifecycle<T> || true);  // Static lifecycle methods are optional

/**
 * @brief Concept that defines the requirements for a logger backend tag
 *
 * @details A logger backend tag must provide both LoggerBackend and
 * ConfigBackend types through BackendTraits. The LoggerBackend type must
 * satisfy the LoggerBackend concept when paired with the ConfigBackend type.
 *
 * This concept is used to constrain the template parameter of the Logger class,
 * ensuring that only valid backend tags can be used.
 *
 * @tparam T The type to check
 *
 * @see LoggerBackend
 * @see HasLoggerBackend
 * @see HasConfigBackend
 */
template <typename T>
concept LoggerBackendType =
    HasLoggerBackend<T> && HasConfigBackend<T> &&
    LoggerBackend<typename traits::BackendTraits<T>::LoggerBackend,
                  typename traits::BackendTraits<T>::ConfigBackend>;

}  // namespace metada::framework