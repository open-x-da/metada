#pragma once

#include <concepts>
#include <string>

#include "LogLevel.hpp"

namespace metada::framework {

/**
 * @brief Concept that checks if a type has a LoggerBackend type member
 */
template <typename T>
concept HasLoggerBackend =
    requires { typename traits::BackendTraits<T>::LoggerBackend; };

/**
 * @brief Concept that checks if a type has a ConfigBackend type member
 */
template <typename T>
concept HasConfigBackend =
    requires { typename traits::BackendTraits<T>::ConfigBackend; };

/**
 * @brief Concept that checks if a type has a LogMessage method
 */
template <typename T>
concept HasLogMessage =
    requires(T t, LogLevel level, const std::string& message) {
      { t.LogMessage(level, message) } -> std::same_as<void>;
    };

/**
 * @brief Concept that checks if a type has static Init and Shutdown methods
 */
template <typename T>
concept HasStaticLifecycle = requires() {
  { T::Init() } -> std::same_as<void>;
  { T::Shutdown() } -> std::same_as<void>;
};

/**
 * @brief Concept that defines the requirements for a logger backend
 *
 * A logger backend must have a LogMessage method that takes a LogLevel and
 * string message. It may optionally have static Init and Shutdown methods for
 * lifecycle management.
 */
template <typename T>
concept LoggerBackend =
    HasLogMessage<T> &&
    (HasStaticLifecycle<T> || true);  // Static lifecycle methods are optional

/**
 * @brief Concept that defines the requirements for a logger backend tag
 *
 * A logger backend tag must provide a LoggerBackend type through BackendTraits
 * that satisfies the LoggerBackend concept.
 */
template <typename T>
concept LoggerBackendTag =
    HasLoggerBackend<T> && HasConfigBackend<T> &&
    LoggerBackend<typename traits::BackendTraits<T>::LoggerBackend>;

}  // namespace metada::framework