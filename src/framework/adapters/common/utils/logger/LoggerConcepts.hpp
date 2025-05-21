/**
 * @file LoggerConcepts.hpp
 * @brief Concept definitions for logger classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the Logger adapter class. These concepts ensure that backend
 * implementations provide all the necessary functionality required by the
 * logging operations.
 */

#pragma once

#include <concepts>
#include <string>

#include "CommonConcepts.hpp"
#include "LogLevel.hpp"

namespace metada::framework {

//-----------------------------------------------------------------------------
// Component concepts
//-----------------------------------------------------------------------------

/**
 * @brief Concept that defines requirements for a logger backend implementation
 *
 * @details A valid logger backend implementation must provide:
 * - A LogMessage method for writing log entries
 * - Static Init/Shutdown methods for lifecycle management
 * - Configuration constructor
 * - Proper resource management with deleted default constructor, copy
 * constructor, and copy assignment operator
 *
 * This concept is used to ensure that backend implementations provide
 * all the necessary functionality required by the Logger class.
 *
 * @tparam T The logger backend implementation type
 * @tparam ConfigBackend The configuration backend type
 *
 * @see HasDeletedDefaultConstructor
 * @see HasDeletedCopyConstructor
 * @see HasDeletedCopyAssignment
 */
template <typename T, typename ConfigBackend>
concept LoggerBackendImpl =
    requires(T& t, const ConfigBackend& config, LogLevel level,
             const std::string& message, const std::string& app_name) {
      // Core logging functionality
      { t.LogMessage(level, message) } -> std::same_as<void>;

      // Lifecycle management
      { T::Init(app_name) } -> std::same_as<void>;
      { T::Shutdown() } -> std::same_as<void>;

      // Construction and resource management
      { T(config) } -> std::same_as<T>;
      requires HasDeletedDefaultConstructor<T>;
      requires HasDeletedCopyConstructor<T>;
      requires HasDeletedCopyAssignment<T>;
    };

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