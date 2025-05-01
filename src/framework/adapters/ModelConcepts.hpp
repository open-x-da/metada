/**
 * @file ModelConcepts.hpp
 * @brief Concept definitions for model classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the Model adapter class. These concepts ensure that backend
 * implementations provide all the necessary functionality required by the
 * model operations.
 */

#pragma once

#include <concepts>
#include <string>

#include "CommonConcepts.hpp"
#include "DateTime.hpp"

namespace metada::framework {

using base::DateTime;

//-----------------------------------------------------------------------------
// Component concepts
//-----------------------------------------------------------------------------

/**
 * @brief Concept that defines requirements for a model backend implementation
 *
 * @details A valid model backend implementation must provide:
 * - Configuration construction and initialization
 * - Methods for reset and finalization
 * - Parameter access (get/set) by name
 * - A run method that executes the model over a time interval
 * - Proper resource management with deleted default constructor, copy
 * constructor, and copy assignment operator
 *
 * This concept is used to ensure that backend implementations provide
 * all the necessary functionality required by the Model class.
 *
 * @tparam T The model backend implementation type
 * @tparam ConfigBackend The configuration backend type
 * @tparam StateBackend The state backend type
 *
 * @see HasDeletedDefaultConstructor
 * @see HasDeletedCopyConstructor
 * @see HasDeletedCopyAssignment
 */
template <typename T, typename ConfigBackend, typename StateBackend>
concept ModelBackendImpl =
    requires(T& model, const T& const_model, const ConfigBackend& config,
             const StateBackend& initialState, StateBackend& finalState,
             DateTime startTime, DateTime endTime, const std::string& name,
             const std::string& value) {
      // Construction and initialization
      { T(config) } -> std::same_as<T>;
      { model.initialize(config) } -> std::same_as<void>;

      // Lifecycle management
      { model.reset() } -> std::same_as<void>;
      { model.finalize() } -> std::same_as<void>;

      // Parameter access
      { const_model.getParameter(name) } -> std::convertible_to<std::string>;
      { model.setParameter(name, value) } -> std::same_as<void>;

      // Model execution
      {
        model.run(initialState, finalState, startTime, endTime)
      } -> std::same_as<void>;

      // Resource management constraints
      requires HasDeletedDefaultConstructor<T>;
      requires HasDeletedCopyConstructor<T>;
      requires HasDeletedCopyAssignment<T>;
    };

/**
 * @brief Concept that defines requirements for a model backend tag type
 *
 * @details A valid backend tag must:
 * - Provide ModelBackend, ConfigBackend, and StateBackend types via
 * BackendTraits
 * - Ensure the ModelBackend type satisfies the ModelBackendImpl concept
 *   when paired with the ConfigBackend and StateBackend types
 *
 * This concept constrains the template parameter of the Model class,
 * ensuring that only valid backend configurations can be used. It provides
 * compile-time validation of backend compatibility.
 *
 * @tparam T The backend tag type to check
 *
 * @see HasModelBackend
 * @see HasConfigBackend
 * @see HasStateBackend
 * @see ModelBackendImpl
 */
template <typename T>
concept ModelBackendType =
    HasModelBackend<T> && HasConfigBackend<T> && HasStateBackend<T> &&
    ModelBackendImpl<typename traits::BackendTraits<T>::ModelBackend,
                     typename traits::BackendTraits<T>::ConfigBackend,
                     typename traits::BackendTraits<T>::StateBackend>;

}  // namespace metada::framework