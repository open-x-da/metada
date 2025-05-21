/**
 * @file ObservationConcepts.hpp
 * @brief Concept definitions for observation classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the Observation adapter class. These concepts ensure that backend
 * implementations provide all the necessary functionality required by the
 * observation operations in data assimilation systems, including data access,
 * variable information, and arithmetic operations.
 */

#pragma once

#include <concepts>
#include <memory>
#include <string>
#include <vector>

#include "CommonConcepts.hpp"

namespace metada::framework {

/**
 * @brief Concept that defines requirements for an observation backend
 * implementation
 *
 * @details A valid observation backend implementation must provide:
 * - Data access methods (getData for both mutable and const access)
 * - Variable name access method
 * - Dimension information for each variable
 * - Support for cloning the observation
 * - Vector arithmetic operations (add, subtract, multiply)
 * - Equality comparison
 * - Initialization method
 * - Quality control application
 * - File I/O operations
 *
 * These requirements enable the Observation adapter to perform common
 * observation-space operations needed in data assimilation algorithms.
 *
 * @tparam T The observation backend implementation type
 */
template <typename T>
concept ObservationBackendImpl =
    requires(T& t, const T& ct, const T& other, const std::string& varName,
             double scalar) {
      // Data access
      { t.getData() } -> std::same_as<void*>;
      { ct.getData() } -> std::same_as<const void*>;

      // Variable information
      {
        ct.getVariableNames()
      } -> std::same_as<const std::vector<std::string>&>;
      { ct.getDimensions(varName) } -> std::same_as<const std::vector<size_t>&>;

      // Cloning
      { ct.clone() } -> std::convertible_to<std::unique_ptr<T>>;

      // Vector arithmetic
      { t.add(other) } -> std::same_as<void>;
      { t.subtract(other) } -> std::same_as<void>;
      { t.multiply(scalar) } -> std::same_as<void>;

      // Comparison
      { ct.equals(other) } -> std::convertible_to<bool>;

      // Lifecycle management
      { t.initialize() } -> std::same_as<void>;
      { t.applyQC() } -> std::same_as<void>;

      // File I/O
      { t.loadFromFile(varName) } -> std::same_as<void>;
      { ct.saveToFile(varName) } -> std::same_as<void>;
    };

/**
 * @brief Concept that defines requirements for an observation backend tag type
 *
 * @details A valid backend tag must:
 * - Provide both ObservationBackend and ConfigBackend types via BackendTraits
 * - Ensure the ObservationBackend type satisfies the ObservationBackendImpl
 *   concept when paired with the ConfigBackend type
 *
 * This concept is used by the Observation adapter to validate that a backend
 * implementation provides all necessary functionality at compile time.
 *
 * @tparam T The backend tag type to check
 *
 * @see ObservationBackendImpl
 * @see HasObservationBackend
 * @see HasConfigBackend
 */
template <typename T>
concept ObservationBackendType =
    HasObservationBackend<T> && HasConfigBackend<T> &&
    ObservationBackendImpl<
        typename traits::BackendTraits<T>::ObservationBackend>;
}  // namespace metada::framework