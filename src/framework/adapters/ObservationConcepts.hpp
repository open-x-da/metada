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
 * - Iteration capabilities (begin, end, size, operator[])
 * - Data access methods (getData for both mutable and const access)
 * - Template data access method for backward compatibility
 * - Variable name access method
 * - Support for cloning the observation
 * - Vector arithmetic operations (add, subtract, multiply)
 * - Equality comparison
 * - Initialization method
 * - Quality control application
 * - File I/O operations
 * - Covariance matrix access
 * - Geographic filtering operations
 *
 * These requirements enable the Observation adapter to perform common
 * observation-space operations needed in data assimilation algorithms.
 *
 * @tparam T The observation backend implementation type
 */
template <typename T>
concept ObservationBackendImpl =
    requires(T& t, const T& ct, const T& other, const std::string& typeName,
             const std::string& varName, const std::string& filename,
             double scalar, size_t index, double min_lat, double max_lat,
             double min_lon, double max_lon, double min_level, double max_level,
             double error, double missing_value) {
      // Iteration capabilities
      { ct.begin() } -> std::convertible_to<typename T::iterator_type>;
      { ct.end() } -> std::convertible_to<typename T::iterator_type>;
      { ct.size() } -> std::convertible_to<size_t>;
      { ct[index] } -> std::convertible_to<const typename T::value_type&>;

      // Data access
      { t.getData() } -> std::same_as<void*>;
      { ct.getData() } -> std::same_as<const void*>;

      // Template data access for backward compatibility
      {
        ct.template getData<std::vector<double>>()
      } -> std::convertible_to<std::vector<double>>;

      // Variable information
      { ct.getTypeNames() } -> std::convertible_to<std::vector<std::string>>;
      {
        ct.getVariableNames(typeName)
      } -> std::convertible_to<std::vector<std::string>>;

      // Covariance matrix
      { ct.getCovariance() } -> std::convertible_to<const std::vector<double>&>;

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
      { t.loadFromFile(filename, error, missing_value) } -> std::same_as<void>;
      { ct.saveToFile(filename) } -> std::same_as<void>;

      // Geographic filtering
      {
        ct.getObservationsInBox(min_lat, max_lat, min_lon, max_lon)
      } -> std::convertible_to<std::vector<typename T::value_type>>;
      {
        ct.getObservationsInVerticalRange(min_level, max_level)
      } -> std::convertible_to<std::vector<typename T::value_type>>;
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