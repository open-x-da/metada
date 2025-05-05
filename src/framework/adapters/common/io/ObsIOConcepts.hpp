/**
 * @file ObsIOConcepts.hpp
 * @brief Concept definitions for observation I/O classes
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains concept definitions that constrain the types that can be
 * used with the ObsIO adapter class. These concepts ensure that backend
 * implementations provide all the necessary functionality required for reading
 * and writing observation data in different file formats like Bufr, NetCDF, and
 * HDF5.
 */

#pragma once

#include <chrono>
#include <concepts>
#include <optional>
#include <string>
#include <vector>

#include "CommonConcepts.hpp"
#include "DateTime.hpp"

namespace metada::framework {

/**
 * @brief Structure representing a single observation record
 *
 * @details Contains all essential metadata and measurement data for a single
 * observation point, including type, value, location, timestamp, and quality
 * information.
 */
struct ObservationRecord {
  std::string type;  ///< Type of observation (e.g., "temperature", "pressure")
  double value;      ///< Observed value
  std::optional<std::string> unit;        ///< Unit of the observation
  std::optional<std::string> identifier;  ///< ID of the observation
  std::string
      location;  ///< Observation location (e.g., station ID or coordinates)
  DateTime datetime;      ///< Date and time of observation
  std::size_t qc_marker;  ///< Quality control marker/flag
  std::optional<std::string>
      error_message;  ///< Optional error message or additional information
};

/**
 * @brief Concept that defines requirements for an observation I/O backend
 * implementation
 *
 * @details A valid observation I/O backend implementation must provide:
 * - Methods for reading observation data from files
 * - Methods for writing observation data to files
 * - Support for querying file format capabilities
 * - Support for checking if a file is of a supported format
 *
 * These requirements enable the ObsIO adapter to perform file I/O
 * operations with different file formats like Bufr, NetCDF, HDF5, etc.
 *
 * @tparam T The observation I/O backend implementation type
 */
template <typename T>
concept ObsIOBackendImpl =
    requires(T& t, const T& ct, const std::string& filename,
             const std::vector<ObservationRecord>& records) {
      // Reading methods
      {
        t.readObservations(filename)
      } -> std::same_as<std::vector<ObservationRecord>>;
      { ct.canRead(filename) } -> std::same_as<bool>;

      // Writing methods
      { t.writeObservations(filename, records) } -> std::same_as<void>;
      { ct.canWrite() } -> std::same_as<bool>;

      // Format information
      { ct.getFormatName() } -> std::same_as<std::string>;
      { ct.getFileExtensions() } -> std::same_as<std::vector<std::string>>;

      // Proper implementation constraints
      requires HasDeletedDefaultConstructor<T>;
      requires HasDeletedCopyConstructor<T>;
      requires HasDeletedCopyAssignment<T>;
    };

/**
 * @brief Concept that defines requirements for an observation I/O backend tag
 * type
 *
 * @details A valid backend tag must:
 * - Provide an ObsIOBackend type via BackendTraits
 * - Ensure the ObsIOBackend type satisfies the ObsIOBackendImpl
 * concept
 *
 * This concept is used to validate that a backend tag provides the necessary
 * type information and that the associated backend implementation meets the
 * requirements defined in ObsIOBackendImpl.
 *
 * @tparam T The backend tag type to check
 */
template <typename T>
concept ObsIOBackendType = requires {
  typename traits::BackendTraits<T>::ObsIOBackend;
} && ObsIOBackendImpl<typename traits::BackendTraits<T>::ObsIOBackend>;

}  // namespace metada::framework