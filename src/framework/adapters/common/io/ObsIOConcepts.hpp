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
#include "ObsRecord.hpp"

namespace metada::framework {

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
 * @tparam ConfigBackend The configuration backend type
 */
template <typename T, typename ConfigBackend>
concept ObsIOBackendImpl =
    requires(T& t, const T& ct, const std::string& filename,
             const std::vector<ObsRecord>& records, ConfigBackend&& config) {
      // Constructor from Config
      { T(std::move(config)) } -> std::same_as<T>;

      // Reading methods
      { t.read(filename) } -> std::same_as<std::vector<ObsRecord>>;
      { ct.canRead(filename) } -> std::same_as<bool>;

      // Writing methods
      { t.write(filename, records) } -> std::same_as<void>;
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
concept ObsIOBackendType =
    HasObsIOBackend<T> && HasConfigBackend<T> &&
    ObsIOBackendImpl<typename traits::BackendTraits<T>::ObsIOBackend,
                     typename traits::BackendTraits<T>::ConfigBackend>;
}  // namespace metada::framework