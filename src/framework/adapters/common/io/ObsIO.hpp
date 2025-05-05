/**
 * @file ObsIO.hpp
 * @brief Adapter class for observation I/O operations
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains the ObsIO adapter class that provides a unified
 * interface for reading and writing observation data in different file formats.
 * It leverages C++20 concepts to ensure backends implement the required
 * interface.
 */

#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "BackendTraits.hpp"
#include "Config.hpp"
#include "NonCopyable.hpp"
#include "ObsIOConcepts.hpp"

namespace metada::framework {

/**
 * @brief Adapter class for observation I/O operations
 *
 * @details This class provides a unified interface for reading and writing
 * observation data in different file formats like Bufr, NetCDF, and HDF5.
 * It wraps a backend implementation that satisfies the ObsIOBackendType
 * concept, providing consistent access patterns across different backends.
 *
 * The ObsIO class inherits from NonCopyable to prevent unintended
 * copying, but supports move semantics for efficient resource management.
 *
 * Features:
 * - File format detection and validation
 * - Reading observations from various file formats
 * - Writing observations to various file formats
 * - Format metadata access (extensions, capabilities)
 * - Move semantics for efficient resource management
 *
 * @tparam BackendTag The backend tag type that must satisfy
 * ObsIOBackendType concept
 *
 * @see ObsIOBackendType
 * @see NonCopyable
 */
template <typename BackendTag>
  requires ObsIOBackendType<BackendTag>
class ObsIO : public NonCopyable {
 public:
  /** @brief Type alias for the backend implementation type */
  using ObsIOBackend = typename traits::BackendTraits<BackendTag>::ObsIOBackend;

  /** @brief Default constructor deleted - must provide initialization params */
  ObsIO() = delete;

  /**
   * @brief Constructor that initializes with format-specific parameters
   *
   * @details Creates and initializes the observation I/O backend with
   * the provided parameters. The exact meaning of the parameters depends
   * on the specific backend implementation.
   *
   * @param params Format-specific initialization parameters
   */
  explicit ObsIO(Config<BackendTag>&& config)
      : backend_(std::move(config.backend())) {}

  /**
   * @brief Move constructor
   *
   * @details Transfers ownership of the backend from another observation I/O
   * object.
   *
   * @param other The observation I/O object to move from
   */
  ObsIO(ObsIO&& other) noexcept : backend_(std::move(other.backend_)) {}

  /**
   * @brief Move assignment operator
   *
   * @details Transfers ownership of the backend from another observation I/O
   * object.
   *
   * @param other The observation I/O object to move from
   * @return Reference to this observation I/O object after assignment
   */
  ObsIO& operator=(ObsIO&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
    }
    return *this;
  }

  /**
   * @brief Get const reference to the backend implementation
   *
   * @return Const reference to the backend implementation
   */
  const ObsIOBackend& backend() const { return backend_; }

  /**
   * @brief Get mutable reference to the backend implementation
   *
   * @return Mutable reference to the backend implementation
   */
  ObsIOBackend& backend() { return backend_; }

  /**
   * @brief Check if the backend can read the specified file
   *
   * @details Delegates to the backend implementation to check if it
   * supports reading the specified file format.
   *
   * @param filename Path to the file to check
   * @return True if the backend can read the file, false otherwise
   */
  bool canRead(const std::string& filename) const {
    return backend_.canRead(filename);
  }

  /**
   * @brief Check if the backend can write observations
   *
   * @details Delegates to the backend implementation to check if it
   * supports writing observations in its format.
   *
   * @return True if the backend can write observations, false otherwise
   */
  bool canWrite() const { return backend_.canWrite(); }

  /**
   * @brief Get the name of the file format supported by this backend
   *
   * @details Delegates to the backend implementation to get the name
   * of the file format it supports.
   *
   * @return The name of the file format (e.g., "Bufr", "NetCDF", "HDF5")
   */
  std::string getFormatName() const { return backend_.getFormatName(); }

  /**
   * @brief Get the file extensions supported by this backend
   *
   * @details Delegates to the backend implementation to get the list
   * of file extensions it supports.
   *
   * @return Vector of supported file extensions (e.g., {".bufr", ".nc", ".h5"})
   */
  std::vector<std::string> getFileExtensions() const {
    return backend_.getFileExtensions();
  }

  /**
   * @brief Read observations from a file
   *
   * @details Delegates to the backend implementation to read observation
   * data from the specified file. Throws an exception if the file cannot
   * be read or is not in a supported format.
   *
   * @param filename Path to the file to read
   * @return Vector of observation records read from the file
   * @throws std::runtime_error If the file cannot be read
   */
  std::vector<ObsRecord> read(const std::string& filename) {
    if (!canRead(filename)) {
      throw std::runtime_error("Cannot read file format: " + filename);
    }

    try {
      return backend_.read(filename);
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to read observations: " +
                               std::string(e.what()));
    }
  }

  /**
   * @brief Write observations to a file
   *
   * @details Delegates to the backend implementation to write observation
   * data to the specified file. Throws an exception if the backend does
   * not support writing or if writing fails.
   *
   * @param filename Path to the file to write
   * @param records Vector of observation records to write
   * @throws std::runtime_error If writing fails or is not supported
   */
  void write(const std::string& filename,
             const std::vector<ObsRecord>& records) {
    if (!canWrite()) {
      throw std::runtime_error(
          "This backend does not support writing observations");
    }

    try {
      backend_.write(filename, records);
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to write observations: " +
                               std::string(e.what()));
    }
  }

 private:
  ObsIOBackend backend_; /**< Backend implementation */
};

}  // namespace metada::framework