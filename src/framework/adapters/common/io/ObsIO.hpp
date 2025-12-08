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
 * - Reading observations from various formats
 * - Writing observations to various formats
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
   * @brief Read observations from the configured data source
   *
   * @details Delegates to the backend implementation to read observation
   * data. The backend will throw appropriate exceptions if the operation
   * is not supported or fails.
   *
   * @return Vector of observation records read from the data source
   * @throws std::runtime_error If the reading operation fails
   */
  std::vector<ObsRecord> read() {
    try {
      return backend_.read();
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to read observations: " +
                               std::string(e.what()));
    }
  }

  /**
   * @brief Write observations to the configured data destination
   *
   * @details Delegates to the backend implementation to write observation
   * data. The backend will throw appropriate exceptions if the operation
   * is not supported or fails.
   *
   * @param records Vector of observation records to write
   * @throws std::runtime_error If the writing operation fails
   */
  void write(const std::vector<ObsRecord>& records) {
    try {
      backend_.write(records);
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to write observations: " +
                               std::string(e.what()));
    }
  }

 private:
  ObsIOBackend backend_; /**< Backend implementation */
};

}  // namespace metada::framework