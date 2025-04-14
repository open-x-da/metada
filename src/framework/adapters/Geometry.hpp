/**
 * @file Geometry.hpp
 * @brief Template class providing a generic interface to geometry
 * implementations
 * @ingroup repr
 * @author Metada Framework Team
 *
 * @details
 * This header provides a template class that wraps geometry backend
 * implementations and provides a unified interface for geometry operations. The
 * Geometry class delegates operations to the backend while providing type
 * safety and a consistent API.
 *
 * The Geometry class template is designed to:
 * - Provide a generic interface to different geometry backend implementations
 * - Support initialization from configuration objects
 * - Enable grid traversal through iterators
 * - Implement proper move semantics (non-copyable)
 * - Support halo exchange operations
 * - Allow geometry information queries
 *
 * @see Config
 * @see State
 * @see GeometryIterator
 * @see GeometryBackendType
 */

#pragma once
#include <stdexcept>  // for exceptions (if needed)
#include <utility>    // for std::move

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "GeometryConcepts.hpp"  // Added include for GeometryBackendType concept
#include "GeometryIterator.hpp"
#include "NonCopyable.hpp"

namespace metada::framework {

// Forward declarations
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;  // Config adapter for backend

/**
 * @brief Main geometry class template providing a generic interface to geometry
 * implementations
 *
 * @details
 * This class template wraps a geometry backend implementation and provides a
 * type-safe interface for all geometry operations. It delegates operations to
 * the backend while adding:
 * - Type safety through templates
 * - Proper move semantics (non-copyable)
 * - Consistent interface across different backends
 * - Support for grid traversal through iterators
 * - Integration with state operations like halo exchange
 *
 * The backend tag must satisfy the GeometryBackendType concept, which ensures
 * it provides valid backend implementation types through BackendTraits.
 *
 * Example usage:
 * @code
 * Config<T> config;
 * Geometry<Backend> geometry(config);
 * for (auto& point : geometry) {
 *     // Process each grid point
 * }
 * @endcode
 *
 * @tparam BackendTag The tag type that defines the geometry backend through
 * BackendTraits
 *
 * @see GeometryBackendType
 * @see Config
 * @see State
 * @see GeometryIterator
 */
template <typename BackendTag>
  requires GeometryBackendType<BackendTag>
class Geometry : private NonCopyable {
 public:
  using GeometryBackend =
      typename traits::BackendTraits<BackendTag>::GeometryBackend;
  using Iterator =
      class GeometryIterator<BackendTag>;  // defined in GeometryIterator.hpp

  // Disable default constructor – a Geometry must be created with a config or
  // backend.
  Geometry() = delete;

  // Destructor (default ok as backend will clean up)
  ~Geometry() = default;

  /**
   * @brief Main constructor: build Geometry from a configuration
   *
   * @param[in] config Configuration object containing initialization parameters
   * @throws std::runtime_error If backend initialization fails
   */
  explicit Geometry(const Config<BackendTag>& config)
      : config_(config),  // store a copy of configuration
        backend_(
            config.backend()),  // initialize backend with config's backend data
        initialized_(true) {
    // If backend initialization fails, consider throwing an exception
    if (!initialized_) {
      throw std::runtime_error("Geometry backend initialization failed");
    }
  }

  /**
   * @brief Constructor with config and existing backend
   *
   * @param[in] config Configuration object
   * @param[in] backend Existing geometry backend to use
   */
  Geometry(const Config<BackendTag>& config, GeometryBackend& backend)
      : config_(config), backend_(backend), initialized_(true) {}

  /**
   * @brief Move constructor
   * @param other Geometry instance to move from
   */
  Geometry(Geometry&& other) noexcept
      : config_(std::move(other.config_)),
        backend_(std::move(other.backend_)),
        initialized_(other.initialized_) {
    other.initialized_ = false;  // mark moved-from as uninitialized
  }

  /**
   * @brief Move assignment operator
   * @param other Geometry instance to move from
   * @return Reference to this instance
   */
  Geometry& operator=(Geometry&& other) noexcept {
    if (this != &other) {
      config_ = std::move(other.config_);
      backend_ = std::move(other.backend_);
      initialized_ = other.initialized_;
      other.initialized_ = false;
    }
    return *this;
  }

  /**
   * @brief Clone method: create a new Geometry with the same configuration and
   * backend state
   * @return A new Geometry instance with identical configuration and backend
   * state
   */
  Geometry clone() const {
    // Use backend's clone functionality to duplicate the geometry backend
    GeometryBackend clonedBackend = *backend_.clone();
    // Construct a new Geometry with the cloned backend and same config
    return Geometry(std::move(clonedBackend), config_);
  }

  /**
   * @brief Get iterator to the beginning of the grid
   * @return Iterator pointing to the first grid point
   */
  Iterator begin() { return Iterator(backend_.begin()); }

  /**
   * @brief Get iterator to the end of the grid
   * @return Iterator pointing past the last grid point
   */
  Iterator end() { return Iterator(backend_.end()); }

  /**
   * @brief Get const iterator to the beginning of the grid
   * @return Const iterator pointing to the first grid point
   */
  Iterator begin() const { return Iterator(backend_.begin()); }

  /**
   * @brief Get const iterator to the end of the grid
   * @return Const iterator pointing past the last grid point
   */
  Iterator end() const { return Iterator(backend_.end()); }

  /**
   * @brief Get the total number of grid points
   * @return Total number of grid points in the geometry
   */
  size_t totalGridSize() const {
    return backend_
        .size();  // assume backend_ provides total number of grid points
  }
  // Example: query dimension sizes (if applicable)
  // int nx() const { return backend_.nx(); }
  // int ny() const { return backend_.ny(); }
  // int nz() const { return backend_.nz(); }

  /**
   * @brief Check if the geometry is periodic in X dimension
   * @return True if periodic in X, false otherwise
   */
  bool isPeriodicX() const { return backend_.isPeriodicX(); }

  /**
   * @brief Check if the geometry is periodic in Y dimension
   * @return True if periodic in Y, false otherwise
   */
  bool isPeriodicY() const { return backend_.isPeriodicY(); }

  /**
   * @brief Check if the geometry is periodic in Z dimension
   * @return True if periodic in Z, false otherwise
   */
  bool isPeriodicZ() const { return backend_.isPeriodicZ(); }

  /**
   * @brief Check if the geometry is periodic in the given dimension
   * @param dimension Dimension index (0=X, 1=Y, 2=Z)
   * @return True if periodic in the specified dimension, false otherwise
   */
  bool isPeriodic(int dimension) const {
    return backend_.isPeriodic(dimension);
  }

  /**
   * @brief Get the number of dimensions in the geometry
   * @return Number of dimensions
   */
  int getDimensions() const { return backend_.getDimensions(); }

  /**
   * @brief Get the size along a specific dimension
   * @param dimension Dimension index (0=X, 1=Y, 2=Z)
   * @return Size along the specified dimension
   */
  int getSize(int dimension) const { return backend_.getSize(dimension); }

  /**
   * @brief Get the total size (number of grid points)
   * @return Total number of grid points
   */
  int getTotalSize() const { return backend_.getTotalSize(); }

  /**
   * @brief Perform halo exchange on a State using this geometry
   * @param state The state on which to perform halo exchange
   */
  template <typename StateType>
    requires requires(StateType state) { state.backend(); }
  void haloExchange(StateType& state) const {
    backend_.haloExchange(state.backend());
  }

  // (Additional geometry-related methods can be added as needed, e.g.,
  // coordinate transformations)

  /**
   * @brief Check if geometry is properly initialized
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Access the underlying backend
   * @return Reference to the backend implementation
   */
  GeometryBackend& backend() { return backend_; }

  /**
   * @brief Access the underlying backend (const version)
   * @return Const reference to the backend implementation
   */
  const GeometryBackend& backend() const { return backend_; }

  /**
   * @brief Access the stored configuration
   * @return Const reference to the configuration
   */
  const Config<BackendTag>& config() const { return config_; }

 private:
  /**
   * @brief Private constructor used by clone(): create Geometry from an
   * existing backend instance
   * @param backend Backend instance to use
   * @param config Configuration to associate with this geometry
   */
  Geometry(GeometryBackend&& backend, const Config<BackendTag>& config)
      : config_(config), backend_(std::move(backend)), initialized_(true) {}

  const Config<BackendTag>&
      config_;  // Holds reference to externally managed config
  GeometryBackend
      backend_;  // The actual geometry implementation (grid data/operations)
  bool initialized_{false};  // Flag to indicate if geometry is initialized
};
}  // namespace metada::framework
