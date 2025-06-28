/**
 * @file Geometry.hpp
 * @brief Template class providing a generic interface to geometry
 * implementations
 * @ingroup adapters
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
 * - Support cloning for creating identical geometry instances
 *
 * @see Config
 * @see State
 * @see GeometryIterator
 * @see GeometryBackendType
 */

#pragma once

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "GeometryConcepts.hpp"
#include "GeometryIterator.hpp"
#include "Logger.hpp"
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
 * - Periodicity queries for boundary conditions
 * - Size information about the underlying grid
 *
 * The backend tag must satisfy the GeometryBackendType concept, which ensures
 * it provides valid backend implementation types through BackendTraits.
 *
 * Example usage:
 * @code
 * Config<BackendTag> config;
 * Geometry<BackendTag> geometry(config);
 *
 * // Iterate through grid points
 * for (auto& point : geometry) {
 *     // Process each grid point
 * }
 *
 * // Perform halo exchange on a state
 * State<BackendTag> state(config);
 * geometry.haloExchange(state);
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
  using value_type = typename GeometryBackend::value_type;
  using reference = typename GeometryBackend::reference;
  using const_reference = typename GeometryBackend::const_reference;
  using pointer = typename GeometryBackend::pointer;
  using const_pointer = typename GeometryBackend::const_pointer;
  using size_type = typename GeometryBackend::size_type;
  using difference_type = typename GeometryBackend::difference_type;
  using iterator = typename GeometryBackend::iterator;
  using const_iterator = typename GeometryBackend::const_iterator;

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
      : backend_(config.backend()) {
    logger_.Debug() << "Geometry constructed";
  }

  /**
   * @brief Move constructor
   * @param other Geometry instance to move from
   */
  Geometry(Geometry&& other) noexcept : backend_(std::move(other.backend_)) {}

  /**
   * @brief Move assignment operator
   * @param other Geometry instance to move from
   * @return Reference to this instance
   */
  Geometry& operator=(Geometry&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
    }
    return *this;
  }

  /**
   * @brief Clone method: create a new Geometry with the same configuration and
   * backend state
   * @return A new Geometry instance with identical configuration and backend
   * state
   */
  Geometry clone() const { return Geometry(std::move(backend_.clone())); }

  /**
   * @brief Get iterator to the beginning of the grid
   * @return Iterator pointing to the first grid point
   */
  iterator begin() { return backend_.begin(); }

  /**
   * @brief Get iterator to the end of the grid
   * @return Iterator pointing past the last grid point
   */
  iterator end() { return backend_.end(); }

  /**
   * @brief Get const iterator to the beginning of the grid
   * @return Const iterator pointing to the first grid point
   */
  const_iterator begin() const { return backend_.begin(); }

  /**
   * @brief Get const iterator to the end of the grid
   * @return Const iterator pointing past the last grid point
   */
  const_iterator end() const { return backend_.end(); }

  /**
   * @brief Get const iterator to the beginning of the grid
   * @return Const iterator pointing to the first grid point
   */
  const_iterator cbegin() const { return backend_.cbegin(); }

  /**
   * @brief Get const iterator to the end of the grid
   * @return Const iterator pointing past the last grid point
   */
  const_iterator cend() const { return backend_.cend(); }

  /**
   * @brief Get the total number of grid points
   * @return Total number of grid points in the geometry
   */
  size_type size() const { return backend_.size(); }

  /**
   * @brief Check if the geometry is empty
   * @return True if empty, false otherwise
   */
  bool empty() const { return backend_.empty(); }

  /**
   * @brief Get the maximum number of grid points
   * @return Maximum number of grid points in the geometry
   */
  size_type max_size() const { return backend_.max_size(); }

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

  // STL-compliant element access
  reference operator[](size_type idx) { return backend_[idx]; }
  const_reference operator[](size_type idx) const { return backend_[idx]; }
  reference at(size_type idx) { return backend_.at(idx); }
  const_reference at(size_type idx) const { return backend_.at(idx); }
  reference front() { return backend_.front(); }
  const_reference front() const { return backend_.front(); }
  reference back() { return backend_.back(); }
  const_reference back() const { return backend_.back(); }

 private:
  /**
   * @brief Private constructor used by clone(): create Geometry from an
   * existing backend instance
   * @param backend Backend instance to use
   */
  Geometry(GeometryBackend&& backend) : backend_(std::move(backend)) {}

  GeometryBackend backend_;
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};
}  // namespace metada::framework
