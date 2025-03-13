/**
 * @file Geometry.hpp
 * @brief Implementation of a Geometry adapter
 * @ingroup adapters
 *
 * @details
 * This header provides a concrete implementation of a Geometry adapter.
 * It represents the grid structure and physical domain used in numerical
 * models and supports traversal of all grid points via iterators.
 */

#pragma once

// Include dependencies
#include "GeometryPointIterator.hpp"
#include "IGeometry.hpp"

// Standard library includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace metada::framework {

// Forward declarations
class IConfig;

/**
 * @brief Concrete implementation of a Geometry adapter
 *
 * @details
 * This class implements a geometry adapter that delegates all operations
 * to a backend implementation of IGeometry. It uses composition rather than
 * inheritance, following the adapter pattern.
 *
 * @tparam Backend Type that implements IGeometry interface
 */
template <typename Backend>
class Geometry {
 public:
  /**
   * @brief Default constructor creates geometry with a default backend
   * @note Backend must be default constructible
   */
  Geometry() : backend_ptr_(std::make_unique<Backend>()) {}

  /**
   * @brief Construct a geometry with a specific backend instance
   *
   * @param backend Backend implementation
   */
  explicit Geometry(const Backend& backend)
      : backend_ptr_(std::make_unique<Backend>()) {
    // Instead of copying, use backend methods to get properties and set them
    if (backend.getDimensions() > 0) {
      // Copy dimensions
      std::vector<size_t> resolution = backend.getResolution();
      backend_ptr_->setDimensions(backend.getDimensions());
      backend_ptr_->setResolution(resolution);

      // Copy domain bounds
      for (size_t i = 0; i < backend.getDimensions(); ++i) {
        auto bounds = backend.getDomainBounds(i);
        backend_ptr_->setDomainBounds(i, bounds[0], bounds[1]);
      }
    }
  }

  /**
   * @brief Construct a geometry with a specific backend pointer
   *
   * @param backend_ptr Raw pointer to backend implementation (ownership is
   * taken)
   */
  explicit Geometry(Backend* backend_ptr)
      : backend_ptr_(backend_ptr ? backend_ptr : new Backend()) {}

  /**
   * @brief Construct a geometry with specified dimensions and resolution
   *
   * @param dimensions Number of spatial dimensions
   * @param resolution Vector of grid points along each dimension
   */
  Geometry(size_t dimensions, const std::vector<size_t>& resolution)
      : backend_ptr_(std::make_unique<Backend>()) {
    backend_ptr_->setDimensions(dimensions);
    backend_ptr_->setResolution(resolution);
  }

  /**
   * @brief Construct a geometry with dimensions, resolution, and domain bounds
   *
   * @param dimensions Number of spatial dimensions
   * @param resolution Vector of grid points along each dimension
   * @param min_bounds Minimum value for each dimension
   * @param max_bounds Maximum value for each dimension
   */
  Geometry(size_t dimensions, const std::vector<size_t>& resolution,
           const std::vector<double>& min_bounds,
           const std::vector<double>& max_bounds)
      : backend_ptr_(std::make_unique<Backend>()) {
    backend_ptr_->setDimensions(dimensions);
    backend_ptr_->setResolution(resolution);

    // Set domain bounds for each dimension
    for (size_t i = 0; i < dimensions; ++i) {
      backend_ptr_->setDomainBounds(i, min_bounds[i], max_bounds[i]);
    }
  }

  /**
   * @brief Construct from configuration
   *
   * @param config Configuration object
   */
  explicit Geometry(const IConfig& config)
      : backend_ptr_(std::make_unique<Backend>()) {
    initialize(config);
  }

  /**
   * @brief Destructor
   */
  ~Geometry() = default;

  /**
   * @brief Initialize geometry from configuration
   *
   * @param config Configuration object containing geometry parameters
   */
  void initialize(const IConfig& config) { backend_ptr_->initialize(config); }

  /**
   * @brief Get the number of dimensions
   *
   * @return Number of spatial dimensions
   */
  size_t getDimensions() const { return backend_ptr_->getDimensions(); }

  /**
   * @brief Get resolution along each dimension
   *
   * @return Vector containing number of grid points along each dimension
   */
  std::vector<size_t> getResolution() const {
    return backend_ptr_->getResolution();
  }

  /**
   * @brief Get total number of grid points
   *
   * @return Total grid point count
   */
  size_t getTotalPoints() const { return backend_ptr_->getTotalPoints(); }

  /**
   * @brief Check if a point is within the domain
   *
   * @param coordinates Vector of coordinates to check
   * @return True if point is within domain, false otherwise
   */
  bool containsPoint(const std::vector<double>& coordinates) const {
    return backend_ptr_->containsPoint(coordinates);
  }

  /**
   * @brief Get domain bounds for a specific dimension
   *
   * @param dimension Dimension index (0-based)
   * @return Array containing min and max values for the dimension
   */
  std::array<double, 2> getDomainBounds(size_t dimension) const {
    return backend_ptr_->getDomainBounds(dimension);
  }

  /**
   * @brief Set domain bounds for a specific dimension
   *
   * @param dimension Dimension index (0-based)
   * @param min_value Minimum value
   * @param max_value Maximum value
   */
  void setDomainBounds(size_t dimension, double min_value, double max_value) {
    backend_ptr_->setDomainBounds(dimension, min_value, max_value);
  }

  /**
   * @brief Get grid spacing for each dimension
   *
   * @return Vector containing grid spacing for each dimension
   */
  std::vector<double> getGridSpacing() const {
    return backend_ptr_->getGridSpacing();
  }

  /**
   * @brief Set the resolution along each dimension
   *
   * @param resolution Vector containing new resolution for each dimension
   */
  void setResolution(const std::vector<size_t>& resolution) {
    backend_ptr_->setResolution(resolution);
  }

  /**
   * @brief Get coordinates at a specific grid index
   *
   * @param indices Grid indices
   * @return Vector of physical coordinates
   */
  std::vector<double> getCoordinates(const std::vector<size_t>& indices) const {
    return backend_ptr_->getCoordinates(indices);
  }

  /**
   * @brief Get grid indices for a physical coordinate
   *
   * @param coordinates Physical coordinates
   * @return Vector of grid indices
   */
  std::vector<size_t> getIndices(const std::vector<double>& coordinates) const {
    return backend_ptr_->getIndices(coordinates);
  }

  /**
   * @brief Get iterator to the beginning of the geometry grid
   *
   * @return Iterator pointing to the first grid point
   */
  std::unique_ptr<IGeometryIterator<double>> begin() const {
    auto dimensions = getDimensions();
    std::vector<size_t> start_position(dimensions, 0);
    std::vector<double> start_coordinates = getCoordinates(start_position);

    return std::make_unique<GeometryPointIterator<double>>(
        start_position, getResolution(), start_coordinates, this);
  }

  /**
   * @brief Get iterator to the end of the geometry grid
   *
   * @return Iterator pointing past the last grid point
   */
  std::unique_ptr<IGeometryIterator<double>> end() const {
    return std::make_unique<GeometryPointIterator<double>>();
  }

  /**
   * @brief Create a new geometry with the same configuration
   *
   * @return Unique pointer to a new geometry instance
   */
  std::unique_ptr<Geometry<Backend>> clone() const {
    auto backend_clone = backend_ptr_->clone();
    Backend* typed_backend = static_cast<Backend*>(backend_clone.release());
    return std::make_unique<Geometry<Backend>>(typed_backend);
  }

  /**
   * @brief Check if this geometry is compatible with another geometry
   *
   * @param other Geometry to compare with
   * @return True if geometries are compatible, false otherwise
   */
  bool isCompatible(const Geometry<Backend>& other) const {
    return backend_ptr_->isCompatible(*other.backend_ptr_);
  }

 private:
  /**
   * @brief Backend implementation
   */
  std::unique_ptr<Backend> backend_ptr_;
};

}  // namespace metada::framework