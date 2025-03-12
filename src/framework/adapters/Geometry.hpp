/**
 * @file Geometry.hpp
 * @brief Implementation of a Geometry adapter
 * @ingroup adapters
 *
 * @details
 * This file provides a concrete implementation of a Geometry adapter.
 * It represents the grid structure and physical domain used in numerical
 * models and supports traversal of all grid points via iterators.
 */

#pragma once

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

#include "GeometryPointIterator.hpp"
#include "IGeometry.hpp"
#include "utils/config/IConfig.hpp"

namespace metada::framework {

// Forward declarations
template <typename Backend>
class Geometry;

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
  Geometry() : backend_() {}

  /**
   * @brief Construct a geometry with a specific backend instance
   *
   * @param backend Backend implementation
   */
  explicit Geometry(const Backend& backend) : backend_(backend) {}

  /**
   * @brief Construct a geometry with specified dimensions and resolution
   *
   * @param dimensions Number of spatial dimensions
   * @param resolution Vector of grid points along each dimension
   */
  Geometry(size_t dimensions, const std::vector<size_t>& resolution) {
    backend_.setDimensions(dimensions);
    backend_.setResolution(resolution);
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
           const std::vector<double>& max_bounds) {
    backend_.setDimensions(dimensions);
    backend_.setResolution(resolution);

    // Set domain bounds for each dimension
    for (size_t i = 0; i < dimensions; ++i) {
      backend_.setDomainBounds(i, min_bounds[i], max_bounds[i]);
    }
  }

  /**
   * @brief Construct from configuration
   *
   * @param config Configuration object
   */
  explicit Geometry(const IConfig& config) { initialize(config); }

  /**
   * @brief Destructor
   */
  ~Geometry() = default;

  /**
   * @brief Initialize geometry from configuration
   *
   * @param config Configuration object containing geometry parameters
   */
  void initialize(const IConfig& config) { backend_.initialize(config); }

  /**
   * @brief Get the number of dimensions
   *
   * @return Number of spatial dimensions
   */
  size_t getDimensions() const { return backend_.getDimensions(); }

  /**
   * @brief Get resolution along each dimension
   *
   * @return Vector containing number of grid points along each dimension
   */
  std::vector<size_t> getResolution() const { return backend_.getResolution(); }

  /**
   * @brief Get total number of grid points
   *
   * @return Total grid point count
   */
  size_t getTotalPoints() const { return backend_.getTotalPoints(); }

  /**
   * @brief Check if a point is within the domain
   *
   * @param coordinates Vector of coordinates to check
   * @return True if point is within domain, false otherwise
   */
  bool containsPoint(const std::vector<double>& coordinates) const {
    return backend_.containsPoint(coordinates);
  }

  /**
   * @brief Get domain bounds for a specific dimension
   *
   * @param dimension Dimension index (0-based)
   * @return Array containing min and max values for the dimension
   */
  std::array<double, 2> getDomainBounds(size_t dimension) const {
    return backend_.getDomainBounds(dimension);
  }

  /**
   * @brief Set domain bounds for a specific dimension
   *
   * @param dimension Dimension index (0-based)
   * @param min_value Minimum value
   * @param max_value Maximum value
   */
  void setDomainBounds(size_t dimension, double min_value, double max_value) {
    backend_.setDomainBounds(dimension, min_value, max_value);
  }

  /**
   * @brief Get grid spacing for each dimension
   *
   * @return Vector containing grid spacing for each dimension
   */
  std::vector<double> getGridSpacing() const {
    return backend_.getGridSpacing();
  }

  /**
   * @brief Set the resolution along each dimension
   *
   * @param resolution Vector containing new resolution for each dimension
   */
  void setResolution(const std::vector<size_t>& resolution) {
    backend_.setResolution(resolution);
  }

  /**
   * @brief Get coordinates at a specific grid index
   *
   * @param indices Grid indices
   * @return Vector of physical coordinates
   */
  std::vector<double> getCoordinates(const std::vector<size_t>& indices) const {
    return backend_.getCoordinates(indices);
  }

  /**
   * @brief Get grid indices for a physical coordinate
   *
   * @param coordinates Physical coordinates
   * @return Vector of grid indices
   */
  std::vector<size_t> getIndices(const std::vector<double>& coordinates) const {
    return backend_.getIndices(coordinates);
  }

  /**
   * @brief Get iterator to the beginning of the geometry grid
   *
   * @return Iterator pointing to the first grid point
   */
  std::unique_ptr<IGeometryIterator<double>> begin() const {
    if (getDimensions() == 0 || getTotalPoints() == 0) {
      return end();
    }

    std::vector<size_t> position(getDimensions(), 0);
    std::vector<double> coordinates = getCoordinates(position);

    return std::make_unique<GeometryPointIterator<double>>(
        position, getResolution(), coordinates, this);
  }

  /**
   * @brief Get iterator to the end of the geometry grid
   *
   * @return Iterator pointing past the last grid point
   */
  std::unique_ptr<IGeometryIterator<double>> end() const {
    return std::make_unique<GeometryIterator<double>>();
  }

  /**
   * @brief Create a new geometry with the same configuration
   *
   * @return Unique pointer to a new geometry instance
   */
  std::unique_ptr<IGeometry> clone() const {
    return std::make_unique<Geometry<Backend>>(backend_);
  }

  /**
   * @brief Check if this geometry is compatible with another geometry adapter
   *
   * @param other Geometry adapter to compare with
   * @return True if geometries are compatible, false otherwise
   */
  template <typename OtherBackend>
  bool isCompatible(const Geometry<OtherBackend>& other) const {
    return backend_.isCompatible(other.backend());
  }

  /**
   * @brief Check if this geometry is compatible with another geometry
   *
   * @param other Reference to IGeometry to compare with
   * @return True if geometries are compatible, false otherwise
   */
  bool isCompatible(const IGeometry& other) const {
    // Delegate to the backend to handle compatibility check
    return backend_.isCompatible(other);
  }

  /**
   * @brief Access to the underlying backend implementation
   *
   * @return Reference to the backend implementation
   */
  Backend& backend() { return backend_; }

  /**
   * @brief Const access to the underlying backend implementation
   *
   * @return Const reference to the backend implementation
   */
  const Backend& backend() const { return backend_; }

 private:
  Backend backend_{};  ///< Backend implementation
};

}  // namespace metada::framework