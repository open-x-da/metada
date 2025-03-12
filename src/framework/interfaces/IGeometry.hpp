/**
 * @file IGeometry.hpp
 * @brief Interface defining the contract for geometry implementations
 * @ingroup interfaces
 *
 * @details
 * This header provides the abstract interface that all geometry implementations
 * must follow. It defines a unified API for handling grid/geometry/resolution
 * settings of models and states.
 *
 * The interface includes:
 * - Grid dimension and resolution management
 * - Coordinate system definition
 * - Iterator for traversing grid points
 * - Boundary conditions specification
 * - Domain decomposition support for parallel computing
 */

#pragma once

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace metada::framework {

class IConfig;  // Forward declaration

template <typename T>
class IGeometryIterator;

/**
 * @brief Abstract interface for geometry implementations
 *
 * @details
 * This interface defines the contract that all geometry implementations must
 * follow. It provides a unified API for handling grid/geometry/resolution
 * settings of models and states.
 */
class IGeometry {
 protected:
  /**
   * @brief Protected constructor for abstract class
   */
  IGeometry() = default;

 public:
  /**
   * @brief Virtual destructor
   */
  virtual ~IGeometry() = default;

  /**
   * @brief Copy constructor
   */
  IGeometry(const IGeometry&) = default;

  /**
   * @brief Copy assignment operator
   */
  IGeometry& operator=(const IGeometry&) = default;

  /**
   * @brief Move constructor
   */
  IGeometry(IGeometry&&) = default;

  /**
   * @brief Move assignment operator
   */
  IGeometry& operator=(IGeometry&&) = default;

  /**
   * @brief Initialize geometry from configuration
   *
   * @param config Configuration object containing geometry parameters
   */
  virtual void initialize(const IConfig& config) = 0;

  /**
   * @brief Get the number of dimensions
   *
   * @return Number of spatial dimensions
   */
  virtual size_t getDimensions() const = 0;

  /**
   * @brief Get resolution along each dimension
   *
   * @return Vector containing number of grid points along each dimension
   */
  virtual std::vector<size_t> getResolution() const = 0;

  /**
   * @brief Get total number of grid points
   *
   * @return Total grid point count
   */
  virtual size_t getTotalPoints() const = 0;

  /**
   * @brief Check if a point is within the domain
   *
   * @param coordinates Vector of coordinates to check
   * @return True if point is within domain, false otherwise
   */
  virtual bool containsPoint(const std::vector<double>& coordinates) const = 0;

  /**
   * @brief Get domain bounds for a specific dimension
   *
   * @param dimension Dimension index (0-based)
   * @return Array containing min and max values for the dimension
   */
  virtual std::array<double, 2> getDomainBounds(size_t dimension) const = 0;

  /**
   * @brief Set domain bounds for a specific dimension
   *
   * @param dimension Dimension index (0-based)
   * @param min_value Minimum value
   * @param max_value Maximum value
   */
  virtual void setDomainBounds(size_t dimension, double min_value,
                               double max_value) = 0;

  /**
   * @brief Get grid spacing for each dimension
   *
   * @return Vector containing grid spacing for each dimension
   */
  virtual std::vector<double> getGridSpacing() const = 0;

  /**
   * @brief Set the resolution along each dimension
   *
   * @param resolution Vector containing new resolution for each dimension
   */
  virtual void setResolution(const std::vector<size_t>& resolution) = 0;

  /**
   * @brief Get coordinates at a specific grid index
   *
   * @param indices Grid indices
   * @return Vector of physical coordinates
   */
  virtual std::vector<double> getCoordinates(
      const std::vector<size_t>& indices) const = 0;

  /**
   * @brief Get grid indices for a physical coordinate
   *
   * @param coordinates Physical coordinates
   * @return Vector of grid indices
   */
  virtual std::vector<size_t> getIndices(
      const std::vector<double>& coordinates) const = 0;

  /**
   * @brief Get iterator to the beginning of the geometry grid
   *
   * @return Iterator pointing to the first grid point
   */
  virtual std::unique_ptr<IGeometryIterator<double>> begin() const = 0;

  /**
   * @brief Get iterator to the end of the geometry grid
   *
   * @return Iterator pointing past the last grid point
   */
  virtual std::unique_ptr<IGeometryIterator<double>> end() const = 0;

  /**
   * @brief Create a new geometry with the same configuration
   *
   * @return Unique pointer to a new geometry instance
   */
  virtual std::unique_ptr<IGeometry> clone() const = 0;

  /**
   * @brief Check if this geometry is compatible with another geometry
   *
   * @param other Geometry to compare with
   * @return True if geometries are compatible, false otherwise
   */
  virtual bool isCompatible(const IGeometry& other) const = 0;
};

}  // namespace metada::framework