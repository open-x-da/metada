/**
 * @file MockGeometry.hpp
 * @brief Mock implementation of IGeometry interface for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the IGeometry interface using
 * Google Mock. It allows testing code that depends on IGeometry by providing
 * mock implementations of all interface methods that can be configured with
 * expectations and behaviors.
 *
 * The mock implementation supports:
 * - Setting expectations on method calls
 * - Configuring return values and behaviors
 * - Verifying interaction patterns
 * - Testing error conditions
 *
 * @see IGeometry
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include "IGeometry.hpp"
#include "IGeometryIterator.hpp"
#include "utils/config/Config.hpp"

namespace metada::tests {

using framework::Config;
using framework::GeometryIterator;
using framework::IConfig;
using framework::IGeometry;
using framework::IGeometryIterator;

/**
 * @brief Mock implementation of IGeometry for testing
 *
 * @details
 * Provides mock methods for all IGeometry interface operations, organized into
 * the following categories:
 *
 * @par Core Operations
 * - initialize() - Initialize geometry from configuration
 * - getDimensions() - Get the number of dimensions
 * - getResolution() - Get resolution along each dimension
 * - containsPoint() - Check if a point is within the domain
 *
 * @par Coordinate System
 * - getDomainBounds() - Get domain bounds for a specific dimension
 * - setDomainBounds() - Set domain bounds for a specific dimension
 * - getGridSpacing() - Get grid spacing for each dimension
 * - setResolution() - Set the resolution along each dimension
 *
 * @par Coordinate Transformations
 * - getCoordinates() - Get coordinates at a specific grid index
 * - getIndices() - Get grid indices for a physical coordinate
 *
 * @par Grid Iteration
 * - begin() - Get iterator to the beginning of the geometry grid
 * - end() - Get iterator to the end of the geometry grid
 *
 * @par Other
 * - clone() - Create a new geometry with the same configuration
 * - isCompatible() - Check if this geometry is compatible with another
 */

/**
 * @brief Mock implementation of IGeometry using Google Mock
 */
class MockGeometry : public IGeometry {
 private:
  // Store dimensions and resolution for iterator returns
  size_t dimensions_ = 2;
  std::vector<size_t> resolution_ = {10, 10};
  std::vector<double> min_bounds_ = {0.0, 0.0};
  std::vector<double> max_bounds_ = {1.0, 1.0};
  std::vector<double> grid_spacing_ = {0.1, 0.1};
  size_t total_points_ = 100;

  // Iterator state for tests
  mutable std::vector<size_t> current_position_ = {0, 0};
  mutable std::vector<double> current_coordinates_ = {0.0, 0.0};

 public:
  /**
   * @brief Default constructor
   */
  MockGeometry() = default;

  /**
   * @brief Constructor with config
   * @param config Configuration object
   */
  explicit MockGeometry(const IConfig& config) { initialize(config); }

  // Core methods
  MOCK_METHOD(void, initialize, (const IConfig& config), (override));
  MOCK_METHOD(size_t, getDimensions, (), (const, override));
  MOCK_METHOD(std::vector<size_t>, getResolution, (), (const, override));
  MOCK_METHOD(size_t, getTotalPoints, (), (const, override));
  MOCK_METHOD(bool, containsPoint, (const std::vector<double>& coordinates),
              (const, override));

  // Domain methods
  MOCK_METHOD((std::array<double, 2>), getDomainBounds, (size_t dimension),
              (const, override));
  MOCK_METHOD(void, setDomainBounds,
              (size_t dimension, double min_value, double max_value),
              (override));
  MOCK_METHOD(std::vector<double>, getGridSpacing, (), (const, override));
  MOCK_METHOD(void, setResolution, (const std::vector<size_t>& resolution),
              (override));

  // Coordinate transformation methods
  MOCK_METHOD(std::vector<double>, getCoordinates,
              (const std::vector<size_t>& indices), (const, override));
  MOCK_METHOD(std::vector<size_t>, getIndices,
              (const std::vector<double>& coordinates), (const, override));

  /**
   * @brief Get iterator to beginning of grid
   * @return Iterator pointing to first grid point
   */
  std::unique_ptr<IGeometryIterator<double>> begin() const override {
    current_position_ = std::vector<size_t>(dimensions_, 0);
    current_coordinates_ = std::vector<double>(dimensions_, 0.0);
    return std::make_unique<GeometryIterator<double>>(
        current_position_, resolution_, current_coordinates_);
  }

  /**
   * @brief Get iterator to end of grid
   * @return Iterator pointing past last grid point
   */
  std::unique_ptr<IGeometryIterator<double>> end() const override {
    return std::make_unique<GeometryIterator<double>>();
  }

  MOCK_METHOD(std::unique_ptr<IGeometry>, clone, (), (const, override));
  MOCK_METHOD(bool, isCompatible, (const IGeometry& other), (const, override));

  /**
   * @brief Set dimension count for test
   * @param dimensions Number of dimensions
   */
  void setTestDimensions(size_t dimensions) {
    dimensions_ = dimensions;
    ON_CALL(*this, getDimensions()).WillByDefault(testing::Return(dimensions));

    // Make sure iterators have the right size
    current_position_.resize(dimensions_, 0);
    current_coordinates_.resize(dimensions_, 0.0);
  }

  /**
   * @brief Set grid resolution for test
   * @param resolution Vector of grid points per dimension
   */
  void setTestResolution(const std::vector<size_t>& resolution) {
    resolution_ = resolution;
    ON_CALL(*this, getResolution()).WillByDefault(testing::Return(resolution));

    // Set the total points based on dimensions and resolution
    total_points_ = 1;
    for (auto r : resolution) {
      total_points_ *= r;
    }
    ON_CALL(*this, getTotalPoints())
        .WillByDefault(testing::Return(total_points_));
  }

  /**
   * @brief Set domain bounds for test
   * @param min_bounds Minimum bounds for each dimension
   * @param max_bounds Maximum bounds for each dimension
   */
  void setTestBounds(const std::vector<double>& min_bounds,
                     const std::vector<double>& max_bounds) {
    min_bounds_ = min_bounds;
    max_bounds_ = max_bounds;

    // Calculate grid spacing
    grid_spacing_.resize(min_bounds.size());
    for (size_t i = 0; i < min_bounds.size(); ++i) {
      if (resolution_[i] <= 1) {
        grid_spacing_[i] = max_bounds[i] - min_bounds[i];
      } else {
        grid_spacing_[i] =
            (max_bounds[i] - min_bounds[i]) / (resolution_[i] - 1);
      }
    }

    ON_CALL(*this, getGridSpacing())
        .WillByDefault(testing::Return(grid_spacing_));
  }

  /**
   * @brief Set grid spacing for test
   * @param spacing Grid spacing for each dimension
   */
  void setTestGridSpacing(const std::vector<double>& spacing) {
    grid_spacing_ = spacing;
    ON_CALL(*this, getGridSpacing()).WillByDefault(testing::Return(spacing));
  }

  /**
   * @brief Set domain bounds for specific dimension
   * @param dimension Dimension index
   * @param min_val Minimum bound
   * @param max_val Maximum bound
   */
  void setTestDomainBounds(size_t dimension, double min_val, double max_val) {
    if (dimension < min_bounds_.size() && dimension < max_bounds_.size()) {
      min_bounds_[dimension] = min_val;
      max_bounds_[dimension] = max_val;
    }

    std::array<double, 2> bounds = {min_val, max_val};
    ON_CALL(*this, getDomainBounds(dimension))
        .WillByDefault(testing::Return(bounds));

    // Update grid spacing for this dimension
    if (dimension < grid_spacing_.size()) {
      if (resolution_[dimension] <= 1) {
        grid_spacing_[dimension] = max_val - min_val;
      } else {
        grid_spacing_[dimension] =
            (max_val - min_val) / (resolution_[dimension] - 1);
      }
    }
  }

  /**
   * @brief Helper to create a default implementation for clone
   * @return A clone of this geometry
   */
  std::unique_ptr<IGeometry> defaultClone() const {
    auto clone = std::make_unique<MockGeometry>();
    clone->setTestDimensions(dimensions_);
    clone->setTestResolution(resolution_);
    clone->setTestBounds(min_bounds_, max_bounds_);
    return clone;
  }

  /**
   * @brief Helper to implement default compatibility checking
   * @param other The geometry to compare with
   * @return True if geometries are compatible
   */
  bool defaultIsCompatible(const IGeometry& other) const {
    // Basic compatibility checks
    if (getDimensions() != other.getDimensions()) {
      return false;
    }

    if (getResolution() != other.getResolution()) {
      return false;
    }

    return true;
  }
};

}  // namespace metada::tests