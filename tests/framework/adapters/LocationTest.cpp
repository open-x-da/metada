/**
 * @file LocationTest.cpp
 * @brief Unit tests for the new Location class
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This test suite verifies the functionality of the new Location class,
 * which provides a unified way to represent locations across the framework.
 * The tests cover:
 *
 * - Construction with different coordinate systems
 * - Coordinate extraction methods
 * - Type safety and error handling
 * - Conversion between coordinate systems
 * - Equality comparisons
 */

#include <gtest/gtest.h>

#include "PointObservation.hpp"

namespace metada::tests {

using namespace metada::framework;

class LocationTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(LocationTest, GridCoordinates3D) {
  // Test 3D grid coordinates
  Location loc(10, 20, 5);

  EXPECT_EQ(loc.getCoordinateSystem(), CoordinateSystem::GRID);

  auto [i, j, k] = loc.getGridCoords();
  EXPECT_EQ(i, 10);
  EXPECT_EQ(j, 20);
  EXPECT_EQ(k, 5);

  auto [i2, j2] = loc.getGridCoords2D();
  EXPECT_EQ(i2, 10);
  EXPECT_EQ(j2, 20);
}

TEST_F(LocationTest, GridCoordinates2D) {
  // Test 2D grid coordinates
  Location loc(15, 25);

  EXPECT_EQ(loc.getCoordinateSystem(), CoordinateSystem::GRID);

  auto [i, j] = loc.getGridCoords2D();
  EXPECT_EQ(i, 15);
  EXPECT_EQ(j, 25);

  auto [i2, j2, k2] = loc.getGridCoords();
  EXPECT_EQ(i2, 15);
  EXPECT_EQ(j2, 25);
  EXPECT_EQ(k2, 0);  // Default k value
}

TEST_F(LocationTest, GeographicCoordinates) {
  // Test geographic coordinates
  Location loc(45.5, -120.3, 1000.0);

  EXPECT_EQ(loc.getCoordinateSystem(), CoordinateSystem::GEOGRAPHIC);

  auto [lat, lon, level] = loc.getGeographicCoords();
  EXPECT_DOUBLE_EQ(lat, 45.5);
  EXPECT_DOUBLE_EQ(lon, -120.3);
  EXPECT_DOUBLE_EQ(level, 1000.0);
}

TEST_F(LocationTest, CartesianCoordinates) {
  // Test Cartesian coordinates
  Location loc(100.0, 200.0, 300.0, CoordinateSystem::CARTESIAN);

  EXPECT_EQ(loc.getCoordinateSystem(), CoordinateSystem::CARTESIAN);

  auto [x, y, z] = loc.getCartesianCoords();
  EXPECT_DOUBLE_EQ(x, 100.0);
  EXPECT_DOUBLE_EQ(y, 200.0);
  EXPECT_DOUBLE_EQ(z, 300.0);
}

TEST_F(LocationTest, EqualityComparison) {
  Location loc1(10, 20, 5);
  Location loc2(10, 20, 5);
  Location loc3(15, 20, 5);

  EXPECT_EQ(loc1, loc2);
  EXPECT_NE(loc1, loc3);

  Location geo1(45.0, -120.0, 1000.0);
  Location geo2(45.0, -120.0, 1000.0);
  Location geo3(46.0, -120.0, 1000.0);

  EXPECT_EQ(geo1, geo2);
  EXPECT_NE(geo1, geo3);
}

TEST_F(LocationTest, ErrorHandling) {
  // Test error handling for wrong coordinate system access
  Location grid_loc(10, 20);
  Location geo_loc(45.0, -120.0, 1000.0);

  // Should not throw for correct access
  EXPECT_NO_THROW(grid_loc.getGridCoords());
  EXPECT_NO_THROW(geo_loc.getGeographicCoords());

  // Should throw for wrong coordinate system access
  EXPECT_THROW(grid_loc.getGeographicCoords(), std::runtime_error);
  EXPECT_THROW(geo_loc.getGridCoords(), std::runtime_error);
  EXPECT_THROW(geo_loc.getCartesianCoords(), std::runtime_error);
}

TEST_F(LocationTest, ObservationPointCompatibility) {
  // Test that ObservationPoint works with Location
  Location loc(45.0, -120.0, 1000.0);

  // Create ObservationPoint with Location
  ObservationPoint point1(loc, 25.5, 0.5);
  EXPECT_EQ(point1.location, loc);
  EXPECT_DOUBLE_EQ(point1.value, 25.5);
  EXPECT_DOUBLE_EQ(point1.error, 0.5);
  EXPECT_TRUE(point1.is_valid);

  // Test that location coordinates can be accessed
  auto [lat, lon, level] = point1.location.getGeographicCoords();
  EXPECT_DOUBLE_EQ(lat, 45.0);
  EXPECT_DOUBLE_EQ(lon, -120.0);
  EXPECT_DOUBLE_EQ(level, 1000.0);
}

TEST_F(LocationTest, MixedCoordinateSystems) {
  // Test that different coordinate systems are handled correctly
  Location grid_loc(10, 20);
  Location geo_loc(45.0, -120.0, 1000.0);
  Location cart_loc(100.0, 200.0, 300.0, CoordinateSystem::CARTESIAN);

  // Each should have the correct coordinate system
  EXPECT_EQ(grid_loc.getCoordinateSystem(), CoordinateSystem::GRID);
  EXPECT_EQ(geo_loc.getCoordinateSystem(), CoordinateSystem::GEOGRAPHIC);
  EXPECT_EQ(cart_loc.getCoordinateSystem(), CoordinateSystem::CARTESIAN);

  // They should not be equal even with similar values
  EXPECT_NE(grid_loc, geo_loc);
  EXPECT_NE(geo_loc, cart_loc);
  EXPECT_NE(grid_loc, cart_loc);
}

}  // namespace metada::tests