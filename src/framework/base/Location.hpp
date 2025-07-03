#pragma once

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <variant>

namespace metada::framework {

/**
 * @brief Enumeration for different coordinate systems
 *
 * Defines the supported coordinate systems for spatial locations:
 * - GRID: Integer-based grid coordinates for discrete spatial indexing
 * - GEOGRAPHIC: Latitude/longitude/level coordinates for Earth-based
 * positioning
 * - CARTESIAN: 3D Cartesian coordinates for mathematical calculations
 */
enum class CoordinateSystem {
  GRID,        ///< Integer grid coordinates (i, j, k)
  GEOGRAPHIC,  ///< Geographic coordinates (lat, lon, level)
  CARTESIAN    ///< Cartesian coordinates (x, y, z)
};

/**
 * @brief Unified spatial location representation supporting multiple coordinate
 * systems
 *
 * The Location class provides a type-safe, unified interface for representing
 * spatial positions across different coordinate systems. It uses std::variant
 * for efficient storage and compile-time type safety, supporting:
 *
 * - 2D/3D grid coordinates (integer-based)
 * - Geographic coordinates (latitude, longitude, level)
 * - Cartesian coordinates (x, y, z)
 *
 * @note All coordinate transformations must be handled externally. This class
 *       only stores and provides access to coordinates in their native systems.
 *
 * @example
 * @code
 * // Create grid-based location
 * Location gridLoc(10, 20, 5);
 *
 * // Create geographic location
 * Location geoLoc(45.0, -122.0, 100.0, CoordinateSystem::GEOGRAPHIC);
 *
 * // Calculate distance between compatible locations
 * double dist = gridLoc1.distance_to(gridLoc2);
 * @endcode
 */
class Location {
 public:
  // ============================================================================
  // CONSTRUCTORS
  // ============================================================================

  /**
   * @brief Construct a 3D grid location
   * @param i Grid index in the i-direction
   * @param j Grid index in the j-direction
   * @param k Grid index in the k-direction (vertical/depth)
   */
  Location(int i, int j, int k)
      : system(CoordinateSystem::GRID), coordinates(std::make_tuple(i, j, k)) {}

  /**
   * @brief Construct a 2D grid location
   * @param i Grid index in the i-direction
   * @param j Grid index in the j-direction
   */
  Location(int i, int j)
      : system(CoordinateSystem::GRID), coordinates(std::make_pair(i, j)) {}

  /**
   * @brief Construct a geographic or Cartesian location
   * @param x First coordinate (latitude for geographic, x for Cartesian)
   * @param y Second coordinate (longitude for geographic, y for Cartesian)
   * @param z Third coordinate (level/altitude for geographic, z for Cartesian)
   * @param coord_system Coordinate system type (defaults to GEOGRAPHIC)
   */
  Location(double x, double y, double z,
           CoordinateSystem coord_system = CoordinateSystem::GEOGRAPHIC)
      : system(coord_system), coordinates(std::make_tuple(x, y, z)) {}

  // ============================================================================
  // COORDINATE ACCESSORS
  // ============================================================================

  /**
   * @brief Get 3D grid coordinates
   * @return Tuple containing (i, j, k) grid indices
   * @throws std::runtime_error If location is not in grid coordinate system
   * @throws std::runtime_error If coordinate storage is invalid for grid system
   *
   * @note For 2D grid locations, k coordinate is returned as 0
   */
  std::tuple<int, int, int> getGridCoords() const {
    if (system != CoordinateSystem::GRID) {
      throw std::runtime_error("Location is not in grid coordinate system");
    }
    if (std::holds_alternative<std::tuple<int, int, int>>(coordinates)) {
      return std::get<std::tuple<int, int, int>>(coordinates);
    } else if (std::holds_alternative<std::pair<int, int>>(coordinates)) {
      auto [i, j] = std::get<std::pair<int, int>>(coordinates);
      return std::make_tuple(i, j, 0);
    }
    throw std::runtime_error("Invalid coordinate storage for grid system");
  }

  /**
   * @brief Get 2D grid coordinates
   * @return Pair containing (i, j) grid indices
   * @throws std::runtime_error If location is not in grid coordinate system
   * @throws std::runtime_error If coordinate storage is invalid for grid system
   *
   * @note For 3D grid locations, only i and j coordinates are returned
   */
  std::pair<int, int> getGridCoords2D() const {
    if (system != CoordinateSystem::GRID) {
      throw std::runtime_error("Location is not in grid coordinate system");
    }
    if (std::holds_alternative<std::pair<int, int>>(coordinates)) {
      return std::get<std::pair<int, int>>(coordinates);
    } else if (std::holds_alternative<std::tuple<int, int, int>>(coordinates)) {
      auto [i, j, k] = std::get<std::tuple<int, int, int>>(coordinates);
      return std::make_pair(i, j);
    }
    throw std::runtime_error("Invalid coordinate storage for grid system");
  }

  /**
   * @brief Get geographic coordinates
   * @return Tuple containing (latitude, longitude, level) in degrees and meters
   * @throws std::runtime_error If location is not in geographic coordinate
   * system
   * @throws std::runtime_error If coordinate storage is invalid for geographic
   * system
   */
  std::tuple<double, double, double> getGeographicCoords() const {
    if (system != CoordinateSystem::GEOGRAPHIC) {
      throw std::runtime_error(
          "Location is not in geographic coordinate system");
    }
    if (!std::holds_alternative<std::tuple<double, double, double>>(
            coordinates)) {
      throw std::runtime_error(
          "Invalid coordinate storage for geographic system");
    }
    return std::get<std::tuple<double, double, double>>(coordinates);
  }

  /**
   * @brief Get Cartesian coordinates
   * @return Tuple containing (x, y, z) coordinates
   * @throws std::runtime_error If location is not in Cartesian coordinate
   * system
   * @throws std::runtime_error If coordinate storage is invalid for Cartesian
   * system
   */
  std::tuple<double, double, double> getCartesianCoords() const {
    if (system != CoordinateSystem::CARTESIAN) {
      throw std::runtime_error(
          "Location is not in Cartesian coordinate system");
    }
    if (!std::holds_alternative<std::tuple<double, double, double>>(
            coordinates)) {
      throw std::runtime_error(
          "Invalid coordinate storage for Cartesian system");
    }
    return std::get<std::tuple<double, double, double>>(coordinates);
  }

  /**
   * @brief Get the active coordinate system
   * @return The coordinate system type for this location
   */
  CoordinateSystem getCoordinateSystem() const { return system; }

  // ============================================================================
  // UTILITY METHODS
  // ============================================================================

  /**
   * @brief Compute the distance to another Location
   *
   * Calculates distance based on the coordinate system:
   * - Grid coordinates: Euclidean distance in (i,j) space
   * - Geographic coordinates: Great circle distance using Haversine formula
   * - Cartesian coordinates: 3D Euclidean distance
   *
   * @param other The target location to measure distance to
   * @return Distance value (units depend on coordinate system)
   * @throws std::runtime_error If coordinate systems are incompatible
   *
   * @note Grid distances are dimensionless. Geographic distances are in
   * kilometers. Cartesian distances use the same units as the input
   * coordinates.
   */
  double distance_to(const Location& other) const {
    if (system == CoordinateSystem::GRID &&
        other.system == CoordinateSystem::GRID) {
      auto [i1, j1] = getGridCoords2D();
      auto [i2, j2] = other.getGridCoords2D();
      double dx = static_cast<double>(i1 - i2);
      double dy = static_cast<double>(j1 - j2);
      return std::sqrt(dx * dx + dy * dy);
    } else if (system == CoordinateSystem::GEOGRAPHIC &&
               other.system == CoordinateSystem::GEOGRAPHIC) {
      auto [lat1, lon1, level1] = getGeographicCoords();
      auto [lat2, lon2, level2] = other.getGeographicCoords();
      return haversine(lat1, lon1, lat2, lon2);
    } else if (system == CoordinateSystem::CARTESIAN &&
               other.system == CoordinateSystem::CARTESIAN) {
      auto [x1, y1, z1] = getCartesianCoords();
      auto [x2, y2, z2] = other.getCartesianCoords();
      double dx = x1 - x2;
      double dy = y1 - y2;
      double dz = z1 - z2;
      return std::sqrt(dx * dx + dy * dy + dz * dz);
    } else {
      throw std::runtime_error(
          "Cannot compute distance between incompatible Location coordinate "
          "systems");
    }
  }

  // ============================================================================
  // COMPARISON OPERATORS
  // ============================================================================

  /**
   * @brief Equality comparison operator
   * @param other Location to compare with
   * @return true if both coordinate system and coordinates are identical
   */
  bool operator==(const Location& other) const {
    return system == other.system && coordinates == other.coordinates;
  }

  /**
   * @brief Inequality comparison operator
   * @param other Location to compare with
   * @return true if coordinate system or coordinates differ
   */
  bool operator!=(const Location& other) const { return !(*this == other); }

  // ============================================================================
  // STREAM OUTPUT
  // ============================================================================

  /**
   * @brief Stream output operator for Location class
   *
   * Provides human-readable string representation of the location:
   * - Grid: "Grid(i,j,k)" or "Grid(i,j,0)" for 2D
   * - Geographic: "Geo(lat,lon,level)"
   * - Cartesian: "Cart(x,y,z)"
   *
   * @param os Output stream
   * @param loc Location to output
   * @return Reference to the output stream for chaining
   */
  friend std::ostream& operator<<(std::ostream& os, const Location& loc) {
    switch (loc.system) {
      case CoordinateSystem::GRID: {
        if (std::holds_alternative<std::tuple<int, int, int>>(
                loc.coordinates)) {
          auto [i, j, k] = std::get<std::tuple<int, int, int>>(loc.coordinates);
          os << "Grid(" << i << "," << j << "," << k << ")";
        } else if (std::holds_alternative<std::pair<int, int>>(
                       loc.coordinates)) {
          auto [i, j] = std::get<std::pair<int, int>>(loc.coordinates);
          os << "Grid(" << i << "," << j << ",0)";
        }
        break;
      }
      case CoordinateSystem::GEOGRAPHIC: {
        auto [lat, lon, level] =
            std::get<std::tuple<double, double, double>>(loc.coordinates);
        os << "Geo(" << lat << "," << lon << "," << level << ")";
        break;
      }
      case CoordinateSystem::CARTESIAN: {
        auto [x, y, z] =
            std::get<std::tuple<double, double, double>>(loc.coordinates);
        os << "Cart(" << x << "," << y << "," << z << ")";
        break;
      }
    }
    return os;
  }

 private:
  // ============================================================================
  // MEMBER VARIABLES
  // ============================================================================

  CoordinateSystem system;  ///< Active coordinate system for this location

  /**
   * @brief Storage for coordinate values using type-safe variant
   *
   * Supports three coordinate storage types:
   * - std::tuple<int, int, int>: 3D grid coordinates (i, j, k)
   * - std::tuple<double, double, double>: Geographic (lat, lon, level) or
   * Cartesian (x, y, z)
   * - std::pair<int, int>: 2D grid coordinates (i, j)
   */
  std::variant<std::tuple<int, int, int>,  // 3D grid coordinates (i, j, k)
               std::tuple<double, double, double>,  // Geographic/Cartesian (x,
                                                    // y, z)
               std::pair<int, int>  // 2D grid coordinates (i, j)
               >
      coordinates;

  // ============================================================================
  // CONSTANTS AND HELPER METHODS
  // ============================================================================

  static constexpr double kEarthRadiusKm =
      6371.0;  ///< Earth radius in kilometers
  static constexpr double pi =
      3.14159265358979323846;  ///< Mathematical constant π

  /**
   * @brief Convert degrees to radians
   * @param deg Angle in degrees
   * @return Angle in radians
   */
  static double deg2rad(double deg) { return deg * pi / 180.0; }

  /**
   * @brief Calculate great circle distance using Haversine formula
   *
   * Computes the shortest distance between two points on the surface of a
   * sphere (Earth) given their latitude and longitude coordinates.
   *
   * @param lat1 Latitude of first point in degrees
   * @param lon1 Longitude of first point in degrees
   * @param lat2 Latitude of second point in degrees
   * @param lon2 Longitude of second point in degrees
   * @return Distance in kilometers
   */
  static double haversine(double lat1, double lon1, double lat2, double lon2) {
    double dlat = deg2rad(lat2 - lat1);
    double dlon = deg2rad(lon2 - lon1);
    double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
               std::cos(deg2rad(lat1)) * std::cos(deg2rad(lat2)) *
                   std::sin(dlon / 2) * std::sin(dlon / 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
    return kEarthRadiusKm * c;
  }
};

}  // namespace metada::framework