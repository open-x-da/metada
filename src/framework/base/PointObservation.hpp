/**
 * @file PointObservation.hpp
 * @brief Generic data structures for point-based observations
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains generic data structures for representing point-based
 * observations, including their location, value, and error. These structures
 * are intended to be reused across different observation-related backend
 * implementations to promote consistency and code reuse.
 */
#pragma once

#include <stdexcept>
#include <tuple>
#include <utility>
#include <variant>

namespace metada::framework {

/**
 * @brief Enumeration for different coordinate systems
 */
enum class CoordinateSystem {
  GRID,        ///< Integer grid coordinates (i, j, k)
  GEOGRAPHIC,  ///< Geographic coordinates (lat, lon, level)
  CARTESIAN    ///< Cartesian coordinates (x, y, z)
};

/**
 * @brief Generic location class that can handle different coordinate systems
 *
 * This class provides a unified way to represent locations across the
 * framework, supporting grid coordinates, geographic coordinates, and Cartesian
 * coordinates. It uses std::variant to store different coordinate types
 * efficiently.
 */
struct Location {
  CoordinateSystem system;
  std::variant<std::tuple<int, int, int>,  // Grid coordinates (i, j, k)
               std::tuple<double, double, double>,  // Geographic/Cartesian
                                                    // (lat, lon, level) or (x,
                                                    // y, z)
               std::pair<int, int>  // 2D grid coordinates (i, j)
               >
      coordinates;

  /**
   * @brief Constructor for 3D grid coordinates
   * @param i X coordinate (grid index)
   * @param j Y coordinate (grid index)
   * @param k Z coordinate (grid index)
   */
  Location(int i, int j, int k)
      : system(CoordinateSystem::GRID), coordinates(std::make_tuple(i, j, k)) {}

  /**
   * @brief Constructor for 2D grid coordinates
   * @param i X coordinate (grid index)
   * @param j Y coordinate (grid index)
   */
  Location(int i, int j)
      : system(CoordinateSystem::GRID), coordinates(std::make_pair(i, j)) {}

  /**
   * @brief Constructor for geographic or Cartesian coordinates
   * @param x First coordinate (latitude for geographic, x for Cartesian)
   * @param y Second coordinate (longitude for geographic, y for Cartesian)
   * @param z Third coordinate (level for geographic, z for Cartesian)
   * @param coord_system Coordinate system (defaults to GEOGRAPHIC)
   */
  Location(double x, double y, double z,
           CoordinateSystem coord_system = CoordinateSystem::GEOGRAPHIC)
      : system(coord_system), coordinates(std::make_tuple(x, y, z)) {}

  /**
   * @brief Get 3D grid coordinates
   * @return Tuple of (i, j, k) grid coordinates
   * @throws std::runtime_error if coordinates are not in grid format
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
   * @return Pair of (i, j) grid coordinates
   * @throws std::runtime_error if coordinates are not in grid format
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
   * @return Tuple of (latitude, longitude, level)
   * @throws std::runtime_error if coordinates are not in geographic format
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
   * @return Tuple of (x, y, z)
   * @throws std::runtime_error if coordinates are not in Cartesian format
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
   * @brief Get coordinate system
   * @return The coordinate system of this location
   */
  CoordinateSystem getCoordinateSystem() const { return system; }

  /**
   * @brief Equality comparison
   */
  bool operator==(const Location& other) const {
    return system == other.system && coordinates == other.coordinates;
  }

  /**
   * @brief Inequality comparison
   */
  bool operator!=(const Location& other) const { return !(*this == other); }
};

/**
 * @brief Structure to hold a single observation point
 */
struct ObservationPoint {
  Location location;  ///< Location of the observation
  double value;       ///< Observed value
  double error;       ///< Observation error
  bool is_valid;      ///< Whether the observation is valid

  /**
   * @brief Constructor with Location object
   * @param loc Location of the observation
   * @param val Observed value
   * @param err Observation error
   */
  ObservationPoint(const Location& loc, double val, double err)
      : location(loc), value(val), error(err), is_valid(true) {}

  /**
   * @brief Constructor with Location object (invalid observation)
   * @param loc Location of the observation
   */
  ObservationPoint(const Location& loc)
      : location(loc), value(0.0), error(0.0), is_valid(false) {}

  bool operator==(const ObservationPoint& other) const {
    return location == other.location && value == other.value &&
           error == other.error && is_valid == other.is_valid;
  }
};

}  // namespace metada::framework