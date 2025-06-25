#pragma once

#include <cmath>
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
 * @brief Class representing a spatial location in various coordinate systems
 *
 * This class provides a unified way to represent locations across the
 * framework, supporting grid coordinates, geographic coordinates, and Cartesian
 * coordinates. It uses std::variant to store different coordinate types
 * efficiently.
 */
class Location {
 public:
  CoordinateSystem system;
  std::variant<std::tuple<int, int, int>,  // Grid coordinates (i, j, k)
               std::tuple<double, double, double>,  // Geographic/Cartesian
               std::pair<int, int>  // 2D grid coordinates (i, j)
               >
      coordinates;

  // Constructors
  Location(int i, int j, int k)
      : system(CoordinateSystem::GRID), coordinates(std::make_tuple(i, j, k)) {}
  Location(int i, int j)
      : system(CoordinateSystem::GRID), coordinates(std::make_pair(i, j)) {}
  Location(double x, double y, double z,
           CoordinateSystem coord_system = CoordinateSystem::GEOGRAPHIC)
      : system(coord_system), coordinates(std::make_tuple(x, y, z)) {}

  // Coordinate accessors
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
  CoordinateSystem getCoordinateSystem() const { return system; }
  bool operator==(const Location& other) const {
    return system == other.system && coordinates == other.coordinates;
  }
  bool operator!=(const Location& other) const { return !(*this == other); }

  /**
   * @brief Compute the distance to another Location
   *
   * For grid coordinates, computes Euclidean distance in (i,j).
   * For geographic coordinates, computes Euclidean distance in (lat,lon).
   * (Extend to Haversine for more accuracy if needed.)
   * Throws if coordinate systems are incompatible.
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

  // Add this method for 2D grid index calculation
  int getGridIndex(int x_dim) const {
    auto [x, y] = getGridCoords2D();
    return y * x_dim + x;
  }

 private:
  static constexpr double kEarthRadiusKm = 6371.0;
  static constexpr double pi = 3.14159265358979323846;
  static double deg2rad(double deg) { return deg * pi / 180.0; }
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