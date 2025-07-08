#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "GridObservationIterator.hpp"
#include "PointObservation.hpp"

namespace metada::backends::common::observation {

using framework::CoordinateSystem;
using framework::Location;
using framework::ObservationPoint;

/**
 * @brief Grid-based observation backend implementation
 *
 * @details
 * This class implements a grid-based observation backend that:
 * - Stores observations as point observations with location information
 * - Each observation has latitude, longitude, and vertical level
 * - Supports iteration over observations
 * - Handles missing values and quality control
 * - Provides basic arithmetic operations
 * - Organizes data by observation type and variable
 */
class GridObservation {
 public:
  // Type aliases for concept compliance
  using iterator_type = GridObservationIterator;
  using value_type = ObservationPoint;

  // Delete default constructor
  GridObservation() = delete;

  // Delete copy constructor and assignment
  GridObservation(const GridObservation&) = delete;
  GridObservation& operator=(const GridObservation&) = delete;

  /**
   * @brief Constructor that initializes from configuration
   * @tparam ConfigBackend Type of configuration backend
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  explicit GridObservation(const ConfigBackend& config) {
    // Get observation files from config
    const auto type_configs = config.Get("types").asVectorMap();

    // Load each observation type
    for (const auto& type_map : type_configs) {
      const auto& [type_name, type_config] = *type_map.begin();
      const auto& type_backend = ConfigBackend(type_config.asMap());
      bool type_if_use = type_backend.Get("if_use").asBool();
      if (!type_if_use) {
        continue;
      }
      std::string filename = type_backend.Get("file").asString();

      // Get coordinate system from config (default to "geographic" for backward
      // compatibility)
      std::string coordinate_system = "geographic";
      try {
        coordinate_system = type_backend.Get("coordinate").asString();
      } catch (...) {
        // Use default if not specified
      }

      // Loop over variables for this type
      const auto& variables_configs =
          type_backend.Get("variables").asVectorMap();
      for (const auto& var_map : variables_configs) {
        for (const auto& [var_name, var_config] : var_map) {
          const auto& var_backend = ConfigBackend(var_config.asMap());
          bool var_if_use = var_backend.Get("if_use").asBool();
          if (!var_if_use) {
            continue;
          }
          double error = var_backend.Get("error").asFloat();
          double missing_value = var_backend.Get("missing_value").asFloat();

          // Load observation file with coordinate system
          loadFromFile(filename, error, missing_value, coordinate_system);
        }
      }
    }
  }

  /**
   * @brief Move constructor
   * @param other Observation to move from
   */
  GridObservation(GridObservation&& other) noexcept
      : observations_(std::move(other.observations_)),
        type_variable_map_(std::move(other.type_variable_map_)),
        covariance_(std::move(other.covariance_)) {}

  /**
   * @brief Move assignment operator
   * @param other Observation to move from
   * @return Reference to this observation
   */
  GridObservation& operator=(GridObservation&& other) noexcept {
    if (this != &other) {
      observations_ = std::move(other.observations_);
      type_variable_map_ = std::move(other.type_variable_map_);
      covariance_ = std::move(other.covariance_);
    }
    return *this;
  }

  /**
   * @brief Get iterator to beginning of observations
   * @return Iterator to first observation
   */
  GridObservationIterator begin() const {
    return GridObservationIterator(&observations_, 0);
  }

  /**
   * @brief Get iterator to end of observations
   * @return Iterator past last observation
   */
  GridObservationIterator end() const {
    return GridObservationIterator(&observations_, observations_.size());
  }

  /**
   * @brief Get number of observations
   * @return Total number of observations
   */
  size_t size() const { return observations_.size(); }

  /**
   * @brief Get observation at specific index
   * @param index Index of the observation
   * @return Reference to the observation
   */
  const ObservationPoint& operator[](size_t index) const {
    return observations_[index];
  }

  /**
   * @brief Get raw pointer to observation data (for backward compatibility)
   * @return Pointer to data
   */
  void* getData() {
    if (observations_.empty()) {
      return nullptr;
    }
    return const_cast<void*>(static_cast<const void*>(observations_.data()));
  }

  /**
   * @brief Get const raw pointer to observation data (for backward
   * compatibility)
   * @return Const pointer to data
   */
  const void* getData() const {
    if (observations_.empty()) {
      return nullptr;
    }
    return static_cast<const void*>(observations_.data());
  }

  /**
   * @brief Get template-specific data (for backward compatibility)
   * @tparam T Type of data to return
   * @return Vector of observation values
   */
  template <typename T>
  T getData() const {
    if constexpr (std::is_same_v<T, std::vector<double>>) {
      std::vector<double> values;
      values.reserve(observations_.size());
      for (const auto& obs : observations_) {
        values.push_back(obs.value);
      }
      return values;
    }
    return T{};
  }

  /**
   * @brief Get names of observation types
   * @return Vector of observation type names
   */
  std::vector<std::string> getTypeNames() const {
    std::vector<std::string> types;
    for (const auto& [type_name, _] : type_variable_map_) {
      types.push_back(type_name);
    }
    return types;
  }

  /**
   * @brief Get names of variables for a specific type
   * @param type_name Name of the observation type
   * @return Vector of variable names for the type
   */
  std::vector<std::string> getVariableNames(
      const std::string& type_name) const {
    std::vector<std::string> variables;
    auto it = type_variable_map_.find(type_name);
    if (it != type_variable_map_.end()) {
      for (const auto& [var_name, _] : it->second) {
        variables.push_back(var_name);
      }
    }
    return variables;
  }

  /**
   * @brief Get the observation error variances
   * @return Vector of observation error variances (one per observation)
   */
  std::vector<double> getCovariance() const {
    std::vector<double> errors;
    errors.reserve(observations_.size());
    for (const auto& obs : observations_) {
      if (obs.is_valid) {
        errors.push_back(obs.error * obs.error);  // Convert error to variance
      } else {
        errors.push_back(
            std::numeric_limits<double>::infinity());  // Invalid observations
                                                       // get infinite variance
      }
    }
    return errors;
  }

  /**
   * @brief Clone this observation
   * @return Unique pointer to cloned observation
   */
  std::unique_ptr<GridObservation> clone() const {
    return std::unique_ptr<GridObservation>(new GridObservation(*this, true));
  }

  /**
   * @brief Add another observation to this one
   * @param other Observation to add
   */
  void add(const GridObservation& other) {
    if (observations_.size() != other.observations_.size()) {
      throw std::runtime_error("Cannot add observations of different sizes");
    }

    for (size_t i = 0; i < observations_.size(); ++i) {
      if (observations_[i].is_valid && other.observations_[i].is_valid) {
        observations_[i].value += other.observations_[i].value;
      }
    }
  }

  /**
   * @brief Subtract another observation from this one
   * @param other Observation to subtract
   */
  void subtract(const GridObservation& other) {
    if (observations_.size() != other.observations_.size()) {
      throw std::runtime_error(
          "Cannot subtract observations of different sizes");
    }

    for (size_t i = 0; i < observations_.size(); ++i) {
      if (observations_[i].is_valid && other.observations_[i].is_valid) {
        observations_[i].value -= other.observations_[i].value;
      }
    }
  }

  /**
   * @brief Multiply this observation by a scalar
   * @param scalar Value to multiply by
   */
  void multiply(double scalar) {
    for (auto& obs : observations_) {
      if (obs.is_valid) {
        obs.value *= scalar;
      }
    }
  }

  /**
   * @brief Compare equality with another observation
   * @param other Observation to compare with
   * @return true if equal, false otherwise
   */
  bool equals(const GridObservation& other) const {
    return observations_ == other.observations_ &&
           type_variable_map_ == other.type_variable_map_ &&
           covariance_ == other.covariance_;
  }

  /**
   * @brief Initialize observation
   */
  void initialize() {
    // Initialization is handled in constructor
  }

  /**
   * @brief Apply quality control to observations
   */
  void applyQC() {
    for (auto& obs : observations_) {
      if (obs.is_valid) {
        // Simple QC: mark values outside reasonable range as invalid
        if (obs.value < -100.0 || obs.value > 100.0) {
          obs.is_valid = false;
        }

        // Check for reasonable latitude/longitude
        if (obs.location.getCoordinateSystem() ==
            CoordinateSystem::GEOGRAPHIC) {
          auto [lat, lon, level] = obs.location.getGeographicCoords();
          if (lat < -90.0 || lat > 90.0 || lon < -180.0 || lon > 180.0) {
            obs.is_valid = false;
          }
        }
      }
    }
  }

  /**
   * @brief Load observation data from file
   * @param filename Path to observation file
   * @param error The observation error to assign to all valid points
   * @param missing_value The value that indicates a missing observation
   * @param coordinate_system Coordinate system to use ("grid" or "geographic")
   */
  void loadFromFile(const std::string& filename, double error,
                    double missing_value,
                    const std::string& coordinate_system = "geographic") {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open observation file: " + filename);
    }

    // Clear existing observations
    observations_.clear();

    std::string line;
    bool header_found = false;
    bool data_started = false;

    while (std::getline(file, line)) {
      // Skip empty lines
      if (line.empty()) {
        continue;
      }

      // Look for the header line that contains column names
      if (line.find("Time") != std::string::npos &&
          line.find("XLONG_U") != std::string::npos) {
        header_found = true;
        continue;
      }

      // Skip lines with only dashes (separators)
      if (line.find("---") != std::string::npos &&
          line.find_first_not_of("- \t") == std::string::npos) {
        if (header_found) {
          data_started = true;  // After header and separator, data begins
        }
        continue;
      }

      // Skip lines before data starts
      if (!data_started) {
        continue;
      }

      // Parse data line: Time Z Y X ZNU XLONG_U XLAT_U Value
      std::istringstream iss(line);
      int time, z, y, x;
      double znu, xlong_u, xlat_u, obs_value;

      if (iss >> time >> z >> y >> x >> znu >> xlong_u >> xlat_u >> obs_value) {
        // Skip missing values
        if (obs_value != missing_value) {
          Location location(0, 0, 0);  // Default initialization

          if (coordinate_system == "grid") {
            // For grid coordinates: use Z, Y, X and ignore geographic
            // coordinates
            location = Location(x, y, z);
          } else if (coordinate_system == "geographic") {
            // For geographic coordinates: use ZNU (level), XLAT_U (latitude),
            // XLONG_U (longitude)
            location =
                Location(xlat_u, xlong_u, znu, CoordinateSystem::GEOGRAPHIC);
          } else {
            throw std::runtime_error("Unsupported coordinate system: " +
                                     coordinate_system);
          }

          observations_.emplace_back(location, obs_value, error);
        }
      }
    }

    if (observations_.empty()) {
      throw std::runtime_error("No valid observations loaded from file: " +
                               filename);
    }
  }

 private:
  /**
   * @brief Save observation data to file
   * @param filename Path to save observation file
   */
  void saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open file for writing: " + filename);
    }

    for (const auto& obs : observations_) {
      if (obs.is_valid) {
        if (obs.location.getCoordinateSystem() ==
            CoordinateSystem::GEOGRAPHIC) {
          auto [lat, lon, level] = obs.location.getGeographicCoords();
          file << std::fixed << std::setprecision(6) << std::setw(10) << lat
               << std::setw(10) << lon << std::setw(10) << level
               << std::setw(12) << obs.value << std::setw(12) << obs.error
               << "\n";
        }
      }
    }
  }

  /**
   * @brief Get observations within a geographic bounding box
   * @param min_lat Minimum latitude
   * @param max_lat Maximum latitude
   * @param min_lon Minimum longitude
   * @param max_lon Maximum longitude
   * @return Vector of observations within the bounding box
   */
  std::vector<ObservationPoint> getObservationsInBox(double min_lat,
                                                     double max_lat,
                                                     double min_lon,
                                                     double max_lon) const {
    std::vector<ObservationPoint> result;
    for (const auto& obs : observations_) {
      if (obs.is_valid &&
          obs.location.getCoordinateSystem() == CoordinateSystem::GEOGRAPHIC) {
        auto [lat, lon, level] = obs.location.getGeographicCoords();
        if (lat >= min_lat && lat <= max_lat && lon >= min_lon &&
            lon <= max_lon) {
          result.push_back(obs);
        }
      }
    }
    return result;
  }

  /**
   * @brief Get observations within a vertical range
   * @param min_level Minimum vertical level
   * @param max_level Maximum vertical level
   * @return Vector of observations within the vertical range
   */
  std::vector<ObservationPoint> getObservationsInVerticalRange(
      double min_level, double max_level) const {
    std::vector<ObservationPoint> result;
    for (const auto& obs : observations_) {
      if (obs.is_valid &&
          obs.location.getCoordinateSystem() == CoordinateSystem::GEOGRAPHIC) {
        auto [lat, lon, level] = obs.location.getGeographicCoords();
        if (level >= min_level && level <= max_level) {
          result.push_back(obs);
        }
      }
    }
    return result;
  }

 private:
  /**
   * @brief Private constructor for cloning
   * @param other Observation to clone from
   * @param is_clone Flag to indicate this is a clone operation
   */
  GridObservation(const GridObservation& other, bool)
      : observations_(other.observations_),
        type_variable_map_(other.type_variable_map_),
        covariance_(other.covariance_) {}

  std::vector<ObservationPoint>
      observations_;  ///< Vector of observation points

  // Mapping from type/variable to observation indices (for backward
  // compatibility)
  std::unordered_map<std::string,
                     std::unordered_map<std::string, std::vector<size_t>>>
      type_variable_map_;

  std::vector<double> covariance_;  ///< Observation error covariance matrix
};

}  // namespace metada::backends::common::observation