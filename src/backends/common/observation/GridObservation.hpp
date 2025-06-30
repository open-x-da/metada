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

#include "PointObservation.hpp"

namespace metada::backends::common::observation {

using framework::CoordinateSystem;
using framework::Location;
using framework::ObservationPoint;

/**
 * @brief Iterator for observations
 */
class ObservationIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = ObservationPoint;
  using difference_type = std::ptrdiff_t;
  using pointer = ObservationPoint*;
  using reference = ObservationPoint&;

  ObservationIterator() = default;

  ObservationIterator(const std::vector<ObservationPoint>* data, size_t index)
      : data_(data), index_(index) {}

  // Dereference operators
  reference operator*() { return const_cast<reference>((*data_)[index_]); }
  pointer operator->() { return &const_cast<reference>((*data_)[index_]); }

  // Increment operators
  ObservationIterator& operator++() {
    ++index_;
    return *this;
  }

  ObservationIterator operator++(int) {
    ObservationIterator tmp = *this;
    ++index_;
    return tmp;
  }

  // Comparison operators
  bool operator==(const ObservationIterator& other) const {
    return data_ == other.data_ && index_ == other.index_;
  }

  bool operator!=(const ObservationIterator& other) const {
    return !(*this == other);
  }

 private:
  const std::vector<ObservationPoint>* data_;
  size_t index_;
};

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
  using iterator_type = ObservationIterator;
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

          // Load WRF observation file
          loadFromFile(filename, error, missing_value);
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
  ObservationIterator begin() const {
    return ObservationIterator(&observations_, 0);
  }

  /**
   * @brief Get iterator to end of observations
   * @return Iterator past last observation
   */
  ObservationIterator end() const {
    return ObservationIterator(&observations_, observations_.size());
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
   */
  void loadFromFile(const std::string& filename, double error,
                    double missing_value) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open observation file: " + filename);
    }

    // Check if this is a WRF observation file format
    std::string first_line;
    std::getline(file, first_line);
    file.seekg(0);  // Reset to beginning

    if (first_line.find("Time") != std::string::npos &&
        first_line.find("XLONG_U") != std::string::npos) {
      // This is a WRF observation format file
      loadFromWRFObsFile(filename, error, missing_value);
      return;
    }

    // For non-WRF formats, throw an error since we only support WRF observation
    // files now
    throw std::runtime_error(
        "Only WRF observation format files are supported. "
        "Expected file with 'Time' and 'XLONG_U' headers.");
  }

 private:
  /**
   * @brief Load observation data from WRF observation format file
   * @param filename Path to WRF observation file
   * @param error The observation error to assign to all valid points
   * @param missing_value The value that indicates a missing observation
   */
  void loadFromWRFObsFile(const std::string& filename, double error,
                          double missing_value) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open WRF observation file: " +
                               filename);
    }

    // Clear existing observations
    observations_.clear();

    std::string line;
    bool header_found = false;

    while (std::getline(file, line)) {
      // Skip header lines
      if (line.find("Time") != std::string::npos ||
          line.find("---") != std::string::npos) {
        header_found = true;
        continue;
      }

      if (!header_found || line.empty()) {
        continue;
      }

      // Parse data line: Time Z Y X ZNU XLONG_U XLAT_U U
      std::istringstream iss(line);
      int time, z, y, x;
      double znu, xlong_u, xlat_u, u_value;

      if (iss >> time >> z >> y >> x >> znu >> xlong_u >> xlat_u >> u_value) {
        // Skip missing values
        if (u_value != missing_value) {
          // Create geographic location using XLAT_U (latitude), XLONG_U
          // (longitude), ZNU (level)
          Location location(xlat_u, xlong_u, znu, CoordinateSystem::GEOGRAPHIC);
          observations_.emplace_back(location, u_value, error);
        }
      }
    }

    if (observations_.empty()) {
      throw std::runtime_error("No valid observations loaded from WRF file: " +
                               filename);
    }
  }

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