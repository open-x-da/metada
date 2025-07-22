/**
 * @file MACOMObservation.hpp
 * @brief MACOM observation backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This class provides an observation backend specifically designed for MACOM's
 * requirements. It supports geographic coordinate format observations which are
 * common in oceanographic data assimilation.
 */

#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "Location.hpp"
#include "PointObservation.hpp"

namespace metada::backends::macom {

using framework::CoordinateSystem;
using framework::Location;
using framework::ObservationPoint;

/**
 * @brief Iterator for MACOM observations
 */
class MACOMObservationIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = ObservationPoint;
  using difference_type = std::ptrdiff_t;
  using pointer = ObservationPoint*;
  using reference = ObservationPoint&;

  MACOMObservationIterator() = default;

  MACOMObservationIterator(const std::vector<ObservationPoint>* data,
                           size_t index)
      : data_(data), index_(index) {}

  // Dereference operators
  reference operator*() { return const_cast<reference>((*data_)[index_]); }
  pointer operator->() { return &const_cast<reference>((*data_)[index_]); }

  // Increment operators
  MACOMObservationIterator& operator++() {
    ++index_;
    return *this;
  }

  MACOMObservationIterator operator++(int) {
    MACOMObservationIterator tmp = *this;
    ++index_;
    return tmp;
  }

  // Comparison operators
  bool operator==(const MACOMObservationIterator& other) const {
    return data_ == other.data_ && index_ == other.index_;
  }

  bool operator!=(const MACOMObservationIterator& other) const {
    return !(*this == other);
  }

 private:
  const std::vector<ObservationPoint>* data_;
  size_t index_;
};

/**
 * @brief MACOM observation backend implementation
 *
 * @details
 * This class implements an observation backend specifically for MACOM that:
 * - Supports geographic coordinate format: lat lon level value error
 * - Stores observations as point observations with location information
 * - Handles missing values and quality control for oceanographic data
 * - Provides efficient spatial filtering operations
 * - Organizes data by observation type and variable
 */
class MACOMObservation {
 public:
  // =============================================================================
  // FRAMEWORK CONCEPTS REQUIRED INTERFACES
  // Required by ObservationBackendImpl concept
  // =============================================================================

  // Type aliases for concept compliance
  using iterator_type = MACOMObservationIterator;
  using value_type = ObservationPoint;

  // --- Resource management (required by framework) ---

  /**
   * @brief Default constructor (required for clone)
   */
  MACOMObservation() = default;

  /**
   * @brief Copy constructor is deleted (required by framework)
   */
  MACOMObservation(const MACOMObservation&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  MACOMObservation& operator=(const MACOMObservation&) = delete;

  /**
   * @brief Constructor that initializes from configuration (required by
   * framework)
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  explicit MACOMObservation(const ConfigBackend& config) {
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

      // Check if format is specified, default to "geographic" for MACOM
      std::string format = "geographic";
      if (type_backend.HasKey("format")) {
        format = type_backend.Get("format").asString();
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

          double missing_value = var_backend.Get("missing_value").asFloat();

          if (format == "geographic") {
            loadFromGeographicFile(filename, missing_value);
          } else {
            throw std::runtime_error(
                "MACOMObservation only supports geographic format");
          }
        }
      }
    }
  }

  /**
   * @brief Move constructor (required by framework)
   */
  MACOMObservation(MACOMObservation&& other) noexcept
      : observations_(std::move(other.observations_)),
        type_variable_map_(std::move(other.type_variable_map_)),
        covariance_(std::move(other.covariance_)) {}

  /**
   * @brief Move assignment operator (required by framework)
   */
  MACOMObservation& operator=(MACOMObservation&& other) noexcept {
    if (this != &other) {
      observations_ = std::move(other.observations_);
      type_variable_map_ = std::move(other.type_variable_map_);
      covariance_ = std::move(other.covariance_);
    }
    return *this;
  }

  /**
   * @brief Clone this observation (required by framework)
   * @return A new instance that is a copy of this one
   */
  std::unique_ptr<MACOMObservation> clone() const {
    auto new_obs = std::make_unique<MACOMObservation>();
    new_obs->observations_ = observations_;
    new_obs->type_variable_map_ = type_variable_map_;
    new_obs->covariance_ = covariance_;
    return new_obs;
  }

  // --- Iterator interface (required by framework) ---

  /**
   * @brief Get iterator to beginning of observations
   * @return Iterator to first observation
   */
  MACOMObservationIterator begin() const {
    return MACOMObservationIterator(&observations_, 0);
  }

  /**
   * @brief Get iterator to end of observations
   * @return Iterator past last observation
   */
  MACOMObservationIterator end() const {
    return MACOMObservationIterator(&observations_, observations_.size());
  }

  // --- Size and access interface (required by framework) ---

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

  // --- Data access interface (required by framework) ---

  /**
   * @brief Get raw pointer to observation data (required by framework)
   * @return Pointer to data
   */
  void* getData() {
    if (observations_.empty()) {
      return nullptr;
    }
    return const_cast<void*>(static_cast<const void*>(observations_.data()));
  }

  /**
   * @brief Get const raw pointer to observation data (required by framework)
   * @return Const pointer to data
   */
  const void* getData() const {
    if (observations_.empty()) {
      return nullptr;
    }
    return static_cast<const void*>(observations_.data());
  }

  /**
   * @brief Get template-specific data (required by framework)
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

  // --- Covariance interface (required by framework) ---

  /**
   * @brief Get the observation error variances (required by framework)
   * @return Vector of observation error variances (one per observation)
   */
  std::vector<double> getCovariance() const {
    std::vector<double> variances;
    variances.reserve(observations_.size());
    for (const auto& obs : observations_) {
      variances.push_back(obs.error * obs.error);  // Variance = error^2
    }
    return variances;
  }

  // --- Type/variable interface (required by framework) ---

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

  // --- Arithmetic operations (required by framework) ---

  /**
   * @brief Add another observation to this one
   * @param other Observation to add
   */
  void add(const MACOMObservation& other) {
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
  void subtract(const MACOMObservation& other) {
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

  // --- Comparison operations (required by framework) ---

  /**
   * @brief Compare equality with another observation
   * @param other Observation to compare with
   * @return true if equal, false otherwise
   */
  bool equals(const MACOMObservation& other) const {
    return observations_ == other.observations_ &&
           type_variable_map_ == other.type_variable_map_ &&
           covariance_ == other.covariance_;
  }

  // --- Initialization interface (required by framework) ---

  /**
   * @brief Initialize observation
   */
  void initialize() {
    // Initialization is handled in constructor
  }

  /**
   * @brief Load observation data from file (required by framework)
   * @param filename Path to observation file
   * @param error The observation error to assign to all valid points
   * @param missing_value The value that indicates a missing observation
   */
  void loadFromFile(const std::string& filename, double error,
                    double missing_value) {
    loadFromGeographicFile(filename, missing_value);

    // Update error values for all loaded observations
    for (auto& obs : observations_) {
      if (obs.is_valid) {
        obs.error = error;
      }
    }
  }

  // =============================================================================
  // MACOM SPECIFIC FUNCTIONALITY
  // These are MACOM-specific methods beyond framework requirements
  // =============================================================================

  /**
   * @brief Apply quality control to oceanographic observations
   */
  void applyQC() {
    for (auto& obs : observations_) {
      if (obs.is_valid) {
        // Oceanographic QC: Check for reasonable temperature/salinity ranges
        if (obs.value < -5.0 || obs.value > 40.0) {  // Temperature range
          obs.is_valid = false;
        }

        // Check for reasonable latitude/longitude
        if (obs.location.getCoordinateSystem() ==
            CoordinateSystem::GEOGRAPHIC) {
          auto [lat, lon, level] = obs.location.getGeographicCoords();
          if (lat < -90.0 || lat > 90.0 || lon < -180.0 || lon > 360.0) {
            obs.is_valid = false;
          }
        }
      }
    }
  }

  /**
   * @brief Load observation data from geographic coordinate file
   * @param filename Path to observation file with geographic coordinates
   * @param missing_value The value that indicates a missing observation.
   *
   * Expected format: lat lon level value error (one observation per line)
   * Lines starting with # are treated as comments and ignored
   */
  void loadFromGeographicFile(const std::string& filename,
                              double missing_value) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open observation file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
      // Skip comment lines
      if (line.empty() || line[0] == '#') {
        continue;
      }

      std::istringstream iss(line);
      double lat, lon, level, value, error;

      if (iss >> lat >> lon >> level >> value >> error) {
        if (value != missing_value) {
          // Create geographic location
          Location location(lat, lon, level, CoordinateSystem::GEOGRAPHIC);
          observations_.emplace_back(location, value, error);
        }
      }
      // If parsing fails, skip the line (could log warning in future)
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

    file << "# lat lon level value error\n";
    for (const auto& obs : observations_) {
      if (obs.is_valid) {
        if (obs.location.getCoordinateSystem() ==
            CoordinateSystem::GEOGRAPHIC) {
          auto [lat, lon, level] = obs.location.getGeographicCoords();
          file << std::fixed << std::setprecision(6) << lat << " " << lon << " "
               << level << " " << obs.value << " " << obs.error << "\n";
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

  /**
   * @brief Get the number of valid observations
   * @return Number of valid observations
   */
  size_t getValidCount() const {
    size_t count = 0;
    for (const auto& obs : observations_) {
      if (obs.is_valid) {
        count++;
      }
    }
    return count;
  }

  /**
   * @brief Get observation statistics
   * @return String with observation statistics
   */
  std::string getStatistics() const {
    std::stringstream ss;
    ss << "Observation Statistics:\n";
    ss << "  Total observations: " << observations_.size() << "\n";
    ss << "  Valid observations: " << getValidCount() << "\n";
    ss << "  Invalid observations: " << (observations_.size() - getValidCount())
       << "\n";

    if (!observations_.empty()) {
      double min_val = std::numeric_limits<double>::max();
      double max_val = std::numeric_limits<double>::lowest();
      double sum = 0.0;
      size_t valid_count = 0;

      for (const auto& obs : observations_) {
        if (obs.is_valid) {
          min_val = std::min(min_val, obs.value);
          max_val = std::max(max_val, obs.value);
          sum += obs.value;
          valid_count++;
        }
      }

      if (valid_count > 0) {
        ss << "  Value range: [" << min_val << ", " << max_val << "]\n";
        ss << "  Mean value: " << (sum / valid_count) << "\n";
      }
    }

    return ss.str();
  }

  /**
   * @brief Compute quadratic form with observation error covariance (required
   * by framework)
   * @param obs_increment Observation increment
   * @return Quadratic form value
   */
  double quadraticForm(const std::vector<double>& obs_increment) const {
    if (obs_increment.size() != observations_.size()) {
      throw std::invalid_argument("Observation increment size mismatch");
    }

    // For diagonal covariance: x^T R^-1 x = sum(x_i^2 / sigma_i^2)
    double result = 0.0;
    for (size_t i = 0; i < observations_.size(); ++i) {
      if (observations_[i].is_valid) {
        double error = observations_[i].error;
        if (error > 0.0) {
          result += (obs_increment[i] * obs_increment[i]) / (error * error);
        }
      }
    }
    return result;
  }

  /**
   * @brief Apply inverse observation error covariance (required by framework)
   * @param obs_increment Input observation increment
   * @return Result observation increment
   */
  std::vector<double> applyInverseCovariance(
      const std::vector<double>& obs_increment) const {
    if (obs_increment.size() != observations_.size()) {
      throw std::invalid_argument("Observation increment size mismatch");
    }

    std::vector<double> result(obs_increment.size());
    for (size_t i = 0; i < observations_.size(); ++i) {
      if (observations_[i].is_valid) {
        double error = observations_[i].error;
        if (error > 0.0) {
          result[i] = obs_increment[i] / (error * error);
        } else {
          result[i] = 0.0;
        }
      } else {
        result[i] = 0.0;
      }
    }
    return result;
  }

 private:
  std::vector<ObservationPoint>
      observations_;  ///< Vector of observation points

  // Mapping from type/variable to observation indices (for backward
  // compatibility)
  std::unordered_map<std::string,
                     std::unordered_map<std::string, std::vector<size_t>>>
      type_variable_map_;

  std::vector<double> covariance_;  ///< Observation error covariance matrix
};

}  // namespace metada::backends::macom