#pragma once

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "ObsRecord.hpp"
#include "WRFDAObsOperator_c_api.h"

namespace metada::backends::wrf {
using ObsRecord = framework::ObsRecord;

// Simple iterator implementation
template <typename GeometryBackend>
class WRFObservationIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = framework::ObsRecord;
  using difference_type = std::ptrdiff_t;
  using pointer = const value_type*;
  using reference = const value_type&;

  WRFObservationIterator() : index_(0) {}
  WRFObservationIterator(size_t index) : index_(index) {}

  bool operator==(const WRFObservationIterator& other) const {
    return index_ == other.index_;
  }
  bool operator!=(const WRFObservationIterator& other) const {
    return !(*this == other);
  }

  WRFObservationIterator& operator++() {
    ++index_;
    return *this;
  }

  reference operator*() const {
    static value_type dummy;
    return dummy;
  }

 private:
  size_t index_;
};

/**
 * @brief WRF observation backend that bridges framework adapters to observation
 * database
 *
 * This class serves as the WRF observation backend, handling the conversion
 * between METADA's observation representation and WRF-compatible formats
 * (including WRFDA iv_type/y_type structures). It provides methods to:
 * 1. Convert PrepBUFR observations to WRF format
 * 2. Handle different observation types (SYNOP, SOUND, etc.)
 * 3. Manage observation errors and QC flags
 * 4. Support WRF-specific coordinate transformations (Arakawa-C grid)
 * 5. Bridge between framework adapters and observation database
 */
template <typename GeometryBackend>
class WRFObservation {
 public:
  // Required type definitions for ObservationBackendImpl concept
  using iterator_type = WRFObservationIterator<GeometryBackend>;
  using value_type = ObsRecord;

  // Level-specific observation data structure
  struct LevelData {
    double level;     // Vertical level (pressure, height, or model level)
    double pressure;  // Pressure at this level (Pa)
    double height;    // Height at this level (m)

    // Variable-specific data
    struct VariableData {
      double value;    // Observed value
      double error;    // Observation error
      double bias;     // Bias correction
      int qc;          // Quality control flag
      bool available;  // Whether this variable is available
    };

    // All possible variables at this level
    VariableData u, v, t, q, p, slp, pw, ref;  // Standard variables
    VariableData td2m, rh2m, t2m, u10m, v10m;  // 2m/10m variables
  };

  // WRFDA-based observation operator data structure
  struct ObsOperatorData {
    // Number of stations (not levels)
    size_t num_obs;

    // Station-level data (one per station - no duplication)
    std::vector<double> lats;                // Station latitudes
    std::vector<double> lons;                // Station longitudes
    std::vector<double> elevations;          // Station elevations
    std::vector<std::string> station_ids;    // Station identifiers
    std::vector<std::string> station_names;  // Station names
    std::vector<std::string> types;          // Station observation types

    // Station-level metadata
    std::vector<int> qc_flags;         // Station quality flags
    std::vector<int> usage_flags;      // Usage flags
    std::vector<double> time_offsets;  // Time offsets

    // Level-specific data (per station)
    std::vector<std::vector<LevelData>> levels;  // levels[station][level]

    // Compatibility methods for WRFDAObsOperator
    std::vector<double> getValues() const {
      std::vector<double> result;
      for (const auto& station_levels : levels) {
        for (const auto& level : station_levels) {
          if (level.u.available) result.push_back(level.u.value);
          if (level.v.available) result.push_back(level.v.value);
          if (level.t.available) result.push_back(level.t.value);
          if (level.p.available) result.push_back(level.p.value);
          if (level.q.available) result.push_back(level.q.value);
        }
      }
      return result;
    }

    std::vector<double> getErrors() const {
      std::vector<double> result;
      for (const auto& station_levels : levels) {
        for (const auto& level : station_levels) {
          if (level.u.available) result.push_back(level.u.error);
          if (level.v.available) result.push_back(level.v.error);
          if (level.t.available) result.push_back(level.t.error);
          if (level.p.available) result.push_back(level.p.error);
          if (level.q.available) result.push_back(level.q.error);
        }
      }
      return result;
    }

    // Legacy compatibility members (computed on demand)
    std::vector<double> values() const { return getValues(); }
    std::vector<double> errors() const { return getErrors(); }
  };

  // Delete default constructor and copying
  WRFObservation() = delete;
  WRFObservation(const WRFObservation&) = delete;
  WRFObservation& operator=(const WRFObservation&) = delete;

  /**
   * @brief Construct from config and geometry (with geometry filtering)
   * @param config Configuration object containing observation settings
   * @param geometry WRF geometry for domain filtering
   *
   * @details Uses WRFDA's standard observation reading pipeline to load
   * observations from BUFR file. Calls da_setup_obs_structures_bufr which
   * reads the file and allocates both iv_type and y_type structures.
   * The iv_type pointer is passed to WRFObsOperator for innovation computation.
   */
  WRFObservation(const backends::config::YamlConfig& config,
                 const GeometryBackend& geometry)
      : geometry_(std::shared_ptr<GeometryBackend>(), &geometry) {
    // Initialize WRFDA data structures using WRFDA's standard pipeline
    initializeWithWRFDA(config, geometry);
  }

  // Move semantics
  WRFObservation(WRFObservation&& other) noexcept = default;
  WRFObservation& operator=(WRFObservation&& other) noexcept = default;

  /**
   * @brief Initialize observation data structures (framework requirement)
   * @details This is called by the framework. For WRFDA-based initialization,
   * the work is already done in the constructor via initializeWithWRFDA.
   */
  void initialize() {
    // For WRFDA-based initialization, all work is done in constructor
    // This method is kept for framework compatibility
    if (!iv_type_data_ || !y_type_data_) {
      throw std::runtime_error(
          "WRFObservation not properly initialized. "
          "Use constructor with config and geometry.");
    }
  }

  /**
   * @brief Initialize observation data structures using WRFDA pipeline
   * @param config Configuration containing observation file path
   * @param geometry WRF geometry for domain bounds checking
   */
  void initializeWithWRFDA(const backends::config::YamlConfig& config,
                           const GeometryBackend& geometry) {
    // Get observation file path from config
    // For now, skip WRFDA reading if obs_file not configured
    // This allows gradual migration from PrepBUFRObservation to WRFDA
    std::string obs_file;
    bool use_wrfda_reading = false;

    try {
      auto obs_file_value = config.Get("obs_file");
      if (obs_file_value.isString()) {
        obs_file = obs_file_value.asString();
        use_wrfda_reading = true;
      }
    } catch (...) {
      // obs_file not found - fallback to PrepBUFRObservation
      std::cout << "Note: obs_file not configured, using PrepBUFRObservation "
                   "(legacy mode)"
                << std::endl;
      std::cout << "To use WRFDA observation reading, add 'obs_file: <path>' "
                   "to observation config"
                << std::endl;
      use_wrfda_reading = false;
    }

    if (!use_wrfda_reading) {
      // Legacy mode: PrepBUFRObservation already initialized in constructor
      return;
    }

    std::cout << "Loading observations with WRFDA from: " << obs_file
              << std::endl;

    // Get WRFDA grid pointer from geometry
    void* grid_ptr = geometry.getGridPtr();
    if (!grid_ptr) {
      throw std::runtime_error("Invalid grid pointer from geometry");
    }

    // Call WRFDA's standard observation reading pipeline
    // This allocates and populates both iv_type and y_type structures
    int rc = wrfda_read_and_allocate_observations(
        grid_ptr, obs_file.c_str(), static_cast<int>(obs_file.length()),
        &iv_type_data_, &y_type_data_);

    if (rc != 0) {
      throw std::runtime_error("Failed to read observations with WRFDA, code " +
                               std::to_string(rc));
    }

    // Extract observation metadata for framework adapter queries
    extractObservationMetadata();

    std::cout << "WRFDA observation structures initialized successfully"
              << std::endl;
  }

  /**
   * @brief Get pointer to WRFDA iv_type structure
   * @return Raw pointer to iv_type data
   */
  void* getIVTypeData() const { return iv_type_data_; }

  /**
   * @brief Get pointer to WRFDA y_type structure
   * @return Raw pointer to y_type data
   */
  void* getYTypeData() const { return y_type_data_; }

  /**
   * @brief Get iterator to beginning of observations
   * @return Iterator to first observation
   */
  iterator_type begin() const { return iterator_type(0); }

  /**
   * @brief Get iterator to end of observations
   * @return Iterator past last observation
   */
  iterator_type end() const { return iterator_type(size()); }

  /**
   * @brief Get number of observations
   * @return Total number of observations
   */
  size_t size() const { return getTotalObservationCount(); }

  /**
   * @brief Get observation at specific index
   * @param index Index of the observation
   * @return Reference to the observation record
   */
  const value_type& operator[](size_t index) const {
    if (index >= size()) {
      throw std::out_of_range("Observation index out of range");
    }
    // TODO: Implement access to specific observation record
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Get raw data pointer
   * @return Void pointer to underlying data
   */
  void* getData() { return iv_type_data_; }

  /**
   * @brief Get raw data pointer (const)
   * @return Const void pointer to underlying data
   */
  const void* getData() const { return iv_type_data_; }

  /**
   * @brief Get all observation values as a flat vector
   * @return Vector containing all observation values in order
   *
   * @details Extracts all individual observation values (u, v, t, p, q)
   * from all stations and levels into a single flat vector. Values are
   * ordered consistently across all observations.
   */
  std::vector<double> getObservationValues() const {
    std::vector<double> values;

    if (!iv_type_data_ || !y_type_data_) {
      return values;  // Return empty if not initialized
    }

    // Extract observation values from WRFDA y_type and iv_type structures
    // First, get the count from wrfda_count_innovations
    int num_obs = 0;
    int rc = wrfda_count_innovations(iv_type_data_, &num_obs);
    if (rc != 0 || num_obs <= 0) {
      return values;  // Return empty if no observations
    }

    // Allocate and extract observations
    values.resize(num_obs);
    int actual_obs = 0;
    rc = wrfda_extract_observations(iv_type_data_, y_type_data_, values.data(),
                                    &actual_obs);
    if (rc != 0 || actual_obs != num_obs) {
      values.clear();  // Return empty if extraction failed
      return values;
    }

    return values;
  }

  /**
   * @brief Get the size of observation values vector
   * @return Number of individual observation values
   *
   * @details Returns the total number of individual observation values
   * (across all variables and levels) that would be returned by
   * getObservationValues(). This is useful for pre-allocating
   * vectors or checking sizes without extracting the full data.
   */
  size_t getObservationValuesSize() const { return getTotalObservationCount(); }

  /**
   * @brief Get number of observations by type
   * @param obs_type Observation type (e.g., "ADPSFC", "ADPSFC_u", "ADPSFC_t")
   * @return Number of observations of that type
   */
  size_t getObservationCount(const std::string& obs_type) const {
    auto it = obs_counts_.find(obs_type);
    return it != obs_counts_.end() ? it->second : 0;
  }

  /**
   * @brief Get number of observations by base type (sum of all variables)
   * @param base_type Base observation type (e.g., "ADPSFC", "ADPUPA")
   * @return Total number of observations of that base type across all variables
   */
  size_t getObservationCountByBaseType(const std::string& base_type) const {
    size_t total = 0;
    for (const auto& [type, count] : obs_counts_) {
      if (type.substr(0, base_type.length()) == base_type) {
        total += count;
      }
    }
    return total;
  }

  /**
   * @brief Get total number of observations across all types
   * @return Total observation count
   */
  size_t getTotalObservationCount() const {
    size_t total = 0;
    for (const auto& [type, count] : obs_counts_) {
      total += count;
    }
    return total;
  }

  /**
   * @brief Apply quality control checks
   * This performs WRFDA-specific QC in addition to basic checks
   */
  void applyQC() {
    // TODO: Implement WRFDA QC checks including:
    // - Domain boundary checks
    // - Variable range checks
    // - Temporal checks
    // - Consistency checks between variables
  }

  /**
   * @brief Get observation error statistics
   * @return Map of observation type to error statistics
   */
  std::unordered_map<std::string, double> getErrorStats() const {
    return error_stats_;
  }

  /**
   * @brief Get observation error covariance matrix
   * @return Vector containing covariance matrix diagonal elements
   *
   * @details Returns the diagonal elements of the observation error covariance
   * matrix R. WRFDA handles error covariance internally via iv_type structure.
   * For now, return unit errors as placeholder.
   */
  std::vector<double> getCovariance() const {
    // WRFDA handles error covariance internally
    // Return unit variance as placeholder
    return std::vector<double>(getTotalObservationCount(), 1.0);
  }

  /**
   * @brief Get size of observations for a specific type and variable
   * @param typeName Observation type name
   * @param varName Variable name
   * @return Number of observations
   */
  size_t getSize(const std::string& typeName,
                 const std::string& varName) const {
    // Extract base type name (remove variable suffix if present)
    std::string base_type = typeName;
    size_t last_underscore = typeName.find_last_of('_');
    if (last_underscore != std::string::npos) {
      base_type = typeName.substr(0, last_underscore);
    }

    // Convert to lowercase for WRFDA lookup
    std::transform(base_type.begin(), base_type.end(), base_type.begin(),
                   ::tolower);

    // Get size from WRFDA structures
    if (!iv_type_data_) {
      return 0;  // Return 0 if not initialized
    }

    // Look up the count for this type and variable combination
    std::string type_var_key = base_type + "_" + varName;
    if (obs_counts_.find(type_var_key) != obs_counts_.end()) {
      return obs_counts_.at(type_var_key);
    }

    // Fallback: return count for base type only
    if (obs_counts_.find(base_type) != obs_counts_.end()) {
      return obs_counts_.at(base_type);
    }

    return 0;
  }

  /**
   * @brief Get observation operator data
   * @return ObsOperatorData structure
   */
  ObsOperatorData getObsOperatorData() const {
    ObsOperatorData data;

    if (!iv_type_data_) {
      return data;  // Return empty if not initialized
    }

    // TODO: Extract observation operator data from WRFDA iv_type structure
    // This requires accessing WRFDA Fortran structures directly via C API
    // For now, return empty structure - this needs to be implemented
    // using WRFDA's native data access patterns

    return data;
  }

  // Required by ObservationBackendImpl concept
  double quadraticForm(const std::vector<double>& innovation) const {
    // Compute quadratic form: innovation^T * R^-1 * innovation
    // WRFDA handles error covariance internally via iv_type structure
    // For now, use simple diagonal with unit errors
    double result = 0.0;
    for (size_t i = 0; i < innovation.size(); ++i) {
      result += innovation[i] * innovation[i];
    }
    return result;
  }

  std::vector<double> applyInverseCovariance(
      const std::vector<double>& vector) const {
    // Apply R^-1 element-wise for diagonal covariance matrix
    // WRFDA handles error covariance internally via iv_type structure
    // For now, use identity (unit errors)
    return vector;
  }

  std::vector<double> applyCovariance(const std::vector<double>& vector) const {
    // Apply R element-wise for diagonal covariance matrix
    // WRFDA handles error covariance internally via iv_type structure
    // For now, use identity (unit errors)
    return vector;
  }

  std::vector<double> getInverseCovarianceDiagonal() const {
    // WRFDA handles error covariance internally
    // Return unit inverse variance as placeholder
    return std::vector<double>(getTotalObservationCount(), 1.0);
  }

  bool isDiagonalCovariance() const {
    return true;  // WRFDA uses diagonal R by default
  }

  /**
   * @brief Add another observation to this one
   * @param other Observation to add
   */
  void add(const WRFObservation& other) {
    // TODO: Implement observation addition
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Subtract another observation from this one
   * @param other Observation to subtract
   */
  void subtract(const WRFObservation& other) {
    // TODO: Implement observation subtraction
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Multiply observation by scalar
   * @param scalar Value to multiply by
   */
  void multiply(double scalar) {
    // TODO: Implement scalar multiplication
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Check equality with another observation
   * @param other Observation to compare with
   * @return True if equal, false otherwise
   */
  bool equals(const WRFObservation& other) const {
    // TODO: Implement equality comparison
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Get observation type names
   * @return Vector of observation type names (includes variable-specific types)
   */
  std::vector<std::string> getTypeNames() const {
    std::vector<std::string> types;
    for (const auto& [type, count] : obs_counts_) {
      types.push_back(type);
    }
    return types;
  }

  /**
   * @brief Get base observation type names (without variable suffixes)
   * @return Vector of base observation type names
   */
  std::vector<std::string> getBaseTypeNames() const {
    std::set<std::string> base_types;
    for (const auto& [type, count] : obs_counts_) {
      // Extract base type by finding the last underscore and taking everything
      // before it
      size_t last_underscore = type.find_last_of('_');
      if (last_underscore != std::string::npos) {
        base_types.insert(type.substr(0, last_underscore));
      } else {
        base_types.insert(type);
      }
    }
    return std::vector<std::string>(base_types.begin(), base_types.end());
  }

  /**
   * @brief Get variable names for observation type
   * @param typeName Observation type name (may include variable suffix like
   * "ADPSFC_u")
   * @return Vector of variable names
   */
  std::vector<std::string> getVariableNames(const std::string& typeName) const {
    // Extract base type name (remove variable suffix if present)
    std::string base_type = typeName;
    std::string var_suffix;
    size_t last_underscore = typeName.find_last_of('_');
    if (last_underscore != std::string::npos) {
      base_type = typeName.substr(0, last_underscore);
      var_suffix = typeName.substr(last_underscore + 1);
    }

    // Convert to lowercase for WRFDA lookup
    std::transform(base_type.begin(), base_type.end(), base_type.begin(),
                   ::tolower);

    // Get variables from WRFDA type definitions
    std::vector<std::string> variables;

    if (!iv_type_data_) {
      return variables;  // Return empty if not initialized
    }

    // Define standard variables for each WRFDA observation type
    // Based on WRFDA's iv_type structure and error factors (da_control.f90)
    if (base_type == "synop" || base_type == "metar" || base_type == "ships" ||
        base_type == "buoy" || base_type == "sonde_sfc" ||
        base_type == "tamdar_sfc") {
      variables = {"U", "V", "T", "P", "Q"};
    } else if (base_type == "sound" || base_type == "airep" ||
               base_type == "tamdar" || base_type == "mtgirs") {
      variables = {"U", "V", "T", "Q"};
    } else if (base_type == "pilot" || base_type == "profiler") {
      variables = {"U", "V"};
    } else if (base_type == "gpspw") {
      variables = {"TPW"};
    } else if (base_type == "gpsrf" || base_type == "gpsref") {
      variables = {"REF", "P", "T", "Q"};
    } else if (base_type == "gpseph") {
      variables = {"EPH"};
    } else if (base_type == "geoamv" || base_type == "polaramv") {
      variables = {"U", "V"};
    } else if (base_type == "qscat") {
      variables = {"U", "V"};
    } else if (base_type == "radar") {
      variables = {"RV", "RF", "RR"};
    } else if (base_type == "bogus") {
      variables = {"U", "V", "T", "P", "Q", "SLP"};
    } else if (base_type == "airsr") {
      variables = {"T", "Q"};
    } else if (base_type == "rain") {
      variables = {"R"};
    } else if (base_type == "lightning") {
      variables = {"W", "DIV", "QV"};
    } else if (base_type == "ssmi_rv") {
      variables = {"SPEED", "TPW"};
    } else if (base_type == "satem") {
      variables = {"THICKNESS"};
    } else if (base_type == "ssmt1") {
      variables = {"T"};
    } else if (base_type == "ssmt2") {
      variables = {"RH"};
    } else {
      // Default: return empty for unknown types
      return std::vector<std::string>();
    }

    // If type had a variable suffix, filter to that variable only
    if (!var_suffix.empty()) {
      std::transform(var_suffix.begin(), var_suffix.end(), var_suffix.begin(),
                     ::toupper);
      if (std::find(variables.begin(), variables.end(), var_suffix) !=
          variables.end()) {
        return {var_suffix};
      } else {
        return {};
      }
    }

    return variables;
  }

  /**
   * @brief Create a clone of this observation
   * @return Unique pointer to cloned observation
   */
  std::unique_ptr<WRFObservation> clone() const {
    // TODO: Implement deep copy
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Get observations in geographic bounding box
   * @param min_lat Minimum latitude
   * @param max_lat Maximum latitude
   * @param min_lon Minimum longitude
   * @param max_lon Maximum longitude
   * @return Vector of observations in box
   */
  std::vector<value_type> getObservationsInBox(double min_lat, double max_lat,
                                               double min_lon,
                                               double max_lon) const {
    // TODO: Implement geographic filtering
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Get observations in vertical range
   * @param min_level Minimum level
   * @param max_level Maximum level
   * @return Vector of observations in range
   */
  std::vector<value_type> getObservationsInVerticalRange(
      double min_level, double max_level) const {
    // TODO: Implement vertical filtering
    throw std::runtime_error("Not implemented");
  }

 private:
  /**
   * @brief Allocate and populate WRFDA y_type structure
   * @param obs_data Organized observation data
   */
  void allocateYType(const ObsOperatorData& obs_data) {
    // Delegate to WRFDA bridge for y_type construction
    // This uses the same wrfda_construct_y_type function as WRFDAObsOperator

    if (obs_data.num_obs == 0) {
      return;
    }

    // Determine the operator family (default to "synop" for surface
    // observations)
    std::string family = "synop";

    // Extract structured data for WRFDA y_type construction
    std::vector<double> u_values, v_values, t_values, p_values, q_values;
    std::vector<double> u_errors, v_errors, t_errors, p_errors, q_errors;
    std::vector<double> lats, lons;
    std::vector<int> u_available, v_available, t_available, p_available,
        q_available;

    // Prepare observation types as a flat string
    std::string obs_types_flat;
    for (size_t i = 0; i < obs_data.num_obs; ++i) {
      std::string obs_type =
          (i < obs_data.types.size()) ? obs_data.types[i] : "ADPSFC";
      obs_types_flat += obs_type + std::string(20 - obs_type.length(), '\0');
    }

    // Extract data from the hierarchical structure
    for (size_t station = 0; station < obs_data.num_obs; ++station) {
      lats.push_back(obs_data.lats[station]);
      lons.push_back(obs_data.lons[station]);

      // Process each level for this station
      for (const auto& level : obs_data.levels[station]) {
        // Use WRFDA standard missing value (-888888.0) for unavailable
        // observations
        const double missing_r = -888888.0;
        u_values.push_back(level.u.available ? level.u.value : missing_r);
        v_values.push_back(level.v.available ? level.v.value : missing_r);
        t_values.push_back(level.t.available ? level.t.value : missing_r);
        p_values.push_back(level.p.available ? level.p.value : missing_r);
        q_values.push_back(level.q.available ? level.q.value : missing_r);

        u_errors.push_back(level.u.available ? level.u.error : 1.0);
        v_errors.push_back(level.v.available ? level.v.error : 1.0);
        t_errors.push_back(level.t.available ? level.t.error : 1.0);
        p_errors.push_back(level.p.available ? level.p.error : 1.0);
        q_errors.push_back(level.q.available ? level.q.error : 1.0);

        // Extract availability flags
        u_available.push_back(level.u.available ? 1 : 0);
        v_available.push_back(level.v.available ? 1 : 0);
        t_available.push_back(level.t.available ? 1 : 0);
        p_available.push_back(level.p.available ? 1 : 0);
        q_available.push_back(level.q.available ? 1 : 0);
      }
    }

    // Call the WRFDA bridge function to construct y_type
    int num_obs_int = static_cast<int>(obs_data.num_obs);
    int num_levels_int =
        static_cast<int>(u_values.size());  // Total number of level records

    void* y_ptr = wrfda_construct_y_type(
        &num_obs_int, &num_levels_int, const_cast<double*>(u_values.data()),
        const_cast<double*>(v_values.data()),
        const_cast<double*>(t_values.data()),
        const_cast<double*>(p_values.data()),
        const_cast<double*>(q_values.data()),
        const_cast<double*>(u_errors.data()),
        const_cast<double*>(v_errors.data()),
        const_cast<double*>(t_errors.data()),
        const_cast<double*>(p_errors.data()),
        const_cast<double*>(q_errors.data()),
        const_cast<int*>(u_available.data()),
        const_cast<int*>(v_available.data()),
        const_cast<int*>(t_available.data()),
        const_cast<int*>(p_available.data()),
        const_cast<int*>(q_available.data()), const_cast<double*>(lats.data()),
        const_cast<double*>(lons.data()),
        const_cast<char*>(obs_types_flat.c_str()),
        const_cast<char*>(family.c_str()));

    // Store the pointer (memory is managed by Fortran persistent structures)
    if (y_ptr) {
      y_type_data_ = y_ptr;
    }

    // Populate observation counts by type
    updateObservationCounts(obs_data);

    // Set up observation error covariance matrix R
    setupCovarianceMatrix(obs_data);
  }

  /**
   * @brief Set up observation error covariance matrix R
   * @param obs_data Observation data to extract error information
   *
   * @details Sets up the diagonal elements of R (variances) and R^-1 (inverse
   * variances) once during initialization. This avoids recomputing them in
   * every covariance operation.
   */
  void setupCovarianceMatrix(const ObsOperatorData& obs_data) {
    // Clear existing matrices
    r_matrix_.clear();
    r_inverse_matrix_.clear();

    // Reserve space for efficiency
    size_t total_obs = getTotalObservationCount();
    r_matrix_.reserve(total_obs);
    r_inverse_matrix_.reserve(total_obs);

    // Extract diagonal elements of R (variances) for diagonal covariance matrix
    // R[i,i] = error[i]^2, R^-1[i,i] = 1/error[i]^2
    for (size_t station = 0; station < obs_data.num_obs; ++station) {
      for (const auto& level : obs_data.levels[station]) {
        // Process each available variable in consistent order (u, v, t, p, q)
        if (level.u.available) {
          double variance = level.u.error * level.u.error;
          r_matrix_.push_back(variance);
          r_inverse_matrix_.push_back((variance > 0.0) ? 1.0 / variance : 0.0);
        }

        if (level.v.available) {
          double variance = level.v.error * level.v.error;
          r_matrix_.push_back(variance);
          r_inverse_matrix_.push_back((variance > 0.0) ? 1.0 / variance : 0.0);
        }

        if (level.t.available) {
          double variance = level.t.error * level.t.error;
          r_matrix_.push_back(variance);
          r_inverse_matrix_.push_back((variance > 0.0) ? 1.0 / variance : 0.0);
        }

        if (level.p.available) {
          double variance = level.p.error * level.p.error;
          r_matrix_.push_back(variance);
          r_inverse_matrix_.push_back((variance > 0.0) ? 1.0 / variance : 0.0);
        }

        if (level.q.available) {
          double variance = level.q.error * level.q.error;
          r_matrix_.push_back(variance);
          r_inverse_matrix_.push_back((variance > 0.0) ? 1.0 / variance : 0.0);
        }
      }
    }
  }

  /**
   * @brief Update observation counts by type
   * @param obs_data Observation data to count
   *
   * @details Populates obs_counts_ map with the count of each observation type
   * found in the data. Counts all individual observation values including
   * all variables (u, v, t, p, q) and all levels per station.
   */
  void updateObservationCounts(const ObsOperatorData& obs_data) {
    obs_counts_.clear();

    // Count all individual observation values by type
    for (size_t station = 0; station < obs_data.num_obs; ++station) {
      std::string obs_type = (station < obs_data.types.size())
                                 ? obs_data.types[station]
                                 : "ADPSFC";

      // Count all variables and levels for this station
      for (const auto& level : obs_data.levels[station]) {
        // Count each available variable as a separate observation
        if (level.u.available) obs_counts_[obs_type + "_u"]++;
        if (level.v.available) obs_counts_[obs_type + "_v"]++;
        if (level.t.available) obs_counts_[obs_type + "_t"]++;
        if (level.p.available) obs_counts_[obs_type + "_p"]++;
        if (level.q.available) obs_counts_[obs_type + "_q"]++;
      }
    }
  }

  /**
   * @brief Extract observation metadata from WRFDA iv_type structure
   * @details Extracts observation counts and other metadata from the
   * WRFDA-allocated iv_type structure for framework adapter queries
   */
  void extractObservationMetadata() {
    if (!iv_type_data_) {
      return;
    }

    // Get total observation count
    int total_count = 0;
    int rc = wrfda_get_total_obs_count(iv_type_data_, &total_count);
    if (rc != 0) {
      throw std::runtime_error("Failed to extract total observation count");
    }

    std::cout << "Total observations loaded: " << total_count << std::endl;

    // Get observation counts by type (num_ob_indexes = 31)
    constexpr int num_ob_indexes = 31;
    std::vector<int> type_counts(num_ob_indexes, 0);
    rc = wrfda_get_obs_type_counts(iv_type_data_, type_counts.data());
    if (rc != 0) {
      throw std::runtime_error("Failed to extract observation type counts");
    }

    // Map WRFDA type indices to names and populate obs_counts_
    const std::vector<std::string> type_names = {
        "sound",      "synop",      "pilot",   "satem",     "geoamv",
        "polaramv",   "airep",      "gpspw",   "gpsref",    "metar",
        "ships",      "ssmi_rv",    "ssmi_tb", "ssmt1",     "ssmt2",
        "qscat",      "profiler",   "buoy",    "bogus",     "pseudo",
        "radar",      "radiance",   "airsr",   "sonde_sfc", "mtgirs",
        "tamdar",     "tamdar_sfc", "rain",    "gpseph",    "lightning",
        "chemic_surf"};

    obs_counts_.clear();
    for (int i = 0;
         i < num_ob_indexes && i < static_cast<int>(type_names.size()); ++i) {
      if (type_counts[i] > 0) {
        obs_counts_[type_names[i]] = type_counts[i];
        std::cout << "  " << type_names[i] << ": " << type_counts[i]
                  << " observations" << std::endl;
      }
    }
  }

  /**
   * @brief Apply geometry filtering to observations
   * @param obs_data Observation data to filter
   *
   * @details Geometry filtering is handled by WRFDA during observation reading
   * (da_read_obs_bufr checks domain bounds automatically). This method is
   * kept for backward compatibility but is no longer needed.
   */
  void applyGeometryFiltering(ObsOperatorData& obs_data) {
    // Geometry filtering is now handled by WRFDA's da_read_obs_bufr
    // which checks domain bounds and sets outside flags
    // No additional filtering needed here
    (void)obs_data;  // Suppress unused parameter warning
  }

  // Data members
  std::shared_ptr<const GeometryBackend>
      geometry_;  // Geometry reference (for grid pointer access)

  // WRFDA observation metadata (extracted from WRFDA structures)
  std::map<std::string, size_t> obs_counts_;  // Type -> count mapping

  // WRFDA data structures (non-owning pointers to WRFDA-managed memory)
  // These persist in WRFDA Fortran modules and are managed by WRFDA
  void* iv_type_data_ = nullptr;  // Non-owning pointer to WRFDA iv_type
  void* y_type_data_ = nullptr;   // Non-owning pointer to WRFDA y_type

  // Statistics
  std::unordered_map<std::string, double> error_stats_;  // Error stats by type

  // Observation error covariance matrix (diagonal elements)
  std::vector<double> r_matrix_;  // R diagonal elements (variances)
  std::vector<double>
      r_inverse_matrix_;  // R^-1 diagonal elements (1/variances)
};

}  // namespace metada::backends::wrf