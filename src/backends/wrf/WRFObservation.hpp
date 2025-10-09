#pragma once

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "ObsRecord.hpp"
#include "PrepBUFRObservation.hpp"
#include "WRFDAObsOperator_c_api.h"

namespace metada::backends::wrf {
using ObsRecord = framework::ObsRecord;
using PrepBUFRObservation = common::observation::PrepBUFRObservation;

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
  using ObsOperatorData = PrepBUFRObservation::ObsOperatorData;

  // Delete default constructor and copying
  WRFObservation() = delete;
  WRFObservation(const WRFObservation&) = delete;
  WRFObservation& operator=(const WRFObservation&) = delete;

  /**
   * @brief Construct from config only (no geometry filtering)
   * @param config Configuration object containing observation settings
   *
   * @details This constructor does not use geometry.
   * No geometry filtering is applied - all observations are included.
   */
  explicit WRFObservation(const backends::config::YamlConfig& config)
      : obs_(std::make_shared<PrepBUFRObservation>(config)),
        geometry_(nullptr),
        apply_geometry_filtering_(false) {
    // Initialize WRFDA data structures
    initialize();
  }

  /**
   * @brief Construct from config and geometry (with geometry filtering)
   * @param config Configuration object containing observation settings
   * @param geometry WRF geometry for domain filtering
   *
   * @details This constructor uses the provided geometry to filter
   * observations. Only observations within the geometry domain will be
   * included.
   */
  WRFObservation(const backends::config::YamlConfig& config,
                 const GeometryBackend& geometry)
      : obs_(std::make_shared<PrepBUFRObservation>(config)),
        geometry_(std::shared_ptr<GeometryBackend>(), &geometry),
        apply_geometry_filtering_(true) {
    // Initialize WRFDA data structures with geometry filtering
    initialize();
  }

  // Move semantics
  WRFObservation(WRFObservation&& other) noexcept = default;
  WRFObservation& operator=(WRFObservation&& other) noexcept = default;

  /**
   * @brief Initialize observation data structures
   */
  void initialize() {
    // Extract observation data organized by WRFDA types
    auto obs_data = obs_->getObsOperatorData();

    // Apply geometry filtering if requested
    if (apply_geometry_filtering_) {
      applyGeometryFiltering(obs_data);
    }

    // Allocate and populate y_type structure (this also updates obs_counts_)
    allocateYType(obs_data);
  }

  /**
   * @brief Get pointer to WRFDA iv_type structure
   * @return Raw pointer to iv_type data
   */
  void* getIVTypeData() { return iv_type_data_; }

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

    // Get the observation data from the underlying PrepBUFRObservation
    auto obs_data = obs_->getObsOperatorData();

    // Extract all observation values in order
    for (size_t station = 0; station < obs_data.num_obs; ++station) {
      for (const auto& level : obs_data.levels[station]) {
        // Add available observation values in a consistent order (u, v, t, p,
        // q)
        if (level.u.available) values.push_back(level.u.value);
        if (level.v.available) values.push_back(level.v.value);
        if (level.t.available) values.push_back(level.t.value);
        if (level.p.available) values.push_back(level.p.value);
        if (level.q.available) values.push_back(level.q.value);
      }
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
   * matrix R. For diagonal matrices, this is equivalent to the variance of each
   * observation. The elements are returned in the same order as the observation
   * values.
   */
  std::vector<double> getCovariance() const {
    // Return pre-computed R matrix diagonal elements (variances)
    return r_matrix_;
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

    // Convert to lowercase for PrepBUFR lookup
    std::transform(base_type.begin(), base_type.end(), base_type.begin(),
                   ::tolower);

    // Delegate to PrepBUFR observation
    return obs_->getSize(base_type, varName);
  }

  /**
   * @brief Get observation operator data
   * @return ObsOperatorData structure
   */
  ObsOperatorData getObsOperatorData() const {
    return obs_->getObsOperatorData();
  }

  // Required by ObservationBackendImpl concept
  double quadraticForm(const std::vector<double>& innovation) const {
    // Validate input vector size
    if (innovation.size() != r_inverse_matrix_.size()) {
      throw std::invalid_argument(
          "Innovation vector size mismatch in quadraticForm: expected " +
          std::to_string(r_inverse_matrix_.size()) + ", got " +
          std::to_string(innovation.size()));
    }

    // Compute quadratic form: innovation^T * R^-1 * innovation
    // For diagonal R: result = sum(innovation[i]^2 / error[i]^2)
    double result = 0.0;
    for (size_t i = 0; i < innovation.size(); ++i) {
      result += innovation[i] * innovation[i] * r_inverse_matrix_[i];
    }

    return result;
  }

  std::vector<double> applyInverseCovariance(
      const std::vector<double>& vector) const {
    // Validate input vector size
    if (vector.size() != r_inverse_matrix_.size()) {
      throw std::invalid_argument(
          "Vector size mismatch in applyInverseCovariance: expected " +
          std::to_string(r_inverse_matrix_.size()) + ", got " +
          std::to_string(vector.size()));
    }

    // Apply R^-1 element-wise for diagonal covariance matrix
    // R^-1[i,i] = 1/error[i]^2, so R^-1 * vector[i] = vector[i] / error[i]^2
    std::vector<double> result(vector.size());
    for (size_t i = 0; i < vector.size(); ++i) {
      result[i] = vector[i] * r_inverse_matrix_[i];
    }

    return result;
  }

  std::vector<double> applyCovariance(const std::vector<double>& vector) const {
    // Validate input vector size
    if (vector.size() != r_matrix_.size()) {
      throw std::invalid_argument(
          "Vector size mismatch in applyCovariance: expected " +
          std::to_string(r_matrix_.size()) + ", got " +
          std::to_string(vector.size()));
    }

    // Apply R element-wise for diagonal covariance matrix
    // R[i,i] = error[i]^2, so R * vector[i] = error[i]^2 * vector[i]
    std::vector<double> result(vector.size());
    for (size_t i = 0; i < vector.size(); ++i) {
      result[i] = r_matrix_[i] * vector[i];
    }

    return result;
  }

  std::vector<double> getInverseCovarianceDiagonal() const {
    // Return pre-computed R^-1 matrix diagonal elements
    return r_inverse_matrix_;
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

    // Convert to lowercase for PrepBUFR lookup
    std::transform(base_type.begin(), base_type.end(), base_type.begin(),
                   ::tolower);

    // Get variables from PrepBUFR
    auto variables = obs_->getVariableNames(base_type);

    // If type had a variable suffix, filter to that variable only
    if (!var_suffix.empty()) {
      std::transform(var_suffix.begin(), var_suffix.end(), var_suffix.begin(),
                     ::toupper);
      variables.erase(std::remove_if(variables.begin(), variables.end(),
                                     [&var_suffix](const std::string& v) {
                                       return v != var_suffix;
                                     }),
                      variables.end());
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
   * @brief Apply geometry filtering to observations
   * @param obs_data Observation data to filter
   *
   * @details Filters observations based on geometry domain:
   * - Checks if observation location is within model domain
   * - Removes out-of-domain observations
   * - Updates observation counts accordingly
   */
  void applyGeometryFiltering(ObsOperatorData& obs_data) {
    // TODO: Implement geometry filtering:
    // 1. For each observation, check if location is within geometry domain
    // 2. Use geometry_.isInDomain(lat, lon) or similar method
    // 3. Remove observations outside domain
    // 4. Update observation metadata
    (void)obs_data;  // Suppress unused parameter warning
  }

  // Data members
  std::shared_ptr<const PrepBUFRObservation>
      obs_;  // Observation (owned or aliased)
  std::shared_ptr<const GeometryBackend>
      geometry_;                   // Geometry (owned or aliased)
  bool apply_geometry_filtering_;  // Whether to apply geometry filtering

  // WRFDA data structures (non-owning pointers to Fortran-managed memory)
  void* iv_type_data_ = nullptr;  // Non-owning pointer to Fortran iv_type
  void* y_type_data_ = nullptr;   // Non-owning pointer to Fortran y_type

  // Statistics
  std::unordered_map<std::string, size_t> obs_counts_;   // Counts by type
  std::unordered_map<std::string, double> error_stats_;  // Error stats by type

  // Observation error covariance matrix (diagonal elements)
  std::vector<double> r_matrix_;  // R diagonal elements (variances)
  std::vector<double>
      r_inverse_matrix_;  // R^-1 diagonal elements (1/variances)
};

}  // namespace metada::backends::wrf