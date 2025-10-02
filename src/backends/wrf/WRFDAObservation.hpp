#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "Location.hpp"
#include "ObsRecord.hpp"
#include "ObservationIterator.hpp"
#include "PrepBUFRObservation.hpp"
#include "WRFGeometry.hpp"

namespace metada::backends::wrf {
using ObsRecord = framework::ObsRecord;
using ObsLevelRecord = framework::ObsLevelRecord;
using Location = framework::Location;
using PrepBUFRObservation = backends::common::observation::PrepBUFRObservation;

// Simple iterator implementation
template <typename GeometryBackend>
class WRFDAObservationIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = framework::ObsRecord;
  using difference_type = std::ptrdiff_t;
  using pointer = const value_type*;
  using reference = const value_type&;

  WRFDAObservationIterator() : index_(0) {}
  WRFDAObservationIterator(size_t index) : index_(index) {}

  bool operator==(const WRFDAObservationIterator& other) const {
    return index_ == other.index_;
  }
  bool operator!=(const WRFDAObservationIterator& other) const {
    return !(*this == other);
  }

  WRFDAObservationIterator& operator++() {
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
 * @brief Adapter class to convert METADA observations to WRFDA iv_type/y_type
 *
 * This class handles the conversion between METADA's observation
 * representation and WRFDA's internal iv_type/y_type structures. It provides
 * methods to:
 * 1. Convert PrepBUFR observations to WRFDA format
 * 2. Handle different observation types (SYNOP, SOUND, etc.)
 * 3. Manage observation errors and QC flags
 * 4. Support WRF-specific coordinate transformations (Arakawa-C grid)
 */
template <typename GeometryBackend>
class WRFDAObservation {
 public:
  // Required type definitions for ObservationBackendImpl concept
  using iterator_type = WRFDAObservationIterator<GeometryBackend>;
  using value_type = ObsRecord;
  using ObsOperatorData = PrepBUFRObservation::ObsOperatorData;

  // Delete default constructor and copying
  WRFDAObservation() = delete;
  WRFDAObservation(const WRFDAObservation&) = delete;
  WRFDAObservation& operator=(const WRFDAObservation&) = delete;

  /**
   * @brief Construct from config only (no geometry filtering)
   * @param config Configuration object containing observation settings
   *
   * @details This constructor creates its own geometry from config.
   * No geometry filtering is applied - all observations are included.
   */
  explicit WRFDAObservation(const backends::config::YamlConfig& config)
      : obs_ptr_(std::make_shared<PrepBUFRObservation>(config)),
        geometry_ptr_(std::make_shared<GeometryBackend>(config)),
        obs_(*obs_ptr_),
        geometry_(*geometry_ptr_),
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
  WRFDAObservation(const backends::config::YamlConfig& config,
                   const GeometryBackend& geometry)
      : obs_ptr_(std::make_shared<PrepBUFRObservation>(config)),
        geometry_ptr_(nullptr),
        obs_(*obs_ptr_),
        geometry_(geometry),
        apply_geometry_filtering_(true) {
    // Initialize WRFDA data structures with geometry filtering
    initialize();
  }

  /**
   * @brief Construct from PrepBUFR observation and WRF geometry
   * @param obs Source PrepBUFR observation
   * @param geometry WRF geometry for domain information and Arakawa-C grid
   * transforms
   *
   * @details This constructor is for direct construction with existing objects.
   * Geometry filtering is applied.
   */
  WRFDAObservation(const PrepBUFRObservation& obs,
                   const GeometryBackend& geometry)
      : obs_ptr_(nullptr),
        geometry_ptr_(nullptr),
        obs_(obs),
        geometry_(geometry),
        apply_geometry_filtering_(true) {
    // Initialize WRFDA data structures
    initialize();
  }

  // Move semantics
  WRFDAObservation(WRFDAObservation&& other) noexcept = default;
  WRFDAObservation& operator=(WRFDAObservation&& other) noexcept = default;

  /**
   * @brief Initialize observation data structures
   */
  void initialize() {
    // Extract observation data organized by WRFDA types
    auto obs_data = obs_.getObsOperatorData();

    // Apply geometry filtering if requested
    if (apply_geometry_filtering_) {
      applyGeometryFiltering(obs_data);
    }

    // Allocate and populate iv_type structure
    allocateIVType(obs_data);

    // Allocate and populate y_type structure
    allocateYType(obs_data);

    // Compute observation counts and error statistics
    computeStats(obs_data);
  }

  /**
   * @brief Get pointer to WRFDA iv_type structure
   * @return Raw pointer to iv_type data
   */
  void* getIVTypeData() { return iv_type_data_.get(); }

  /**
   * @brief Get pointer to WRFDA y_type structure
   * @return Raw pointer to y_type data
   */
  void* getYTypeData() { return y_type_data_.get(); }

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
  void* getData() { return iv_type_data_.get(); }

  /**
   * @brief Get raw data pointer (const)
   * @return Const void pointer to underlying data
   */
  const void* getData() const { return iv_type_data_.get(); }

  /**
   * @brief Get typed data access
   * @tparam T Expected data type
   * @return Data cast to requested type
   */
  template <typename T>
  T getData() const {
    // TODO: Implement type-safe data access
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Get number of observations by type
   * @param obs_type Observation type (e.g., "synop", "sound")
   * @return Number of observations of that type
   */
  size_t getObservationCount(const std::string& obs_type) const {
    auto it = obs_counts_.find(obs_type);
    return it != obs_counts_.end() ? it->second : 0;
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
   * @return Vector containing covariance matrix elements
   */
  std::vector<double> getCovariance() const {
    // TODO: Implement covariance matrix access from WRFDA structures
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Get size of observations for a specific type and variable
   * @param typeName Observation type name
   * @param varName Variable name
   * @return Number of observations
   */
  size_t getSize(const std::string& typeName,
                 const std::string& varName) const {
    // TODO: Implement size calculation for specific type/variable
    return getTotalObservationCount();
  }

  /**
   * @brief Get observation operator data
   * @return ObsOperatorData structure
   */
  ObsOperatorData getObsOperatorData() const {
    return obs_.getObsOperatorData();
  }

  // Required by ObservationBackendImpl concept
  double quadraticForm(const std::vector<double>& innovation) const {
    // TODO: Implement R^-1 quadratic form
    throw std::runtime_error("Not implemented");
  }

  std::vector<double> applyInverseCovariance(
      const std::vector<double>& vector) const {
    // TODO: Implement R^-1 application
    throw std::runtime_error("Not implemented");
  }

  std::vector<double> applyCovariance(const std::vector<double>& vector) const {
    // TODO: Implement R application
    throw std::runtime_error("Not implemented");
  }

  std::vector<double> getInverseCovarianceDiagonal() const {
    // TODO: Implement R^-1 diagonal access
    throw std::runtime_error("Not implemented");
  }

  bool isDiagonalCovariance() const {
    return true;  // WRFDA uses diagonal R by default
  }

  /**
   * @brief Add another observation to this one
   * @param other Observation to add
   */
  void add(const WRFDAObservation& other) {
    // TODO: Implement observation addition
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Subtract another observation from this one
   * @param other Observation to subtract
   */
  void subtract(const WRFDAObservation& other) {
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
  bool equals(const WRFDAObservation& other) const {
    // TODO: Implement equality comparison
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Get observation type names
   * @return Vector of observation type names
   */
  std::vector<std::string> getTypeNames() const {
    std::vector<std::string> types;
    for (const auto& [type, count] : obs_counts_) {
      types.push_back(type);
    }
    return types;
  }

  /**
   * @brief Get variable names for observation type
   * @param typeName Observation type name
   * @return Vector of variable names
   */
  std::vector<std::string> getVariableNames(const std::string& typeName) const {
    // TODO: Implement variable name retrieval from WRFDA structures
    throw std::runtime_error("Not implemented");
  }

  /**
   * @brief Create a clone of this observation
   * @return Unique pointer to cloned observation
   */
  std::unique_ptr<WRFDAObservation> clone() const {
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
   * @brief Allocate and populate WRFDA iv_type structure
   * @param obs_data Organized observation data
   */
  void allocateIVType(const ObsOperatorData& obs_data) {
    // TODO: Implement iv_type allocation and population following WRFDA
    // conventions:
    // - Organize by observation type (synop, sound, etc.)
    // - Handle multiple variables per station
    // - Set proper QC flags and metadata
    // - Convert units to WRFDA convention
  }

  /**
   * @brief Allocate and populate WRFDA y_type structure
   * @param obs_data Organized observation data
   */
  void allocateYType(const ObsOperatorData& obs_data) {
    // TODO: Implement y_type allocation and population:
    // - Match iv_type organization
    // - Initialize innovation vectors
    // - Set proper metadata
  }

  /**
   * @brief Compute observation statistics
   * @param obs_data Organized observation data
   */
  void computeStats(const ObsOperatorData& obs_data) {
    // TODO: Implement statistics computation:
    // - Count by observation type
    // - Error statistics by variable
    // - Domain coverage statistics
    // - QC flag distribution
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
  }

  // Data members
  std::shared_ptr<PrepBUFRObservation>
      obs_ptr_;  // Owned observation (when constructed from config)
  std::shared_ptr<GeometryBackend>
      geometry_ptr_;  // Owned geometry (when constructed from config)
  const PrepBUFRObservation& obs_;   // Reference to source observation
  const GeometryBackend& geometry_;  // Reference to WRF geometry
  bool apply_geometry_filtering_;    // Whether to apply geometry filtering

  // WRFDA data structures
  std::unique_ptr<char[]> iv_type_data_;  // Raw iv_type data
  std::unique_ptr<char[]> y_type_data_;   // Raw y_type data

  // Statistics
  std::unordered_map<std::string, size_t> obs_counts_;   // Counts by type
  std::unordered_map<std::string, double> error_stats_;  // Error stats by type
};

}  // namespace metada::backends::wrf