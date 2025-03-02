#pragma once

#include <map>
#include <string>
#include <vector>

namespace metada::framework {

// Forward declarations
class IConfig;

/**
 * @brief Interface for observation implementations in data assimilation systems
 *
 * @details
 * This interface defines the contract for handling observational data in data
 * assimilation systems. Observations represent measurements of the real system
 * and are a fundamental component in data assimilation algorithms.
 *
 * Unlike the previous design, observations are now copyable to allow for
 * operations like ensemble generation, observation perturbation, and
 * observation-space calculations in data assimilation algorithms.
 *
 * Key features:
 * - Type-safe data access
 * - Comprehensive metadata for spatiotemporal context
 * - Quality control and uncertainty quantification
 * - Arithmetic operations for observation-space calculations
 * - Support for observation operators and observation error covariance
 */
class IObservation {
 public:
  /**
   * @brief Virtual destructor for proper cleanup
   */
  virtual ~IObservation() = default;

  // Lifecycle management
  /**
   * @brief Initialize the observation from configuration
   * @param config Configuration object containing initialization parameters
   * @throws std::runtime_error If initialization fails
   */
  virtual void initialize(const IConfig& config) = 0;

  /**
   * @brief Reset the observation to its initial state
   */
  virtual void reset() = 0;

  /**
   * @brief Validate the observation data
   * @throws std::runtime_error If validation fails
   */
  virtual void validate() const = 0;

  /**
   * @brief Check if the observation is valid
   * @return true if valid, false otherwise
   */
  virtual bool isValid() const = 0;

  /**
   * @brief Check if the observation is initialized
   * @return true if initialized, false otherwise
   */
  virtual bool isInitialized() const = 0;

  // Copy/Move operations
  /**
   * @brief Copy data from another observation
   * @param other Source observation
   */
  virtual void copyFrom(const IObservation& other) = 0;

  /**
   * @brief Move data from another observation
   * @param other Source observation (will be in a valid but unspecified state
   * after the move)
   */
  virtual void moveFrom(IObservation&& other) = 0;

  /**
   * @brief Check if this observation equals another
   * @param other Observation to compare with
   * @return true if equal, false otherwise
   */
  virtual bool equals(const IObservation& other) const = 0;

  // Data access
  /**
   * @brief Get mutable access to observation data
   * @return Void pointer to the underlying data
   */
  virtual void* getData() = 0;

  /**
   * @brief Get immutable access to observation data
   * @return Const void pointer to the underlying data
   */
  virtual const void* getData() const = 0;

  /**
   * @brief Get mutable access to observation error/uncertainty
   * @return Void pointer to the uncertainty information
   */
  virtual void* getUncertainty() = 0;

  /**
   * @brief Get immutable access to observation error/uncertainty
   * @return Const void pointer to the uncertainty information
   */
  virtual const void* getUncertainty() const = 0;

  /**
   * @brief Get the size of the observation vector
   * @return Number of observation elements
   */
  virtual size_t getSize() const = 0;

  // Metadata operations
  /**
   * @brief Set a metadata key-value pair
   * @param key Metadata key
   * @param value Metadata value
   */
  virtual void setMetadata(const std::string& key,
                           const std::string& value) = 0;

  /**
   * @brief Get a metadata value
   * @param key Metadata key
   * @return Metadata value
   * @throws std::runtime_error If key doesn't exist
   */
  virtual std::string getMetadata(const std::string& key) const = 0;

  /**
   * @brief Check if a metadata key exists
   * @param key Metadata key
   * @return true if key exists, false otherwise
   */
  virtual bool hasMetadata(const std::string& key) const = 0;

  // Observation information
  /**
   * @brief Get the names of observation variables
   * @return Vector of variable names
   */
  virtual const std::vector<std::string>& getVariableNames() const = 0;

  /**
   * @brief Check if a variable exists in the observation
   * @param name Variable name
   * @return true if variable exists, false otherwise
   */
  virtual bool hasVariable(const std::string& name) const = 0;

  /**
   * @brief Get the dimensions of the observation space
   * @return Vector of dimensions
   */
  virtual const std::vector<size_t>& getDimensions() const = 0;

  // Spatiotemporal metadata
  /**
   * @brief Set the geospatial locations of the observations
   * @param locations Vector of [lat, lon, elevation] triplets for each
   * observation
   */
  virtual void setLocations(
      const std::vector<std::vector<double>>& locations) = 0;

  /**
   * @brief Set the timestamps of the observations
   * @param timestamps Vector of timestamps for each observation
   */
  virtual void setTimes(const std::vector<double>& timestamps) = 0;

  /**
   * @brief Get the geospatial locations of the observations
   * @return Vector of [lat, lon, elevation] triplets for each observation
   */
  virtual const std::vector<std::vector<double>>& getLocations() const = 0;

  /**
   * @brief Get the timestamps of the observations
   * @return Vector of timestamps for each observation
   */
  virtual const std::vector<double>& getTimes() const = 0;

  // Arithmetic operations
  /**
   * @brief Add another observation to this one
   * @param other Observation to add
   */
  virtual void add(const IObservation& other) = 0;

  /**
   * @brief Subtract another observation from this one
   * @param other Observation to subtract
   */
  virtual void subtract(const IObservation& other) = 0;

  /**
   * @brief Multiply this observation by a scalar
   * @param scalar Scalar value
   */
  virtual void multiply(double scalar) = 0;

  // Quality control
  /**
   * @brief Set the quality control flags
   * @param flags Vector of quality flags for each observation
   */
  virtual void setQualityFlags(const std::vector<int>& flags) = 0;

  /**
   * @brief Get the quality control flags
   * @return Vector of quality flags for each observation
   */
  virtual const std::vector<int>& getQualityFlags() const = 0;

  /**
   * @brief Set the confidence values
   * @param values Vector of confidence values for each observation
   */
  virtual void setConfidenceValues(const std::vector<double>& values) = 0;

  /**
   * @brief Get the confidence values
   * @return Vector of confidence values for each observation
   */
  virtual const std::vector<double>& getConfidenceValues() const = 0;
};

}  // namespace metada::framework