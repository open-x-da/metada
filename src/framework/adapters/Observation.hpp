#pragma once

#include <iomanip>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "Logger.hpp"
#include "NonCopyable.hpp"
#include "ObservationConcepts.hpp"
#include "ObservationIterator.hpp"

namespace metada::framework {

/**
 * @brief Forward declaration of Config class
 *
 * @tparam BackendTag The backend tag type that must satisfy ConfigBackendType
 * concept
 */
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

/**
 * @brief Adapter class for observation implementations in data assimilation
 * systems
 *
 * @details This class provides a type-safe interface for handling observational
 * data in data assimilation workflows. It wraps a backend implementation that
 * satisfies the ObservationBackendType concept, providing consistent access
 * patterns and operations across different backend implementations.
 *
 * The Observation class inherits from NonCopyable to prevent unintended
 * copying, but provides explicit clone functionality for when copies are needed
 * in operations like ensemble generation, observation perturbation, and
 * observation-space calculations in data assimilation algorithms.
 *
 * Features:
 * - Type-safe data access through templates
 * - Delegation to backend implementation
 * - Comprehensive error handling
 * - Arithmetic operations for observation-space calculations
 * - Explicit cloning mechanism for controlled copying
 * - Move semantics for efficient resource management
 * - Hierarchical data organization: type → variable → data
 *
 * @tparam BackendTag The backend tag type that must satisfy
 * ObservationBackendType concept
 *
 * @see ObservationBackendType
 * @see NonCopyable
 */
template <typename BackendTag>
  requires ObservationBackendType<BackendTag>
class Observation : private NonCopyable {
 public:
  /** @brief Type alias for the backend implementation type */
  using ObservationBackend =
      typename traits::BackendTraits<BackendTag>::ObservationBackend;

  /** @brief Default constructor is deleted since we need a backend */
  Observation() = delete;

  /** @brief Destructor */
  ~Observation() = default;

  /**
   * @brief Constructor that initializes observation with configuration only
   *
   * @details This constructor creates an observation without geometry
   * filtering. The backend will create its own geometry instance if needed.
   *
   * @param config Configuration object containing initialization parameters
   */
  explicit Observation(const Config<BackendTag>& config)
      : backend_(config.backend()), initialized_(true) {
    logger_.Info() << "Observation constructed (no geometry filtering)";
  }

  /**
   * @brief Constructor that initializes observation with configuration and
   * geometry
   *
   * @details This constructor creates an observation with geometry filtering.
   * Only observations within the geometry domain will be included.
   *
   * @param config Configuration object containing initialization parameters
   * @param geometry Geometry object for domain filtering
   */
  template <typename GeometryType>
  Observation(const Config<BackendTag>& config, const GeometryType& geometry)
      : backend_(config.backend(), geometry.backend()), initialized_(true) {
    logger_.Info() << "Observation constructed (with geometry filtering)";
  }

  /**
   * @brief Move constructor
   *
   * @details Transfers ownership of the backend from another observation
   * and updates initialization state.
   *
   * @param other The observation to move from
   */
  Observation(Observation&& other) noexcept
      : backend_(std::move(other.backend_)), initialized_(other.initialized_) {
    other.initialized_ = false;
    logger_.Info() << "Observation moved";
  }

  /**
   * @brief Move assignment operator
   *
   * @details Transfers ownership of the backend from another observation
   * and updates initialization state.
   *
   * @param other The observation to move from
   * @return Reference to this observation after assignment
   */
  Observation& operator=(Observation&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
      initialized_ = other.initialized_;
      other.initialized_ = false;
    }
    return *this;
  }

  /**
   * @brief Clone the observation
   *
   * @details Creates a deep copy of this observation, including all data and
   * metadata. This is the preferred way to create copies of observations when
   * needed.
   *
   * @return A new observation with the same data and configuration
   */
  Observation clone() const {
    return Observation(std::move(*backend_.clone()));
  }

  /**
   * @brief Check if the observation is properly initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Apply quality control procedures to the observation
   *
   * @details Delegates to the backend implementation to perform quality control
   * checks and filtering on the observation data.
   */
  void applyQC() { backend_.applyQC(); }

  /**
   * @brief Load observation data from a file
   *
   * @details Delegates to the backend implementation to read observation data
   * from the specified file and initialize this observation object.
   *
   * @param filename Path to the file containing observation data
   * @param error Error value to use for missing data
   * @param missing_value Missing value to use for missing data
   */
  void loadFromFile(const std::string& filename, double error,
                    double missing_value) {
    backend_.loadFromFile(filename, error, missing_value);
    initialized_ = true;
  }

  /**
   * @brief Save observation data to a file
   *
   * @details Delegates to the backend implementation to write the current
   * observation data to the specified file.
   *
   * @param filename Path where observation data will be saved
   */
  void saveToFile(const std::string& filename) const {
    backend_.saveToFile(filename);
  }

  /**
   * @brief Get information about applied filters
   *
   * @details Returns a summary of what filtering was applied during observation
   * loading, including data types, variables, and geographic bounds.
   * Delegates to the concrete backend implementation when available.
   *
   * @return String containing filtering information
   */
  std::string getFilteringInfo() const {
    // Use constexpr if to check if the backend has getFilteringInfo method
    if constexpr (requires { backend_.getFilteringInfo(); }) {
      return backend_.getFilteringInfo();
    } else {
      return "Filtering information not available from this backend";
    }
  }

  /**
   * @brief Equality operator
   *
   * @details Compares this observation with another for equality by delegating
   * to the backend implementation and checking initialization status.
   *
   * @param other The observation to compare with
   * @return True if the observations are equal, false otherwise
   */
  bool operator==(const Observation& other) const {
    return backend_.equals(other.backend_) &&
           initialized_ == other.initialized_;
  }

  /**
   * @brief Inequality operator
   *
   * @details Compares this observation with another for inequality by negating
   * the result of the equality operator.
   *
   * @param other The observation to compare with
   * @return True if the observations are not equal, false otherwise
   */
  bool operator!=(const Observation& other) const { return !(*this == other); }

  /**
   * @brief Get typed access to the underlying data
   *
   * @details Provides access to the observation data cast to the specified
   * type. The caller is responsible for ensuring the type is correct.
   *
   * @tparam T The expected data type
   * @return The data of the specified type
   */
  template <typename T>
  T getData() {
    return backend_.template getData<T>();
  }

  /**
   * @brief Get typed access to the underlying data (const version)
   *
   * @details Provides const access to the observation data cast to the
   * specified type. The caller is responsible for ensuring the type is correct.
   *
   * @tparam T The expected data type
   * @return The data of the specified type
   */
  template <typename T>
  T getData() const {
    return backend_.template getData<T>();
  }

  /**
   * @brief Get the observation error variances
   *
   * @details Returns the observation error variances as a vector of doubles.
   * Each element corresponds to the variance of one observation.
   *
   * @return Vector of observation error variances
   */
  std::vector<double> getCovariance() const { return backend_.getCovariance(); }

  /**
   * @brief Get the names of observation types
   *
   * @details Retrieves the list of observation type names contained in this
   * observation.
   *
   * @return Vector of observation type names
   */
  std::vector<std::string> getTypeNames() const {
    return backend_.getTypeNames();
  }

  /**
   * @brief Get the names of variables for a specific observation type
   *
   * @details Retrieves the list of variable names for the specified observation
   * type.
   *
   * @param typeName Name of the observation type
   * @return Vector of variable names for the type
   */
  std::vector<std::string> getVariableNames(const std::string& typeName) const {
    return backend_.getVariableNames(typeName);
  }

  /**
   * @brief Check if a specific observation type exists
   *
   * @details Searches for the specified observation type name in the list of
   * types contained in this observation.
   *
   * @param typeName Name of the observation type to check
   * @return True if the observation type exists, false otherwise
   */
  bool hasType(const std::string& typeName) const {
    const auto& types = getTypeNames();
    return std::find(types.begin(), types.end(), typeName) != types.end();
  }

  /**
   * @brief Check if a specific variable exists in a specific observation type
   *
   * @details Searches for the specified variable name in the list of variables
   * for the specified observation type.
   *
   * @param typeName Name of the observation type
   * @param varName Name of the variable to check
   * @return True if the variable exists in the type, false otherwise
   */
  bool hasVariable(const std::string& typeName,
                   const std::string& varName) const {
    const auto& variables = getVariableNames(typeName);
    return std::find(variables.begin(), variables.end(), varName) !=
           variables.end();
  }

  /**
   * @brief Get the size of observation data for a specific type/variable
   *
   * @details Retrieves the size of the observation data vector for the
   * specified type and variable.
   *
   * @param typeName Name of the observation type
   * @param varName Name of the variable
   * @return Size of the observation data vector
   */
  size_t getSize(const std::string& typeName,
                 const std::string& varName) const {
    return backend_.getSize(typeName, varName);
  }

  /**
   * @brief Addition operator
   *
   * @details Creates a new observation that is the sum of this observation and
   * another.
   *
   * @param other The observation to add to this one
   * @return A new observation containing the result of the addition
   */
  Observation operator+(const Observation& other) const {
    Observation result(std::move(*backend_.clone()));
    result.backend_.add(other.backend_);
    return result;
  }

  /**
   * @brief Compound addition operator
   *
   * @details Adds another observation to this one in-place.
   *
   * @param other The observation to add to this one
   * @return Reference to this observation after addition
   */
  Observation& operator+=(const Observation& other) {
    backend_.add(other.backend_);
    return *this;
  }

  /**
   * @brief Subtraction operator
   *
   * @details Creates a new observation that is the difference between this
   * observation and another.
   *
   * @param other The observation to subtract from this one
   * @return A new observation containing the result of the subtraction
   */
  Observation operator-(const Observation& other) const {
    Observation result(std::move(*backend_.clone()));
    result.backend_.subtract(other.backend_);
    return result;
  }

  /**
   * @brief Compound subtraction operator
   *
   * @details Subtracts another observation from this one in-place.
   *
   * @param other The observation to subtract from this one
   * @return Reference to this observation after subtraction
   */
  Observation& operator-=(const Observation& other) {
    backend_.subtract(other.backend_);
    return *this;
  }

  /**
   * @brief Multiplication operator
   *
   * @details Creates a new observation that is this observation multiplied by a
   * scalar.
   *
   * @param scalar The scalar value to multiply this observation by
   * @return A new observation containing the result of the multiplication
   */
  Observation operator*(double scalar) const {
    Observation result(std::move(*backend_.clone()));
    result.backend_.multiply(scalar);
    return result;
  }

  /**
   * @brief Compound multiplication operator
   *
   * @details Multiplies this observation by a scalar in-place.
   *
   * @param scalar The scalar value to multiply this observation by
   * @return Reference to this observation after multiplication
   */
  Observation& operator*=(double scalar) {
    backend_.multiply(scalar);
    return *this;
  }

  /**
   * @brief Multiplication operator with scalar
   *
   * @details Creates a new observation that is the result of multiplying
   * a scalar by an observation. This enables the commutative property
   * of scalar multiplication.
   *
   * @param scalar The scalar value to multiply the observation by
   * @param obs The observation to multiply
   * @return A new observation containing the result of the multiplication
   */
  friend Observation operator*(double scalar, const Observation& obs) {
    return obs * scalar;
  }

  /**
   * @brief Access to backend instance
   *
   * @details Provides direct access to the underlying backend implementation.
   * This should be used with caution as it bypasses the adapter's interface.
   *
   * @return Reference to the backend instance
   */
  ObservationBackend& backend() { return backend_; }

  /**
   * @brief Const access to backend instance
   *
   * @details Provides read-only access to the underlying backend
   * implementation. This should be used with caution as it bypasses the
   * adapter's interface.
   *
   * @return Const reference to the backend instance
   */
  const ObservationBackend& backend() const { return backend_; }

  /**
   * @brief Get iterator to beginning of observations
   * @return Iterator to first observation
   */
  auto begin() const { return backend_.begin(); }

  /**
   * @brief Get iterator to end of observations
   * @return Iterator past last observation
   */
  auto end() const { return backend_.end(); }

  /**
   * @brief Get number of observations
   * @return Total number of observations
   */
  size_t size() const { return backend_.size(); }

  /**
   * @brief Get observation at specific index
   * @param index Index of the observation
   * @return Reference to the observation
   */
  auto operator[](size_t index) const { return backend_[index]; }

  /**
   * @brief Get observations within a geographic bounding box
   * @param min_lat Minimum latitude
   * @param max_lat Maximum latitude
   * @param min_lon Minimum longitude
   * @param max_lon Maximum longitude
   * @return Vector of observations within the bounding box
   */
  auto getObservationsInBox(double min_lat, double max_lat, double min_lon,
                            double max_lon) const {
    return backend_.getObservationsInBox(min_lat, max_lat, min_lon, max_lon);
  }

  /**
   * @brief Get observations within a vertical range
   * @param min_level Minimum vertical level
   * @param max_level Maximum vertical level
   * @return Vector of observations within the vertical range
   */
  auto getObservationsInVerticalRange(double min_level,
                                      double max_level) const {
    return backend_.getObservationsInVerticalRange(min_level, max_level);
  }

  /**
   * @brief Compute quadratic form with observation error covariance: v^T R^-1 v
   *
   * @details Computes the quadratic form of a vector with the inverse
   * observation error covariance matrix. This is used in variational data
   * assimilation to compute the observation term of the cost function.
   *
   * @param innovation Innovation vector (y - H(x))
   * @return Quadratic form value
   */
  double quadraticForm(const std::vector<double>& innovation) const {
    return backend_.quadraticForm(innovation);
  }

  /**
   * @brief Apply inverse observation error covariance: R^-1 v
   *
   * @details Applies the inverse observation error covariance matrix to a
   * vector. This is used in variational data assimilation gradient
   * computations.
   *
   * @param innovation Innovation vector
   * @return Result of R^-1 applied to the innovation
   */
  std::vector<double> applyInverseCovariance(
      const std::vector<double>& innovation) const {
    return backend_.applyInverseCovariance(innovation);
  }

  /**
   * @brief Apply observation error covariance: R v
   *
   * @details Applies the observation error covariance matrix to a vector.
   * This can be used for preconditioning or error propagation.
   *
   * @param vector Input vector
   * @return Result of R applied to the vector
   */
  std::vector<double> applyCovariance(const std::vector<double>& vector) const {
    return backend_.applyCovariance(vector);
  }

  /**
   * @brief Get inverse observation error covariance matrix diagonal
   *
   * @details For diagonal observation error covariance matrices, returns the
   * diagonal elements of R^-1.
   *
   * @return Vector containing diagonal elements of R^-1
   */
  std::vector<double> getInverseCovarianceDiagonal() const {
    return backend_.getInverseCovarianceDiagonal();
  }

  /**
   * @brief Check if observation error covariance is diagonal
   *
   * @return True if R is diagonal, false for full covariance matrix
   */
  bool isDiagonalCovariance() const { return backend_.isDiagonalCovariance(); }

  /**
   * @brief Stream insertion operator for observation summary
   *
   * @details Outputs a comprehensive summary of the observation data including:
   * - Total number of observations
   * - Summary for each observation type
   * - Number of variables per type
   * - Total observations per type
   * - Coordinate system information
   *
   * @param os Output stream to write to
   * @param obs Observation object to summarize
   * @return Reference to the output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const Observation& obs) {
    os << "=== Observation Summary ===\n";
    os << "Total observations: " << obs.size() << "\n";
    os << "Initialized: " << (obs.isInitialized() ? "Yes" : "No") << "\n";
    os << "Coordinate system: Geographic\n\n";

    const auto& type_names = obs.getTypeNames();
    if (type_names.empty()) {
      os << "No observation types defined\n";
    } else {
      os << "Observation Types (" << type_names.size() << "):\n";
      os << std::string(50, '-') << "\n";

      for (const auto& type_name : type_names) {
        const auto& var_names = obs.getVariableNames(type_name);
        size_t total_obs_for_type = 0;

        // Count total observations for this type across all variables
        for (const auto& var_name : var_names) {
          total_obs_for_type += obs.getSize(type_name, var_name);
        }

        os << "Type: " << std::setw(15) << std::left << type_name;
        os << " | Variables: " << std::setw(3) << std::right
           << var_names.size();
        os << " | Total obs: " << std::setw(6) << std::right
           << total_obs_for_type << "\n";

        // Show variable details if there are any
        if (!var_names.empty()) {
          os << "  Variables: ";
          for (size_t i = 0; i < var_names.size(); ++i) {
            if (i > 0) os << ", ";
            os << var_names[i] << "(" << obs.getSize(type_name, var_names[i])
               << ")";
          }
          os << "\n";
        }
        os << "\n";
      }
    }

    // Show covariance information
    os << "Covariance Information:\n";
    os << std::string(50, '-') << "\n";
    os << "Diagonal: " << (obs.isDiagonalCovariance() ? "Yes" : "No") << "\n";

    const auto& cov = obs.getCovariance();
    if (!cov.empty()) {
      double min_error = *std::min_element(cov.begin(), cov.end());
      double max_error = *std::max_element(cov.begin(), cov.end());
      os << "Error variance range: [" << min_error << ", " << max_error
         << "]\n";
    }

    return os;
  }

 private:
  /**
   * @brief Private constructor to be used internally by clone
   *
   * @details Creates an observation by taking ownership of an existing backend
   * instance. This is primarily used by the clone method.
   *
   * @param backend Backend instance to move from
   */
  explicit Observation(ObservationBackend&& backend)
      : backend_(std::move(backend)), initialized_(true) {}

  ObservationBackend backend_;  ///< Backend implementation instance
  bool initialized_{false};     ///< Initialization flag
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework