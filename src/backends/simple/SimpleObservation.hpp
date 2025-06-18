/**
 * @file SimpleObservation.hpp
 * @brief Simple implementation of observation backend
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This class provides a simple implementation of the observation backend
 * interface that reads observation data from a text file. The observations are
 * stored in a 2D grid format with missing values represented by -999.0.
 */

#pragma once

#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace metada::backends::simple {

/**
 * @brief Simple observation backend implementation
 *
 * @details
 * This class implements a basic observation backend that:
 * - Reads observation data from multiple text files
 * - Each file contains one observation variable
 * - Stores observations in a 2D grid format
 * - Handles missing values (-999.0)
 * - Provides basic arithmetic operations
 * - Supports quality control operations
 * - Organizes data hierarchically: type → variable → actual data
 */
class SimpleObservation {
 public:
  // Delete default constructor
  SimpleObservation() = delete;

  // Delete copy constructor and assignment
  SimpleObservation(const SimpleObservation&) = delete;
  SimpleObservation& operator=(const SimpleObservation&) = delete;

  /**
   * @brief Constructor that initializes from configuration
   * @tparam ConfigBackend Type of configuration backend
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  explicit SimpleObservation(const ConfigBackend& config) {
    // Get observation files from config
    const auto type_configs = config.Get("types").asVectorMap();

    // Load each observation type
    for (const auto& type_map : type_configs) {
      for (const auto& [type_name, type_config] : type_map) {
        const auto& type_backend = ConfigBackend(type_config.asMap());
        std::string filename = type_backend.Get("file").asString();

        // Get variables for this type
        const auto& variables_config =
            type_backend.Get("variables").asVectorMap();

        for (const auto& var_map : variables_config) {
          for (const auto& [var_name, var_config] : var_map) {
            const auto& var_backend = ConfigBackend(var_config.asMap());
            float error = var_backend.Get("error").asFloat();
            float missing_value = var_backend.Get("missing_value").asFloat();

            // Load data for this type/variable combination
            loadFromFile(type_name, var_name, filename, error, missing_value);
          }
        }
      }
    }
  }

  /**
   * @brief Move constructor
   * @param other Observation to move from
   */
  SimpleObservation(SimpleObservation&& other) noexcept
      : data_(std::move(other.data_)),
        sizes_(std::move(other.sizes_)),
        errors_(std::move(other.errors_)),
        missing_values_(std::move(other.missing_values_)),
        covariance_(std::move(other.covariance_)) {}

  /**
   * @brief Move assignment operator
   * @param other Observation to move from
   * @return Reference to this observation
   */
  SimpleObservation& operator=(SimpleObservation&& other) noexcept {
    if (this != &other) {
      data_ = std::move(other.data_);
      sizes_ = std::move(other.sizes_);
      errors_ = std::move(other.errors_);
      missing_values_ = std::move(other.missing_values_);
      covariance_ = std::move(other.covariance_);
    }
    return *this;
  }

  /**
   * @brief Get raw pointer to observation data
   * @return Pointer to data
   */
  void* getData() {
    if (data_.empty()) {
      return nullptr;
    }
    // Return pointer to the first type's first variable's data
    return data_.begin()->second.begin()->second.data();
  }

  /**
   * @brief Get const raw pointer to observation data
   * @return Const pointer to data
   */
  const void* getData() const {
    if (data_.empty()) {
      return nullptr;
    }
    // Return pointer to the first type's first variable's data
    return data_.begin()->second.begin()->second.data();
  }

  /**
   * @brief Get names of observation types
   * @return Vector of observation type names
   */
  std::vector<std::string> getTypeNames() const {
    std::vector<std::string> types;
    for (const auto& [type_name, _] : data_) {
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
    auto it = data_.find(type_name);
    if (it != data_.end()) {
      for (const auto& [var_name, _] : it->second) {
        variables.push_back(var_name);
      }
    }
    return variables;
  }

  /**
   * @brief Get size of observation data for a specific type/variable
   * @param type_name Name of the observation type
   * @param var_name Name of the variable
   * @return Size of the observation data vector
   */
  size_t getSize(const std::string& type_name,
                 const std::string& var_name) const {
    return sizes_.at(type_name + "." + var_name);
  }

  /**
   * @brief Get error value for a specific type/variable
   * @param type_name Name of the observation type
   * @param var_name Name of the variable
   * @return Error value
   */
  float getError(const std::string& type_name,
                 const std::string& var_name) const {
    return errors_.at(type_name + "." + var_name);
  }

  /**
   * @brief Get missing value for a specific type/variable
   * @param type_name Name of the observation type
   * @param var_name Name of the variable
   * @return Missing value
   */
  float getMissingValue(const std::string& type_name,
                        const std::string& var_name) const {
    return missing_values_.at(type_name + "." + var_name);
  }

  /**
   * @brief Get the observation error covariance matrix
   * @return Const reference to the covariance matrix data
   */
  const std::vector<double>& getCovariance() const { return covariance_; }

  /**
   * @brief Clone this observation
   * @return Unique pointer to cloned observation
   */
  std::unique_ptr<SimpleObservation> clone() const {
    return std::unique_ptr<SimpleObservation>(
        new SimpleObservation(*this, true));
  }

  /**
   * @brief Add another observation to this one
   * @param other Observation to add
   */
  void add(const SimpleObservation& other) {
    for (const auto& [type_name, type_data] : data_) {
      for (const auto& [var_name, var_data] : type_data) {
        std::string key = type_name + "." + var_name;
        if (var_data.size() != other.data_.at(type_name).at(var_name).size()) {
          throw std::runtime_error(
              "Cannot add observations of different sizes");
        }
        for (size_t i = 0; i < var_data.size(); ++i) {
          float missing_val = missing_values_.at(key);
          if (var_data[i] != missing_val &&
              other.data_.at(type_name).at(var_name)[i] != missing_val) {
            data_[type_name][var_name][i] +=
                other.data_.at(type_name).at(var_name)[i];
          }
        }
      }
    }
  }

  /**
   * @brief Subtract another observation from this one
   * @param other Observation to subtract
   */
  void subtract(const SimpleObservation& other) {
    for (const auto& [type_name, type_data] : data_) {
      for (const auto& [var_name, var_data] : type_data) {
        std::string key = type_name + "." + var_name;
        if (var_data.size() != other.data_.at(type_name).at(var_name).size()) {
          throw std::runtime_error(
              "Cannot subtract observations of different sizes");
        }
        for (size_t i = 0; i < var_data.size(); ++i) {
          float missing_val = missing_values_.at(key);
          if (var_data[i] != missing_val &&
              other.data_.at(type_name).at(var_name)[i] != missing_val) {
            data_[type_name][var_name][i] -=
                other.data_.at(type_name).at(var_name)[i];
          }
        }
      }
    }
  }

  /**
   * @brief Multiply this observation by a scalar
   * @param scalar Value to multiply by
   */
  void multiply(double scalar) {
    for (auto& [type_name, type_data] : data_) {
      for (auto& [var_name, var_data] : type_data) {
        std::string key = type_name + "." + var_name;
        float missing_val = missing_values_.at(key);
        for (auto& value : var_data) {
          if (value != missing_val) {
            value *= scalar;
          }
        }
      }
    }
  }

  /**
   * @brief Compare equality with another observation
   * @param other Observation to compare with
   * @return true if equal, false otherwise
   */
  bool equals(const SimpleObservation& other) const {
    return data_ == other.data_ && sizes_ == other.sizes_ &&
           errors_ == other.errors_ &&
           missing_values_ == other.missing_values_ &&
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
    // Simple QC: mark values outside reasonable range as missing
    for (auto& [type_name, type_data] : data_) {
      for (auto& [var_name, var_data] : type_data) {
        std::string key = type_name + "." + var_name;
        float missing_val = missing_values_.at(key);

        for (auto& value : var_data) {
          if (value < -100.0 || value > 100.0) {
            value = missing_val;
          }
        }
      }
    }
  }

  /**
   * @brief Load observation data from file (concept-compliant version)
   * @param filename Path to observation file
   */
  void loadFromFile(const std::string& filename) {
    // For concept compliance, this loads the first type/variable from the file
    // In a real implementation, this would parse the file to determine
    // types/variables
    if (!data_.empty()) {
      auto& first_type = data_.begin()->second;
      if (!first_type.empty()) {
        std::string type_name = data_.begin()->first;
        std::string var_name = first_type.begin()->first;
        std::string key = type_name + "." + var_name;
        float error = errors_.at(key);
        float missing_value = missing_values_.at(key);
        loadFromFile(type_name, var_name, filename, error, missing_value);
      }
    }
  }

  /**
   * @brief Save observation data to file (concept-compliant version)
   * @param filename Path to save observation file
   */
  void saveToFile(const std::string& filename) const {
    // For concept compliance, this saves the first type/variable to the file
    if (!data_.empty()) {
      auto& first_type = data_.begin()->second;
      if (!first_type.empty()) {
        std::string type_name = data_.begin()->first;
        std::string var_name = first_type.begin()->first;
        saveToFile(type_name, var_name, filename);
      }
    }
  }

  /**
   * @brief Load observation data from file
   * @param type_name Name of the observation type
   * @param var_name Name of the variable
   * @param filename Path to observation file
   * @param error Error value for this variable
   * @param missing_value Missing value for this variable
   */
  void loadFromFile(const std::string& type_name, const std::string& var_name,
                    const std::string& filename, float error,
                    float missing_value) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open observation file: " + filename);
    }

    std::vector<double> var_data;
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      double value;
      while (iss >> value) {
        var_data.push_back(value);
      }
    }

    // Set size based on data
    std::string key = type_name + "." + var_name;
    sizes_[key] = var_data.size();
    data_[type_name][var_name] = std::move(var_data);
    errors_[key] = error;
    missing_values_[key] = missing_value;
  }

  /**
   * @brief Save observation data to file
   * @param type_name Name of the observation type
   * @param var_name Name of the variable
   * @param filename Path to save observation file
   */
  void saveToFile(const std::string& type_name, const std::string& var_name,
                  const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open file for writing: " + filename);
    }

    const auto& var_data = data_.at(type_name).at(var_name);
    for (size_t i = 0; i < var_data.size(); i += 36) {
      for (size_t j = 0; j < 36 && (i + j) < var_data.size(); ++j) {
        file << std::fixed << std::setprecision(6) << std::setw(12)
             << var_data[i + j];
      }
      file << "\n";
    }
  }

 private:
  /**
   * @brief Private constructor for cloning
   * @param other Observation to clone from
   * @param is_clone Flag to indicate this is a clone operation
   */
  SimpleObservation(const SimpleObservation& other, bool)
      : data_(other.data_),
        sizes_(other.sizes_),
        errors_(other.errors_),
        missing_values_(other.missing_values_),
        covariance_(other.covariance_) {}

  // Hierarchical data structure: type -> variable -> data
  std::unordered_map<std::string,
                     std::unordered_map<std::string, std::vector<double>>>
      data_;  ///< Observation data organized as type -> variable -> data
  std::unordered_map<std::string, size_t>
      sizes_;  ///< Sizes for each type.variable combination
  std::unordered_map<std::string, float>
      errors_;  ///< Error values for each type.variable combination
  std::unordered_map<std::string, float>
      missing_values_;  ///< Missing values for each type.variable combination
  std::vector<double> covariance_;  ///< Observation error covariance matrix
};

}  // namespace metada::backends::simple