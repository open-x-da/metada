/**
 * @file MACOMState.hpp
 * @brief MACOM state backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <algorithm>  // For std::equal, std::find
#include <cmath>
#include <memory>
#include <numeric>    // For std::inner_product (potentially for dot product)
#include <stdexcept>  // For std::runtime_error
#include <string>
#include <unordered_map>
#include <vector>

namespace metada::backends::macom {

/**
 * @brief MACOM state backend implementation
 *
 * @details
 * This class implements a state backend for the MACOM model. It manages
 * oceanographic state variables and provides operations required by the State
 * adapter.
 */
template <typename ConfigBackend>
class MACOMState {
 public:
  /**
   * @brief Default constructor is deleted
   */
  MACOMState() = delete;

  /**
   * @brief Copy constructor is deleted
   */
  MACOMState(const MACOMState&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  MACOMState& operator=(const MACOMState&) = delete;

  /**
   * @brief Constructor that takes a configuration backend
   *
   * @param config Configuration containing MACOM file path and variables
   */
  explicit MACOMState(const ConfigBackend& config);

  /**
   * @brief Move constructor
   *
   * @param other MACOM state backend to move from
   */
  MACOMState(MACOMState&& other) noexcept;

  /**
   * @brief Move assignment operator
   *
   * @param other MACOM state backend to move from
   * @return Reference to this state after assignment
   */
  MACOMState& operator=(MACOMState&& other) noexcept;

  /**
   * @brief Destructor
   */
  ~MACOMState() = default;

  /**
   * @brief Clone this state
   *
   * @return Unique pointer to a new identical MACOM state backend
   */
  std::unique_ptr<MACOMState> clone() const;

  /**
   * @brief Get mutable access to the underlying data
   *
   * @return Void pointer to the active data array
   */
  void* getData();

  /**
   * @brief Get const access to the underlying data
   *
   * @return Const void pointer to the active data array
   */
  const void* getData() const;

  /**
   * @brief Get the active variable name
   *
   * @return Name of the currently active variable
   */
  const std::string& getActiveVariable() const;

  /**
   * @brief Set the active variable
   *
   * @param name Name of the variable to set as active
   * @throws std::out_of_range If variable doesn't exist
   */
  void setActiveVariable([[maybe_unused]] const std::string& name);

  /**
   * @brief Get all variable names in this state
   *
   * @return Vector of variable names
   */
  const std::vector<std::string>& getVariableNames() const;

  /**
   * @brief Get the dimensions of a variable
   *
   * @param name Name of the variable
   * @return Vector of dimension sizes
   * @throws std::out_of_range If variable doesn't exist
   */
  const std::vector<size_t>& getDimensions(const std::string& name) const;

  /**
   * @brief Set all values to zero
   */
  void zero();

  /**
   * @brief Add another state to this one
   *
   * @param other State to add
   * @throws std::runtime_error If states are incompatible
   */
  void add([[maybe_unused]] const MACOMState& other);

  /**
   * @brief Subtract another state from this one
   *
   * @param other State to subtract
   * @throws std::runtime_error If states are incompatible
   */
  void subtract([[maybe_unused]] const MACOMState& other);

  /**
   * @brief Multiply this state by a scalar
   *
   * @param scalar Value to multiply by
   */
  void multiply([[maybe_unused]] double scalar);

  /**
   * @brief Check if state is properly initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Calculate the norm of this state
   *
   * @return Norm of the state
   */
  double norm() const {
    double sum_sq = 0.0;
    for (const auto& pair : variables_) {
      for (double val : pair.second) {
        sum_sq += val * val;
      }
    }
    return std::sqrt(sum_sq);
  }

  /**
   * @brief 计算与另一个状态的点积
   *
   * @param other 另一个状态
   * @return double 点积结果
   */
  double dot(const MACOMState& other) const {
    if (variableNames_.size() != other.variableNames_.size()) {
      // Or handle more gracefully if some variables might be missing
      throw std::runtime_error(
          "MACOMState::dot: States have different number of variables.");
    }
    double total_dot_product = 0.0;
    for (const auto& name : variableNames_) {
      const auto& vec_this = this->variables_.at(name);
      const auto& vec_other = other.variables_.at(name);
      if (vec_this.size() != vec_other.size()) {
        throw std::runtime_error("MACOMState::dot: Variable '" + name +
                                 "' has different sizes.");
      }
      total_dot_product += std::inner_product(vec_this.begin(), vec_this.end(),
                                              vec_other.begin(), 0.0);
    }
    return total_dot_product;
  }

  /**
   * @brief 检查两个状态是否相等
   *
   * @param other 要比较的另一个状态
   * @return bool 如果相等返回true，否则返回false
   */
  bool equals(const MACOMState& other) const {
    if (initialized_ != other.initialized_ ||
        activeVariable_ != other.activeVariable_ ||
        variableNames_.size() != other.variableNames_.size() ||
        variables_.size() != other.variables_.size() ||
        dimensions_.size() != other.dimensions_.size()) {
      return false;
    }

    // Check variable names
    std::vector<std::string> vn1 = variableNames_;
    std::vector<std::string> vn2 = other.variableNames_;
    std::sort(vn1.begin(), vn1.end());
    std::sort(vn2.begin(), vn2.end());
    if (vn1 != vn2) return false;

    for (const auto& name : variableNames_) {
      if (variables_.at(name).size() != other.variables_.at(name).size() ||
          dimensions_.at(name) != other.dimensions_.at(name)) {
        return false;
      }
      // Approx equals for floating point data
      const auto& data_this = variables_.at(name);
      const auto& data_other = other.variables_.at(name);
      for (size_t i = 0; i < data_this.size(); ++i) {
        if (std::abs(data_this[i] - data_other[i]) >
            1e-9) {  // Tolerance for float comparison
          return false;
        }
      }
    }
    return true;
  }

 private:
  // Reference to configuration
  const ConfigBackend& config_;

  // State information
  bool initialized_ = false;
  std::string inputFile_;                   // Input data file path
  std::vector<std::string> variableNames_;  // Available variables
  std::string activeVariable_;              // Currently active variable

  // Data storage (this would be filled by the actual implementation)
  std::unordered_map<std::string, std::vector<double>> variables_;
  std::unordered_map<std::string, std::vector<size_t>> dimensions_;
};

template <typename ConfigBackend>
MACOMState<ConfigBackend>::MACOMState(const ConfigBackend& config)
    : config_(config), initialized_(false) {}

// Constructor implementation with ConfigBackend
template <typename ConfigBackend>
MACOMState<ConfigBackend>::MACOMState(MACOMState&& other) noexcept
    : config_(other.config_),
      initialized_(other.initialized_),
      inputFile_(std::move(other.inputFile_)),
      variableNames_(std::move(other.variableNames_)),
      activeVariable_(std::move(other.activeVariable_)),
      variables_(std::move(other.variables_)),
      dimensions_(std::move(other.dimensions_)) {
  other.initialized_ = false;
}

template <typename ConfigBackend>
MACOMState<ConfigBackend>& MACOMState<ConfigBackend>::operator=(
    MACOMState&& other) noexcept {
  if (this != &other) {
    initialized_ = other.initialized_;
    inputFile_ = std::move(other.inputFile_);
    variableNames_ = std::move(other.variableNames_);
    activeVariable_ = std::move(other.activeVariable_);
    variables_ = std::move(other.variables_);
    dimensions_ = std::move(other.dimensions_);

    other.initialized_ = false;
  }
  return *this;
}

template <typename ConfigBackend>
std::unique_ptr<MACOMState<ConfigBackend>> MACOMState<ConfigBackend>::clone()
    const {
  auto cloned = std::make_unique<MACOMState>(config_);
  cloned->initialized_ = this->initialized_;
  cloned->inputFile_ = this->inputFile_;
  cloned->variableNames_ = this->variableNames_;
  cloned->activeVariable_ = this->activeVariable_;
  cloned->variables_ = this->variables_;
  cloned->dimensions_ = this->dimensions_;
  return cloned;
}

template <typename ConfigBackend>
void* MACOMState<ConfigBackend>::getData() {
  // 基本实现框架
  return nullptr;
}

template <typename ConfigBackend>
const void* MACOMState<ConfigBackend>::getData() const {
  // 基本实现框架
  return nullptr;
}

template <typename ConfigBackend>
const std::string& MACOMState<ConfigBackend>::getActiveVariable() const {
  return activeVariable_;
}

template <typename ConfigBackend>
void MACOMState<ConfigBackend>::setActiveVariable(
    [[maybe_unused]] const std::string& name) {
  // 实现框架
}

template <typename ConfigBackend>
const std::vector<std::string>& MACOMState<ConfigBackend>::getVariableNames()
    const {
  return variableNames_;
}

template <typename ConfigBackend>
const std::vector<size_t>& MACOMState<ConfigBackend>::getDimensions(
    const std::string& name) const {
  try {
    return dimensions_.at(name);
  } catch (const std::out_of_range&) {
    throw std::out_of_range("Variable not found: " + name);
  }
}

template <typename ConfigBackend>
void MACOMState<ConfigBackend>::zero() {
  // 基本实现框架
}

template <typename ConfigBackend>
void MACOMState<ConfigBackend>::add([[maybe_unused]] const MACOMState& other) {
  // 实现框架
}

template <typename ConfigBackend>
void MACOMState<ConfigBackend>::subtract(
    [[maybe_unused]] const MACOMState& other) {
  // 实现框架
}

template <typename ConfigBackend>
void MACOMState<ConfigBackend>::multiply([[maybe_unused]] double scalar) {
  // 实现框架
}

}  // namespace metada::backends::macom