/**
 * @file WRFState.hpp
 * @brief WRF state backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <memory>
#include <netcdf>
#include <string>
#include <unordered_map>
#include <vector>
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>

namespace metada::backends::wrf {

/**
 * @brief WRF state backend implementation
 *
 * @details
 * This class implements a state backend for the WRF model. It manages
 * meteorological state variables stored in NetCDF files and provides
 * operations required by the State adapter.
 */
template <typename ConfigBackend>
class WRFState {
 public:
  /**
   * @brief Default constructor is deleted
   */
  WRFState() = delete;

  /**
   * @brief Copy constructor is deleted
   */
  WRFState(const WRFState&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  WRFState& operator=(const WRFState&) = delete;

  /**
   * @brief Constructor that takes a configuration backend
   *
   * @param config Configuration containing WRF file path and variables
   */
  explicit WRFState(const ConfigBackend& config);

  /**
   * @brief Move constructor
   *
   * @param other WRF state backend to move from
   */
  WRFState(WRFState&& other) noexcept;

  /**
   * @brief Move assignment operator
   *
   * @param other WRF state backend to move from
   * @return Reference to this state after assignment
   */
  WRFState& operator=(WRFState&& other) noexcept;

  /**
   * @brief Destructor
   */
  ~WRFState() = default;

  /**
   * @brief Clone this state
   *
   * @return Unique pointer to a new identical WRF state backend
   */
  std::unique_ptr<WRFState> clone() const;

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
  void setActiveVariable(const std::string& name);

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
   * @brief Calculate dot product with another state
   *
   * @param other State to calculate dot product with
   * @return Scalar dot product result
   * @throws std::runtime_error If states are incompatible
   */
  double dot(const WRFState& other) const;

  /**
   * @brief Calculate L2 norm of this state
   *
   * @return L2 norm value
   */
  double norm() const;

  /**
   * @brief Check if this state equals another
   *
   * @param other State to compare with
   * @return True if states are equal, false otherwise
   */
  bool equals(const WRFState& other) const;

  /**
   * @brief Add another state to this one
   *
   * @param other State to add
   * @throws std::runtime_error If states are incompatible
   */
  void add(const WRFState& other);

  /**
   * @brief Subtract another state from this one
   *
   * @param other State to subtract
   * @throws std::runtime_error If states are incompatible
   */
  void subtract(const WRFState& other);

  /**
   * @brief Multiply this state by a scalar
   *
   * @param scalar Value to multiply by
   */
  void multiply(double scalar);

  /**
   * @brief Check if state is properly initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const;

  /**
   * @brief Set data for a specific variable
   *
   * @param values Vector of values to set
   * @param variableName Optional variable name (uses active variable if not
   * specified)
   * @throws std::runtime_error If variable doesn't exist or dimensions don't
   * match
   */
  void setData(const std::vector<double>& values,
               const std::string& variableName = "");

  /**
   * @brief Get data for a specific variable
   *
   * @param variableName Optional variable name (uses active variable if not
   * specified)
   * @return const xt::xarray<double>& Reference to the data array
   * @throws std::runtime_error If variable doesn't exist
   */
  const xt::xarray<double>& getData(const std::string& variableName = "") const;

 private:
  /**
   * @brief Load state data from the WRF NetCDF file
   *
   * @param filename Path to the WRF NetCDF file
   * @param variables List of variables to load
   * @param timestamp Timestamp to read from the file
   */
  void loadStateData(const std::string& filename,
                     const std::vector<std::string>& variables,
                     const std::string& timestamp);

  /**
   * @brief Verify that two states have compatible variables and dimensions
   *
   * @param other State to check compatibility with
   * @return True if compatible, false otherwise
   */
  bool isCompatible(const WRFState& other) const;

  const ConfigBackend& config_;

  // WRF NetCDF file information
  std::string wrfFilename_;
  std::string timestamp_;
  bool initialized_ = false;

  // State data
  std::unordered_map<std::string, xt::xarray<double>> variables_;
  std::unordered_map<std::string, std::vector<size_t>> dimensions_;
  std::vector<std::string> variableNames_;
  std::string activeVariable_;  // Currently active variable for data access
};

}  // namespace metada::backends::wrf