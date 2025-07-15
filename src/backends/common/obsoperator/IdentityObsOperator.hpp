/**
 * @file IdentityObsOperator.hpp
 * @brief Identity observation operator backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This class provides an identity observation operator backend interface that maps 
 * any state backend to any observation backend. The observation operator performs nearest-neighbor 
 * interpolation from the model grid to observation locations.
 */

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "Location.hpp"

namespace metada::backends::common::obsoperator {

using framework::CoordinateSystem;

/**
 * @brief Identity observation operator backend implementation
 *
 * @details
 * This class implements an identity observation operator that:
 * - Maps state variables from model grid to observation locations
 * - Uses nearest-neighbor interpolation for grid-to-point mapping
 * - Supports multiple observation types and variables
 * - Handles missing values and quality control
 * - Provides configuration-based initialization
 * - Works with any state and observation backends that provide required interfaces
 *
 * @tparam StateBackend Type of state backend to operate on
 * @tparam ObsBackend Type of observation backend to operate on
 */
template <typename StateBackend, typename ObsBackend>
class IdentityObsOperator {
 public:
  // Delete default constructor
  IdentityObsOperator() = delete;

  // Delete copy constructor and assignment
  IdentityObsOperator(const IdentityObsOperator&) = delete;
  IdentityObsOperator& operator=(const IdentityObsOperator&) = delete;

  /**
   * @brief Constructor that initializes from configuration
   * @tparam ConfigBackend Type of configuration backend
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  explicit IdentityObsOperator(const ConfigBackend& config) {
    initialize(config);
  }

  /**
   * @brief Move constructor
   * @param other Observation operator to move from
   */
  IdentityObsOperator(IdentityObsOperator&& other) noexcept
      : initialized_(other.initialized_),
        required_state_vars_(std::move(other.required_state_vars_)),
        required_obs_vars_(std::move(other.required_obs_vars_)) {
    other.initialized_ = false;
  }

  /**
   * @brief Move assignment operator
   * @param other Observation operator to move from
   * @return Reference to this observation operator
   */
  IdentityObsOperator& operator=(IdentityObsOperator&& other) noexcept {
    if (this != &other) {
      initialized_ = other.initialized_;
      required_state_vars_ = std::move(other.required_state_vars_);
      required_obs_vars_ = std::move(other.required_obs_vars_);
      other.initialized_ = false;
    }
    return *this;
  }

  /**
   * @brief Initialize with configuration
   * @tparam ConfigBackend Type of configuration backend
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  void initialize(const ConfigBackend& config) {
    if (isInitialized()) {
      throw std::runtime_error("IdentityObsOperator already initialized");
    }

    // Get required state variables from config
    try {
      required_state_vars_ = config.Get("required_state_vars").asVectorString();
    } catch (...) {
      // Default to "state" if not specified
      required_state_vars_ = {"state"};
    }

    // Get required observation variables from config
    try {
      required_obs_vars_ = config.Get("required_obs_vars").asVectorString();
    } catch (...) {
      // Default to empty if not specified
      required_obs_vars_ = {};
    }

    initialized_ = true;
  }

  /**
   * @brief Check if the observation operator is initialized
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Apply forward observation operator: H(x)
   *
   * @details Maps state to observation space by interpolating from model grid
   * to observation locations using nearest-neighbor interpolation.
   * For multi-variable states, interpolates from the appropriate variable
   * based on observation type and required state variables.
   *
   * @param state Model state to transform
   * @param obs Output observation to store the result
   * @return Vector of interpolated values at observation locations
   */
  std::vector<double> apply(const StateBackend& state,
                            const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("IdentityObsOperator not initialized");
    }

    std::vector<double> result;
    result.reserve(obs.size());

    // Apply observation operator to each observation point
    for (const auto& obs_point : obs) {
      if (!obs_point.is_valid) {
        result.push_back(0.0);  // Invalid observations get zero
        continue;
      }

      // Determine which state variable to interpolate from
      std::string state_var_name = determineStateVariable(obs_point);
      
      // Perform variable-specific interpolation
      double interpolated_value = interpolateFromVariable(state, obs_point, state_var_name);
      result.push_back(interpolated_value);
    }

    return result;
  }

  /**
   * @brief Determine which state variable to interpolate from based on observation
   *
   * @param obs_point The observation point
   * @return Name of the state variable to interpolate from
   */
  std::string determineStateVariable(const typename ObsBackend::value_type& obs_point) const {
    // If we have required state variables specified, use the first one
    if (!required_state_vars_.empty()) {
      return required_state_vars_[0];
    }
    
    // For WRFState-like backends, try to determine variable from observation type
    if constexpr (requires { obs_point.variable_name; }) {
      // If observation has a variable name, try to map it to state variable
      std::string obs_var = obs_point.variable_name;
      
      // Simple mapping from observation variables to state variables
      if (obs_var == "T" || obs_var == "T2" || obs_var == "temperature") {
        return "T";  // Temperature
      } else if (obs_var == "Q" || obs_var == "Q2" || obs_var == "humidity") {
        return "QVAPOR";  // Water vapor mixing ratio
      } else if (obs_var == "U" || obs_var == "U10" || obs_var == "wind_u") {
        return "U";  // U-wind component
      } else if (obs_var == "V" || obs_var == "V10" || obs_var == "wind_v") {
        return "V";  // V-wind component
      } else if (obs_var == "P" || obs_var == "PSFC" || obs_var == "pressure") {
        return "P";  // Pressure
      } else {
        // Default to first available variable
        return obs_var;
      }
    }
    
    // Default: use first available variable from state
    // Note: We can't access state here, so we'll use a fallback
    // The actual state variable selection will be handled in interpolateFromVariable
    
    // Fallback: assume single variable state
    return "state";
  }

  /**
   * @brief Interpolate from a specific state variable at observation location
   *
   * @param state Model state
   * @param obs_point Observation point
   * @param state_var_name Name of the state variable to interpolate from
   * @return Interpolated value
   */
  double interpolateFromVariable(const StateBackend& state,
                                const typename ObsBackend::value_type& obs_point,
                                const std::string& state_var_name) const {
    // Convert observation location to grid coordinates
    double x, y, z;
    if (obs_point.location.getCoordinateSystem() == CoordinateSystem::GEOGRAPHIC) {
      auto [lat, lon, level] = obs_point.location.getGeographicCoords();
      
      // Convert geographic coordinates to grid coordinates
      auto [grid_i, grid_j, grid_k] = convertGeographicToGrid(lat, lon, level, state.geometry());
      x = static_cast<double>(grid_i);
      y = static_cast<double>(grid_j);
      z = static_cast<double>(grid_k);
    } else {
      // For non-geographic coordinates, use grid coordinates directly
      auto [i, j, k] = obs_point.location.getGridCoords();
      x = static_cast<double>(i);
      y = static_cast<double>(j);
      z = static_cast<double>(k);
    }

    // For multi-variable states, access the specific variable
    if constexpr (requires { state.at(state_var_name, std::declval<size_t>()); }) {
      // WRFState-like backends with variable-specific access
      return nearestNeighborInterpolationVariable(state, x, y, z, state_var_name);
    } else if constexpr (requires { std::declval<StateBackend>().getVariable(state_var_name); }) {
      // Backends with getVariable method
      auto var_data = state.getVariable(state_var_name);
      
      // Get dimensions for the variable
      size_t nx = static_cast<size_t>(state.geometry().x_dim());
      size_t ny = static_cast<size_t>(state.geometry().y_dim());
      size_t nz = 1; // Default to 1 for 2D geometries
      if constexpr (requires { state.geometry().z_dim(); }) {
        nz = state.geometry().z_dim();
      }
      
      return nearestNeighborInterpolationArray(var_data, x, y, z, nx, ny, nz);
    } else {
      // Fallback to single variable state
      // Get dimensions for the state
      size_t nx = static_cast<size_t>(state.geometry().x_dim());
      size_t ny = static_cast<size_t>(state.geometry().y_dim());
      size_t nz = 1; // Default to 1 for 2D geometries
      if constexpr (requires { state.geometry().z_dim(); }) {
        nz = state.geometry().z_dim();
      }
      
      return nearestNeighborInterpolation(state, x, y, z, nx, ny, nz);
    }
  }

  /**
   * @brief Get the state variables required by this observation operator
   * @return Const reference to vector of required state variable names
   */
  const std::vector<std::string>& getRequiredStateVars() const {
    return required_state_vars_;
  }

  /**
   * @brief Get the observation variables required by this observation operator
   * @return Const reference to vector of required observation variable names
   */
  const std::vector<std::string>& getRequiredObsVars() const {
    return required_obs_vars_;
  }

  /**
   * @brief Apply tangent linear observation operator: H dx
   *
   * @details For the identity operator, the tangent linear is the same as the
   * forward operator since H is linear.
   *
   * @param state_increment State increment to transform
   * @param reference_state Reference state (not used for linear operator)
   * @param obs Reference observations for context
   * @return Vector containing the transformed increment in observation space
   */
  std::vector<double> applyTangentLinear(const StateBackend& state_increment,
                                         [[maybe_unused]] const StateBackend& reference_state,
                                         const ObsBackend& obs) const {
    // For identity operator, tangent linear is the same as forward operator
    return apply(state_increment, obs);
  }

  /**
   * @brief Apply adjoint observation operator: H^T delta_y
   *
   * @details Maps from observation space back to state space. For the identity
   * operator, this maps observation increments back to the corresponding grid
   * points using the adjoint of nearest-neighbor interpolation.
   * For multi-variable states, maps to the appropriate variable based on
   * observation type and required state variables.
   *
   * @param obs_increment Observation space increment
   * @param reference_state Reference state to get structure
   * @param result_state State to store the adjoint result
   * @param obs Observations to determine grid coordinates for adjoint mapping
   */
  void applyAdjoint(const std::vector<double>& obs_increment,
                    const StateBackend& reference_state,
                    StateBackend& result_state,
                    const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("IdentityObsOperator not initialized");
    }

    // Initialize result state to zero
    result_state = std::move(*(reference_state.clone()));
    result_state.zero();

    // For the adjoint of nearest-neighbor interpolation:
    // 1. For each observation, determine its grid coordinates (same as in forward pass)
    // 2. Map each observation increment back to its corresponding grid point
    
    size_t obs_size = obs_increment.size();
    if (obs_size != obs.size()) {
      throw std::runtime_error("Observation increment size does not match observation count");
    }

    // Process each observation
    for (size_t obs_idx = 0; obs_idx < obs_size; ++obs_idx) {
      const auto& obs_point = obs[obs_idx];
      if (!obs_point.is_valid) {
        continue; // Skip invalid observations
      }

      // Determine which state variable to update
      std::string state_var_name = determineStateVariable(obs_point);
      
      // Convert observation location to grid coordinates (same as in forward pass)
      double x, y, z;
      if (obs_point.location.getCoordinateSystem() == CoordinateSystem::GEOGRAPHIC) {
        auto [lat, lon, level] = obs_point.location.getGeographicCoords();
        
        // Convert geographic coordinates to grid coordinates
        auto [grid_i, grid_j, grid_k] = convertGeographicToGrid(lat, lon, level, reference_state.geometry());
        x = static_cast<double>(grid_i);
        y = static_cast<double>(grid_j);
        z = static_cast<double>(grid_k);
      } else {
        // For non-geographic coordinates, use grid coordinates directly
        auto [i, j, k] = obs_point.location.getGridCoords();
        x = static_cast<double>(i);
        y = static_cast<double>(j);
        z = static_cast<double>(k);
      }

      // Adjoint of nearest-neighbor interpolation: add the observation increment
      // to the corresponding grid point in the appropriate variable
      adjointNearestNeighborInterpolationVariable(result_state, x, y, z, state_var_name, obs_increment[obs_idx]);
    }
  }

  /**
   * @brief Check if tangent linear and adjoint operators are available
   * @return True (identity operator always supports linearization)
   */
  bool supportsLinearization() const { return true; }

  /**
   * @brief Check if the observation operator is linear
   * @return True (identity operator is linear)
   */
  bool isLinear() const { return true; }

 private:
  /**
   * @brief Convert geographic coordinates to grid coordinates
   *
   * @details This method performs a nearest-neighbor search to find the closest
   * grid point to the given geographic coordinates. It requires the geometry
   * backend to provide coordinate lookup capabilities.
   *
   * @param lat Latitude in degrees
   * @param lon Longitude in degrees
   * @param level Vertical level/pressure
   * @param geometry Reference to the geometry backend
   * @return Tuple of (i, j, k) grid coordinates
   */
  template <typename GeometryBackend>
  std::tuple<size_t, size_t, size_t> convertGeographicToGrid(double lat, double lon, 
                                                   double level,
                                                   const GeometryBackend& geometry) const {
    // Try to use geometry-specific coordinate conversion if available
    if constexpr (requires { geometry.unstaggered_info(); }) {
      // For WRF-type geometry with coordinate arrays
      const auto& grid_info = geometry.unstaggered_info();
      
      if (grid_info.has_2d_coords()) {
        // Find closest grid point by searching through coordinate arrays
        double min_dist = std::numeric_limits<double>::max();
        size_t min_idx = 0;

        for (size_t idx = 0; idx < grid_info.longitude_2d.size(); ++idx) {
          double grid_lon = grid_info.longitude_2d[idx];
          double grid_lat = grid_info.latitude_2d[idx];

          // Calculate Euclidean distance in lat/lon space
          double dist = std::sqrt((lon - grid_lon) * (lon - grid_lon) +
                                 (lat - grid_lat) * (lat - grid_lat));

          if (dist < min_dist) {
            min_dist = dist;
            min_idx = idx;
          }
        }

        // Convert linear index to 2D grid coordinates
        size_t j = min_idx / grid_info.nx;
        size_t i = min_idx % grid_info.nx;
        
        // Find closest vertical level
        size_t k = 0;
        if (grid_info.has_vertical_coords()) {
          double min_vert_dist = std::numeric_limits<double>::max();
          for (size_t z = 0; z < grid_info.vertical_coords.size(); ++z) {
            double dist = std::abs(level - grid_info.vertical_coords[z]);
            if (dist < min_vert_dist) {
              min_vert_dist = dist;
              k = z;
            }
          }
        }
        
        return std::make_tuple(i, j, k);
      }
    }
    
    // Fallback: Search through geometry locations (slower but more general)
    double min_dist = std::numeric_limits<double>::max();
    size_t best_i = 0, best_j = 0, best_k = 0;
    
    size_t z_dim = 1; // Default to 1 for 2D geometries
    if constexpr (requires { geometry.z_dim(); }) {
      z_dim = geometry.z_dim();
    }
    
    for (size_t k = 0; k < z_dim; ++k) {
      for (size_t j = 0; j < static_cast<size_t>(geometry.y_dim()); ++j) {
        for (size_t i = 0; i < static_cast<size_t>(geometry.x_dim()); ++i) {
          try {
            // Try to get geographic coordinates for this grid point
            auto grid_location = [&]() {
              if constexpr (requires { geometry.getLocation(i, j, k); }) {
                return geometry.getLocation(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k));
              } else {
                return geometry.getLocation(static_cast<int>(i), static_cast<int>(j));
              }
            }();
            if (grid_location.getCoordinateSystem() == CoordinateSystem::GEOGRAPHIC) {
              auto [grid_lat, grid_lon, grid_level] = grid_location.getGeographicCoords();
              
              // Calculate distance in 3D space (lat, lon, level)
              double horizontal_dist = std::sqrt((lon - grid_lon) * (lon - grid_lon) +
                                               (lat - grid_lat) * (lat - grid_lat));
              double vertical_dist = std::abs(level - grid_level);
              
              // Combined distance (you might want to weight these differently)
              double dist = horizontal_dist + vertical_dist * 0.001; // Scale vertical distance
              
              if (dist < min_dist) {
                min_dist = dist;
                best_i = i;
                best_j = j;
                best_k = k;
              }
            }
          } catch (...) {
            // Skip points where geographic coordinates are not available
            continue;
          }
        }
      }
    }
    
    return std::make_tuple(best_i, best_j, best_k);
  }

  /**
   * @brief Perform nearest-neighbor interpolation from grid to point
   *
   * @param state Model state
   * @param x Grid i-coordinate (already converted from geographic if needed)
   * @param y Grid j-coordinate (already converted from geographic if needed)
   * @param z Grid k-coordinate (already converted from geographic if needed)
   * @param nx Number of grid points in x direction
   * @param ny Number of grid points in y direction
   * @param nz Number of grid points in z direction
   * @return Value at the nearest grid point
   */
  double nearestNeighborInterpolation(const StateBackend& state, double x, double y, double z,
                            size_t nx, size_t ny, size_t nz) const {
    // Clamp coordinates to grid bounds
    x = std::max(0.0, std::min(static_cast<double>(nx - 1), x));
    y = std::max(0.0, std::min(static_cast<double>(ny - 1), y));
    z = std::max(0.0, std::min(static_cast<double>(nz - 1), z));

    // Find nearest grid point (round to nearest integer)
    size_t i = static_cast<size_t>(std::round(x));
    size_t j = static_cast<size_t>(std::round(y));
    size_t k = static_cast<size_t>(std::round(z));

    // Ensure indices are within bounds
    if (i >= nx) i = nx - 1;
    if (j >= ny) j = ny - 1;
    if (k >= nz) k = nz - 1;

    // Try different access patterns depending on state backend interface
    
    // Pattern 1: Check for linear indexing with operator[]
    if constexpr (requires { state[std::declval<size_t>()]; }) {
      // Convert 3D grid coordinates to linear index using column-major order (like WRFState)
      size_t linear_index = k * (ny * nx) + j * nx + i;
      return state[linear_index];
    }
    // Pattern 2: Check for linear indexing with at(size_t)
    else if constexpr (requires { state.at(std::declval<size_t>()); }) {
      // Convert 3D grid coordinates to linear index using column-major order (like WRFState)
      size_t linear_index = k * (ny * nx) + j * nx + i;
      return state.at(linear_index);
    }
    // Pattern 3: Check for 3D coordinate access
    else if constexpr (requires { state.at(std::declval<int>(), std::declval<int>(), std::declval<int>()); }) {
      return state.at(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k));
    }
    // Pattern 4: Check for 2D coordinate access
    else if constexpr (requires { state.at(std::declval<int>(), std::declval<int>()); }) {
      return state.at(static_cast<int>(i), static_cast<int>(j));
    }
    else {
      throw std::runtime_error("State backend does not support any recognized access methods. "
                               "Backend must provide one of: "
                               "operator[](size_t), at(size_t), "
                               "at(int,int,int), or at(int,int)");
    }
  }

  /**
   * @brief Perform nearest-neighbor interpolation from grid to point for a specific variable
   *
   * @param state Model state
   * @param x Grid i-coordinate (already converted from geographic if needed)
   * @param y Grid j-coordinate (already converted from geographic if needed)
   * @param z Grid k-coordinate (already converted from geographic if needed)
   * @param state_var_name Name of the state variable to interpolate from
   * @return Value at the nearest grid point for the specific variable
   */
  double nearestNeighborInterpolationVariable(const StateBackend& state, double x, double y, double z,
                                            const std::string& state_var_name) const {
    // Get state dimensions
    size_t nx = static_cast<size_t>(state.geometry().x_dim());
    size_t ny = static_cast<size_t>(state.geometry().y_dim());
    size_t nz = 1; // Default to 1 for 2D geometries
    if constexpr (requires { state.geometry().z_dim(); }) {
      nz = state.geometry().z_dim();
    }

    // Clamp coordinates to grid bounds
    x = std::max(0.0, std::min(static_cast<double>(nx - 1), x));
    y = std::max(0.0, std::min(static_cast<double>(ny - 1), y));
    z = std::max(0.0, std::min(static_cast<double>(nz - 1), z));

    // Find nearest grid point (round to nearest integer)
    size_t i = static_cast<size_t>(std::round(x));
    size_t j = static_cast<size_t>(std::round(y));
    size_t k = static_cast<size_t>(std::round(z));

    // Ensure indices are within bounds
    if (i >= nx) i = nx - 1;
    if (j >= ny) j = ny - 1;
    if (k >= nz) k = nz - 1;

    // Try different access patterns depending on state backend interface
    
    // Pattern 1: WRFState-like backends with variable-specific linear indexing
    if constexpr (requires { state.at(state_var_name, std::declval<size_t>()); }) {
      // WRFState uses column-major (Fortran-style) indexing: [Z, Y, X]
      // The linear index should match WRFState's indexing scheme
      // WRFState: k * dims[1] * dims[2] + j * dims[2] + i
      // where dims[1] = Y dimension, dims[2] = X dimension
      
      // For multi-variable states, each variable can have different dimensions
      // We need to access the variable's actual dimensions, not geometry dimensions
      if constexpr (requires { std::declval<StateBackend>().getVariableDimensions(state_var_name); }) {
        // Get variable-specific dimensions
        const auto& var_dims = state.getVariableDimensions(state_var_name);
        
        // Calculate linear index using variable's actual dimensions
        size_t linear_index;
        if (var_dims.size() == 3) {
          // 3D variable: [Z, Y, X] - use column-major indexing
          linear_index = k * var_dims[1] * var_dims[2] + j * var_dims[2] + i;
        } else if (var_dims.size() == 2) {
          // 2D variable: [Y, X] - use column-major indexing
          linear_index = j * var_dims[1] + i;
        } else {
          // Fallback to geometry dimensions
          linear_index = k * (ny * nx) + j * nx + i;
        }
        return state.at(state_var_name, linear_index);
      } else {
        // Fallback to geometry dimensions for other state backends
        size_t linear_index = k * (ny * nx) + j * nx + i;  // Column-major: [Z, Y, X]
        return state.at(state_var_name, linear_index);
      }
    }
    // Pattern 2: Check for 3D coordinate access
    else if constexpr (requires { state.at(state_var_name, std::declval<int>(), std::declval<int>(), std::declval<int>()); }) {
      return state.at(state_var_name, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k));
    }
    // Pattern 4: Check for 2D coordinate access
    else if constexpr (requires { state.at(state_var_name, std::declval<int>(), std::declval<int>()); }) {
      return state.at(state_var_name, static_cast<int>(i), static_cast<int>(j));
    }
    // Pattern 5: Fallback to single variable state with linear indexing
    else if constexpr (requires { state[std::declval<size_t>()]; }) {
      // Convert 3D grid coordinates to linear index using row-major order
      size_t linear_index = k * (nx * ny) + j * nx + i;
      return state[linear_index];
    }
    // Pattern 6: Fallback to single variable state with at(size_t)
    else if constexpr (requires { state.at(std::declval<size_t>()); }) {
      // Convert 3D grid coordinates to linear index using row-major order
      size_t linear_index = k * (nx * ny) + j * nx + i;
      return state.at(linear_index);
    }
    else {
      throw std::runtime_error("State backend does not support any recognized access methods for variable interpolation. "
                               "Backend must provide one of: "
                               "at(variable_name, size_t), at(variable_name, int,int,int), "
                               "at(variable_name, int,int), operator[](size_t), or at(size_t)");
    }
  }

  /**
   * @brief Perform nearest-neighbor interpolation from grid to point for a specific variable (array-based)
   *
   * @param var_data Data array for the specific variable
   * @param x Grid i-coordinate (already converted from geographic if needed)
   * @param y Grid j-coordinate (already converted from geographic if needed)
   * @param z Grid k-coordinate (already converted from geographic if needed)
   * @param nx Number of grid points in x direction
   * @param ny Number of grid points in y direction
   * @param nz Number of grid points in z direction
   * @return Value at the nearest grid point for the specific variable
   */
  double nearestNeighborInterpolationArray(const std::vector<double>& var_data, double x, double y, double z,
                                         size_t nx, size_t ny, size_t nz) const {
    // Clamp coordinates to grid bounds
    x = std::max(0.0, std::min(static_cast<double>(nx - 1), x));
    y = std::max(0.0, std::min(static_cast<double>(ny - 1), y));
    z = std::max(0.0, std::min(static_cast<double>(nz - 1), z));

    // Find nearest grid point (round to nearest integer)
    size_t i = static_cast<size_t>(std::round(x));
    size_t j = static_cast<size_t>(std::round(y));
    size_t k = static_cast<size_t>(std::round(z));

    // Ensure indices are within bounds
    if (i >= nx) i = nx - 1;
    if (j >= ny) j = ny - 1;
    if (k >= nz) k = nz - 1;

    // Convert 3D grid coordinates to linear index using column-major order (like WRFState)
    size_t linear_index = k * (ny * nx) + j * nx + i;
    return var_data[linear_index];
  }

  /**
   * @brief Perform adjoint of nearest-neighbor interpolation from point to grid for a specific variable
   *
   * @param state Model state to update
   * @param x Grid x-coordinate (continuous)
   * @param y Grid y-coordinate (continuous)
   * @param z Grid z-coordinate (continuous)
   * @param state_var_name Name of the state variable to update
   * @param increment Value to add to the grid point
   */
  void adjointNearestNeighborInterpolationVariable(StateBackend& state, double x, double y, double z,
                                                  const std::string& state_var_name, double increment) const {
    // Get state dimensions
    size_t nx = static_cast<size_t>(state.geometry().x_dim());
    size_t ny = static_cast<size_t>(state.geometry().y_dim());
    size_t nz = 1; // Default to 1 for 2D geometries
    if constexpr (requires { state.geometry().z_dim(); }) {
      nz = state.geometry().z_dim();
    }

    // Clamp coordinates to grid bounds (same as in forward pass)
    x = std::max(0.0, std::min(static_cast<double>(nx - 1), x));
    y = std::max(0.0, std::min(static_cast<double>(ny - 1), y));
    z = std::max(0.0, std::min(static_cast<double>(nz - 1), z));

    // Find nearest grid point (same as in forward pass)
    size_t i = static_cast<size_t>(std::round(x));
    size_t j = static_cast<size_t>(std::round(y));
    size_t k = static_cast<size_t>(std::round(z));

    // Ensure indices are within bounds
    if (i >= nx) i = nx - 1;
    if (j >= ny) j = ny - 1;
    if (k >= nz) k = nz - 1;

    // Try different access patterns for variable-specific updates
    if constexpr (requires { state.at(state_var_name, std::declval<size_t>()); }) {
      // WRFState-like backends with variable-specific linear indexing
      // Use variable-specific dimensions for proper indexing
      if constexpr (requires { std::declval<StateBackend>().getVariableDimensions(state_var_name); }) {
        // Get variable-specific dimensions
        const auto& var_dims = state.getVariableDimensions(state_var_name);
        
        // Calculate linear index using variable's actual dimensions
        size_t linear_index;
        if (var_dims.size() == 3) {
          // 3D variable: [Z, Y, X] - use column-major indexing
          linear_index = k * var_dims[1] * var_dims[2] + j * var_dims[2] + i;
        } else if (var_dims.size() == 2) {
          // 2D variable: [Y, X] - use column-major indexing
          linear_index = j * var_dims[1] + i;
        } else {
          // Fallback to geometry dimensions
          linear_index = k * (ny * nx) + j * nx + i;
        }
        state.at(state_var_name, linear_index) += increment;
      } else {
        // Fallback to geometry dimensions for other state backends
        size_t linear_index = k * (ny * nx) + j * nx + i;
        state.at(state_var_name, linear_index) += increment;
      }
    } else if constexpr (requires { state.at(state_var_name, std::declval<int>(), std::declval<int>(), std::declval<int>()); }) {
      // Variable-specific 3D coordinate access
      state.at(state_var_name, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k)) += increment;
    } else if constexpr (requires { state.at(state_var_name, std::declval<int>(), std::declval<int>()); }) {
      // Variable-specific 2D coordinate access
      state.at(state_var_name, static_cast<int>(i), static_cast<int>(j)) += increment;
    } else if constexpr (requires { state[std::declval<size_t>()]; }) {
      // Fallback to single variable state
      // Use column-major order (like WRFState): k * (ny * nx) + j * nx + i
      size_t linear_index = k * (ny * nx) + j * nx + i;
      state[linear_index] += increment;
    } else if constexpr (requires { state.at(std::declval<size_t>()); }) {
      // Fallback to single variable state
      // Use column-major order (like WRFState): k * (ny * nx) + j * nx + i
      size_t linear_index = k * (ny * nx) + j * nx + i;
      state.at(linear_index) += increment;
    } else {
      throw std::runtime_error(
          "State backend does not support required access methods for variable-specific adjoint");
    }
  }

  /**
   * @brief Perform adjoint of nearest-neighbor interpolation from point to grid
   *
   * @param state Model state to update
   * @param x Grid x-coordinate (continuous)
   * @param y Grid y-coordinate (continuous)
   * @param z Grid z-coordinate (continuous)
   * @param nx Number of grid points in x direction
   * @param ny Number of grid points in y direction
   * @param nz Number of grid points in z direction
   * @param increment Value to add to the grid point
   */
  void adjointNearestNeighborInterpolation(StateBackend& state, double x, double y, double z,
                                          size_t nx, size_t ny, size_t nz, double increment) const {
    // Clamp coordinates to grid bounds (same as in forward pass)
    x = std::max(0.0, std::min(static_cast<double>(nx - 1), x));
    y = std::max(0.0, std::min(static_cast<double>(ny - 1), y));
    z = std::max(0.0, std::min(static_cast<double>(nz - 1), z));

    // Find nearest grid point (same as in forward pass)
    size_t i = static_cast<size_t>(std::round(x));
    size_t j = static_cast<size_t>(std::round(y));
    size_t k = static_cast<size_t>(std::round(z));

    // Ensure indices are within bounds
    if (i >= nx) i = nx - 1;
    if (j >= ny) j = ny - 1;
    if (k >= nz) k = nz - 1;

    if constexpr (requires { state[std::declval<size_t>()]; }) {
      // Convert 3D grid coordinates to linear index using column-major order (like WRFState)
      size_t linear_index = k * (ny * nx) + j * nx + i;
      state[linear_index] += increment;
    } else if constexpr (requires { state.at(std::declval<size_t>()); }) {
      // Convert 3D grid coordinates to linear index using column-major order (like WRFState)
      size_t linear_index = k * (ny * nx) + j * nx + i;
      state.at(linear_index) += increment;
    } else if constexpr (requires { state.at(std::declval<int>(), std::declval<int>(), std::declval<int>()); }) {
      state.at(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k)) += increment;
    } else if constexpr (requires { state.at(std::declval<int>(), std::declval<int>()); }) {
      state.at(static_cast<int>(i), static_cast<int>(j)) += increment;
    } else {
      throw std::runtime_error(
          "State backend does not support required access methods for adjoint");
    }
  }

  bool initialized_ = false;                      ///< Initialization status
  std::vector<std::string> required_state_vars_;  ///< Required state variables
  std::vector<std::string>
      required_obs_vars_;  ///< Required observation variables
};

}  // namespace metada::backends::common::obsoperator 