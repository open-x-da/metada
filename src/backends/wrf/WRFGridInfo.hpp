/**
 * @file WRFGridInfo.hpp
 * @brief WRF grid information and staggering configuration structures
 * @ingroup wrf_backend
 * @author Metada Framework Team
 *
 * @details
 * This file contains data structures and utilities for managing WRF grid
 * information, including dimension names, staggering configurations, and
 * grid type determinations. These structures support the Arakawa-C staggered
 * grid system used by the Weather Research and Forecasting (WRF) model.
 *
 * The grid information is essential for proper coordinate transformations,
 * interpolation operations, and variable associations within the WRF backend
 * implementation of the metada framework.
 */
#pragma once

#include <string>

namespace metada::backends::wrf {

/**
 * @brief Structure containing grid information for WRF variables
 *
 * @details This structure stores dimension names and staggering information
 * for WRF variables, allowing proper association with corresponding geometry
 * grids. It supports both 2D and 3D variables with various staggering
 * configurations used in the WRF model's Arakawa-C grid system.
 *
 * The WRF model uses staggered grids to improve numerical accuracy:
 * - Mass points (unstaggered): Temperature, pressure, humidity
 * - U-points (x-staggered): Zonal wind components
 * - V-points (y-staggered): Meridional wind components
 * - W-points (z-staggered): Vertical wind components
 *
 * @example
 * @code
 * // Configure grid info for temperature variable (mass points)
 * VariableGridInfo tempGrid;
 * tempGrid.x_dim_name = "west_east";
 * tempGrid.y_dim_name = "south_north";
 * tempGrid.z_dim_name = "bottom_top";
 * tempGrid.grid_type = VariableGridInfo::GridType::UNSTAGGERED;
 *
 * // Configure grid info for U-wind component (x-staggered)
 * VariableGridInfo uWindGrid;
 * uWindGrid.x_dim_name = "west_east_stag";
 * uWindGrid.y_dim_name = "south_north";
 * uWindGrid.z_dim_name = "bottom_top";
 * uWindGrid.x_staggered = true;
 * uWindGrid.grid_type = uWindGrid.determineGridType();
 * @endcode
 */
struct VariableGridInfo {
  // ============================================================================
  // DIMENSION INFORMATION
  // ============================================================================

  ///@{ @name Dimension Information
  std::string
      x_dim_name;  ///< Name of X dimension in NetCDF file (e.g., "west_east")
  std::string
      y_dim_name;  ///< Name of Y dimension in NetCDF file (e.g., "south_north")
  std::string
      z_dim_name;  ///< Name of Z dimension in NetCDF file (e.g., "bottom_top")
  ///@}

  // ============================================================================
  // STAGGERING CONFIGURATION
  // ============================================================================

  ///@{ @name Staggering Configuration
  bool x_staggered =
      false;  ///< Whether variable is staggered in X direction (U-grid)
  bool y_staggered =
      false;  ///< Whether variable is staggered in Y direction (V-grid)
  bool z_staggered =
      false;  ///< Whether variable is staggered in Z direction (W-grid)
  ///@}

  // ============================================================================
  // GRID TYPE DEFINITIONS
  // ============================================================================

  /**
   * @brief Grid type enumeration for WRF variable association
   *
   * @details Defines the type of staggered grid used by a variable,
   * which determines which geometry grid should be used for coordinate
   * transformations and interpolation operations. Based on the Arakawa-C
   * grid staggering convention used in WRF.
   */
  enum class GridType {
    UNSTAGGERED,  ///< Mass points grid (unstaggered in all directions) - T, P,
                  ///< RH
    U_STAGGERED,  ///< U-wind grid (staggered in X direction) - Zonal wind
                  ///< components
    V_STAGGERED,  ///< V-wind grid (staggered in Y direction) - Meridional wind
                  ///< components
    W_STAGGERED   ///< W-wind grid (staggered in Z direction) - Vertical wind
                  ///< components
  };

  GridType grid_type = GridType::UNSTAGGERED;  ///< Associated grid type

  // ============================================================================
  // UTILITY METHODS
  // ============================================================================

  /**
   * @brief Automatically determine grid type from staggering configuration
   *
   * @details Analyzes the staggering flags to determine the appropriate
   * grid type following WRF conventions where typically only one dimension
   * should be staggered per variable type. This method provides automatic
   * grid type detection based on the dimension staggering pattern.
   *
   * @return GridType The determined grid type based on staggering pattern
   *
   * @note If multiple dimensions are staggered or no staggering is detected,
   *       the method defaults to UNSTAGGERED grid type.
   */
  GridType determineGridType() const {
    if (x_staggered && !y_staggered && !z_staggered) {
      return GridType::U_STAGGERED;  // Staggered in X (U wind component)
    } else if (!x_staggered && y_staggered && !z_staggered) {
      return GridType::V_STAGGERED;  // Staggered in Y (V wind component)
    } else if (!x_staggered && !y_staggered && z_staggered) {
      return GridType::W_STAGGERED;  // Staggered in Z (W wind component)
    } else {
      return GridType::UNSTAGGERED;  // Not staggered or multiple staggering
    }
  }

  /**
   * @brief Get grid type as string representation
   *
   * @details Converts the grid type enumeration to a human-readable
   * string for debugging, logging, and configuration file output.
   * Useful for diagnostics and troubleshooting grid associations.
   *
   * @return std::string String representation of the grid type
   */
  std::string getGridTypeString() const {
    switch (grid_type) {
      case GridType::UNSTAGGERED:
        return "unstaggered";
      case GridType::U_STAGGERED:
        return "u_staggered";
      case GridType::V_STAGGERED:
        return "v_staggered";
      case GridType::W_STAGGERED:
        return "w_staggered";
      default:
        return "unknown";
    }
  }
};

}  // namespace metada::backends::wrf