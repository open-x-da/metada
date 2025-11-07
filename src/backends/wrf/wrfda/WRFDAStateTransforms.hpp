/**
 * @file WRFDAStateTransforms.hpp
 * @brief C++ wrapper for WRFDA state transformation routines
 * @details Provides type-safe C++ interface to WRFDA's proven Fortran routines
 *          for state transformations without modifying WRFDA source code
 */

#pragma once

#include <stdexcept>
#include <string>

namespace metada::backends::wrf::wrfda {

// Forward declarations of C bridge functions from Fortran module
extern "C" {

/**
 * @brief Fortran bridge to da_transform_xtoxa
 * @param grid_ptr Pointer to WRF grid structure
 * @param error_code Output error code (0 = success)
 */
void wrfda_transform_xtoxa(void* grid_ptr, int* error_code);

/**
 * @brief Fortran bridge to da_transfer_xatowrf
 * @param grid_ptr Pointer to WRF grid structure
 * @param config_flags_ptr Pointer to WRF config_flags structure
 * @param error_code Output error code (0 = success)
 */
void wrfda_transfer_xatowrf(void* grid_ptr, void* config_flags_ptr,
                            int* error_code);

/**
 * @brief Fortran bridge to da_update_firstguess
 * @param grid_ptr Pointer to WRF grid structure
 * @param filename Output filename buffer (C-string)
 * @param filename_len Length of filename buffer
 * @param error_code Output error code (0 = success)
 */
void wrfda_update_firstguess(void* grid_ptr, const char* filename,
                             int filename_len, int* error_code);
}

/**
 * @brief C++ wrapper for WRFDA state transformation routines
 * @details Provides type-safe, exception-based interface to proven WRFDA
 *          routines for state transformations. This class delegates
 * WRF-specific physics and grid transformations to WRFDA's battle-tested code.
 */
class WRFDAStateTransforms {
 public:
  /**
   * @brief Compute diagnostic increments from prognostic increments
   * @details Calls WRFDA's da_transform_xtoxa to compute diagnostic variables
   *          needed by observation operators. This includes:
   *          - Pressure increments (xa%p) from surface pressure and moisture
   *          - Density increments (xa%rho) via hydrostatic balance
   *          - Geopotential height increments (xa%geoh)
   *          - Total precipitable water (for GPS observations)
   *          - GPS refractivity (for GPS RO observations)
   *          - Surface wind speed (for satellite wind observations)
   *          - Radar vertical velocity (for radar observations)
   *          - Cloud hydrometeors (for cloud/precipitation observations)
   *
   * @param grid_ptr Pointer to WRF grid structure containing increments in xa
   * @throws std::runtime_error if grid_ptr is null or WRFDA routine fails
   *
   * @note This should be called after control variable transformation (vtox)
   *       and before observation operators are applied
   *
   * @note This is a const method because it only computes derived quantities
   *       from existing increments
   */
  static void transformXToXa(void* grid_ptr) {
    if (grid_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDAStateTransforms::transformXToXa: grid_ptr is null");
    }

    int error_code = 0;
    wrfda_transform_xtoxa(grid_ptr, &error_code);

    if (error_code != 0) {
      throw std::runtime_error(
          "WRFDA da_transform_xtoxa failed with error code: " +
          std::to_string(error_code) + getErrorMessage(error_code));
    }
  }

  /**
   * @brief Add analysis increments to background state
   * @details Calls WRFDA's da_transfer_xatowrf to perform the complete
   *          transformation from analysis increments to full WRF analysis
   * state. This handles all WRF-specific physics and grid transformations:
   *
   *          Physical transformations:
   *          - Converts specific humidity increments to mixing ratio increments
   *          - Computes dry air mass increments from surface pressure and
   * moisture
   *          - Converts temperature increments to potential temperature
   *          - Computes full pressure field from dry air mass and moisture
   *          - Computes geopotential height using WRF's vertical coordinate
   *          - Recomputes 2m/10m diagnostic fields (T2, Q2, U10, V10, TH2)
   *
   *          Grid transformations:
   *          - Handles Arakawa-C grid staggering (A-grid to C-grid conversion)
   *          - Properly interpolates wind components to staggered locations
   *          - Handles boundary conditions for regional domains
   *
   *          Quality control:
   *          - Applies positivity constraints (e.g., moisture >= 0)
   *          - Ensures physical consistency of analysis state
   *
   * @param grid_ptr Pointer to WRF grid structure
   * @param config_flags_ptr Pointer to WRF configuration structure
   * @throws std::runtime_error if pointers are null or WRFDA routine fails
   *
   * @note After this routine:
   *       - grid%u_2, grid%v_2, grid%w_2 contain full analysis wind
   *       - grid%t_2 contains full analysis potential temperature perturbation
   *       - grid%mu_2 contains full analysis dry air mass
   *       - grid%ph_2 contains full analysis geopotential
   *       - grid%moist contains full analysis moisture species
   *       - grid%psfc contains full analysis surface pressure
   *       - All fields are ready for WRF model integration
   *
   * @note This routine implements the fundamental DA equation:
   *       x_analysis = x_background + x_increment
   *       with all WRF-specific transformations applied
   */
  static void transferXaToWRF(void* grid_ptr, void* config_flags_ptr) {
    if (grid_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDAStateTransforms::transferXaToWRF: grid_ptr is null");
    }

    if (config_flags_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDAStateTransforms::transferXaToWRF: config_flags_ptr is null");
    }

    int error_code = 0;
    wrfda_transfer_xatowrf(grid_ptr, config_flags_ptr, &error_code);

    if (error_code != 0) {
      throw std::runtime_error(
          "WRFDA da_transfer_xatowrf failed with error code: " +
          std::to_string(error_code) + getErrorMessage(error_code));
    }
  }

  /**
   * @brief Update first-guess file using WRFDA da_update_firstguess
   * @details Copies the background (fg) file to the requested output and
   *          updates only the analysis variables touched by WRFDA, preserving
   *          all metadata and non-assimilated fields.
   *
   * @param grid_ptr Pointer to WRF grid structure
   * @param filename Desired output filename; empty string uses WRFDA default
   * @throws std::runtime_error if pointers are null or WRFDA routine fails
   */
  static void updateFirstGuess(void* grid_ptr, const std::string& filename) {
    if (grid_ptr == nullptr) {
      throw std::runtime_error(
          "WRFDAStateTransforms::updateFirstGuess: grid_ptr is null");
    }

    int error_code = 0;
    const char* fname = filename.c_str();
    int length = static_cast<int>(filename.size());

    wrfda_update_firstguess(grid_ptr, fname, length, &error_code);

    if (error_code != 0) {
      throw std::runtime_error(
          "WRFDA da_update_firstguess failed with error code: " +
          std::to_string(error_code) + getErrorMessage(error_code));
    }
  }

 private:
  /**
   * @brief Get human-readable error message for error code
   * @param error_code Error code from Fortran routine
   * @return Error message string
   */
  static std::string getErrorMessage(int error_code) {
    switch (error_code) {
      case 0:
        return " (Success)";
      case 1:
        return " (Null grid pointer)";
      case 2:
        return " (Null config_flags pointer)";
      case 3:
        return " (Invalid filename length)";
      default:
        return " (Unknown error in WRFDA routine)";
    }
  }
};

}  // namespace metada::backends::wrf::wrfda
