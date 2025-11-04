/**
 * @file wrfda_types.h
 * @brief Common type definitions for WRFDA C API
 *
 * @details This header contains only type definitions used across
 *          WRFDA integration, avoiding circular dependencies and
 *          minimizing compilation dependencies.
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Field information structure for WRFDA control variables
 * @details Describes which fields are active based on cv_options configuration
 */
typedef struct {
  int cv_options;    /**< WRFDA cv_options value */
  int num_3d_fields; /**< Number of active 3D fields */
  int num_2d_fields; /**< Number of active 2D fields */

  // Field flags (cv_options=5: u, v, t, q, psfc)
  bool has_u;    /**< U-wind component (3D) */
  bool has_v;    /**< V-wind component (3D) */
  bool has_t;    /**< Temperature (3D) */
  bool has_q;    /**< Specific humidity (3D) */
  bool has_psfc; /**< Surface pressure (2D) */

  // Additional fields for cv_options > 5
  bool has_w;        /**< Vertical velocity (3D) */
  bool has_rho;      /**< Density (3D) */
  bool has_p;        /**< Pressure (3D) */
  bool has_qcloud;   /**< Cloud water mixing ratio (3D) */
  bool has_qrain;    /**< Rain water mixing ratio (3D) */
  bool has_qice;     /**< Ice mixing ratio (3D) */
  bool has_qsnow;    /**< Snow mixing ratio (3D) */
  bool has_qgraupel; /**< Graupel mixing ratio (3D) */
} wrfda_field_info_t;

/**
 * @brief Query WRFDA's control variable configuration
 * @details Returns information about which fields are active based on
 * cv_options
 * @param[out] field_info Structure to populate with field information
 * @return 0 on success, non-zero on error
 */
int wrfda_get_field_info(wrfda_field_info_t* field_info);

#ifdef __cplusplus
}
#endif

