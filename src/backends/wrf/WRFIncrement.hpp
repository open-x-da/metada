/**
 * @file WRFIncrement.hpp
 * @brief WRF increment backend implementation operating on grid%xa
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

// WRFDA C bindings for increment operations on grid%xa
extern "C" {
/// Zero analysis increments (grid%xa) - uses WRFDA's da_zero_x
int wrfda_zero_xa(void* grid_ptr);

/// Copy analysis increments: dst = src (ALL 35+ fields) - uses WRFDA's
/// da_copy_xa
void wrfda_copy_xa(void* grid_dst_ptr, void* grid_src_ptr);

/// Randomize ALL analysis increment fields (35+ fields) for testing
void wrfda_randomize_xa(void* grid_ptr);

/// Extract analysis increments to arrays
int wrfda_extract_analysis_increments(void* grid_ptr, double* u, double* v,
                                      double* t, double* q, double* psfc,
                                      int* nx, int* ny, int* nz);

/// Compute norm of analysis increments
double wrfda_xa_norm(void* grid_ptr);

/// Compute dot product of two increment structures
double wrfda_xa_dot(void* grid_ptr1, void* grid_ptr2);

/// Perform AXPY on analysis increments: xa1 = xa1 + alpha * xa2
void wrfda_xa_axpy(double alpha, void* grid_ptr1, void* grid_ptr2);

/// Scale analysis increments: xa = alpha * xa
void wrfda_xa_scale(double alpha, void* grid_ptr);

/// Update analysis increments from 5-variable arrays (control variables only)
/// @note For complete x_type alignment, use wrfda_copy_xa instead
void wrfda_update_analysis_increments(const double* u, const double* v,
                                      const double* t, const double* q,
                                      const double* psfc, const void* grid_ptr);
}

namespace metada::backends::wrf {

/**
 * @brief WRF increment backend implementation for analysis increments
 *
 * @details This class implements an increment backend for the Weather Research
 * and Forecasting (WRF) model, specifically managing analysis increments
 * stored in grid%xa (NOT grid%xb which is the background state).
 *
 * ## Architecture: Full WRFDA x_type Alignment
 *
 * WRFIncrement operates directly on WRFDA's grid%xa which contains **35+
 * fields**:
 *
 * ### 3D Fields (17 fields):
 * - **Control variables**: u, v, t, q (wind, temperature, moisture)
 * - **Diagnostic variables**: w, p, geoh, rh, wh, rho, ref
 * - **Hydrometeors**: qcw, qrn, qci, qsn, qgr, qt
 *
 * ### 2D Fields (18+ fields):
 * - **Surface**: psfc, mu, tgrn, u10, v10, t2, q2
 * - **Derived**: ztd, tpw, speed
 * - **Brightness temps**: tb19v, tb19h, tb22v, tb37v, tb37h, tb85v, tb85h
 *
 * ### WRFDA Native Functions Used:
 * - `da_zero_x`: Zero ALL increment fields properly
 * - `da_copy_xa`: Copy complete increment structure (ALL 35+ fields)
 * - `da_add_xa`: Add increments with proper field handling
 *
 * This ensures **complete alignment** with WRFDA, not just control variables!
 *
 * The increment represents δx = x - xb in the incremental formulation:
 * J(δx) = 1/2 * δx^T B^-1 δx + 1/2 * (d - H(xb + δx))^T R^-1 (d - H(xb + δx))
 *
 * @tparam ConfigBackend Configuration backend type
 * @tparam GeometryBackend Geometry backend type providing grid information
 *
 * @see WRFState for background state (grid%xb)
 * @see framework::Increment
 * @see WRFDA's da_define_structures.f90 for complete x_type definition
 */
template <typename ConfigBackend, typename GeometryBackend>
class WRFIncrement {
 public:
  /**
   * @brief Construct increment from geometry (creates zero increment)
   *
   * @param geometry Grid geometry defining the increment space
   */
  explicit WRFIncrement(const GeometryBackend& geometry)
      : geometry_(geometry),
        nx_(static_cast<int>(geometry.x_dim())),
        ny_(static_cast<int>(geometry.y_dim())),
        nz_(static_cast<int>(geometry.z_dim())),
        size_3d_(nx_ * ny_ * nz_),
        size_2d_(nx_ * ny_),
        total_size_(3 * size_3d_ + size_3d_ +
                    size_2d_) {  // u,v,t,q (3D) + psfc (2D)

    // Zero the analysis increments in grid%xa
    zero();
  }

  /**
   * @brief Copy constructor
   *
   * @details Uses WRFDA's da_copy_xa to copy ALL 35+ increment fields
   */
  WRFIncrement(const WRFIncrement& other)
      : geometry_(other.geometry_),
        nx_(other.nx_),
        ny_(other.ny_),
        nz_(other.nz_),
        size_3d_(other.size_3d_),
        size_2d_(other.size_2d_),
        total_size_(other.total_size_) {
    // Use WRFDA's proven da_copy_xa function for complete increment copy
    copyFrom(other);
  }

  /**
   * @brief Copy assignment operator
   */
  WRFIncrement& operator=(const WRFIncrement& other) {
    if (this != &other) {
      if (nx_ != other.nx_ || ny_ != other.ny_ || nz_ != other.nz_) {
        throw std::runtime_error(
            "Cannot assign increments with different dimensions");
      }
      copyFrom(other);
    }
    return *this;
  }

  /**
   * @brief Move constructor
   */
  WRFIncrement(WRFIncrement&& other) noexcept = default;

  /**
   * @brief Move assignment operator
   */
  WRFIncrement& operator=(WRFIncrement&& other) noexcept = default;

  /**
   * @brief Destructor
   */
  ~WRFIncrement() = default;

  /**
   * @brief Zero all analysis increments (grid%xa = 0)
   */
  void zero() {
    void* grid_ptr = geometry_.getGridPtr();
    int rc = wrfda_zero_xa(grid_ptr);
    if (rc != 0) {
      throw std::runtime_error(
          "Failed to zero analysis increments: error code " +
          std::to_string(rc));
    }
  }

  /**
   * @brief Scale increment by scalar: xa *= alpha
   *
   * @param alpha Scaling factor
   */
  void scale(double alpha) {
    void* grid_ptr = geometry_.getGridPtr();
    wrfda_xa_scale(alpha, grid_ptr);
  }

  /**
   * @brief AXPY operation: this += alpha * other
   *
   * @param alpha Scaling factor
   * @param other Increment to add
   */
  void axpy(double alpha, const WRFIncrement& other) {
    void* grid_ptr = geometry_.getGridPtr();
    void* other_grid_ptr = other.geometry_.getGridPtr();
    wrfda_xa_axpy(alpha, grid_ptr, other_grid_ptr);
  }

  /**
   * @brief Compute dot product with another increment
   *
   * @param other Increment to compute dot product with
   * @return Dot product value
   */
  double dot(const WRFIncrement& other) const {
    void* grid_ptr = geometry_.getGridPtr();
    void* other_grid_ptr = other.geometry_.getGridPtr();
    return wrfda_xa_dot(grid_ptr, other_grid_ptr);
  }

  /**
   * @brief Compute L2 norm of increment
   *
   * @return L2 norm value
   */
  double norm() const {
    void* grid_ptr = geometry_.getGridPtr();
    return wrfda_xa_norm(grid_ptr);
  }

  /**
   * @brief Addition assignment operator
   */
  WRFIncrement& operator+=(const WRFIncrement& other) {
    axpy(1.0, other);
    return *this;
  }

  /**
   * @brief Subtraction assignment operator
   */
  WRFIncrement& operator-=(const WRFIncrement& other) {
    axpy(-1.0, other);
    return *this;
  }

  /**
   * @brief Scalar multiplication assignment operator
   */
  WRFIncrement& operator*=(double alpha) {
    scale(alpha);
    return *this;
  }

  /**
   * @brief Scalar division assignment operator
   */
  WRFIncrement& operator/=(double alpha) {
    scale(1.0 / alpha);
    return *this;
  }

  /**
   * @brief Get total number of elements in increment
   */
  size_t size() const { return total_size_; }

  /**
   * @brief Get 3D grid dimensions
   */
  int getNx() const { return nx_; }
  int getNy() const { return ny_; }
  int getNz() const { return nz_; }

  /**
   * @brief Get geometry
   */
  const GeometryBackend& geometry() const { return geometry_; }

  /**
   * @brief Extract increment data to arrays (for testing/debugging)
   *
   * @param[out] u U-wind increment (size: nx*ny*nz)
   * @param[out] v V-wind increment (size: nx*ny*nz)
   * @param[out] t Temperature increment (size: nx*ny*nz)
   * @param[out] q Moisture increment (size: nx*ny*nz)
   * @param[out] psfc Surface pressure increment (size: nx*ny)
   */
  void extract(double* u, double* v, double* t, double* q, double* psfc) const {
    void* grid_ptr = geometry_.getGridPtr();
    int nx = nx_, ny = ny_, nz = nz_;
    int rc = wrfda_extract_analysis_increments(grid_ptr, u, v, t, q, psfc, &nx,
                                               &ny, &nz);
    if (rc != 0) {
      throw std::runtime_error(
          "Failed to extract analysis increments: error code " +
          std::to_string(rc));
    }
  }

  /**
   * @brief Get all increment data as a single vector
   * @return Vector containing all increment fields concatenated
   */
  std::vector<double> getData() const {
    std::vector<double> u(nx_ * ny_ * nz_);
    std::vector<double> v(nx_ * ny_ * nz_);
    std::vector<double> t(nx_ * ny_ * nz_);
    std::vector<double> q(nx_ * ny_ * nz_);
    std::vector<double> psfc(nx_ * ny_);

    extract(u.data(), v.data(), t.data(), q.data(), psfc.data());

    std::vector<double> result;
    result.reserve(u.size() + v.size() + t.size() + q.size() + psfc.size());
    result.insert(result.end(), u.begin(), u.end());
    result.insert(result.end(), v.begin(), v.end());
    result.insert(result.end(), t.begin(), t.end());
    result.insert(result.end(), q.begin(), q.end());
    result.insert(result.end(), psfc.begin(), psfc.end());

    return result;
  }

  /**
   * @brief Randomize ALL increment fields for testing
   *
   * @details Fills ALL 35+ increment fields with random values in [-0.5, 0.5].
   * This includes:
   * - 3D fields: u, v, w, t, q, p, geoh, rh, wh, rho, ref, qcw, qrn, qci, qsn,
   * qgr, qt
   * - 2D fields: psfc, mu, tgrn, u10, v10, t2, q2, ztd, tpw, speed, brightness
   * temps
   *
   * Uses WRFDA's da_zero_x to initialize, then randomizes following da_copy_xa
   * structure.
   */
  void randomize() {
    void* grid_ptr = geometry_.getGridPtr();
    wrfda_randomize_xa(grid_ptr);
  }

  /**
   * @brief Copy increment from another WRFIncrement
   *
   * @details Uses WRFDA's proven da_copy_xa function to copy ALL fields
   * in the increment structure (35+ fields: u, v, w, t, q, p, psfc, mu,
   * hydrometeors, diagnostics, surface vars, brightness temps, etc.)
   *
   * This is the RECOMMENDED way to update increments - it ensures complete
   * alignment with WRFDA's increment structure.
   *
   * @param[in] other Source increment to copy from
   */
  void copyFrom(const WRFIncrement& other) {
    void* dst_grid_ptr = geometry_.getGridPtr();
    void* src_grid_ptr = other.geometry_.getGridPtr();
    wrfda_copy_xa(dst_grid_ptr, src_grid_ptr);
  }

  /**
   * @brief Transfer state difference to increment
   *
   * @details Extracts data from state backend and updates this increment.
   * This is used when creating increments from state differences.
   *
   * @note Currently only transfers 5 control variables (u,v,t,q,psfc).
   * This is sufficient for standard 3DVAR where these are the control
   * variables.
   *
   * @param state_backend State backend containing difference data
   */
  template <typename StateBackend>
  void transferFromState(const StateBackend& state_backend) {
    auto* u = static_cast<const double*>(state_backend.getData("U"));
    auto* v = static_cast<const double*>(state_backend.getData("V"));
    auto* t = static_cast<const double*>(state_backend.getData("T"));
    auto* q = static_cast<const double*>(state_backend.getData("QVAPOR"));
    auto* psfc = static_cast<const double*>(state_backend.getData("PSFC"));

    void* grid_ptr = geometry_.getGridPtr();
    wrfda_update_analysis_increments(u, v, t, q, psfc, grid_ptr);
  }

  /**
   * @brief Get data pointer for direct access (for advanced users)
   *
   * @warning This provides direct access to grid%xa. Use with caution!
   */
  template <typename T>
  T* getDataPtr() {
    static_assert(std::is_same_v<T, double>,
                  "WRFIncrement only supports double precision");
    // For now, we don't expose direct pointer access since grid%xa
    // is managed by WRFDA. Users should use extract() instead.
    throw std::runtime_error(
        "Direct pointer access not supported for WRFIncrement. Use extract() "
        "instead.");
  }

 private:
  const GeometryBackend& geometry_;  ///< Grid geometry reference
  int nx_;                           ///< Number of grid points in X
  int ny_;                           ///< Number of grid points in Y
  int nz_;                           ///< Number of vertical levels
  size_t size_3d_;                   ///< Size of 3D fields
  size_t size_2d_;                   ///< Size of 2D fields
  size_t total_size_;                ///< Total number of elements
};

}  // namespace metada::backends::wrf
