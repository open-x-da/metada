/**
 * @file WRFIncrement.hpp
 * @brief WRF increment backend implementation operating on grid%xa
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
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

/// Extract all analysis increment fields to arrays (18 fields)
int wrfda_extract_all_analysis_increments(
    void* grid_ptr, double* u, double* v, double* w, double* t, double* p,
    double* q, double* qt, double* rh, double* rho, double* geoh, double* wh,
    double* qcw, double* qrn, double* qci, double* qsn, double* qgr,
    double* psfc, double* mu, int* nx, int* ny, int* nz);

/// Compute norm of analysis increments
double wrfda_xa_norm(void* grid_ptr);

/// Compute dot product of two increment structures
double wrfda_xa_dot(void* grid_ptr1, void* grid_ptr2);

/// Perform AXPY on analysis increments: xa1 = xa1 + alpha * xa2
void wrfda_xa_axpy(double alpha, void* grid_ptr1, void* grid_ptr2);

/// Scale analysis increments: xa = alpha * xa
void wrfda_xa_scale(double alpha, void* grid_ptr);

/// Update analysis increments from all 18 grid%xa field arrays
/// @note Includes 16 3D fields and 2 2D fields
void wrfda_update_analysis_increments(
    const double* u, const double* v, const double* w, const double* t,
    const double* p, const double* q, const double* qt, const double* rh,
    const double* rho, const double* geoh, const double* wh, const double* qcw,
    const double* qrn, const double* qci, const double* qsn, const double* qgr,
    const double* psfc, const double* mu, const void* grid_ptr);
}

namespace metada::backends::wrf {

/**
 * @brief WRF increment backend implementation for analysis increments
 *
 * @details This class implements an increment backend for the Weather Research
 * and Forecasting (WRF) model, using **Control Variable (CV) space storage**
 * for true increment independence.
 *
 * ## Architecture: CV Space Storage with WRFDA Integration
 *
 * Each WRFIncrement instance:
 * 1. Stores data internally as `std::vector<double>` (18 fields, CV space)
 * 2. Uses `grid%xa` as a **staging area** for WRFDA operations only
 * 3. Provides true independence between increment objects
 *
 * This design aligns with WRFDA's minimization approach where increments are
 * stored as control variable vectors, and `grid%xa` is working memory.
 *
 * ### Stored Fields (18 fields):
 * **3D Fields (16)**:
 * - **Core**: u, v, w, t, p, q
 * - **Moisture/Thermo**: qt, rh
 * - **Derived**: rho, geoh, wh
 * - **Hydrometeors**: qcw, qrn, qci, qsn, qgr
 *
 * **2D Fields (2)**:
 * - **Surface**: psfc, mu
 *
 * ### WRFDA Integration:
 * - `syncToGrid()`: Write internal data → `grid%xa` before WRFDA operations
 * - `syncFromGrid()`: Extract `grid%xa` → internal data after WRFDA operations
 * - WRFDA operations (norm, dot, axpy) use `grid%xa` as staging area
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
        total_size_(16 * size_3d_ +
                    2 * size_2d_),  // 16 3D fields + 2 2D fields
        data_(total_size_, 0.0) {   // Initialize internal CV storage with zeros

    // Note: Internal data is already zeroed via initialization
    // grid%xa will be synced when needed for WRFDA operations
  }

  /**
   * @brief Copy constructor
   *
   * @details Copies internal CV storage - provides true increment independence
   */
  WRFIncrement(const WRFIncrement& other)
      : geometry_(other.geometry_),
        nx_(other.nx_),
        ny_(other.ny_),
        nz_(other.nz_),
        size_3d_(other.size_3d_),
        size_2d_(other.size_2d_),
        total_size_(other.total_size_),
        data_(other.data_) {  // Direct copy of internal CV storage
    // Each increment now has independent data storage
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
      data_ = other.data_;  // Direct copy of internal CV storage
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
   * @brief Zero all increments in internal storage
   */
  void zero() { std::fill(data_.begin(), data_.end(), 0.0); }

  /**
   * @brief Scale increment by scalar: xa *= alpha
   *
   * @param alpha Scaling factor
   */
  void scale(double alpha) {
    for (auto& val : data_) {
      val *= alpha;
    }
  }

  /**
   * @brief AXPY operation: this += alpha * other
   *
   * @param alpha Scaling factor
   * @param other Increment to add
   */
  void axpy(double alpha, const WRFIncrement& other) {
    if (data_.size() != other.data_.size()) {
      throw std::runtime_error(
          "Cannot perform AXPY with increments of different sizes");
    }
    for (size_t i = 0; i < data_.size(); ++i) {
      data_[i] += alpha * other.data_[i];
    }
  }

  /**
   * @brief Compute dot product with another increment
   *
   * @param other Increment to compute dot product with
   * @return Dot product value
   */
  double dot(const WRFIncrement& other) const {
    if (data_.size() != other.data_.size()) {
      throw std::runtime_error(
          "Cannot compute dot product with increments of different sizes");
    }
    double result = 0.0;
    for (size_t i = 0; i < data_.size(); ++i) {
      result += data_[i] * other.data_[i];
    }
    return result;
  }

  /**
   * @brief Compute L2 norm of increment
   *
   * @return L2 norm value
   */
  double norm() const {
    double sum_sq = 0.0;
    for (const auto& val : data_) {
      sum_sq += val * val;
    }
    return std::sqrt(sum_sq);
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
   * @brief Extract all increment fields to arrays (18 fields)
   *
   * @details Extracts from internal CV storage:
   *          - 3D fields: u, v, w, t, p, q, qt, rh, rho, geoh, wh, qcw, qrn,
   * qci, qsn, qgr
   *          - 2D fields: psfc, mu
   *
   * @param[out] u, v, w, t, p, q, qt, rh, rho, geoh, wh, qcw, qrn, qci, qsn,
   * qgr 3D arrays (size: nx*ny*nz)
   * @param[out] psfc, mu 2D arrays (size: nx*ny)
   */
  void extract(double* u, double* v, double* w, double* t, double* p, double* q,
               double* qt, double* rh, double* rho, double* geoh, double* wh,
               double* qcw, double* qrn, double* qci, double* qsn, double* qgr,
               double* psfc, double* mu) const {
    const size_t size_3d = size_3d_;
    const size_t size_2d = size_2d_;

    std::copy_n(&data_[0 * size_3d], size_3d, u);
    std::copy_n(&data_[1 * size_3d], size_3d, v);
    std::copy_n(&data_[2 * size_3d], size_3d, w);
    std::copy_n(&data_[3 * size_3d], size_3d, t);
    std::copy_n(&data_[4 * size_3d], size_3d, p);
    std::copy_n(&data_[5 * size_3d], size_3d, q);
    std::copy_n(&data_[6 * size_3d], size_3d, qt);
    std::copy_n(&data_[7 * size_3d], size_3d, rh);
    std::copy_n(&data_[8 * size_3d], size_3d, rho);
    std::copy_n(&data_[9 * size_3d], size_3d, geoh);
    std::copy_n(&data_[10 * size_3d], size_3d, wh);
    std::copy_n(&data_[11 * size_3d], size_3d, qcw);
    std::copy_n(&data_[12 * size_3d], size_3d, qrn);
    std::copy_n(&data_[13 * size_3d], size_3d, qci);
    std::copy_n(&data_[14 * size_3d], size_3d, qsn);
    std::copy_n(&data_[15 * size_3d], size_3d, qgr);
    std::copy_n(&data_[16 * size_3d], size_2d, psfc);
    std::copy_n(&data_[16 * size_3d + size_2d], size_2d, mu);
  }

  /**
   * @brief Get all increment data as a single vector
   *
   * @details Returns internal CV storage containing all 18 fields:
   *          1. 3D fields (16): u, v, w, t, p, q, qt, rh, rho, geoh, wh,
   *             qcw, qrn, qci, qsn, qgr
   *          2. 2D fields (2): psfc, mu
   *
   * @return Vector containing all increment fields (direct copy of internal
   * data) Total size: 16 * (nx*ny*nz) + 2 * (nx*ny)
   */
  std::vector<double> getData() const {
    return data_;  // Return copy of internal CV storage
  }

  /**
   * @brief Randomize all increment fields for testing
   *
   * @details Fills internal CV storage with random values in [-0.5, 0.5]
   */
  void randomize() {
    for (auto& val : data_) {
      val = (double(std::rand()) / RAND_MAX) - 0.5;
    }
  }

  /**
   * @brief Copy increment from another WRFIncrement
   *
   * @details Copies internal CV storage directly
   *
   * @param[in] other Source increment to copy from
   */
  void copyFrom(const WRFIncrement& other) {
    if (data_.size() != other.data_.size()) {
      throw std::runtime_error("Cannot copy increments with different sizes");
    }
    data_ = other.data_;
  }

  /**
   * @brief Sync internal CV storage to grid%xa
   *
   * @details Writes internal data to WRFDA's grid%xa structure
   *          This is called before WRFDA operations that read from grid%xa
   */
  void syncToGrid() const {
    // Decompose internal CV vector into 18 field arrays
    const size_t size_3d = size_3d_;
    const size_t size_2d = size_2d_;

    const double* u = &data_[0 * size_3d];
    const double* v = &data_[1 * size_3d];
    const double* w = &data_[2 * size_3d];
    const double* t = &data_[3 * size_3d];
    const double* p = &data_[4 * size_3d];
    const double* q = &data_[5 * size_3d];
    const double* qt = &data_[6 * size_3d];
    const double* rh = &data_[7 * size_3d];
    const double* rho = &data_[8 * size_3d];
    const double* geoh = &data_[9 * size_3d];
    const double* wh = &data_[10 * size_3d];
    const double* qcw = &data_[11 * size_3d];
    const double* qrn = &data_[12 * size_3d];
    const double* qci = &data_[13 * size_3d];
    const double* qsn = &data_[14 * size_3d];
    const double* qgr = &data_[15 * size_3d];
    const double* psfc = &data_[16 * size_3d];
    const double* mu = &data_[16 * size_3d + size_2d];

    void* grid_ptr = geometry_.getGridPtr();
    wrfda_update_analysis_increments(u, v, w, t, p, q, qt, rh, rho, geoh, wh,
                                     qcw, qrn, qci, qsn, qgr, psfc, mu,
                                     grid_ptr);
  }

  /**
   * @brief Sync grid%xa to internal CV storage
   *
   * @details Extracts WRFDA's grid%xa structure to internal data
   *          This is called after WRFDA operations that write to grid%xa
   */
  void syncFromGrid() {
    // Decompose internal CV vector into 18 field arrays
    const size_t size_3d = size_3d_;
    const size_t size_2d = size_2d_;

    double* u = &data_[0 * size_3d];
    double* v = &data_[1 * size_3d];
    double* w = &data_[2 * size_3d];
    double* t = &data_[3 * size_3d];
    double* p = &data_[4 * size_3d];
    double* q = &data_[5 * size_3d];
    double* qt = &data_[6 * size_3d];
    double* rh = &data_[7 * size_3d];
    double* rho = &data_[8 * size_3d];
    double* geoh = &data_[9 * size_3d];
    double* wh = &data_[10 * size_3d];
    double* qcw = &data_[11 * size_3d];
    double* qrn = &data_[12 * size_3d];
    double* qci = &data_[13 * size_3d];
    double* qsn = &data_[14 * size_3d];
    double* qgr = &data_[15 * size_3d];
    double* psfc = &data_[16 * size_3d];
    double* mu = &data_[16 * size_3d + size_2d];

    void* grid_ptr = geometry_.getGridPtr();
    int nx = nx_, ny = ny_, nz = nz_;
    int rc = wrfda_extract_all_analysis_increments(
        grid_ptr, u, v, w, t, p, q, qt, rh, rho, geoh, wh, qcw, qrn, qci, qsn,
        qgr, psfc, mu, &nx, &ny, &nz);

    if (rc != 0) {
      throw std::runtime_error("Failed to sync from grid%xa: error code " +
                               std::to_string(rc));
    }
  }

  /**
   * @brief Transfer state difference to increment
   *
   * @details Extracts all available fields from state backend to internal CV
   * storage. Maps WRF state variables to increment fields:
   *   - Core meteorological: U→u, V→v, W→w, T→t, QVAPOR→q
   *   - Hydrometeors: QCLOUD→qcw, QRAIN→qrn, QICE→qci, QSNOW→qsn, QGRAUPEL→qgr
   *   - Surface: PSFC→psfc, MU→mu
   *   - Diagnostic fields (p, qt, rh, rho, geoh, wh) are set to zero
   *
   * @param state_backend State backend containing difference data
   */
  template <typename StateBackend>
  void transferFromState(const StateBackend& state_backend) {
    // Get pointers to state fields
    auto* u_state = static_cast<const double*>(state_backend.getData("U"));
    auto* v_state = static_cast<const double*>(state_backend.getData("V"));
    auto* w_state = static_cast<const double*>(state_backend.getData("W"));
    auto* t_state = static_cast<const double*>(state_backend.getData("T"));
    auto* q_state = static_cast<const double*>(state_backend.getData("QVAPOR"));
    auto* psfc_state =
        static_cast<const double*>(state_backend.getData("PSFC"));
    auto* mu_state = static_cast<const double*>(state_backend.getData("MU"));
    auto* qcloud_state =
        static_cast<const double*>(state_backend.getData("QCLOUD"));
    auto* qrain_state =
        static_cast<const double*>(state_backend.getData("QRAIN"));
    auto* qice_state =
        static_cast<const double*>(state_backend.getData("QICE"));
    auto* qsnow_state =
        static_cast<const double*>(state_backend.getData("QSNOW"));
    auto* qgraupel_state =
        static_cast<const double*>(state_backend.getData("QGRAUPEL"));

    // Copy to internal CV storage (18 fields in order)
    const size_t size_3d = size_3d_;
    const size_t size_2d = size_2d_;

    std::copy_n(u_state, size_3d, &data_[0 * size_3d]);          // u
    std::copy_n(v_state, size_3d, &data_[1 * size_3d]);          // v
    std::copy_n(w_state, size_3d, &data_[2 * size_3d]);          // w
    std::copy_n(t_state, size_3d, &data_[3 * size_3d]);          // t
    std::fill_n(&data_[4 * size_3d], size_3d, 0.0);              // p (zero)
    std::copy_n(q_state, size_3d, &data_[5 * size_3d]);          // q
    std::fill_n(&data_[6 * size_3d], size_3d, 0.0);              // qt (zero)
    std::fill_n(&data_[7 * size_3d], size_3d, 0.0);              // rh (zero)
    std::fill_n(&data_[8 * size_3d], size_3d, 0.0);              // rho (zero)
    std::fill_n(&data_[9 * size_3d], size_3d, 0.0);              // geoh (zero)
    std::fill_n(&data_[10 * size_3d], size_3d, 0.0);             // wh (zero)
    std::copy_n(qcloud_state, size_3d, &data_[11 * size_3d]);    // qcw
    std::copy_n(qrain_state, size_3d, &data_[12 * size_3d]);     // qrn
    std::copy_n(qice_state, size_3d, &data_[13 * size_3d]);      // qci
    std::copy_n(qsnow_state, size_3d, &data_[14 * size_3d]);     // qsn
    std::copy_n(qgraupel_state, size_3d, &data_[15 * size_3d]);  // qgr
    std::copy_n(psfc_state, size_2d, &data_[16 * size_3d]);      // psfc
    std::copy_n(mu_state, size_2d, &data_[16 * size_3d + size_2d]);  // mu
  }

  /**
   * @brief Get data pointer for direct access (for advanced users)
   *
   * @warning This provides direct access to internal CV storage. Use with
   * caution!
   */
  template <typename T>
  T* getDataPtr() {
    static_assert(std::is_same_v<T, double>,
                  "WRFIncrement only supports double precision");
    return data_.data();
  }

  /**
   * @brief Get const data pointer for direct access
   */
  template <typename T>
  const T* getDataPtr() const {
    static_assert(std::is_same_v<T, double>,
                  "WRFIncrement only supports double precision");
    return data_.data();
  }

 private:
  const GeometryBackend& geometry_;  ///< Grid geometry reference
  int nx_;                           ///< Number of grid points in X
  int ny_;                           ///< Number of grid points in Y
  int nz_;                           ///< Number of vertical levels
  size_t size_3d_;                   ///< Size of 3D fields
  size_t size_2d_;                   ///< Size of 2D fields
  size_t total_size_;                ///< Total number of elements

  /// Internal CV storage: 18 fields (16 3D + 2 2D) concatenated
  /// Layout: [u, v, w, t, p, q, qt, rh, rho, geoh, wh, qcw, qrn, qci, qsn, qgr,
  /// psfc, mu]
  std::vector<double> data_;
};

}  // namespace metada::backends::wrf
