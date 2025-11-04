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
#include <unordered_map>
#include <vector>

// WRFDA C bindings for increment operations on grid%xa
extern "C" {
// Field configuration types and query function
#include "wrfda_types.h"  // For wrfda_field_info_t

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
 * @brief Field descriptor for WRF analysis increments
 * @details Describes a single field in the control variable space
 */
struct FieldDescriptor {
  std::string name;  ///< Field name (e.g., "u", "v", "t")
  size_t offset;     ///< Offset in data_ array
  size_t size;       ///< Number of elements in field
  bool is_3d;        ///< True for 3D fields, false for 2D
  bool is_active;    ///< True if field is active based on cv_options
};

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
 * 1. Stores data internally as `std::vector<double>` (CV space, only active
 * fields)
 * 2. Uses `grid%xa` as a **staging area** for WRFDA operations only
 * 3. Provides true independence between increment objects
 * 4. Dynamically adapts to cv_options configuration
 *
 * This design aligns with WRFDA's minimization approach where increments are
 * stored as control variable vectors, and `grid%xa` is working memory.
 *
 * ### Stored Fields (depends on cv_options):
 * **cv_options=5 (Standard 3D-Var)**:
 * - **3D**: u, v, t, q (4 fields)
 * - **2D**: psfc (1 field)
 *
 * **cv_options=6 (With Hydrometeors)**:
 * - **3D**: u, v, t, q, qcloud, qrain, qice, qsnow, qgraupel (9 fields)
 * - **2D**: psfc (1 field)
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
        size_2d_(nx_ * ny_) {
    // Query WRFDA's field configuration
    initializeFieldDescriptors();

    // Allocate storage for active fields only
    data_.resize(total_size_, 0.0);

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
        fields_(other.fields_),
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
   * @note Only active fields (based on cv_options) contribute to the norm
   */
  double norm() const {
    // Note: data_ now only contains active fields, so this automatically
    // respects cv_options configuration
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
   * @details Extracts active fields from internal storage to output arrays.
   *          Inactive fields are filled with zeros based on cv_options.
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

    // Initialize all output arrays to zero (for inactive fields)
    std::fill_n(u, size_3d, 0.0);
    std::fill_n(v, size_3d, 0.0);
    std::fill_n(w, size_3d, 0.0);
    std::fill_n(t, size_3d, 0.0);
    std::fill_n(p, size_3d, 0.0);
    std::fill_n(q, size_3d, 0.0);
    std::fill_n(qt, size_3d, 0.0);
    std::fill_n(rh, size_3d, 0.0);
    std::fill_n(rho, size_3d, 0.0);
    std::fill_n(geoh, size_3d, 0.0);
    std::fill_n(wh, size_3d, 0.0);
    std::fill_n(qcw, size_3d, 0.0);
    std::fill_n(qrn, size_3d, 0.0);
    std::fill_n(qci, size_3d, 0.0);
    std::fill_n(qsn, size_3d, 0.0);
    std::fill_n(qgr, size_3d, 0.0);
    std::fill_n(psfc, size_2d, 0.0);
    std::fill_n(mu, size_2d, 0.0);

    // Map field names to output pointers
    std::unordered_map<std::string, double*> outputs = {
        {"u", u},       {"v", v},       {"w", w},       {"t", t},
        {"p", p},       {"q", q},       {"qt", qt},     {"rh", rh},
        {"rho", rho},   {"geoh", geoh}, {"wh", wh},     {"qcloud", qcw},
        {"qrain", qrn}, {"qice", qci},  {"qsnow", qsn}, {"qgraupel", qgr},
        {"psfc", psfc}, {"mu", mu}};

    // Copy only active fields from data_ to output arrays
    for (const auto& field : fields_) {
      if (!field.is_active) continue;

      auto it = outputs.find(field.name);
      if (it != outputs.end()) {
        std::copy_n(&data_[field.offset], field.size, it->second);
      }
    }
  }

  /**
   * @brief Get all increment data as a single vector
   *
   * @details Returns internal CV storage containing only active fields based on
   * cv_options. For cv_options=5: u, v, t, q, psfc (4 3D fields + 1 2D field)
   *          For cv_options=6: adds qcloud, qrain, qice, qsnow, qgraupel
   *
   * @return Vector containing active increment fields (direct copy of internal
   * data) Total size depends on cv_options configuration
   */
  std::vector<double> getData() const {
    return data_;  // Return copy of internal CV storage (active fields only)
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
   * @details Writes active fields from internal data to WRFDA's grid%xa
   * structure. Inactive fields are set to zero. This is called before WRFDA
   * operations that read from grid%xa.
   */
  void syncToGrid() const {
    const size_t size_3d = size_3d_;
    const size_t size_2d = size_2d_;

    // Allocate temporary buffers for all 18 fields (WRFDA expects all fields)
    // Active fields get data from data_, inactive fields get zeros
    std::vector<double> u_buf(size_3d, 0.0), v_buf(size_3d, 0.0),
        w_buf(size_3d, 0.0);
    std::vector<double> t_buf(size_3d, 0.0), p_buf(size_3d, 0.0),
        q_buf(size_3d, 0.0);
    std::vector<double> qt_buf(size_3d, 0.0), rh_buf(size_3d, 0.0),
        rho_buf(size_3d, 0.0);
    std::vector<double> geoh_buf(size_3d, 0.0), wh_buf(size_3d, 0.0);
    std::vector<double> qcw_buf(size_3d, 0.0), qrn_buf(size_3d, 0.0),
        qci_buf(size_3d, 0.0);
    std::vector<double> qsn_buf(size_3d, 0.0), qgr_buf(size_3d, 0.0);
    std::vector<double> psfc_buf(size_2d, 0.0), mu_buf(size_2d, 0.0);

    // Map field names to buffer pointers
    std::unordered_map<std::string, double*> buffers = {
        {"u", u_buf.data()},       {"v", v_buf.data()},
        {"w", w_buf.data()},       {"t", t_buf.data()},
        {"p", p_buf.data()},       {"q", q_buf.data()},
        {"qt", qt_buf.data()},     {"rh", rh_buf.data()},
        {"rho", rho_buf.data()},   {"geoh", geoh_buf.data()},
        {"wh", wh_buf.data()},     {"qcloud", qcw_buf.data()},
        {"qrain", qrn_buf.data()}, {"qice", qci_buf.data()},
        {"qsnow", qsn_buf.data()}, {"qgraupel", qgr_buf.data()},
        {"psfc", psfc_buf.data()}, {"mu", mu_buf.data()}};

    // Copy active fields from data_ to appropriate buffers
    for (const auto& field : fields_) {
      if (!field.is_active) continue;

      auto it = buffers.find(field.name);
      if (it != buffers.end()) {
        std::copy_n(&data_[field.offset], field.size, it->second);
      }
    }

    // Call WRFDA with all 18 field buffers
    void* grid_ptr = geometry_.getGridPtr();
    wrfda_update_analysis_increments(
        u_buf.data(), v_buf.data(), w_buf.data(), t_buf.data(), p_buf.data(),
        q_buf.data(), qt_buf.data(), rh_buf.data(), rho_buf.data(),
        geoh_buf.data(), wh_buf.data(), qcw_buf.data(), qrn_buf.data(),
        qci_buf.data(), qsn_buf.data(), qgr_buf.data(), psfc_buf.data(),
        mu_buf.data(), grid_ptr);
  }

  /**
   * @brief Sync grid%xa to internal CV storage
   *
   * @details Extracts active fields from WRFDA's grid%xa structure to internal
   * data. Only active fields are extracted; inactive fields are ignored. This
   * is called after WRFDA operations that write to grid%xa.
   */
  void syncFromGrid() {
    const size_t size_3d = size_3d_;
    const size_t size_2d = size_2d_;

    // Allocate temporary buffers for all 18 fields (WRFDA returns all fields)
    std::vector<double> u_buf(size_3d), v_buf(size_3d), w_buf(size_3d);
    std::vector<double> t_buf(size_3d), p_buf(size_3d), q_buf(size_3d);
    std::vector<double> qt_buf(size_3d), rh_buf(size_3d), rho_buf(size_3d);
    std::vector<double> geoh_buf(size_3d), wh_buf(size_3d);
    std::vector<double> qcw_buf(size_3d), qrn_buf(size_3d), qci_buf(size_3d);
    std::vector<double> qsn_buf(size_3d), qgr_buf(size_3d);
    std::vector<double> psfc_buf(size_2d), mu_buf(size_2d);

    // Extract all 18 fields from WRFDA
    void* grid_ptr = geometry_.getGridPtr();
    int nx = nx_, ny = ny_, nz = nz_;
    int rc = wrfda_extract_all_analysis_increments(
        grid_ptr, u_buf.data(), v_buf.data(), w_buf.data(), t_buf.data(),
        p_buf.data(), q_buf.data(), qt_buf.data(), rh_buf.data(),
        rho_buf.data(), geoh_buf.data(), wh_buf.data(), qcw_buf.data(),
        qrn_buf.data(), qci_buf.data(), qsn_buf.data(), qgr_buf.data(),
        psfc_buf.data(), mu_buf.data(), &nx, &ny, &nz);

    if (rc != 0) {
      throw std::runtime_error("Failed to sync from grid%xa: error code " +
                               std::to_string(rc));
    }

    // Map field names to buffer pointers
    std::unordered_map<std::string, const double*> buffers = {
        {"u", u_buf.data()},       {"v", v_buf.data()},
        {"w", w_buf.data()},       {"t", t_buf.data()},
        {"p", p_buf.data()},       {"q", q_buf.data()},
        {"qt", qt_buf.data()},     {"rh", rh_buf.data()},
        {"rho", rho_buf.data()},   {"geoh", geoh_buf.data()},
        {"wh", wh_buf.data()},     {"qcloud", qcw_buf.data()},
        {"qrain", qrn_buf.data()}, {"qice", qci_buf.data()},
        {"qsnow", qsn_buf.data()}, {"qgraupel", qgr_buf.data()},
        {"psfc", psfc_buf.data()}, {"mu", mu_buf.data()}};

    // Copy only active fields from buffers to data_
    for (const auto& field : fields_) {
      if (!field.is_active) continue;

      auto it = buffers.find(field.name);
      if (it != buffers.end()) {
        std::copy_n(it->second, field.size, &data_[field.offset]);
      }
    }
  }

  /**
   * @brief Transfer state difference to increment
   *
   * @details Extracts active fields from state backend to internal CV storage.
   *          Only fields configured as active (based on cv_options) are copied.
   *          Maps WRF state variables to increment fields:
   *   - Core meteorological: U→u, V→v, T→t, QVAPOR→q
   *   - Hydrometeors: QCLOUD→qcloud, QRAIN→qrain, QICE→qice, QSNOW→qsnow,
   * QGRAUPEL→qgraupel
   *   - Surface: PSFC→psfc
   *
   * @param state_backend State backend containing difference data
   */
  template <typename StateBackend>
  void transferFromState(const StateBackend& state_backend) {
    // Map field names to WRF state variable names
    const std::unordered_map<std::string, std::string> field_to_var = {
        {"u", "U"},       {"v", "V"},           {"w", "W"},
        {"t", "T"},       {"q", "QVAPOR"},      {"p", "P"},
        {"rho", "RHO"},   {"qcloud", "QCLOUD"}, {"qrain", "QRAIN"},
        {"qice", "QICE"}, {"qsnow", "QSNOW"},   {"qgraupel", "QGRAUPEL"},
        {"psfc", "PSFC"}};

    // Copy active fields from state to internal storage
    for (const auto& field : fields_) {
      if (!field.is_active) continue;

      // Find corresponding state variable name
      auto it = field_to_var.find(field.name);
      if (it == field_to_var.end()) {
        // Skip fields not in the mapping (they may be diagnostic)
        continue;
      }

      // Get data from state backend
      const auto* state_data =
          static_cast<const double*>(state_backend.getData(it->second));

      // Copy to appropriate location in data_
      std::copy_n(state_data, field.size, &data_[field.offset]);
    }
  }

  /**
   * @brief Get active field descriptors
   * @return Vector of active field descriptors
   * @note Useful for debugging and introspection of cv_options configuration
   */
  const std::vector<FieldDescriptor>& getActiveFields() const {
    return fields_;
  }

  /**
   * @brief Get total size of increment (number of active field elements)
   * @return Total number of elements in data_ array
   */
  size_t getTotalSize() const { return total_size_; }

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
  /**
   * @brief Initialize field descriptors based on WRFDA cv_options
   * @details Queries WRFDA and builds field map for active fields only
   */
  void initializeFieldDescriptors() {
    // Query WRFDA's field configuration
    wrfda_field_info_t field_info;
    int rc = wrfda_get_field_info(&field_info);
    if (rc != 0) {
      throw std::runtime_error("Failed to query WRFDA field configuration");
    }

    // Build field descriptors based on active fields
    size_t offset = 0;
    fields_.clear();

    // cv_options=5: u, v, t, q, psfc (4 3D + 1 2D)
    if (field_info.has_u) {
      fields_.push_back({"u", offset, size_3d_, true, true});
      offset += size_3d_;
    }
    if (field_info.has_v) {
      fields_.push_back({"v", offset, size_3d_, true, true});
      offset += size_3d_;
    }
    if (field_info.has_t) {
      fields_.push_back({"t", offset, size_3d_, true, true});
      offset += size_3d_;
    }
    if (field_info.has_q) {
      fields_.push_back({"q", offset, size_3d_, true, true});
      offset += size_3d_;
    }
    if (field_info.has_psfc) {
      fields_.push_back({"psfc", offset, size_2d_, false, true});
      offset += size_2d_;
    }

    // cv_options=6: + hydrometeors
    if (field_info.has_qcloud) {
      fields_.push_back({"qcloud", offset, size_3d_, true, true});
      offset += size_3d_;
    }
    if (field_info.has_qrain) {
      fields_.push_back({"qrain", offset, size_3d_, true, true});
      offset += size_3d_;
    }
    if (field_info.has_qice) {
      fields_.push_back({"qice", offset, size_3d_, true, true});
      offset += size_3d_;
    }
    if (field_info.has_qsnow) {
      fields_.push_back({"qsnow", offset, size_3d_, true, true});
      offset += size_3d_;
    }
    if (field_info.has_qgraupel) {
      fields_.push_back({"qgraupel", offset, size_3d_, true, true});
      offset += size_3d_;
    }

    // Additional fields for cv_options > 6
    if (field_info.has_w) {
      fields_.push_back({"w", offset, size_3d_, true, true});
      offset += size_3d_;
    }
    if (field_info.has_p) {
      fields_.push_back({"p", offset, size_3d_, true, true});
      offset += size_3d_;
    }
    if (field_info.has_rho) {
      fields_.push_back({"rho", offset, size_3d_, true, true});
      offset += size_3d_;
    }

    // Store total size
    total_size_ = offset;
  }

  const GeometryBackend& geometry_;  ///< Grid geometry reference
  int nx_;                           ///< Number of grid points in X
  int ny_;                           ///< Number of grid points in Y
  int nz_;                           ///< Number of vertical levels
  size_t size_3d_;                   ///< Size of 3D fields
  size_t size_2d_;                   ///< Size of 2D fields
  size_t
      total_size_;  ///< Total number of elements (computed from active fields)

  /// Field descriptors (populated based on cv_options)
  std::vector<FieldDescriptor> fields_;

  /// Internal CV storage: Only active fields concatenated
  /// Layout depends on cv_options (e.g., cv_options=5: [u, v, t, q, psfc])
  std::vector<double> data_;
};

}  // namespace metada::backends::wrf
