/**
 * @file MACOMGeometry.hpp
 * @brief MACOM geometry backend implementation for hexagonal grid
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <iostream>
// #include <map>
#include <algorithm>  // for std::sort
#include <cmath>
#include <limits>
#include <memory>
#include <netcdf>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>  // for std::pair
#include <vector>

// Include nanoflann library (add to project dependencies)
#include "include/MACOMlogging.hpp"
#include "nanoflann.hpp"

// Forward declarations
namespace metada::backends::macom {
template <typename ConfigBackend>
class MACOMGeometryIterator;
template <typename ConfigBackend>
class MACOMGeometryConstIterator;

// GeoPoint structure for geographical coordinates and 2D KD-tree needs
struct GeoPoint {
  size_t index;     // Grid point index
  double lon;       // Longitude (degrees)
  double lat;       // Latitude (degrees)
  double x;         // Cartesian coordinate X
  double y;         // Cartesian coordinate Y
  double distance;  // Distance from query point (km)
};

// Query result structure
struct GeoQueryResult {
  std::vector<GeoPoint> points;
};

// Point cloud adapter for 2D KD-tree interface
template <typename T>
struct PointCloud {
  std::vector<T> pts;

  // Required interface methods
  inline size_t kdtree_get_point_count() const { return pts.size(); }

  inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
    if (dim == 0)
      return pts[idx].x;
    else
      return pts[idx].y;
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const {
    return false;
  }
};

/**
 * @brief Structure to hold vertical interpolation information
 */
struct VerticalPoint {
  size_t lower_index;   // Index of the lower grid point
  size_t upper_index;   // Index of the upper grid point
  double interp_coef;   // Interpolation coefficient (0-1)
  bool is_outside;      // Whether the point is outside the vertical range
  double lower_depth;   // Depth of the lower grid point
  double upper_depth;   // Depth of the upper grid point
  double target_depth;  // Target depth that was queried
};
}  // namespace metada::backends::macom

namespace metada::backends::macom {

/**
 * @brief MACOM geometry backend implementation for hexagonal grid
 *
 * @details
 * This implementation handles MACOM's hexagonal grid structure and reads
 * grid data from NetCDF files. The grid data includes cell centers, edges,
 * areas and masks.
 */
template <typename ConfigBackend>
class MACOMGeometry {
 public:
  // Iterator type aliases
  using iterator = MACOMGeometryIterator<ConfigBackend>;
  using const_iterator = MACOMGeometryConstIterator<ConfigBackend>;

  // --- Deleted constructors and assignment operators ---
  MACOMGeometry() = delete;
  MACOMGeometry(const MACOMGeometry&) = delete;
  MACOMGeometry& operator=(const MACOMGeometry&) = delete;

  /**
   * @brief Constructor that takes a configuration backend
   *
   * @param config Configuration containing geometry settings
   */
  explicit MACOMGeometry(const ConfigBackend& config);

  /**
   * @brief Move constructor
   *
   * @param other MACOM geometry backend to move from
   */
  MACOMGeometry(MACOMGeometry&& other) noexcept;

  /**
   * @brief Move assignment operator
   *
   * @param other MACOM geometry backend to move from
   * @return Reference to this geometry after assignment
   */
  MACOMGeometry& operator=(MACOMGeometry&& other) noexcept;

  /**
   * @brief Destructor
   */
  ~MACOMGeometry() = default;

  /**
   * @brief Clone this geometry
   *
   * @return A new MACOMGeometry instance (by value)
   */
  MACOMGeometry clone() const;

  /**
   * @brief Get iterator to the beginning of the grid
   *
   * @return Iterator pointing to the first grid point
   */
  iterator begin();

  /**
   * @brief Get iterator to the end of the grid
   *
   * @return Iterator pointing past the last grid point
   */
  iterator end();

  /**
   * @brief Get const iterator to the beginning of the grid
   *
   * @return Const iterator pointing to the first grid point
   */
  const_iterator begin() const;

  /**
   * @brief Get const iterator to the end of the grid
   *
   * @return Const iterator pointing past the last grid point
   */
  const_iterator end() const;

  /**
   * @brief Get the total number of grid points
   *
   * @return Total number of grid points in the geometry
   */
  std::size_t totalGridSize() const;

  /**
   * @brief Check if geometry is properly initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Check if the geometry is periodic in X dimension
   *
   * @return Always false in this implementation
   */
  bool isPeriodicX() const { return false; }

  /**
   * @brief Check if the geometry is periodic in Y dimension
   *
   * @return Always false in this implementation
   */
  bool isPeriodicY() const { return false; }

  /**
   * @brief Check if the geometry is periodic in Z dimension
   *
   * @return Always false in this implementation
   */
  bool isPeriodicZ() const { return false; }

  /**
   * @brief Get the number of grid points (nlpb)
   * @return nlpb_
   */
  std::size_t getNlpb() const { return nlpb_; }

  /**
   * @brief Get the number of vertical levels (nk)
   * @return nk_
   */
  std::size_t getNk() const { return nk_; }

  /**
   * @brief Get the number of vertical levels + 1 (nkp1)
   * @return nkp1_
   */
  std::size_t getNkp1() const { return nkp1_; }

  /**
   * @brief Get the total number of grid points (nl)
   * @return nl_
   */
  std::size_t getNl() const { return nl_; }

  /**
   * @brief Get the number of grid points for vorticity (nlpbz)
   * @return nlpbz_
   */
  std::size_t getNlpbz() const { return nlpbz_; }

  /**
   * @brief Get the total number of grid points for vorticity (nlz)
   * @return nlz_
   */
  std::size_t getNlz() const { return nlz_; }

  /**
   * @brief Get the number of boundary points (nlbdy)
   * @return nlbdy_
   */
  std::size_t getNlbdy() const { return nlbdy_; }

  /**
   * @brief Get the number of iterations (ni)
   * @return ni_
   */
  std::size_t getNi() const { return ni_; }

  template <typename StateBackend>
  void haloExchange([[maybe_unused]] StateBackend& state) {}

  /**
   * @brief Find the nearest grid point to the given latitude and longitude
   *
   * @param lon Longitude (degrees)
   * @param lat Latitude (degrees)
   * @return GeoPoint Containing index and coordinates of the nearest grid point
   */
  GeoPoint findNearestGridPoint(double lon, double lat) const;

  /**
   * @brief Find all grid points within the specified radius of the given
   * coordinates
   *
   * @param lon Longitude (degrees)
   * @param lat Latitude (degrees)
   * @param radius Search radius (kilometers)
   * @return GeoQueryResult Containing indices and coordinates of all points
   * within radius
   */
  GeoQueryResult findGridPointsInRadius(double lon, double lat,
                                        double radius) const;

  /**
   * @brief Print statistics about the horizontal grid data to help with
   * debugging
   */
  void debugGridDataHorizontal();

  /**
   * @brief Print statistics about the vertical grid data to help with debugging
   */
  void debugGridDataVertical();

  /**
   * @brief Alternative implementation for direct nearest point search using
   * Haversine
   */
  GeoPoint findNearestGridPointDirect(double lon, double lat) const;

  /**
   * @brief Find nearest grid points for multiple query locations
   *
   * @param query_lons Vector of query longitudes
   * @param query_lats Vector of query latitudes
   * @param use_check_multiple_neighbors Whether to use multiple neighbors for
   * more accurate results
   * @return std::vector<GeoPoint> Vector of nearest points for each query
   * location
   */
  std::vector<GeoPoint> findNearestGridPointsBatch(
      const std::vector<double>& query_lons,
      const std::vector<double>& query_lats) const;

  /**
   * @brief Find grid points within radius for multiple query locations
   *
   * @param query_lons Vector of query longitudes
   * @param query_lats Vector of query latitudes
   * @param radius Search radius in kilometers
   * @return std::vector<GeoQueryResult> Vector of query results for each
   * location
   */
  std::vector<GeoQueryResult> findGridPointsInRadiusBatch(
      const std::vector<double>& query_lons,
      const std::vector<double>& query_lats, double radius) const;

  /**
   * @brief Get the vertical levels from the geometry
   *
   * @return Vector of vertical level depths (meters)
   */
  const std::vector<double>& getVerticalLevels() const { return rC_z_; }

  /**
   * @brief Find the nearest vertical grid points for a given depth
   *
   * @param depth Target depth (meters, positive downward)
   * @return VerticalPoint Containing indices and interpolation information
   */
  VerticalPoint findNearestVerticalPoints(double depth) const;

  /**
   * @brief Find the nearest vertical grid points for multiple depths
   *
   * @param depths Vector of target depths (meters, positive downward)
   * @return std::vector<VerticalPoint> Vector of vertical interpolation
   * information
   */
  std::vector<VerticalPoint> findNearestVerticalPointsBatch(
      const std::vector<double>& depths) const;

 private:
  /**
   * @brief Load grid dimensions from NetCDF file
   *
   * @param ncFile NetCDF file handle
   */
  void loadGridDimensions(netCDF::NcFile& ncFile);

  /**
   * @brief Load grid arrays from NetCDF file
   *
   * @param ncFile NetCDF file handle
   */
  void loadGridArrays(netCDF::NcFile& ncFile);

  /**
   * @brief Initialize grid from NetCDF file
   *
   * @param filename Path to the NetCDF grid file
   */
  void loadGeometryData(const std::string& filename);

  /**
   * @brief Special constructor for cloning
   *
   * @param nlpb Number of grid points
   * @param nk Number of vertical levels
   * @param init Initialization status
   */

  // Grid dimensions
  std::size_t nlpb_ = 0;   // Number of grid points
  std::size_t nk_ = 0;     // Number of vertical levels
  std::size_t nkp1_ = 0;   // Number of vertical levels + 1
  std::size_t nl_ = 0;     // Total number of grid points
  std::size_t nlpbz_ = 0;  // Number of grid points for vorticity
  std::size_t nlz_ = 0;    // Total number of grid points for vorticity
  std::size_t nlbdy_ = 0;  // Number of boundary points
  std::size_t ni_ = 0;     // Number of iterations
  bool initialized_ = false;
  bool kdtree_initialized_ = false;
  const ConfigBackend* config_ptr_ = nullptr;

  // Grid data
  std::vector<double> latC_;  // Latitude at tracer centers
  std::vector<double> lonC_;  // Longitude at tracer centers
  std::vector<double> latW_;  // Latitude of u point
  std::vector<double> lonW_;  // Longitude of u point
  std::vector<double> latS_;  // Latitude of v point
  std::vector<double> lonS_;  // Longitude of v point
  // std::vector<double> latZ_;  // Latitude of vorticity points
  // std::vector<double> lonZ_;  // Longitude of vorticity points

  // std::vector<double> tw_;  // Latitude at west edges
  // std::vector<double> te_;  // Latitude at east edges
  // std::vector<double> tn_;  // Latitude at north edges
  // std::vector<double> ts_;  // Latitude at south edges

  std::vector<double> rC_z_;  // Cell thickness at tracer centers
  std::vector<double> rF_z_;  // Cell thickness at vorticity points

  // std::vector<double> dxC_;    // Grid spacing in x direction at cell centers
  // std::vector<double> dyC_;    // Grid spacing in y direction at cell centers
  // std::vector<double> dxW_;    // Grid spacing in x direction at west edges
  // std::vector<double> dyW_;    // Grid spacing in y direction at west edges
  // std::vector<double> dxS_;    // Grid spacing in x direction at south edges
  // std::vector<double> dyS_;    // Grid spacing in y direction at south edges
  // std::vector<double> rAc_;    // Cell areas
  // std::vector<double> rAw_;    // Areas at west edges
  // std::vector<double> rAs_;    // Areas at south edges
  // std::vector<double> maskC_;  // Land/sea mask at cell centers
  // std::vector<double> maskW_;  // Land/sea mask at west edges
  // std::vector<double> maskS_;  // Land/sea mask at south edges
  // std::vector<double> hFacC_;  // Cell thickness factors
  // std::vector<double> hFacW_;  // West edge thickness factors
  // std::vector<double> hFacS_;  // South edge thickness factors

  // Friend declaration for iterator
  friend class MACOMGeometryIterator<ConfigBackend>;
  friend class MACOMGeometryConstIterator<ConfigBackend>;

  // KD-tree related members and methods
  // using KDTreeType = nanoflann::KDTreeSingleIndexAdaptor<
  //     nanoflann::L2_Simple_Adaptor<double, PointCloud<GeoPoint>>,
  //     PointCloud<GeoPoint>, 2>;
  // void initializeKDTree();
  // PointCloud<GeoPoint> pointCloud_;
  // std::unique_ptr<KDTreeType> kdtree_;
  // bool kdtree_initialized_ = false;
};

}  // namespace metada::backends::macom

// ===============================================================
// Template Implementation Section
// ===============================================================

#include "MACOMGeometryIterator.hpp"

namespace metada::backends::macom {

// ConfigBackend constructor implementation
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend>::MACOMGeometry(const ConfigBackend& config)
    : initialized_(false), config_ptr_(&config) {
  std::string input_filename = config.Get("input_file").asString();
  if (input_filename.empty()) {
    throw std::runtime_error(
        "MACOM input grid file path not specified in configuration");
  }

  // Load geometry data from MACOM NetCDF file
  loadGeometryData(input_filename);
  initialized_ = true;

  // // Add debug output
  // debugGridDataHorizontal();
  // debugGridDataVertical();

  // Initialize the static KD-tree in iterator (only initializes once)
  iterator::initializeKDTree(lonC_, latC_, nlpb_);

  // double test_lon = 160.0, test_lat = 36.0;
  // GeoPoint kd_result = findNearestGridPoint(test_lon, test_lat);
  // GeoPoint direct_result = findNearestGridPointDirect(test_lon, test_lat);

  // MACOM_LOG_INFO("MACOMGeometry", "Comparison for point (" +
  //                                     std::to_string(test_lon) + ", " +
  //                                     std::to_string(test_lat) + "):");
  // MACOM_LOG_INFO("MACOMGeometry",
  //                "  KD-tree result: index=" + std::to_string(kd_result.index)
  //                +
  //                    ", distance=" + std::to_string(kd_result.distance) +
  //                    " km");
  // MACOM_LOG_INFO(
  //     "MACOMGeometry",
  //     "  Direct result: index=" + std::to_string(direct_result.index) +
  //         ", distance=" + std::to_string(direct_result.distance) + " km");

  // double test_depth = 100.0;  // Test at 1000 meters
  // VerticalPoint v_result = findNearestVerticalPoints(test_depth);

  // MACOM_LOG_INFO("MACOMGeometry", "Vertical interpolation test for depth " +
  //                                     std::to_string(test_depth) + "
  //                                     meters:");
  // MACOM_LOG_INFO("MACOMGeometry",
  //                "  Lower level: " + std::to_string(v_result.lower_index) +
  //                    " at " + std::to_string(v_result.lower_depth) + "
  //                    meters");
  // MACOM_LOG_INFO("MACOMGeometry",
  //                "  Upper level: " + std::to_string(v_result.upper_index) +
  //                    " at " + std::to_string(v_result.upper_depth) + "
  //                    meters");
  // MACOM_LOG_INFO("MACOMGeometry", "  Interpolation coefficient: " +
  //                                     std::to_string(v_result.interp_coef));
  // MACOM_LOG_INFO("MACOMGeometry",
  //                "  Is outside range: " +
  //                    std::string(v_result.is_outside ? "true" : "false"));
}

// Implementation of loadGridDimensions
template <typename ConfigBackend>
void MACOMGeometry<ConfigBackend>::loadGridDimensions(netCDF::NcFile& ncFile) {
  auto readDimension = [&ncFile](const std::string& name, std::size_t& value) {
    // MACOM_LOG_INFO("MACOMGeometry", "Attempting to read dimension: " + name);
    auto dim = ncFile.getDim(name);
    if (dim.isNull()) {
      MACOM_LOG_ERROR("MACOMGeometry",
                      "Dimension '" + name + "' not found in grid file");
      throw std::runtime_error("Dimension '" + name +
                               "' not found in grid file");
    }
    // MACOM_LOG_INFO("MACOMGeometry", "Found dimension: " + name);
    value = dim.getSize();
    // MACOM_LOG_INFO("MACOMGeometry",
    //                "Read size for " + name + ": " + std::to_string(value));
  };

  // Read all grid dimensions
  readDimension("nlpb", nlpb_);
  readDimension("nl", nl_);
  readDimension("nlpbz", nlpbz_);
  readDimension("nlz", nlz_);
  readDimension("nk", nk_);
  readDimension("nkp1", nkp1_);
  readDimension("ni", ni_);

  MACOM_LOG_INFO("MACOMGeometry", "Loaded grid dimensions:");
  MACOM_LOG_INFO("MACOMGeometry", "  nlpb=" + std::to_string(nlpb_) +
                                      ", nk=" + std::to_string(nk_) +
                                      ", nkp1=" + std::to_string(nkp1_));
  // MACOM_LOG_INFO("MACOMGeometry", "  nl=" + std::to_string(nl_) +
  //                                     ", nlpbz=" + std::to_string(nlpbz_) +
  //                                     ", nlz=" + std::to_string(nlz_));
  // MACOM_LOG_INFO("MACOMGeometry", "  ni=" + std::to_string(ni_));
}

// Implementation of loadGridArrays
template <typename ConfigBackend>
void MACOMGeometry<ConfigBackend>::loadGridArrays(netCDF::NcFile& ncFile) {
  // Resize vectors
  latC_.resize(nlpb_);
  lonC_.resize(nlpb_);
  latW_.resize(nlpb_);
  lonW_.resize(nlpb_);
  latS_.resize(nlpb_);
  lonS_.resize(nlpb_);
  // latZ_.resize(nlpb_);
  // lonZ_.resize(nlpb_);

  // tw_.resize(nlpb_);
  // te_.resize(nlpb_);
  // tn_.resize(nlpb_);
  // ts_.resize(nlpb_);

  rC_z_.resize(nk_);
  rF_z_.resize(nkp1_);

  // dxC_.resize(nlpb_);
  // dyC_.resize(nlpb_);
  // dxW_.resize(nlpb_);
  // dyW_.resize(nlpb_);
  // dxS_.resize(nlpb_);
  // dyS_.resize(nlpb_);
  // rAc_.resize(nlpb_);
  // rAw_.resize(nlpb_);
  // rAs_.resize(nlpb_);
  // maskC_.resize(nlpb_ * nk_);
  // maskW_.resize(nlpb_ * nk_);
  // maskS_.resize(nlpb_ * nk_);
  // hFacC_.resize(nlpb_ * nk_);
  // hFacW_.resize(nlpb_ * nk_);
  // hFacS_.resize(nlpb_ * nk_);

  // Read variables
  auto readVar = [&ncFile](const std::string& name, std::vector<double>& data) {
    auto var = ncFile.getVar(name);
    if (var.isNull()) {
      throw std::runtime_error("Variable '" + name +
                               "' not found in grid file");
    }
    var.getVar(data.data());
  };

  // Read grid data
  readVar("lat_c", latC_);
  readVar("lon_c", lonC_);
  readVar("lat_w", latW_);
  readVar("lon_w", lonW_);
  readVar("lat_s", latS_);
  readVar("lon_s", lonS_);
  // readVar("lat_z", latZ_);
  // readVar("lon_z", lonZ_);

  // readVar("tw", tw_);
  // readVar("te", te_);
  // readVar("tn", tn_);
  // readVar("ts", ts_);

  readVar("rC_z", rC_z_);
  readVar("rF_z", rF_z_);

  // readVar("dxC", dxC_);
  // readVar("dyC", dyC_);
  // readVar("dxW", dxW_);
  // readVar("dyW", dyW_);
  // readVar("dxS", dxS_);
  // readVar("dyS", dyS_);
  // readVar("rAc", rAc_);
  // readVar("rAw", rAw_);
  // readVar("rAs", rAs_);
  // readVar("maskC", maskC_);
  // readVar("maskW", maskW_);
  // readVar("maskS", maskS_);
  // readVar("h0FacC", hFacC_);
  // readVar("h0FacW", hFacW_);
  // readVar("h0FacS", hFacS_);

  // MACOM_LOG_INFO("MACOMGeometry", "Loaded grid arrays for " +
  //                                     std::to_string(nlpb_) +
  //                                     " grid points and " +
  //                                     std::to_string(nk_) + " vertical
  //                                     levels");
}

// Implementation of loadGeometryData
template <typename ConfigBackend>
void MACOMGeometry<ConfigBackend>::loadGeometryData(
    const std::string& filename) {
  try {
    // MACOM_LOG_INFO("MACOMGeometry",
    //                "Attempting to open NetCDF file: " + filename);

    // Open the NetCDF file
    netCDF::NcFile ncFile(filename, netCDF::NcFile::read);

    if (ncFile.isNull()) {
      // MACOM_LOG_ERROR("MACOMGeometry",
      //                 "Failed to open NetCDF file: " + filename);
      throw std::runtime_error("Failed to open NetCDF file: " + filename);
    }

    // MACOM_LOG_INFO("MACOMGeometry", "Successfully opened NetCDF file");

    // Step 1: Load grid dimensions
    loadGridDimensions(ncFile);

    // Step 2: Load grid arrays
    loadGridArrays(ncFile);

    // Close the file
    ncFile.close();

    initialized_ = true;

    // MACOM_LOG_INFO("MACOMGeometry",
    //                "Successfully initialized hexagonal grid from " +
    //                filename);

  } catch (const netCDF::exceptions::NcException& e) {
    MACOM_LOG_ERROR("MACOMGeometry", "NetCDF error while reading grid file: " +
                                         std::string(e.what()));
    throw std::runtime_error("NetCDF error while reading grid file: " +
                             std::string(e.what()));
  } catch (const std::exception& e) {
    MACOM_LOG_ERROR("MACOMGeometry",
                    "Error initializing grid: " + std::string(e.what()));
    throw std::runtime_error("Error initializing grid: " +
                             std::string(e.what()));
  }
}

// Find nearest point implementation with improved search
template <typename ConfigBackend>
GeoPoint MACOMGeometry<ConfigBackend>::findNearestGridPoint(double lon,
                                                            double lat) const {
  if (!initialized_) {
    MACOM_LOG_ERROR("MACOMGeometry",
                    "Cannot find nearest point: Geometry not initialized");
    throw std::runtime_error("Geometry not initialized");
  }

  // Call iterator's static method with default parameter (false for multiple
  // neighbors) The user can override this by modifying the call directly as
  // needed
  return iterator::findNearestGridPoint(lonC_, latC_, nlpb_, lon, lat, false);
}

/**
 * @brief Find all grid points within the specified radius of the given
 * coordinates
 *
 * @param lon Longitude (degrees)
 * @param lat Latitude (degrees)
 * @param radius Search radius (kilometers)
 * @return GeoQueryResult Containing indices and coordinates of all points
 * within radius
 */
template <typename ConfigBackend>
GeoQueryResult MACOMGeometry<ConfigBackend>::findGridPointsInRadius(
    double lon, double lat, double radius) const {
  if (!initialized_) {
    MACOM_LOG_ERROR("MACOMGeometry",
                    "Cannot find points in radius: Geometry not initialized");
    throw std::runtime_error("Geometry not initialized");
  }

  // Call iterator's static method
  return iterator::findGridPointsInRadius(lonC_, latC_, nlpb_, lon, lat,
                                          radius);
}

// Alternative implementation for direct nearest point search using Haversine
template <typename ConfigBackend>
GeoPoint MACOMGeometry<ConfigBackend>::findNearestGridPointDirect(
    double lon, double lat) const {
  if (!initialized_) {
    throw std::runtime_error("Geometry not initialized");
  }

  // Call iterator's static method
  return iterator::findNearestGridPointDirect(lonC_, latC_, nlpb_, lon, lat);
}

// Debug grid data implementation
template <typename ConfigBackend>
void MACOMGeometry<ConfigBackend>::debugGridDataHorizontal() {
  if (!initialized_ || nlpb_ == 0) {
    MACOM_LOG_WARNING(
        "MACOMGeometry",
        "Cannot debug horizontal grid data: Grid not initialized or empty");
    return;
  }

  // Find min/max values
  double min_lat = latC_[0], max_lat = latC_[0];
  double min_lon = lonC_[0], max_lon = lonC_[0];

  for (size_t i = 0; i < nlpb_; ++i) {
    min_lat = std::min(min_lat, latC_[i]);
    max_lat = std::max(max_lat, latC_[i]);
    min_lon = std::min(min_lon, lonC_[i]);
    max_lon = std::max(max_lon, lonC_[i]);
  }

  MACOM_LOG_INFO("MACOMGeometry", "Horizontal grid data statistics:");
  MACOM_LOG_INFO("MACOMGeometry",
                 "  Latitude range: " + std::to_string(min_lat) + " to " +
                     std::to_string(max_lat));
  MACOM_LOG_INFO("MACOMGeometry",
                 "  Longitude range: " + std::to_string(min_lon) + " to " +
                     std::to_string(max_lon));

  // Print some sample points
  MACOM_LOG_INFO("MACOMGeometry", "Sample horizontal grid points:");
  size_t step = nlpb_ / 10;  // Print 10 points evenly distributed
  if (step == 0) step = 1;

  for (size_t i = 0; i < nlpb_ && i < 10 * step; i += step) {
    MACOM_LOG_INFO("MACOMGeometry", "  Point " + std::to_string(i) +
                                        ": lon=" + std::to_string(lonC_[i]) +
                                        ", lat=" + std::to_string(latC_[i]));
  }
}

// Implementation of debugGridDataVertical
template <typename ConfigBackend>
void MACOMGeometry<ConfigBackend>::debugGridDataVertical() {
  if (!initialized_ || nk_ == 0) {
    MACOM_LOG_WARNING(
        "MACOMGeometry",
        "Cannot debug vertical grid data: Grid not initialized or empty");
    return;
  }

  // Find min/max values
  double min_depth = rC_z_[0], max_depth = rC_z_[0];
  double total_depth = 0.0;

  for (size_t i = 0; i < nk_; ++i) {
    min_depth = std::min(min_depth, rC_z_[i]);
    max_depth = std::max(max_depth, rC_z_[i]);
    total_depth += rC_z_[i];
  }

  MACOM_LOG_INFO("MACOMGeometry", "Vertical grid data statistics:");
  MACOM_LOG_INFO("MACOMGeometry",
                 "  Number of vertical levels: " + std::to_string(nk_));
  MACOM_LOG_INFO("MACOMGeometry",
                 "  Depth range: " + std::to_string(min_depth) + " to " +
                     std::to_string(max_depth) + " meters");
  MACOM_LOG_INFO("MACOMGeometry", "  Average level thickness: " +
                                      std::to_string(total_depth / nk_) +
                                      " meters");

  // Print all vertical levels
  MACOM_LOG_INFO("MACOMGeometry", "Vertical levels (depths in meters):");
  for (size_t i = 0; i < nk_; ++i) {
    MACOM_LOG_INFO("MACOMGeometry", "  Level " + std::to_string(i) +
                                        ": depth=" + std::to_string(rC_z_[i]));
  }
}

// Move constructor implementation
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend>::MACOMGeometry(MACOMGeometry&& other) noexcept
    : nlpb_(other.nlpb_),
      nk_(other.nk_),
      nkp1_(other.nkp1_),
      nl_(other.nl_),
      nlpbz_(other.nlpbz_),
      nlz_(other.nlz_),
      nlbdy_(other.nlbdy_),
      ni_(other.ni_),
      initialized_(other.initialized_),
      kdtree_initialized_(other.kdtree_initialized_),
      config_ptr_(other.config_ptr_),
      latC_(std::move(other.latC_)),
      lonC_(std::move(other.lonC_)) {
  other.initialized_ = false;
  other.nlpb_ = 0;
  other.nk_ = 0;
  other.config_ptr_ = nullptr;
}

// Move assignment operator implementation
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend>& MACOMGeometry<ConfigBackend>::operator=(
    MACOMGeometry&& other) noexcept {
  if (this != &other) {
    nlpb_ = other.nlpb_;
    nk_ = other.nk_;
    nkp1_ = other.nkp1_;
    nl_ = other.nl_;
    nlpbz_ = other.nlpbz_;
    nlz_ = other.nlz_;
    nlbdy_ = other.nlbdy_;
    ni_ = other.ni_;
    initialized_ = other.initialized_;
    kdtree_initialized_ = other.kdtree_initialized_;
    config_ptr_ = other.config_ptr_;
    latC_ = std::move(other.latC_);
    lonC_ = std::move(other.lonC_);

    other.initialized_ = false;
    other.nlpb_ = 0;
    other.nk_ = 0;
    other.config_ptr_ = nullptr;
  }
  return *this;
}

// Clone implementation
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend> MACOMGeometry<ConfigBackend>::clone() const {
  if (!initialized_) {
    throw std::runtime_error("Cannot clone uninitialized geometry");
  }

  if (!config_ptr_) {
    throw std::runtime_error("Cannot clone: config pointer is null");
  }

  return MACOMGeometry(*config_ptr_);
}

// Batch spatial query implementations
template <typename ConfigBackend>
std::vector<GeoPoint> MACOMGeometry<ConfigBackend>::findNearestGridPointsBatch(
    const std::vector<double>& query_lons,
    const std::vector<double>& query_lats) const {
  if (!initialized_) {
    MACOM_LOG_ERROR("MACOMGeometry",
                    "Cannot find nearest points: Geometry not initialized");
    throw std::runtime_error("Geometry not initialized");
  }

  return iterator::findNearestGridPointsBatch(lonC_, latC_, nlpb_, query_lons,
                                              query_lats);
}

template <typename ConfigBackend>
std::vector<GeoQueryResult>
MACOMGeometry<ConfigBackend>::findGridPointsInRadiusBatch(
    const std::vector<double>& query_lons,
    const std::vector<double>& query_lats, double radius) const {
  if (!initialized_) {
    MACOM_LOG_ERROR("MACOMGeometry",
                    "Cannot find points in radius: Geometry not initialized");
    throw std::runtime_error("Geometry not initialized");
  }

  return iterator::findGridPointsInRadiusBatch(lonC_, latC_, nlpb_, query_lons,
                                               query_lats, radius);
}

template <typename ConfigBackend>
VerticalPoint MACOMGeometry<ConfigBackend>::findNearestVerticalPoints(
    double depth) const {
  if (!initialized_) {
    throw std::runtime_error("Geometry not initialized");
  }
  return iterator::findNearestVerticalPoints(rC_z_, nk_, depth);
}

template <typename ConfigBackend>
std::vector<VerticalPoint>
MACOMGeometry<ConfigBackend>::findNearestVerticalPointsBatch(
    const std::vector<double>& depths) const {
  if (!initialized_) {
    throw std::runtime_error("Geometry not initialized");
  }
  return iterator::findNearestVerticalPointsBatch(rC_z_, nk_, depths);
}

// Iterator methods
template <typename ConfigBackend>
inline typename MACOMGeometry<ConfigBackend>::iterator
MACOMGeometry<ConfigBackend>::begin() {
  return iterator(this, 0);
}

template <typename ConfigBackend>
inline typename MACOMGeometry<ConfigBackend>::iterator
MACOMGeometry<ConfigBackend>::end() {
  return iterator(this, totalGridSize());
}

template <typename ConfigBackend>
inline typename MACOMGeometry<ConfigBackend>::const_iterator
MACOMGeometry<ConfigBackend>::begin() const {
  return const_iterator(this, 0);
}

template <typename ConfigBackend>
inline typename MACOMGeometry<ConfigBackend>::const_iterator
MACOMGeometry<ConfigBackend>::end() const {
  return const_iterator(this, totalGridSize());
}

template <typename ConfigBackend>
inline std::size_t MACOMGeometry<ConfigBackend>::totalGridSize() const {
  return nlpb_ * nk_;
}

}  // namespace metada::backends::macom