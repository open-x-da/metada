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
   * @brief Print statistics about the grid data to help with debugging
   */
  void debugGridData() {
    if (!initialized_ || nlpb_ == 0) {
      MACOM_LOG_WARNING(
          "MACOMGeometry",
          "Cannot debug grid data: Grid not initialized or empty");
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

    MACOM_LOG_INFO("MACOMGeometry", "Grid data statistics:");
    MACOM_LOG_INFO("MACOMGeometry",
                   "  Latitude range: " + std::to_string(min_lat) + " to " +
                       std::to_string(max_lat));
    MACOM_LOG_INFO("MACOMGeometry",
                   "  Longitude range: " + std::to_string(min_lon) + " to " +
                       std::to_string(max_lon));

    // Print some sample points
    MACOM_LOG_INFO("MACOMGeometry", "Sample grid points:");
    size_t step = nlpb_ / 10;  // Print 10 points evenly distributed
    if (step == 0) step = 1;

    for (size_t i = 0; i < nlpb_ && i < 10 * step; i += step) {
      MACOM_LOG_INFO("MACOMGeometry", "  Point " + std::to_string(i) +
                                          ": lon=" + std::to_string(lonC_[i]) +
                                          ", lat=" + std::to_string(latC_[i]));
    }
  }

  // Alternative implementation for direct nearest point search using Haversine
  GeoPoint findNearestGridPointDirect(double lon, double lat) const;

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
  std::vector<double> latC_;  // Latitude at cell centers
  std::vector<double> lonC_;  // Longitude at cell centers

  std::vector<double> tw_;  // Latitude at west edges
  std::vector<double> te_;  // Latitude at east edges
  std::vector<double> tn_;  // Latitude at north edges
  std::vector<double> ts_;  // Latitude at south edges

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
  using KDTreeType = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, PointCloud<GeoPoint>>,
      PointCloud<GeoPoint>, 2 /* dimensions - changed from 3 to 2 */
      >;

  /**
   * @brief Initialize the KD-tree
   */
  void initializeKDTree();

  /**
   * @brief Convert latitude and longitude to 2D Cartesian coordinates
   *
   * @param lon Longitude (degrees)
   * @param lat Latitude (degrees)
   * @param x Output parameter, Cartesian coordinate x
   * @param y Output parameter, Cartesian coordinate y
   */
  static void geoToCartesian(double lon, double lat, double& x, double& y);

  /**
   * @brief Calculate great-circle distance between two points
   *
   * @param lon1 First point longitude
   * @param lat1 First point latitude
   * @param lon2 Second point longitude
   * @param lat2 Second point latitude
   * @return double Distance (kilometers)
   */
  static double haversineDistance(double lon1, double lat1, double lon2,
                                  double lat2);

  // KD-tree members
  PointCloud<GeoPoint> pointCloud_;
  std::unique_ptr<KDTreeType> kdtree_;
  static constexpr double EARTH_RADIUS_KM = 6371.0;
};

// ConfigBackend constructor implementation
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend>::MACOMGeometry(const ConfigBackend& config)
    : initialized_(false), kdtree_initialized_(false) {
  std::string input_filename = config.Get("input_file").asString();
  if (input_filename.empty()) {
    throw std::runtime_error(
        "MACOM input grid file path not specified in configuration");
  }

  // Load geometry data from MACOM NetCDF file
  loadGeometryData(input_filename);
  initialized_ = true;

  // // Add debug output
  // debugGridData();

  // KD-tree can be initialized here, or deferred until the first query
  initializeKDTree();
  // kdtree_initialized_ = true;

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

  tw_.resize(nlpb_);
  te_.resize(nlpb_);
  tn_.resize(nlpb_);
  ts_.resize(nlpb_);

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

  readVar("tw", tw_);
  readVar("te", te_);
  readVar("tn", tn_);
  readVar("ts", ts_);

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

// Better implementation of geoToCartesian for 2D coordinates
template <typename ConfigBackend>
void MACOMGeometry<ConfigBackend>::geoToCartesian(double lon, double lat,
                                                  double& x, double& y) {
  // Use the Mercator projection which better preserves distance relationships
  // Normalize longitude to -180 to 180 for consistency
  if (lon > 180.0) lon -= 360.0;

  // Convert to radians
  double lon_rad = lon * std::numbers::pi / 180.0;
  double lat_rad = std::clamp(lat, -85.0, 85.0) * std::numbers::pi /
                   180.0;  // Clamp latitude to avoid singularities

  // Mercator projection (scaled by Earth's radius)
  x = EARTH_RADIUS_KM * lon_rad;
  y = EARTH_RADIUS_KM * log(tan(std::numbers::pi / 4.0 + lat_rad / 2.0));
}

// Calculate the great-circle distance between two points
template <typename ConfigBackend>
double MACOMGeometry<ConfigBackend>::haversineDistance(double lon1, double lat1,
                                                       double lon2,
                                                       double lat2) {
  // Normalize longitude to 0-360
  lon1 = fmod(lon1 + 360.0, 360.0);
  lon2 = fmod(lon2 + 360.0, 360.0);

  // Convert to radians
  double lon1_rad = lon1 * std::numbers::pi / 180.0;
  double lat1_rad = lat1 * std::numbers::pi / 180.0;
  double lon2_rad = lon2 * std::numbers::pi / 180.0;
  double lat2_rad = lat2 * std::numbers::pi / 180.0;

  // Haversine equation
  double dlon = lon2_rad - lon1_rad;
  double dlat = lat2_rad - lat1_rad;
  double a = sin(dlat / 2) * sin(dlat / 2) +
             cos(lat1_rad) * cos(lat2_rad) * sin(dlon / 2) * sin(dlon / 2);
  double c = 2 * atan2(sqrt(a), sqrt(1 - a));
  return EARTH_RADIUS_KM * c;
}

// KD-tree initialization implementation
template <typename ConfigBackend>
void MACOMGeometry<ConfigBackend>::initializeKDTree() {
  if (kdtree_initialized_) {
    return;
  }

  MACOM_LOG_INFO("MACOMGeometry",
                 "Initializing 2D KD-tree for spatial queries");

  // Clear previous data
  pointCloud_.pts.clear();

  // Fill point cloud data
  pointCloud_.pts.resize(nlpb_);
  for (size_t i = 0; i < nlpb_; ++i) {
    GeoPoint& point = pointCloud_.pts[i];
    point.index = i;
    point.lon = lonC_[i];
    point.lat = latC_[i];

    // Convert to 2D Cartesian coordinates
    geoToCartesian(point.lon, point.lat, point.x, point.y);
  }

  MACOM_LOG_INFO("MACOMGeometry", "Point cloud data filled");

  // Build KD-tree
  const int leaf_size =
      10;  // Adjust this value to affect tree balance and query performance
  kdtree_.reset(new KDTreeType(
      2, pointCloud_, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_size)));
  kdtree_->buildIndex();

  kdtree_initialized_ = true;
  MACOM_LOG_INFO("MACOMGeometry", "2D KD-tree initialized with " +
                                      std::to_string(nlpb_) + " points");
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

  if (!kdtree_initialized_) {
    const_cast<MACOMGeometry*>(this)->initializeKDTree();
  }

  // Normalize longitude to match dataset orientation (if needed)
  if (lon > 180.0) lon -= 360.0;

  // First, try direct search with KD-tree
  double query_x, query_y;
  geoToCartesian(lon, lat, query_x, query_y);

  // We'll consider multiple points and choose the closest by actual great
  // circle distance
  const size_t knn = 5;  // Check 5 nearest neighbors
  std::vector<size_t> indices(knn);
  std::vector<double> distances_sq(knn);

  nanoflann::KNNResultSet<double> resultSet(knn);
  resultSet.init(indices.data(), distances_sq.data());

  double query_pt[2] = {query_x, query_y};
  kdtree_->findNeighbors(resultSet, query_pt, nanoflann::SearchParameters());

  // Find closest point by actual great circle distance
  GeoPoint result;
  double min_distance = std::numeric_limits<double>::max();

  // Log candidate points for debugging
  MACOM_LOG_INFO("MACOMGeometry", "Candidates for nearest point to (" +
                                      std::to_string(lon) + ", " +
                                      std::to_string(lat) + "):");

  for (size_t i = 0; i < resultSet.size(); ++i) {
    size_t idx = indices[i];
    double haversine_dist = haversineDistance(lon, lat, lonC_[idx], latC_[idx]);

    MACOM_LOG_INFO(
        "MACOMGeometry",
        "  Candidate " + std::to_string(i) + ": index=" + std::to_string(idx) +
            " lon=" + std::to_string(lonC_[idx]) +
            " lat=" + std::to_string(latC_[idx]) +
            " Euclidean dist=" + std::to_string(sqrt(distances_sq[i])) +
            " Haversine dist=" + std::to_string(haversine_dist) + " km");

    if (haversine_dist < min_distance) {
      min_distance = haversine_dist;
      result.index = idx;
      result.lon = lonC_[idx];
      result.lat = latC_[idx];
      result.distance = haversine_dist;
    }
  }

  // Check if a point was found
  if (min_distance < std::numeric_limits<double>::max()) {
    std::string logMsg =
        "Best nearest point to (" + std::to_string(lon) + ", " +
        std::to_string(lat) + "): index=" + std::to_string(result.index) +
        ", coords=(" + std::to_string(result.lon) + ", " +
        std::to_string(result.lat) + ")" +
        ", distance=" + std::to_string(result.distance) + " km";
    MACOM_LOG_INFO("MACOMGeometry", logMsg);
  } else {
    // Fall back to brute force search if no point found
    MACOM_LOG_WARNING("MACOMGeometry",
                      "KD-tree search failed, trying brute force search");

    for (size_t i = 0; i < nlpb_; ++i) {
      double dist = haversineDistance(lon, lat, lonC_[i], latC_[i]);
      if (dist < min_distance) {
        min_distance = dist;
        result.index = i;
        result.lon = lonC_[i];
        result.lat = latC_[i];
        result.distance = dist;
      }
    }

    if (min_distance < std::numeric_limits<double>::max()) {
      MACOM_LOG_INFO(
          "MACOMGeometry",
          "Found point by brute force: index=" + std::to_string(result.index) +
              ", distance=" + std::to_string(result.distance) + " km");
    } else {
      MACOM_LOG_WARNING("MACOMGeometry", "No points found (unexpected)");
    }
  }

  return result;
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
  GeoQueryResult result;

  using IndexType = typename KDTreeType::IndexType;
  using DistanceType = typename KDTreeType::DistanceType;
  using MatchItem = nanoflann::ResultItem<IndexType, DistanceType>;

  if (!initialized_) {
    MACOM_LOG_ERROR("MACOMGeometry",
                    "Cannot find points in radius: Geometry not initialized");
    throw std::runtime_error("Geometry not initialized");
  }

  if (!kdtree_initialized_) {
    // KD-tree not initialized yet, initialize it first
    const_cast<MACOMGeometry*>(this)->initializeKDTree();
  }

  // Convert query point to 2D Cartesian coordinates
  double query_x, query_y;
  geoToCartesian(lon, lat, query_x, query_y);

  // In Cartesian space, the search radius needs conversion
  // Since we're using 2D Cartesian coordinates, the search radius needs to be
  // slightly larger to ensure we capture all relevant points
  double search_radius =
      radius *
      1.2;  // Add 20% to ensure all potentially in-range points are included

  // Prepare query point
  double query_pt[2] = {query_x, query_y};  // 2D query point

  // This fixes the vector type error - using std::pair<size_t, double> instead
  // of ResultItem
  std::vector<MatchItem> matches;

  // Create search params
  nanoflann::SearchParameters params;
  params.sorted = true;  // Sort results by distance

  // Perform radius search
  const size_t num_matches = kdtree_->radiusSearch(
      query_pt, search_radius * search_radius, matches, params);

  // Construct the result
  if (num_matches > 0) {
    result.points.reserve(num_matches);

    for (size_t i = 0; i < num_matches; ++i) {
      size_t idx = matches[i].first;
      GeoPoint point = pointCloud_.pts[idx];

      // Calculate precise great circle distance
      double exact_distance = haversineDistance(lon, lat, point.lon, point.lat);

      // Only keep points truly within range
      if (exact_distance <= radius) {
        point.distance = exact_distance;
        result.points.push_back(point);
      }
    }

    // Sort by distance
    std::sort(result.points.begin(), result.points.end(),
              [](const GeoPoint& a, const GeoPoint& b) {
                return a.distance < b.distance;
              });

    std::string logMsg = "Found " + std::to_string(result.points.size()) +
                         " points within " + std::to_string(radius) +
                         " km of (" + std::to_string(lon) + ", " +
                         std::to_string(lat) + ")";
    MACOM_LOG_INFO("MACOMGeometry", logMsg);
  } else {
    std::string logMsg = "No points found within " + std::to_string(radius) +
                         " km of (" + std::to_string(lon) + ", " +
                         std::to_string(lat) + ")";
    MACOM_LOG_INFO("MACOMGeometry", logMsg);
  }

  return result;
}

// Alternative implementation for direct nearest point search using Haversine
template <typename ConfigBackend>
GeoPoint MACOMGeometry<ConfigBackend>::findNearestGridPointDirect(
    double lon, double lat) const {
  if (!initialized_) {
    throw std::runtime_error("Geometry not initialized");
  }

  GeoPoint result;
  double min_distance = std::numeric_limits<double>::max();

  // Brute force search using Haversine distance
  for (size_t i = 0; i < nlpb_; ++i) {
    double dist = haversineDistance(lon, lat, lonC_[i], latC_[i]);
    if (dist < min_distance) {
      min_distance = dist;
      result.index = i;
      result.lon = lonC_[i];
      result.lat = latC_[i];
      result.distance = dist;
    }
  }

  if (min_distance < std::numeric_limits<double>::max()) {
    MACOM_LOG_INFO(
        "MACOMGeometry",
        "Direct search found point: index=" + std::to_string(result.index) +
            ", coords=(" + std::to_string(result.lon) + ", " +
            std::to_string(result.lat) + ")" +
            ", distance=" + std::to_string(result.distance) + " km");
  }

  return result;
}

}  // namespace metada::backends::macom

// Include the iterator implementation first
#include "MACOMGeometryIterator.hpp"

// Then define the iterator methods
namespace metada::backends::macom {

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

}  // namespace metada::backends::macom