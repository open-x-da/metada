/**
 * @file MACOMGeometryIterator.hpp
 * @brief MACOM geometry iterator implementation with KD-tree functionality
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <algorithm>  // For std::clamp support
#include <cmath>      // For math functions
#include <cstddef>
#include <memory>   // For std::unique_ptr
#include <mutex>    // For mutex support
#include <numbers>  // For std::numbers
#include <stdexcept>
#include <tuple>
#include <vector>

#include "include/MACOMlogging.hpp"
#include "nanoflann.hpp"

namespace metada::backends::macom {

// Forward declarations
template <typename ConfigBackend>
class MACOMGeometry;

struct GeoPoint;
struct GeoQueryResult;
template <typename T>
struct PointCloud;

// Define iterator class with KD-tree functionality
template <typename ConfigBackend>
class MACOMGeometryIterator {
 public:
  // Iterator traits
  using iterator_category = std::forward_iterator_tag;
  using value_type = std::tuple<size_t, size_t, size_t>;  // i, j, k indices
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  // Earth radius constant used for spatial calculations
  static constexpr double EARTH_RADIUS_KM = 6371.0;

  // KD-tree type definition
  using KDTreeType = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, PointCloud<GeoPoint>>,
      PointCloud<GeoPoint>, 2>;

  /**
   * @brief Default constructor is deleted
   */
  MACOMGeometryIterator() = delete;

  /**
   * @brief Constructor for creating a valid iterator
   *
   * @param geometry Pointer to the MACOM geometry backend
   * @param index Linear index into the grid (0 for begin, grid size for end)
   */
  MACOMGeometryIterator(const MACOMGeometry<ConfigBackend>* geometry,
                        size_t index)
      : geometry_(geometry), index_(index) {
    MACOM_LOG_DEBUG(
        "MACOMGeometryIterator",
        "Created iterator for geometry with index " + std::to_string(index));
  }

  /**
   * @brief Dereference operator
   *
   * @return Tuple containing (i, j, k) indices of current position
   */
  value_type operator*() const { return std::make_tuple(i_, j_, k_); }

  /**
   * @brief Pre-increment operator
   *
   * @return Reference to this iterator after incrementing
   */
  MACOMGeometryIterator& operator++() {
    // TODO: Implement
    return *this;
  }

  /**
   * @brief Post-increment operator
   *
   * @return Copy of iterator before incrementing
   */
  MACOMGeometryIterator operator++(int) {
    MACOMGeometryIterator temp = *this;
    ++(*this);
    return temp;
  }

  /**
   * @brief Equality comparison
   *
   * @param other Iterator to compare with
   * @return True if iterators point to the same position
   */
  bool operator==(const MACOMGeometryIterator& other) const {
    return (geometry_ == other.geometry_ && index_ == other.index_);
  }

  /**
   * @brief Inequality comparison
   *
   * @param other Iterator to compare with
   * @return True if iterators point to different positions
   */
  bool operator!=(const MACOMGeometryIterator& other) const {
    return !(*this == other);
  }

  /**
   * @brief Get i index
   *
   * @return Current i index
   */
  size_t i() const { return i_; }

  /**
   * @brief Get j index
   *
   * @return Current j index
   */
  size_t j() const { return j_; }

  /**
   * @brief Get k index
   *
   * @return Current k index
   */
  size_t k() const { return k_; }

  /**
   * @brief Convert latitude and longitude to 2D Cartesian coordinates
   *
   * @param lon Longitude (degrees)
   * @param lat Latitude (degrees)
   * @param x Output parameter, Cartesian coordinate x
   * @param y Output parameter, Cartesian coordinate y
   */
  static void geoToCartesian(double lon, double lat, double& x, double& y) {
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

  /**
   * @brief Initialize KD-tree (static method, ensures one-time initialization)
   *
   * @param lonC Longitude array
   * @param latC Latitude array
   * @param nlpb Number of grid points
   * @return bool Initialization success
   */
  static bool initializeKDTree(const std::vector<double>& lonC,
                               const std::vector<double>& latC, size_t nlpb) {
    // Use mutex to ensure thread safety
    std::lock_guard<std::mutex> lock(kdtree_mutex_);

    // If already initialized, return immediately
    if (kdtree_initialized_) {
      return true;
    }

    try {
      MACOM_LOG_INFO("MACOMGeometryIterator",
                     "Initializing static KD-tree for spatial queries");

      // Clear previous data
      static_point_cloud_.pts.clear();
      static_point_cloud_.pts.reserve(nlpb);

      // Fill point cloud data
      for (size_t i = 0; i < nlpb; ++i) {
        GeoPoint point;
        point.index = i;
        point.lon = lonC[i];
        point.lat = latC[i];

        // Calculate Cartesian coordinates
        geoToCartesian(point.lon, point.lat, point.x, point.y);

        static_point_cloud_.pts.push_back(point);
      }

      MACOM_LOG_INFO(
          "MACOMGeometryIterator",
          "Point cloud data filled with " + std::to_string(nlpb) + " points");

      // Build KD-tree
      const int leaf_size = 10;
      static_kdtree_.reset(
          new KDTreeType(2, static_point_cloud_,
                         nanoflann::KDTreeSingleIndexAdaptorParams(leaf_size)));
      static_kdtree_->buildIndex();

      kdtree_initialized_ = true;
      MACOM_LOG_INFO("MACOMGeometryIterator",
                     "Static KD-tree successfully initialized");
      return true;

    } catch (const std::exception& e) {
      MACOM_LOG_ERROR("MACOMGeometryIterator",
                      "Failed to initialize KD-tree: " + std::string(e.what()));
      kdtree_initialized_ = false;
      return false;
    }
  }

  /**
   * @brief Find the nearest grid point (using KD-tree)
   *
   * @param lonC Longitude array (for fallback method)
   * @param latC Latitude array (for fallback method)
   * @param nlpb Number of grid points (for fallback method)
   * @param lon Query point longitude
   * @param lat Query point latitude
   * @param use_multiple_neighbors Whether to check multiple neighbors for more
   * accurate results
   * @return GeoPoint Containing index, coordinates and distance of the nearest
   * point
   */
  static GeoPoint findNearestGridPoint(const std::vector<double>& lonC,
                                       const std::vector<double>& latC,
                                       size_t nlpb, double lon, double lat,
                                       bool use_multiple_neighbors = false) {
    // Normalize longitude
    if (lon > 180.0) lon -= 360.0;

    // Ensure KD-tree is initialized
    if (!kdtree_initialized_) {
      if (!initializeKDTree(lonC, latC, nlpb)) {
        // If initialization fails, use brute force search
        return findNearestGridPointDirect(lonC, latC, nlpb, lon, lat);
      }
    }

    // Use KD-tree search
    double query_x, query_y;
    geoToCartesian(lon, lat, query_x, query_y);

    // Determine number of neighbors to check
    const size_t knn = use_multiple_neighbors ? 5 : 1;

    std::vector<size_t> indices(knn);
    std::vector<double> distances_sq(knn);

    nanoflann::KNNResultSet<double> resultSet(knn);
    resultSet.init(indices.data(), distances_sq.data());

    double query_pt[2] = {query_x, query_y};
    static_kdtree_->findNeighbors(resultSet, query_pt,
                                  nanoflann::SearchParameters());

    // Find closest point
    GeoPoint result;
    double min_distance = std::numeric_limits<double>::max();

    // Log candidate points if using multiple neighbors mode
    if (use_multiple_neighbors) {
      MACOM_LOG_INFO("MACOMGeometryIterator",
                     "Candidates for nearest point to (" + std::to_string(lon) +
                         ", " + std::to_string(lat) + "):");
    }

    for (size_t i = 0; i < resultSet.size(); ++i) {
      size_t idx = indices[i];
      const GeoPoint& point = static_point_cloud_.pts[idx];
      double haversine_dist = haversineDistance(lon, lat, point.lon, point.lat);

      if (use_multiple_neighbors) {
        MACOM_LOG_INFO(
            "MACOMGeometryIterator",
            "  Candidate " + std::to_string(i) + ": index=" +
                std::to_string(idx) + " lon=" + std::to_string(point.lon) +
                " lat=" + std::to_string(point.lat) +
                " Euclidean dist=" + std::to_string(sqrt(distances_sq[i])) +
                " Haversine dist=" + std::to_string(haversine_dist) + " km");
      }

      if (haversine_dist < min_distance) {
        min_distance = haversine_dist;
        result = point;
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
      MACOM_LOG_INFO("MACOMGeometryIterator", logMsg);
      return result;
    }

    // If KD-tree search fails, use brute force search
    MACOM_LOG_WARNING("MACOMGeometryIterator",
                      "KD-tree search failed, using brute force search");
    return findNearestGridPointDirect(lonC, latC, nlpb, lon, lat);
  }

  /**
   * @brief Find the nearest grid point (direct search, not using KD-tree)
   *
   * @param lonC Longitude array
   * @param latC Latitude array
   * @param nlpb Number of grid points
   * @param lon Query point longitude
   * @param lat Query point latitude
   * @return GeoPoint Containing index, coordinates and distance of the nearest
   * point
   */
  static GeoPoint findNearestGridPointDirect(const std::vector<double>& lonC,
                                             const std::vector<double>& latC,
                                             size_t nlpb, double lon,
                                             double lat) {
    GeoPoint result;
    double min_distance = std::numeric_limits<double>::max();

    // Normalize longitude
    if (lon > 180.0) lon -= 360.0;

    // Brute force search
    for (size_t i = 0; i < nlpb; ++i) {
      double dist = haversineDistance(lon, lat, lonC[i], latC[i]);
      if (dist < min_distance) {
        min_distance = dist;
        result.index = i;
        result.lon = lonC[i];
        result.lat = latC[i];
        result.distance = dist;
      }
    }

    return result;
  }

  /**
   * @brief Find all grid points within specified radius
   *
   * @param lonC Longitude array (for fallback method)
   * @param latC Latitude array (for fallback method)
   * @param nlpb Number of grid points (for fallback method)
   * @param lon Query point longitude
   * @param lat Query point latitude
   * @param radius Search radius (kilometers)
   * @return GeoQueryResult Containing all points within radius
   */
  static GeoQueryResult findGridPointsInRadius(const std::vector<double>& lonC,
                                               const std::vector<double>& latC,
                                               size_t nlpb, double lon,
                                               double lat, double radius) {
    // Normalize longitude
    if (lon > 180.0) lon -= 360.0;

    // Ensure KD-tree is initialized
    if (!kdtree_initialized_) {
      if (!initializeKDTree(lonC, latC, nlpb)) {
        // If initialization fails, use brute force search
        return findGridPointsInRadiusDirect(lonC, latC, nlpb, lon, lat, radius);
      }
    }

    // Use KD-tree search
    using IndexType = typename KDTreeType::IndexType;
    using DistanceType = typename KDTreeType::DistanceType;
    using MatchItem = nanoflann::ResultItem<IndexType, DistanceType>;

    double query_x, query_y;
    geoToCartesian(lon, lat, query_x, query_y);

    // In Cartesian space, search radius needs conversion
    // With 2D Cartesian coordinates, search radius needs to be slightly larger
    double search_radius =
        radius *
        1.2;  // Add 20% to ensure all potentially in-range points are included

    // Prepare query point
    double query_pt[2] = {query_x, query_y};  // 2D query point

    // Fix vector type error by using MatchItem instead of ResultItem
    std::vector<MatchItem> matches;

    // Create search parameters
    nanoflann::SearchParameters params;
    params.sorted = true;  // Sort results by distance

    // Perform radius search
    const size_t num_matches = static_kdtree_->radiusSearch(
        query_pt, search_radius * search_radius, matches, params);

    // Build result
    GeoQueryResult result;

    if (num_matches > 0) {
      result.points.reserve(num_matches);

      for (size_t i = 0; i < num_matches; ++i) {
        size_t idx = matches[i].first;
        GeoPoint point = static_point_cloud_.pts[idx];

        // Calculate precise great circle distance
        double exact_distance =
            haversineDistance(lon, lat, point.lon, point.lat);

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
      MACOM_LOG_INFO("MACOMGeometryIterator", logMsg);

      return result;
    } else {
      std::string logMsg = "No points found within " + std::to_string(radius) +
                           " km of (" + std::to_string(lon) + ", " +
                           std::to_string(lat) +
                           ") using KD-tree, trying brute force";
      MACOM_LOG_INFO("MACOMGeometryIterator", logMsg);
    }

    // If KD-tree search finds no points, use brute force search
    return findGridPointsInRadiusDirect(lonC, latC, nlpb, lon, lat, radius);
  }

  /**
   * @brief Find all grid points within specified radius (direct search, not
   * using KD-tree)
   *
   * @param lonC Longitude array
   * @param latC Latitude array
   * @param nlpb Number of grid points
   * @param lon Query point longitude
   * @param lat Query point latitude
   * @param radius Search radius (kilometers)
   * @return GeoQueryResult Containing all points within radius
   */
  static GeoQueryResult findGridPointsInRadiusDirect(
      const std::vector<double>& lonC, const std::vector<double>& latC,
      size_t nlpb, double lon, double lat, double radius) {
    GeoQueryResult result;

    // Normalize longitude
    if (lon > 180.0) lon -= 360.0;

    // Brute force search for all points within radius
    for (size_t i = 0; i < nlpb; ++i) {
      double dist = haversineDistance(lon, lat, lonC[i], latC[i]);
      if (dist <= radius) {
        GeoPoint point;
        point.index = i;
        point.lon = lonC[i];
        point.lat = latC[i];
        point.distance = dist;

        // Calculate Cartesian coordinates (for consistency with original
        // structure)
        geoToCartesian(point.lon, point.lat, point.x, point.y);

        result.points.push_back(point);
      }
    }

    // Sort by distance
    std::sort(result.points.begin(), result.points.end(),
              [](const GeoPoint& a, const GeoPoint& b) {
                return a.distance < b.distance;
              });

    return result;
  }

 private:
  const MACOMGeometry<ConfigBackend>* geometry_;  // Pointer to parent geometry
  size_t i_;                                      // X index
  size_t j_;                                      // Y index
  size_t k_;                                      // Z index
  size_t index_;                                  // Linear index into the grid

  // Static KD-tree related members
  static bool kdtree_initialized_;
  static std::mutex kdtree_mutex_;
  static PointCloud<GeoPoint> static_point_cloud_;
  static std::unique_ptr<KDTreeType> static_kdtree_;
};

// Initialize static members
template <typename ConfigBackend>
bool MACOMGeometryIterator<ConfigBackend>::kdtree_initialized_ = false;

template <typename ConfigBackend>
std::mutex MACOMGeometryIterator<ConfigBackend>::kdtree_mutex_;

template <typename ConfigBackend>
PointCloud<GeoPoint> MACOMGeometryIterator<ConfigBackend>::static_point_cloud_;

template <typename ConfigBackend>
std::unique_ptr<typename MACOMGeometryIterator<ConfigBackend>::KDTreeType>
    MACOMGeometryIterator<ConfigBackend>::static_kdtree_ = nullptr;

template <typename ConfigBackend>
class MACOMGeometryConstIterator : public MACOMGeometryIterator<ConfigBackend> {
 public:
  // Inherit constructors from base class
  using MACOMGeometryIterator<ConfigBackend>::MACOMGeometryIterator;
};

}  // namespace metada::backends::macom