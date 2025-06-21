/**
 * @file PointObservation.hpp
 * @brief Generic data structures for point-based observations
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains generic data structures for representing point-based
 * observations, including their location, value, and error. These structures
 * are intended to be reused across different observation-related backend
 * implementations to promote consistency and code reuse.
 */
#pragma once

namespace metada::framework {

/**
 * @brief Structure to hold observation location information
 */
struct ObservationLocation {
  double latitude;   ///< Latitude in degrees
  double longitude;  ///< Longitude in degrees
  double level;      ///< Vertical level (pressure in hPa, height in m, etc.)

  ObservationLocation(double lat, double lon, double lev)
      : latitude(lat), longitude(lon), level(lev) {}

  bool operator==(const ObservationLocation& other) const {
    return latitude == other.latitude && longitude == other.longitude &&
           level == other.level;
  }
};

/**
 * @brief Structure to hold a single observation point
 */
struct ObservationPoint {
  ObservationLocation location;  ///< Location of the observation
  double value;                  ///< Observed value
  double error;                  ///< Observation error
  bool is_valid;                 ///< Whether the observation is valid

  ObservationPoint(const ObservationLocation& loc, double val, double err)
      : location(loc), value(val), error(err), is_valid(true) {}

  ObservationPoint(const ObservationLocation& loc)
      : location(loc), value(0.0), error(0.0), is_valid(false) {}

  bool operator==(const ObservationPoint& other) const {
    return location == other.location && value == other.value &&
           error == other.error && is_valid == other.is_valid;
  }
};

}  // namespace metada::framework