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

#include "Location.hpp"

namespace metada::framework {

/**
 * @brief Class representing a single observation at a location
 *
 * Contains a Location and observation-specific fields (value, error, validity).
 */
class ObservationPoint {
 public:
  Location location;  ///< Location of the observation
  double value;       ///< Observed value
  double error;       ///< Observation error
  bool is_valid;      ///< Whether the observation is valid

  // Constructors
  ObservationPoint(const Location& loc, double val, double err)
      : location(loc), value(val), error(err), is_valid(true) {}
  ObservationPoint(const Location& loc)
      : location(loc), value(0.0), error(0.0), is_valid(false) {}

  bool operator==(const ObservationPoint& other) const {
    return location == other.location && value == other.value &&
           error == other.error && is_valid == other.is_valid;
  }
};

}  // namespace metada::framework