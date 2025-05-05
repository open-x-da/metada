#pragma once

#include <optional>
#include <string>

#include "DateTime.hpp"

namespace metada::framework {

/**
 * @brief Structure representing a single observation record
 *
 * @details Contains all essential metadata and measurement data for a single
 * observation point, including type, value, location, timestamp, and quality
 * information.
 */
struct ObsRecord {
  std::string type;  ///< Type of observation (e.g., "temperature", "pressure")
  double value;      ///< Observed value
  std::optional<std::string> unit;        ///< Unit of the observation
  std::optional<std::string> identifier;  ///< ID of the observation
  std::string
      location;  ///< Observation location (e.g., station ID or coordinates)
  DateTime datetime;      ///< Date and time of observation
  std::size_t qc_marker;  ///< Quality control marker/flag
  std::optional<std::string>
      error_message;  ///< Optional error message or additional information
};

}  // namespace metada::framework