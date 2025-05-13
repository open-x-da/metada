#pragma once

#include <iomanip>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "DateTime.hpp"

namespace metada::framework {

/**
 * @brief Structure representing a single observation record
 *
 * @details Contains all essential metadata and measurement data for a single
 * observation point, including type, value, location, timestamp, and quality
 * information. Extended to fully support BUFR data structures.
 */
struct ObsRecord {
  // Basic observation data
  std::string type;  ///< Type of observation (e.g., "temperature", "pressure")
  double value;      ///< Observed value
  std::optional<std::string> unit;  ///< Unit of the observation

  // Station information
  std::string station_id;  ///< Station identifier (SID)
  double longitude;        ///< Station longitude (XOB)
  double latitude;         ///< Station latitude (YOB)
  double elevation;        ///< Station elevation (ELV)

  // Time information
  DateTime datetime;  ///< Date and time of observation

  // Report metadata
  std::string report_type;        ///< PREPBUFR report type (TYP)
  std::string input_report_type;  ///< Input report type (T29)
  std::string instrument_type;    ///< Instrument type (ITP)

  // Quality control information
  std::size_t qc_marker;  ///< Quality control marker/flag
  std::optional<std::string>
      error_message;  ///< Optional error message or additional information
};

/**
 * @brief Stream insertion operator for ObsRecord
 *
 * @param os Output stream
 * @param record The observation record to output
 * @return std::ostream& Reference to the output stream
 */
inline std::ostream& operator<<(std::ostream& os, const ObsRecord& record) {
  os << std::left << std::setw(15) << record.type << std::setw(10)
     << record.value << std::setw(15) << record.station_id << std::setw(10)
     << record.longitude << std::setw(10) << record.latitude << std::setw(10)
     << record.elevation << record.datetime.iso8601();
  return os;
}

}  // namespace metada::framework