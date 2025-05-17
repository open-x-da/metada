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

// Shared part (station/report metadata)
struct ObsRecordShared {
  std::string station_id;
  double longitude;
  double latitude;
  double elevation;
  std::string report_type;
  std::string input_report_type;
  std::string instrument_type;
  DateTime datetime;
};

// Per-level part (variables for each level)
struct ObsLevelRecord {
  std::string type;  // e.g., "PRESS", "HUMID", etc.
  double value;
  std::size_t qc_marker;
  std::optional<std::string> unit;
};

// Main record: shared + all levels
struct ObsRecord {
  ObsRecordShared shared;
  std::vector<ObsLevelRecord> levels;
};

/**
 * @brief Stream insertion operator for ObsRecord
 *
 * @param os Output stream
 * @param record The observation record to output
 * @return std::ostream& Reference to the output stream
 */
inline std::ostream& operator<<(std::ostream& os, const ObsRecord& record) {
  const auto& s = record.shared;
  os << std::left << std::setw(15) << s.report_type << std::setw(15)
     << s.station_id << std::setw(10) << s.longitude << std::setw(10)
     << s.latitude << std::setw(10) << s.elevation << std::setw(25)
     << s.datetime.iso8601();
  for (const auto& level : record.levels) {
    os << std::setw(10) << level.type << std::setw(10) << level.value;
  }
  return os;
}

}  // namespace metada::framework