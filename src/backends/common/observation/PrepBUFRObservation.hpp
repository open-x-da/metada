#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "BufrObsIO.hpp"
#include "ConfigValue.hpp"
#include "DateTime.hpp"
#include "Duration.hpp"
#include "ObsRecord.hpp"
#include "ParquetObservation.hpp"
#include "PointObservation.hpp"
#include "PrepBUFRObservationIterator.hpp"

namespace metada::backends::common::observation {

using ObsRecord = framework::ObsRecord;
using ObsLevelRecord = framework::ObsLevelRecord;
using ConfigValue = metada::backends::config::ConfigValue;
using ParquetObservation = framework::ParquetObservation;
using DateTime = metada::DateTime;
using Duration = metada::Duration;

// Lightweight containers describing PrepBUFR observation geometries.
struct SurfaceBatch {
  std::vector<double> lats;    // one per observation
  std::vector<double> lons;    // one per observation
  std::vector<double> levels;  // surface pressure or model-level proxy
};

struct ProfileBatch {
  std::vector<int> counts;          // number of levels per profile
  std::vector<double> lats;         // one per profile
  std::vector<double> lons;         // one per profile
  std::vector<double> levels_flat;  // concatenated per-level pressure/height
};

struct FamilyBatches {
  SurfaceBatch metar, synop, ships, buoy, sonde_sfc;
  ProfileBatch airep, pilot, sound;
};

// Build PrepBUFR batches grouped by WMO family/subset for downstream operators
// and provide a concept-compliant observation backend over flattened points.
class PrepBUFRObservation {
 public:
  // Types for concept compliance
  using iterator_type = PrepBUFRObservationIterator;
  using value_type = framework::ObservationPoint;

  // Delete default constructor and copying
  PrepBUFRObservation() = delete;
  PrepBUFRObservation(const PrepBUFRObservation&) = delete;
  PrepBUFRObservation& operator=(const PrepBUFRObservation&) = delete;

  // Construct from configuration backend reading actual PrepBUFR via BufrObsIO
  template <typename ConfigBackend>
  explicit PrepBUFRObservation(const ConfigBackend& config) {
    const auto type_configs = config.Get("types").asVectorMap();

    // Group observation types by file path to avoid reading the same file
    // multiple times
    std::map<
        std::string,
        std::vector<std::pair<std::string, std::map<std::string, ConfigValue>>>>
        file_groups;

    for (const auto& type_map : type_configs) {
      const auto& [type_name, type_config] = *type_map.begin();
      const auto& type_backend = ConfigBackend(type_config.asMap());
      bool type_if_use = type_backend.Get("if_use").asBool();
      if (!type_if_use) continue;

      // Get file path for this type
      std::string file_path;
      try {
        file_path = type_backend.Get("file").asString();
      } catch (...) {
        try {
          file_path = type_backend.Get("filename").asString();
        } catch (...) {
          continue;  // Skip types without valid file path
        }
      }

      // Store the configuration data as a map that can be copied
      std::map<std::string, ConfigValue> config_data;
      try {
        config_data["if_use"] = type_backend.Get("if_use");
        config_data["coordinate"] = type_backend.Get("coordinate");
        config_data["variables"] = type_backend.Get("variables");
        // Add data_types if present for BUFR subset filtering
        if (type_backend.HasKey("data_types")) {
          config_data["data_types"] = type_backend.Get("data_types");
        }
        // Add station_ids if present for station filtering
        if (type_backend.HasKey("station_ids")) {
          config_data["station_ids"] = type_backend.Get("station_ids");
        }
        // Add time window configuration if present
        if (type_backend.HasKey("analysis_time")) {
          config_data["analysis_time"] = type_backend.Get("analysis_time");
        }
        if (type_backend.HasKey("time_window")) {
          config_data["time_window"] = type_backend.Get("time_window");
        }
      } catch (...) {
        // Some keys might not exist, that's okay
      }

      file_groups[file_path].emplace_back(type_name, config_data);
    }

    // Process each file once and distribute records to appropriate types
    for (const auto& [file_path, type_configs] : file_groups) {
      std::cout << "Processing file: " << file_path << std::endl;

      // Read BUFR file once for all types that share it
      // Pass the type configurations to BufrObsIO for early filtering
      std::vector<ObsRecord> all_records;
      {
        // Create a configuration that includes all enabled types, variables,
        // and filters
        std::map<std::string, ConfigValue> bufr_config;
        bufr_config["file"] = file_path;

        // Collect all enabled data types, variables, and filters for early
        // filtering
        std::vector<std::string> enabled_data_types;
        std::vector<std::string> enabled_variables;
        std::vector<int> enabled_station_ids;

        // Time window filtering
        std::optional<DateTime> analysis_time;
        std::optional<Duration> time_window;

        for (const auto& [type_name, config_data] : type_configs) {
          // Only process types that are enabled
          if (config_data.at("if_use").asBool()) {
            // Add data types if present
            auto data_types_it = config_data.find("data_types");
            if (data_types_it != config_data.end()) {
              try {
                const auto& data_types_array =
                    data_types_it->second.asVectorString();
                enabled_data_types.insert(enabled_data_types.end(),
                                          data_types_array.begin(),
                                          data_types_array.end());
              } catch (...) {
                // Skip if parsing fails
              }
            }

            // Add station IDs if present
            auto station_ids_it = config_data.find("station_ids");
            if (station_ids_it != config_data.end()) {
              try {
                const auto& station_ids_array =
                    station_ids_it->second.asVectorInt();
                enabled_station_ids.insert(enabled_station_ids.end(),
                                           station_ids_array.begin(),
                                           station_ids_array.end());
              } catch (...) {
                // Skip if parsing fails
              }
            }

            // Add variables if present
            auto vars_it = config_data.find("variables");
            if (vars_it != config_data.end()) {
              try {
                const auto& variables_configs = vars_it->second.asVectorMap();
                for (const auto& var_map : variables_configs) {
                  for (const auto& [var_name, var_config] : var_map) {
                    const auto& var_backend = ConfigBackend(var_config.asMap());
                    if (var_backend.Get("if_use").asBool()) {
                      std::string bufr_name = mapUserVarToBufr(var_name);
                      if (!bufr_name.empty()) {
                        enabled_variables.push_back(bufr_name);
                      }
                    }
                  }
                }
              } catch (...) {
                // Skip if parsing fails
              }
            }
          }
        }

        // Collect time window information from the first enabled observation
        // type that has time window configuration
        for (const auto& [type_name, config_data] : type_configs) {
          if (config_data.at("if_use").asBool()) {
            try {
              if (config_data.find("analysis_time") != config_data.end()) {
                analysis_time =
                    DateTime(config_data.at("analysis_time").asString());
              }
              if (config_data.find("time_window") != config_data.end()) {
                std::string time_window_str =
                    config_data.at("time_window").asString();
                // Parse time window string like "±3h", "±30m", etc.
                if (time_window_str.starts_with("±")) {
                  time_window_str = time_window_str.substr(1);  // Remove ±
                  time_window = Duration(
                      time_window_str);  // Let Duration handle the parsing
                }
              }
              // Break after finding the first enabled type with time window
              // config
              if (analysis_time.has_value() || time_window.has_value()) {
                break;
              }
            } catch (...) {
              // Time window parsing failed for this type, continue to next
              continue;
            }
          }
        }

        // Add all filtering configuration to BufrObsIO
        if (!enabled_data_types.empty()) {
          bufr_config["data_types"] = enabled_data_types;
        }
        if (!enabled_variables.empty()) {
          bufr_config["variables"] = enabled_variables;
        }
        if (!enabled_station_ids.empty()) {
          bufr_config["station_ids"] = enabled_station_ids;
        }
        if (analysis_time.has_value()) {
          bufr_config["analysis_time"] = analysis_time->iso8601();
        }
        if (time_window.has_value()) {
          bufr_config["time_window"] = time_window->toString();
        }

        auto bufr_io = io::BufrObsIO<ConfigBackend>(ConfigBackend(bufr_config));
        all_records = bufr_io.read();
        // BufrObsIO destructor will close the file when it goes out of scope
      }

      std::cout << "Read " << all_records.size() << " total records from file"
                << std::endl;

      // Store BUFR subset distribution for later display
      for (const auto& rec : all_records) {
        // Trim trailing whitespace from BUFR subset names
        std::string trimmed_type = rec.shared.report_type;
        trimmed_type.erase(trimmed_type.find_last_not_of(" \t\r\n") + 1);
        bufr_subset_counts_[trimmed_type]++;
      }

      // Process each type that uses this file
      for (const auto& [type_name, config_data] : type_configs) {
        std::cout << "Processing type: " << type_name << std::endl;

        // Desired output coordinate system (BUFR is naturally geographic)
        std::string coordinate_system = "geographic";
        try {
          auto coord_it = config_data.find("coordinate");
          if (coord_it != config_data.end()) {
            coordinate_system = coord_it->second.asString();
          }
        } catch (...) {
          std::cout << "No coordinate section; default geographic" << std::endl;
        }

        // Collect enabled variables for this type with their errors
        struct VarSel {
          std::string bufr_name;
          double error;
        };
        std::vector<std::pair<std::string, VarSel>> enabled_vars;
        try {
          auto vars_it = config_data.find("variables");
          if (vars_it != config_data.end()) {
            const auto& variables_configs = vars_it->second.asVectorMap();
            for (const auto& var_map : variables_configs) {
              for (const auto& [var_name, var_config] : var_map) {
                const auto& var_backend = ConfigBackend(var_config.asMap());
                bool var_if_use = var_backend.Get("if_use").asBool();
                if (!var_if_use) continue;
                double error = var_backend.Get("error").asFloat();
                // Map model variable name to BUFR variable short name
                std::string bufr_name = mapUserVarToBufr(var_name);
                if (!bufr_name.empty()) {
                  enabled_vars.push_back({var_name, VarSel{bufr_name, error}});
                  std::cout << "  Enabled variable: " << var_name << " -> "
                            << bufr_name << " (error: " << error << ")"
                            << std::endl;
                }
              }
            }
          }
        } catch (...) {
          std::cout << "No variables section for " << type_name << std::endl;
        }

        // Initialize the type in the map
        type_variable_map_[type_name] =
            std::unordered_map<std::string, std::vector<size_t>>();
        for (const auto& [user_var, _] : enabled_vars) {
          type_variable_map_[type_name][user_var] = std::vector<size_t>();
        }

        // Get allowed BUFR data types for this observation type
        std::vector<std::string> allowed_data_types;
        try {
          auto data_types_it = config_data.find("data_types");
          if (data_types_it != config_data.end()) {
            std::cout << "  Found data_types in config" << std::endl;
            // data_types is a simple array of strings
            const auto& data_types_array =
                data_types_it->second.asVectorString();
            std::cout << "  Parsed as string array with "
                      << data_types_array.size() << " elements" << std::endl;
            for (const auto& dt_name : data_types_array) {
              allowed_data_types.push_back(dt_name);
              std::cout << "    Added allowed type: " << dt_name << std::endl;
            }
          } else {
            std::cout << "  No data_types key found in config_data"
                      << std::endl;
          }
        } catch (const std::exception& e) {
          std::cout << "  Exception parsing data_types: " << e.what()
                    << std::endl;
        }

        // Show which BUFR subsets will be processed for this type
        if (!allowed_data_types.empty()) {
          std::cout << "  Allowed BUFR subsets: ";
          for (size_t i = 0; i < allowed_data_types.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << allowed_data_types[i];
          }
          std::cout << std::endl;
        } else {
          std::cout << "  No BUFR subset filtering - processing all records"
                    << std::endl;
        }

        // Flatten BUFR records into observation points for selected variables
        size_t type_obs_count = 0;
        size_t records_processed = 0;
        size_t records_filtered = 0;
        size_t levels_processed = 0;
        size_t variables_found = 0;

        for (const auto& rec : all_records) {
          // Filter by BUFR subset if data_types are specified
          if (!allowed_data_types.empty()) {
            bool subset_allowed = false;
            // Trim trailing whitespace from BUFR subset name for comparison
            std::string trimmed_report_type = rec.shared.report_type;
            trimmed_report_type.erase(
                trimmed_report_type.find_last_not_of(" \t\r\n") + 1);

            for (const auto& allowed_type : allowed_data_types) {
              if (trimmed_report_type == allowed_type) {
                subset_allowed = true;
                break;
              }
            }
            if (!subset_allowed) {
              records_filtered++;
              continue;  // Skip this record if subset not allowed
            }
          }

          records_processed++;
          const double lat = rec.shared.latitude;
          const double lon = rec.shared.longitude;

          for (const auto& lvl : rec.levels) {
            levels_processed++;
            // Pick a representative level coordinate (prefer PRES then HEIGHT)
            double level_coord = 0.0;
            for (const auto& v : lvl) {
              if (v.type == "PRES") {
                level_coord = v.value;
                break;
              }
              if (v.type == "HEIGHT") {
                level_coord = v.value;
              }
            }

            for (const auto& [user_var, sel] : enabled_vars) {
              // Find the requested BUFR variable in this level
              for (const auto& v : lvl) {
                if (v.type == sel.bufr_name) {
                  variables_found++;
                  framework::Location location(0, 0, 0);
                  if (coordinate_system == "grid") {
                    // No grid indices in BUFR; fallback to geographic
                    location = framework::Location(
                        lat, lon, level_coord,
                        framework::CoordinateSystem::GEOGRAPHIC);
                  } else {
                    location = framework::Location(
                        lat, lon, level_coord,
                        framework::CoordinateSystem::GEOGRAPHIC);
                  }

                  size_t obs_index = observations_.size();
                  observations_.emplace_back(location, v.value, sel.error);
                  type_variable_map_[type_name][user_var].push_back(obs_index);
                  type_obs_count++;
                  break;
                }
              }
            }
          }
        }

        // Debug output for understanding why no observations were found
        std::cout << "  Debug: Records processed: " << records_processed
                  << ", Records filtered out: " << records_filtered
                  << ", Levels processed: " << levels_processed
                  << ", Variables found: " << variables_found << std::endl;

        std::cout << "Added " << type_obs_count << " observations for type "
                  << type_name << std::endl;

        // Store BUFR subset distribution for this type
        if (!allowed_data_types.empty()) {
          for (const auto& rec : all_records) {
            // Trim trailing whitespace from BUFR subset name for comparison
            std::string trimmed_report_type = rec.shared.report_type;
            trimmed_report_type.erase(
                trimmed_report_type.find_last_not_of(" \t\r\n") + 1);

            for (const auto& allowed_type : allowed_data_types) {
              if (trimmed_report_type == allowed_type) {
                type_bufr_subset_counts_[type_name][trimmed_report_type]++;
                break;
              }
            }
          }
        }
      }

      // Store the configuration for this file's types for later reference
      for (const auto& [type_name, config_data] : type_configs) {
        type_configs_[type_name] = config_data;
      }
    }

    std::cout << "Total observations: " << observations_.size() << std::endl;
    std::cout << "Type map size: " << type_variable_map_.size() << std::endl;
    for (const auto& [type_name, var_map] : type_variable_map_) {
      std::cout << "  Type " << type_name << " has " << var_map.size()
                << " variables" << std::endl;
      for (const auto& [var_name, indices] : var_map) {
        std::cout << "    Variable " << var_name << " has " << indices.size()
                  << " observations" << std::endl;
      }
    }
  }

  // Move semantics
  PrepBUFRObservation(PrepBUFRObservation&& other) noexcept
      : observations_(std::move(other.observations_)),
        type_variable_map_(std::move(other.type_variable_map_)),
        covariance_(std::move(other.covariance_)),
        bufr_subset_counts_(std::move(other.bufr_subset_counts_)),
        type_bufr_subset_counts_(std::move(other.type_bufr_subset_counts_)),
        type_configs_(std::move(other.type_configs_)) {}

  PrepBUFRObservation& operator=(PrepBUFRObservation&& other) noexcept {
    if (this != &other) {
      observations_ = std::move(other.observations_);
      type_variable_map_ = std::move(other.type_variable_map_);
      covariance_ = std::move(other.covariance_);
      bufr_subset_counts_ = std::move(other.bufr_subset_counts_);
      type_bufr_subset_counts_ = std::move(other.type_bufr_subset_counts_);
      type_configs_ = std::move(other.type_configs_);
    }
    return *this;
  }

  // Iteration
  PrepBUFRObservationIterator begin() const {
    return PrepBUFRObservationIterator(&observations_, 0);
  }
  PrepBUFRObservationIterator end() const {
    return PrepBUFRObservationIterator(&observations_, observations_.size());
  }

  // Size and element access
  size_t size() const { return observations_.size(); }
  const framework::ObservationPoint& operator[](size_t index) const {
    return observations_[index];
  }

  // Raw data access (legacy compatibility)
  void* getData() {
    if (observations_.empty()) return nullptr;
    return const_cast<void*>(static_cast<const void*>(observations_.data()));
  }
  const void* getData() const {
    if (observations_.empty()) return nullptr;
    return static_cast<const void*>(observations_.data());
  }

  template <typename T>
  T getData() const {
    if constexpr (std::is_same_v<T, std::vector<double>>) {
      std::vector<double> values;
      values.reserve(observations_.size());
      for (const auto& obs : observations_) values.push_back(obs.value);
      return values;
    }
    return T{};
  }

  // Type/variable names (backward mapping placeholders)
  std::vector<std::string> getTypeNames() const {
    std::vector<std::string> types;
    for (const auto& [type_name, _] : type_variable_map_)
      types.push_back(type_name);
    return types;
  }
  std::vector<std::string> getVariableNames(
      const std::string& type_name) const {
    std::vector<std::string> variables;
    auto it = type_variable_map_.find(type_name);
    if (it != type_variable_map_.end()) {
      for (const auto& [var_name, _] : it->second)
        variables.push_back(var_name);
    }
    return variables;
  }

  // Covariance helpers
  std::vector<double> getCovariance() const {
    std::vector<double> errors;
    errors.reserve(observations_.size());
    for (const auto& obs : observations_) {
      if (obs.is_valid) {
        errors.push_back(obs.error * obs.error);
      } else {
        errors.push_back(std::numeric_limits<double>::infinity());
      }
    }
    return errors;
  }

  double quadraticForm(const std::vector<double>& innovation) const {
    if (innovation.size() != observations_.size()) {
      throw std::invalid_argument("Innovation vector size mismatch");
    }
    double result = 0.0;
    for (size_t i = 0; i < observations_.size(); ++i) {
      if (observations_[i].is_valid) {
        double variance = observations_[i].error * observations_[i].error;
        result += (innovation[i] * innovation[i]) / variance;
      }
    }
    return result;
  }

  std::vector<double> applyInverseCovariance(
      const std::vector<double>& innovation) const {
    if (innovation.size() != observations_.size()) {
      throw std::invalid_argument("Innovation vector size mismatch");
    }
    std::vector<double> result(innovation.size());
    for (size_t i = 0; i < innovation.size(); ++i) {
      if (observations_[i].is_valid) {
        double variance = observations_[i].error * observations_[i].error;
        result[i] = innovation[i] / variance;
      } else {
        result[i] = 0.0;
      }
    }
    return result;
  }

  // Clone
  std::unique_ptr<PrepBUFRObservation> clone() const {
    return std::unique_ptr<PrepBUFRObservation>(
        new PrepBUFRObservation(*this, true));
  }

  // Arithmetic
  void add(const PrepBUFRObservation& other) {
    if (observations_.size() != other.observations_.size()) {
      throw std::runtime_error("Cannot add observations of different sizes");
    }
    for (size_t i = 0; i < observations_.size(); ++i) {
      if (observations_[i].is_valid && other.observations_[i].is_valid) {
        observations_[i].value += other.observations_[i].value;
      }
    }
  }
  void subtract(const PrepBUFRObservation& other) {
    if (observations_.size() != other.observations_.size()) {
      throw std::runtime_error(
          "Cannot subtract observations of different sizes");
    }
    for (size_t i = 0; i < observations_.size(); ++i) {
      if (observations_[i].is_valid && other.observations_[i].is_valid) {
        observations_[i].value -= other.observations_[i].value;
      }
    }
  }
  void multiply(double scalar) {
    for (auto& obs : observations_) {
      if (obs.is_valid) obs.value *= scalar;
    }
  }

  // Comparison
  bool equals(const PrepBUFRObservation& other) const {
    return observations_ == other.observations_ &&
           type_variable_map_ == other.type_variable_map_ &&
           covariance_ == other.covariance_;
  }

  // Lifecycle
  void initialize() {}
  void applyQC() {
    for (auto& obs : observations_) {
      if (obs.is_valid) {
        if (obs.value < -100.0 || obs.value > 100.0) obs.is_valid = false;
        if (obs.location.getCoordinateSystem() ==
            framework::CoordinateSystem::GEOGRAPHIC) {
          auto [lat, lon, level] = obs.location.getGeographicCoords();
          if (lat < -90.0 || lat > 90.0 || lon < -180.0 || lon > 180.0) {
            obs.is_valid = false;
          }
        }
      }
    }
  }

  // File I/O (simple CSV-like loader for testing; same as GridObservation)
  void loadFromFile(const std::string& filename, double error,
                    double missing_value,
                    const std::string& coordinate_system = "geographic") {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open observation file: " + filename);
    }
    observations_.clear();

    std::string line;
    bool header_found = false;
    bool data_started = false;
    while (std::getline(file, line)) {
      if (line.empty()) continue;
      if (line.find("Time") != std::string::npos &&
          line.find("XLONG_U") != std::string::npos) {
        header_found = true;
        continue;
      }
      if (line.find("---") != std::string::npos &&
          line.find_first_not_of("- \t") == std::string::npos) {
        if (header_found) data_started = true;
        continue;
      }
      if (!data_started) continue;

      std::istringstream iss(line);
      int time, z, y, x;
      double znu, xlong_u, xlat_u, obs_value;
      if (iss >> time >> z >> y >> x >> znu >> xlong_u >> xlat_u >> obs_value) {
        if (obs_value != missing_value) {
          framework::Location location(0, 0, 0);
          if (coordinate_system == "grid") {
            location = framework::Location(x, y, z);
          } else if (coordinate_system == "geographic") {
            location = framework::Location(
                xlat_u, xlong_u, znu, framework::CoordinateSystem::GEOGRAPHIC);
          } else {
            throw std::runtime_error("Unsupported coordinate system: " +
                                     coordinate_system);
          }
          observations_.emplace_back(location, obs_value, error);
        }
      }
    }
    if (observations_.empty()) {
      throw std::runtime_error("No valid observations loaded from file: " +
                               filename);
    }
  }

  void saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open file for writing: " + filename);
    }
    for (const auto& obs : observations_) {
      if (obs.is_valid) {
        if (obs.location.getCoordinateSystem() ==
            framework::CoordinateSystem::GEOGRAPHIC) {
          auto [lat, lon, level] = obs.location.getGeographicCoords();
          file << std::fixed << std::setprecision(6) << std::setw(10) << lat
               << std::setw(10) << lon << std::setw(10) << level
               << std::setw(12) << obs.value << std::setw(12) << obs.error
               << "\n";
        }
      }
    }
  }

  // Geographic filtering
  std::vector<framework::ObservationPoint> getObservationsInBox(
      double min_lat, double max_lat, double min_lon, double max_lon) const {
    std::vector<framework::ObservationPoint> result;
    for (const auto& obs : observations_) {
      if (obs.is_valid && obs.location.getCoordinateSystem() ==
                              framework::CoordinateSystem::GEOGRAPHIC) {
        auto [lat, lon, level] = obs.location.getGeographicCoords();
        if (lat >= min_lat && lat <= max_lat && lon >= min_lon &&
            lon <= max_lon) {
          result.push_back(obs);
        }
      }
    }
    return result;
  }

  std::vector<framework::ObservationPoint> getObservationsInVerticalRange(
      double min_level, double max_level) const {
    std::vector<framework::ObservationPoint> result;
    for (const auto& obs : observations_) {
      if (obs.is_valid && obs.location.getCoordinateSystem() ==
                              framework::CoordinateSystem::GEOGRAPHIC) {
        auto [lat, lon, level] = obs.location.getGeographicCoords();
        if (level >= min_level && level <= max_level) {
          result.push_back(obs);
        }
      }
    }
    return result;
  }

  /**
   * @brief Get the size of observation data for a specific type/variable
   * @param type_name Name of the observation type
   * @param var_name Name of the variable
   * @return Size of the observation data vector
   */
  size_t getSize(const std::string& type_name,
                 const std::string& var_name) const {
    auto type_it = type_variable_map_.find(type_name);
    if (type_it != type_variable_map_.end()) {
      auto var_it = type_it->second.find(var_name);
      if (var_it != type_it->second.end()) {
        return var_it->second.size();
      }
    }
    return 0;
  }

  /**
   * @brief Check if observation error covariance is diagonal
   * @return True if R is diagonal, false for full covariance matrix
   */
  bool isDiagonalCovariance() const { return true; }

  /**
   * @brief Get information about applied filters
   *
   * @details Returns a summary of what filtering was applied during observation
   * loading, including data types, variables, and geographic bounds.
   *
   * @return String containing filtering information
   */
  std::string getFilteringInfo() const {
    std::ostringstream oss;
    oss << "PrepBUFR Observation Filtering Summary:\n";
    oss << "=====================================\n";

    // Show enabled types and their variables
    for (const auto& [type_name, var_map] : type_variable_map_) {
      oss << "Type: " << type_name << "\n";
      oss << "  Variables: ";
      for (const auto& [var_name, indices] : var_map) {
        oss << var_name << "(" << indices.size() << " obs) ";
      }
      oss << "\n";

      // Show station filtering if applied
      auto type_config_it = type_configs_.find(type_name);
      if (type_config_it != type_configs_.end()) {
        auto station_ids_it = type_config_it->second.find("station_ids");
        if (station_ids_it != type_config_it->second.end()) {
          try {
            const auto& station_ids = station_ids_it->second.asVectorInt();
            if (!station_ids.empty()) {
              oss << "  Station Filter: Only stations ";
              for (size_t i = 0; i < station_ids.size(); ++i) {
                if (i > 0) oss << ", ";
                oss << station_ids[i];
              }
              oss << "\n";
            }
          } catch (...) {
            // Skip if parsing fails
          }
        }

        // Show time window filtering if applied
        auto analysis_time_it = type_config_it->second.find("analysis_time");
        auto time_window_it = type_config_it->second.find("time_window");
        if (analysis_time_it != type_config_it->second.end() &&
            time_window_it != type_config_it->second.end()) {
          try {
            std::string analysis_time_str = analysis_time_it->second.asString();
            std::string time_window_str = time_window_it->second.asString();
            oss << "  Time Filter: " << analysis_time_str << " "
                << time_window_str << "\n";
          } catch (...) {
            // Skip if parsing fails
          }
        }
      }
    }

    // Show BUFR subset distribution
    if (!bufr_subset_counts_.empty()) {
      oss << "\nBUFR Subset Distribution:\n";
      for (const auto& [subset, count] : bufr_subset_counts_) {
        oss << "  " << subset << ": " << count << " records\n";
      }
    }

    return oss.str();
  }

  // Batch builder from ObsRecords (static utility)
  static FamilyBatches groupByFamily(const std::vector<ObsRecord>& records) {
    FamilyBatches batches;
    for (const auto& rec : records) {
      const double lat = rec.shared.latitude;
      const double lon = rec.shared.longitude;
      const std::string& subset = rec.shared.report_type;
      if (isSurfaceSubset(subset)) {
        auto& b = selectSurface(batches, subset);
        b.lats.push_back(lat);
        b.lons.push_back(lon);
        b.levels.push_back(extractSurfacePressure(rec));
      } else if (isAirep(subset)) {
        appendProfile(rec, lat, lon, batches.airep);
      } else if (isPilot(subset)) {
        appendProfileWinds(rec, lat, lon, batches.pilot);
      } else if (isSound(subset)) {
        appendProfile(rec, lat, lon, batches.sound);
      }
    }
    return batches;
  }

 private:
  // Additional covariance helpers required by Observation adapter API
  std::vector<double> applyCovariance(const std::vector<double>& vec) const {
    if (vec.size() != observations_.size()) {
      throw std::invalid_argument("Vector size mismatch in applyCovariance");
    }
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
      if (observations_[i].is_valid) {
        double variance = observations_[i].error * observations_[i].error;
        result[i] = variance * vec[i];
      } else {
        result[i] = 0.0;
      }
    }
    return result;
  }

  std::vector<double> getInverseCovarianceDiagonal() const {
    std::vector<double> diag;
    diag.reserve(observations_.size());
    for (const auto& obs : observations_) {
      if (obs.is_valid) {
        double variance = obs.error * obs.error;
        diag.push_back(variance > 0.0 ? 1.0 / variance : 0.0);
      } else {
        diag.push_back(0.0);
      }
    }
    return diag;
  }

  static bool isSurfaceSubset(const std::string& subset) {
    return subset == "METAR" || subset == "SYNOP" || subset == "BUOY" ||
           subset == "SHIP" || subset == "SOND";
  }

  static bool isAirep(const std::string& subset) {
    return subset == "AIREP" || subset == "AMDAR" || subset == "TAMDAR";
  }

  static bool isPilot(const std::string& subset) { return subset == "PILOT"; }

  static bool isSound(const std::string& subset) {
    return subset == "SOUND" || subset == "RAOB";
  }

  static SurfaceBatch& selectSurface(FamilyBatches& b,
                                     const std::string& subset) {
    if (subset == "METAR") return b.metar;
    if (subset == "SYNOP") return b.synop;
    if (subset == "BUOY") return b.buoy;
    if (subset == "SHIP") return b.ships;
    return b.sonde_sfc;
  }

  static double extractSurfacePressure(const ObsRecord& rec) {
    // Search first level for PRES; fall back to elevation if missing
    for (const auto& lvl : rec.levels) {
      for (const auto& v : lvl) {
        if (v.type == "PRES" || v.type == "PRESS") return v.value;
      }
      break;
    }
    return rec.shared.elevation;
  }

  static void appendProfile(const ObsRecord& rec, double lat, double lon,
                            ProfileBatch& dst) {
    const int num_levels = static_cast<int>(rec.levels.size());
    if (num_levels <= 0) return;
    dst.lats.push_back(lat);
    dst.lons.push_back(lon);
    dst.counts.push_back(num_levels);
    for (const auto& lvl : rec.levels) {
      double p = findValue(lvl, "PRES");
      dst.levels_flat.push_back(p);
    }
  }

  static void appendProfileWinds(const ObsRecord& rec, double lat, double lon,
                                 ProfileBatch& dst) {
    // Keep counts and levels for winds similar to thermal/moisture profiles
    appendProfile(rec, lat, lon, dst);
  }

  static double findValue(const std::vector<ObsLevelRecord>& lvl,
                          const std::string& key) {
    for (const auto& v : lvl) {
      if (v.type == key) return v.value;
    }
    return 0.0;
  }
  static std::string mapUserVarToBufr(const std::string& user) {
    if (user == "U" || user == "u" || user == "UWIND") return "UWIND";
    if (user == "V" || user == "v" || user == "VWIND") return "VWIND";
    if (user == "T" || user == "t" || user == "TEMP" || user == "Temperature")
      return "TEMP";
    if (user == "Q" || user == "q" || user == "HUMID" || user == "RH")
      return "HUMID";
    if (user == "Z" || user == "z" || user == "HEIGHT") return "HEIGHT";
    if (user == "P" || user == "p" || user == "PRES" || user == "PRESS")
      return "PRES";
    return std::string();
  }

  // Private clone constructor
  PrepBUFRObservation(const PrepBUFRObservation& other, bool)
      : observations_(other.observations_),
        type_variable_map_(other.type_variable_map_),
        covariance_(other.covariance_) {}

  // Implementation of ParquetObservation methods
  ParquetObservation toParquetObservation() const {
    ParquetObservation result;

    // Reserve space for efficiency
    result.reserve(observations_.size());

    for (size_t i = 0; i < observations_.size(); ++i) {
      const auto& obs = observations_[i];

      // Create metadata
      ParquetObservation::ObservationMeta meta;
      meta.station_id = "UNKNOWN";  // BUFR doesn't always have station IDs

      // Get geographic coordinates from Location
      auto [lat, lon, level] = obs.location.getGeographicCoords();
      meta.longitude = lon;
      meta.latitude = lat;
      meta.elevation = level;
      meta.pressure = level;  // Use level as pressure proxy
      meta.datetime =
          DateTime();  // Default datetime, would need to be extracted from BUFR
      meta.obs_type = 0;  // Default, would need to be mapped from BUFR subset
      meta.channel = 0;   // Surface observations don't have channels
      meta.obs_error = obs.error;
      meta.qc_flag = obs.is_valid ? 1 : 0;
      meta.report_type = "BUFR";
      meta.instrument_type = "UNKNOWN";

      // Determine category based on pressure/height
      if (level > 0) {
        meta.category = ParquetObservation::ObsType::PROFILE;
      } else {
        meta.category = ParquetObservation::ObsType::SURFACE;
      }

      // Create field values
      std::unordered_map<std::string, ParquetObservation::FieldValue> fields;
      fields["value"] = obs.value;
      fields["error"] = obs.error;
      fields["is_valid"] = obs.is_valid ? 1 : 0;

      // Add location fields
      fields["longitude"] = lon;
      fields["latitude"] = lat;
      fields["height"] = level;

      result.addObservation(meta, fields);
    }

    return result;
  }

  template <typename ConfigBackend>
  static ParquetObservation createParquetObservation(
      const ConfigBackend& config) {
    ParquetObservation result;

    const auto type_configs = config.Get("types").asVectorMap();

    // Group observation types by file path to avoid reading the same file
    // multiple times
    std::map<
        std::string,
        std::vector<std::pair<std::string, std::map<std::string, ConfigValue>>>>
        file_groups;

    for (const auto& type_map : type_configs) {
      const auto& [type_name, type_config] = *type_map.begin();
      const auto& type_backend = ConfigBackend(type_config.asMap());
      bool type_if_use = type_backend.Get("if_use").asBool();
      if (!type_if_use) continue;

      // Get file path for this type
      std::string file_path;
      try {
        file_path = type_backend.Get("file").asString();
      } catch (...) {
        try {
          file_path = type_backend.Get("filename").asString();
        } catch (...) {
          continue;  // Skip types without valid file path
        }
      }

      // Store the configuration data as a map that can be copied
      std::map<std::string, ConfigValue> config_data;
      try {
        config_data["if_use"] = type_backend.Get("if_use");
        config_data["coordinate"] = type_backend.Get("coordinate");
        config_data["variables"] = type_backend.Get("variables");
      } catch (...) {
        // Some keys might not exist, that's okay
      }

      file_groups[file_path].emplace_back(type_name, config_data);
    }

    // Process each file once and distribute records to appropriate types
    for (const auto& [file_path, type_configs] : file_groups) {
      // Read BUFR file once for all types that share it
      // Apply early filtering based on enabled types and variables
      std::map<std::string, ConfigValue> bufr_config;
      bufr_config["file"] = file_path;

      // Collect all enabled data types and variables for early filtering
      std::vector<std::string> enabled_data_types;
      std::vector<std::string> enabled_variables;

      for (const auto& [type_name, config_data] : type_configs) {
        // Only process types that are enabled
        if (config_data.at("if_use").asBool()) {
          // Add data types if present
          auto data_types_it = config_data.find("data_types");
          if (data_types_it != config_data.end()) {
            try {
              const auto& data_types_array =
                  data_types_it->second.asVectorString();
              enabled_data_types.insert(enabled_data_types.end(),
                                        data_types_array.begin(),
                                        data_types_array.end());
            } catch (...) {
              // Skip if parsing fails
            }
          }

          // Add variables if present
          auto vars_it = config_data.find("variables");
          if (vars_it != config_data.end()) {
            try {
              const auto& variables_configs = vars_it->second.asVectorMap();
              for (const auto& var_map : variables_configs) {
                for (const auto& [var_name, var_config] : var_map) {
                  const auto& var_backend = ConfigBackend(var_config.asMap());
                  if (var_backend.Get("if_use").asBool()) {
                    std::string bufr_name = mapUserVarToBufr(var_name);
                    if (!bufr_name.empty()) {
                      enabled_variables.push_back(bufr_name);
                    }
                  }
                }
              }
            } catch (...) {
              // Skip if parsing fails
            }
          }
        }
      }

      // Add filtering configuration to BufrObsIO
      if (!enabled_data_types.empty()) {
        bufr_config["data_types"] = enabled_data_types;
      }
      if (!enabled_variables.empty()) {
        bufr_config["variables"] = enabled_variables;
      }

      auto bufr_io = io::BufrObsIO<ConfigBackend>(ConfigBackend(bufr_config));
      std::vector<ObsRecord> records = bufr_io.read();

      // Process each type that uses this file
      for (const auto& [type_name, config_data] : type_configs) {
        // Desired output coordinate system (BUFR is naturally geographic)
        std::string coordinate_system = "geographic";
        std::cout << "Processing type: " << type_name << std::endl;
        try {
          auto coord_it = config_data.find("coordinate");
          if (coord_it != config_data.end()) {
            coordinate_system = coord_it->second.asString();
          }
        } catch (...) {
          std::cout << "No coordinate section; default geographic" << std::endl;
        }

        // Collect enabled variables for this type with their errors
        struct VarSel {
          std::string bufr_name;
          double error;
        };
        std::vector<std::pair<std::string, VarSel>>
            enabled_vars;  // (user var, sel)
        try {
          auto vars_it = config_data.find("variables");
          if (vars_it != config_data.end()) {
            const auto& variables_configs = vars_it->second.asVectorMap();
            for (const auto& var_map : variables_configs) {
              for (const auto& [var_name, var_config] : var_map) {
                const auto& var_backend = ConfigBackend(var_config.asMap());
                bool var_if_use = var_backend.Get("if_use").asBool();
                if (!var_if_use) continue;
                double error = var_backend.Get("error").asFloat();
                // Map model variable name to BUFR variable short name
                std::string bufr_name = mapUserVarToBufr(var_name);
                if (!bufr_name.empty()) {
                  enabled_vars.push_back({var_name, VarSel{bufr_name, error}});
                }
              }
            }
          }
        } catch (...) {
          // No variables section; default none
        }

        // Flatten BUFR records into ParquetObservation for selected variables
        for (const auto& rec : records) {
          const double lat = rec.shared.latitude;
          const double lon = rec.shared.longitude;

          for (const auto& lvl : rec.levels) {
            // Pick a representative level coordinate (prefer PRES then HEIGHT)
            double level_coord = 0.0;
            for (const auto& v : lvl) {
              if (v.type == "PRES") {
                level_coord = v.value;
                break;
              }
              if (v.type == "HEIGHT") {
                level_coord = v.value;
              }
            }

            for (const auto& [user_var, sel] : enabled_vars) {
              // Find the requested BUFR variable in this level
              for (const auto& v : lvl) {
                if (v.type == sel.bufr_name) {
                  // Create metadata for this observation
                  ParquetObservation::ObservationMeta meta;
                  meta.station_id = rec.shared.station_id;
                  meta.longitude = lon;
                  meta.latitude = lat;
                  meta.elevation = rec.shared.elevation;
                  meta.pressure = level_coord;
                  meta.datetime = rec.shared.datetime;
                  meta.obs_type = getObsTypeFromSubset(type_name);
                  meta.channel =
                      0;  // Surface/profile observations don't have channels
                  meta.obs_error = sel.error;
                  meta.qc_flag = v.qc_marker;
                  meta.report_type = type_name;
                  meta.instrument_type = rec.shared.instrument_type;

                  // Determine category based on observation type
                  if (isSurfaceSubset(type_name)) {
                    meta.category = ParquetObservation::ObsType::SURFACE;
                  } else if (isAirep(type_name) || isPilot(type_name) ||
                             isSound(type_name)) {
                    meta.category = ParquetObservation::ObsType::PROFILE;
                  } else {
                    meta.category = ParquetObservation::ObsType::UNKNOWN;
                  }

                  // Create field values
                  std::unordered_map<std::string,
                                     ParquetObservation::FieldValue>
                      fields;
                  fields[user_var] = v.value;
                  fields["error"] = sel.error;
                  fields["qc_marker"] = static_cast<int>(v.qc_marker);
                  fields["unit"] = v.unit ? *v.unit : "";

                  // Add location fields
                  fields["longitude"] = lon;
                  fields["latitude"] = lat;
                  fields["height"] = level_coord;
                  fields["elevation"] = rec.shared.elevation;

                  result.addObservation(meta, fields);
                  break;
                }
              }
            }
          }
        }
      }
    }

    return result;
  }

  // Helper method to map subset names to WMO observation type codes
  static int getObsTypeFromSubset(const std::string& subset) {
    if (subset == "METAR") return 1;
    if (subset == "SYNOP") return 2;
    if (subset == "SHIP") return 3;
    if (subset == "BUOY") return 4;
    if (subset == "AIREP") return 5;
    if (subset == "AMDAR") return 6;
    if (subset == "TAMDAR") return 7;
    if (subset == "PILOT") return 8;
    if (subset == "SOUND") return 9;
    if (subset == "RAOB") return 10;
    return 0;  // Unknown
  }

  // Data members
  std::vector<framework::ObservationPoint> observations_;
  std::unordered_map<std::string,
                     std::unordered_map<std::string, std::vector<size_t>>>
      type_variable_map_;
  std::vector<double> covariance_;
  std::map<std::string, size_t>
      bufr_subset_counts_;  // Store overall BUFR subset distribution
  std::map<std::string, std::map<std::string, size_t>>
      type_bufr_subset_counts_;  // Store per-type BUFR subset distribution
  std::map<std::string, std::map<std::string, ConfigValue>>
      type_configs_;  // Store original configuration for filtering info

  /**
   * @brief Stream insertion operator for PrepBUFRObservation summary
   *
   * @details Outputs a detailed summary of the PrepBUFR observation data
   * including:
   * - Total number of observations
   * - Summary for each observation type and variable
   * - Geographic distribution information
   * - Quality control statistics
   *
   * @param os Output stream to write to
   * @param obs PrepBUFRObservation object to summarize
   * @return Reference to the output stream
   */
  friend std::ostream& operator<<(std::ostream& os,
                                  const PrepBUFRObservation& obs) {
    os << "=== PrepBUFR Observation Summary ===\n";
    os << "Total observations: " << obs.observations_.size() << "\n\n";

    // Count valid vs invalid observations
    size_t valid_count = 0;
    size_t invalid_count = 0;
    for (const auto& observation : obs.observations_) {
      if (observation.is_valid) {
        valid_count++;
      } else {
        invalid_count++;
      }
    }
    os << "Quality Control:\n";
    os << std::string(50, '-') << "\n";
    os << "Valid observations: " << valid_count << "\n";
    os << "Invalid observations: " << invalid_count << "\n";
    os << "Total: " << (valid_count + invalid_count) << "\n\n";

    // Show observation types and variables
    const auto& type_names = obs.getTypeNames();
    if (type_names.empty()) {
      os << "No observation types defined\n";
    } else {
      os << "Observation Types (" << type_names.size() << "):\n";
      os << std::string(50, '-') << "\n";

      for (const auto& type_name : type_names) {
        const auto& var_names = obs.getVariableNames(type_name);
        size_t total_obs_for_type = 0;

        // Count total observations for this type across all variables
        for (const auto& var_name : var_names) {
          total_obs_for_type += obs.getSize(type_name, var_name);
        }

        os << "Type: " << std::setw(15) << std::left << type_name;
        os << " | Variables: " << std::setw(3) << std::right
           << var_names.size();
        os << " | Total obs: " << std::setw(6) << std::right
           << total_obs_for_type << "\n";

        // Show variable details if there are any
        if (!var_names.empty()) {
          os << "  Variables: ";
          for (size_t i = 0; i < var_names.size(); ++i) {
            if (i > 0) os << ", ";
            os << var_names[i] << "(" << obs.getSize(type_name, var_names[i])
               << ")";
          }
          os << "\n";
        }
        os << "\n";
      }
    }

    // Show BUFR subset distribution
    if (!obs.bufr_subset_counts_.empty()) {
      os << "Overall BUFR Subset Distribution:\n";
      os << std::string(50, '-') << "\n";
      for (const auto& [subset, count] : obs.bufr_subset_counts_) {
        os << std::setw(15) << std::left << subset << ": " << count << "\n";
      }
      os << "\n";
    }

    // Show per-type BUFR subset distribution
    if (!obs.type_bufr_subset_counts_.empty()) {
      os << "Per-Type BUFR Subset Distribution:\n";
      os << std::string(50, '-') << "\n";
      for (const auto& [type_name, subset_counts] :
           obs.type_bufr_subset_counts_) {
        if (!subset_counts.empty()) {
          os << "Type: " << std::setw(15) << std::left << type_name << "\n";
          for (const auto& [subset, count] : subset_counts) {
            os << "  " << std::setw(15) << std::left << subset << ": " << count
               << "\n";
          }
          os << "\n";
        }
      }
    }

    // Show geographic distribution
    if (!obs.observations_.empty()) {
      os << "Geographic Distribution:\n";
      os << std::string(50, '-') << "\n";

      double min_lat = std::numeric_limits<double>::max();
      double max_lat = std::numeric_limits<double>::lowest();
      double min_lon = std::numeric_limits<double>::max();
      double max_lon = std::numeric_limits<double>::lowest();
      double min_level = std::numeric_limits<double>::max();
      double max_level = std::numeric_limits<double>::lowest();

      for (const auto& observation : obs.observations_) {
        if (observation.is_valid &&
            observation.location.getCoordinateSystem() ==
                framework::CoordinateSystem::GEOGRAPHIC) {
          auto [lat, lon, level] = observation.location.getGeographicCoords();
          min_lat = std::min(min_lat, lat);
          max_lat = std::max(max_lat, lat);
          min_lon = std::min(min_lon, lon);
          max_lon = std::max(max_lon, lon);
          min_level = std::min(min_level, level);
          max_level = std::max(max_level, level);
        }
      }

      if (min_lat != std::numeric_limits<double>::max()) {
        os << "Latitude range: [" << std::fixed << std::setprecision(2)
           << min_lat << ", " << max_lat << "]\n";
        os << "Longitude range: [" << std::fixed << std::setprecision(2)
           << min_lon << ", " << max_lon << "]\n";
        os << "Level range: [" << std::fixed << std::setprecision(2)
           << min_level << ", " << max_level << "]\n";
      }
      os << "\n";
    }

    // Show covariance information
    os << "Covariance Information:\n";
    os << std::string(50, '-') << "\n";
    os << "Diagonal: " << (obs.isDiagonalCovariance() ? "Yes" : "No") << "\n";

    const auto& cov = obs.getCovariance();
    if (!cov.empty()) {
      double min_error = *std::min_element(cov.begin(), cov.end());
      double max_error = *std::max_element(cov.begin(), cov.end());
      os << "Error variance range: [" << std::scientific << std::setprecision(2)
         << min_error << ", " << max_error << "]\n";
    }

    return os;
  }
};

}  // namespace metada::backends::common::observation
