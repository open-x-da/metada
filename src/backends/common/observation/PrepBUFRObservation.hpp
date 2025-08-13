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
#include "ObsRecord.hpp"
#include "PointObservation.hpp"
#include "PrepBUFRObservationIterator.hpp"

namespace metada::backends::common::observation {

using ObsRecord = framework::ObsRecord;
using ObsLevelRecord = framework::ObsLevelRecord;
using ConfigValue = metada::backends::config::ConfigValue;

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
      } catch (...) {
        // Some keys might not exist, that's okay
      }

      file_groups[file_path].emplace_back(type_name, config_data);
    }

    // Process each file once and distribute records to appropriate types
    for (const auto& [file_path, type_configs] : file_groups) {
      // Read BUFR file once for all types that share it
      auto bufr_io =
          io::BufrObsIO<ConfigBackend>(ConfigBackend({{"file", file_path}}));
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

        // Flatten BUFR records into observation points for selected variables
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
                  observations_.emplace_back(location, v.value, sel.error);
                  type_variable_map_[type_name][user_var].push_back(
                      observations_.size() - 1);
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  // Move semantics
  PrepBUFRObservation(PrepBUFRObservation&& other) noexcept
      : observations_(std::move(other.observations_)),
        type_variable_map_(std::move(other.type_variable_map_)),
        covariance_(std::move(other.covariance_)) {}

  PrepBUFRObservation& operator=(PrepBUFRObservation&& other) noexcept {
    if (this != &other) {
      observations_ = std::move(other.observations_);
      type_variable_map_ = std::move(other.type_variable_map_);
      covariance_ = std::move(other.covariance_);
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
  // Size for a given type/variable (fallback: entire flattened vector)
  size_t getSize(const std::string& /*type_name*/,
                 const std::string& /*var_name*/) const {
    return observations_.size();
  }

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

  bool isDiagonalCovariance() const { return true; }
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

  // Data members
  std::vector<framework::ObservationPoint> observations_;
  std::unordered_map<std::string,
                     std::unordered_map<std::string, std::vector<size_t>>>
      type_variable_map_;
  std::vector<double> covariance_;
};

}  // namespace metada::backends::common::observation
