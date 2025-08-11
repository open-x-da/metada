#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ObsRecord.hpp"

namespace metada::backends::common::io {

using ObsRecord = framework::ObsRecord;
using ObsRecordShared = framework::ObsRecordShared;
using ObsLevelRecord = framework::ObsLevelRecord;

struct SurfaceBatch {
  // per obs
  std::vector<double> lats;
  std::vector<double> lons;
  std::vector<double> levels;  // pressure or model level proxy
  // output buffer
  std::vector<double> y;
};

struct ProfileBatch {
  std::vector<int> counts;          // levels per obs
  std::vector<double> lats;         // per obs
  std::vector<double> lons;         // per obs
  std::vector<double> levels_flat;  // flattened levels
  std::vector<double> y_flat;       // flattened outputs
};

struct FamilyBatches {
  SurfaceBatch metar, synop, ships, buoy, sonde_sfc;
  ProfileBatch airep, pilot, sound;
};

// Map PrepBUFR ObsRecords into batches by family and variable tokens.
// Variables tokens per family should align with WRFDA y fields:
// surface: u, v, t, q, p; profiles: u, v, t, q.
class PrepBUFRObsAdapter {
 public:
  static FamilyBatches groupByFamily(const std::vector<ObsRecord>& records) {
    FamilyBatches batches;
    for (const auto& rec : records) {
      const double lat = rec.shared.latitude;
      const double lon = rec.shared.longitude;
      const auto subset = rec.shared.report_type;  // user may remap if needed
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
    int nlev = static_cast<int>(rec.levels.size());
    if (nlev <= 0) return;
    dst.lats.push_back(lat);
    dst.lons.push_back(lon);
    dst.counts.push_back(nlev);
    for (const auto& lvl : rec.levels) {
      double p = findValue(lvl, "PRES");
      dst.levels_flat.push_back(p);
    }
  }

  static void appendProfileWinds(const ObsRecord& rec, double lat, double lon,
                                 ProfileBatch& dst) {
    // Similar to appendProfile; keep counts and levels for winds
    appendProfile(rec, lat, lon, dst);
  }

  static double findValue(const std::vector<ObsLevelRecord>& lvl,
                          const std::string& key) {
    for (const auto& v : lvl) {
      if (v.type == key) return v.value;
    }
    return 0.0;
  }
};

}  // namespace metada::backends::common::io
