#pragma once

#include <string>
#include <vector>

#include "../observation/PrepBUFRObservation.hpp"
#include "ObsRecord.hpp"

namespace metada::backends::common::io {

using ObsRecord = framework::ObsRecord;

// Re-export PrepBUFR batching types in the IO namespace for backward
// compatibility.
using SurfaceBatch = metada::backends::common::observation::SurfaceBatch;
using ProfileBatch = metada::backends::common::observation::ProfileBatch;
using FamilyBatches = metada::backends::common::observation::FamilyBatches;

// Thin wrapper that forwards to the centralized implementation under
// observation/.
class PrepBUFRObsAdapter {
 public:
  static FamilyBatches groupByFamily(const std::vector<ObsRecord>& records) {
    return metada::backends::common::observation::PrepBUFRObservation::
        groupByFamily(records);
  }
};

}  // namespace metada::backends::common::io
