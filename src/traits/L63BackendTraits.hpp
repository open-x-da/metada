#pragma once

#include "backend_traits.hpp"
#include "wrf_state_backend.hpp"
#include "yaml_config_backend.hpp"
#include "google_logger_backend.hpp"
#include "wrf_model_backend.hpp"
#include "wrf_geometry_backend.hpp"
#include "wrf_observation_backend.hpp"

namespace metada::traits {

template<>
struct BackendTraits<L63BackendTag> {
  using ConfigBackend       = YamlConfigBackend;
  using LoggerBackend       = GoogleLoggerBackend;
  using StateBackend        = WrfStateBackend;
  using ModelBackend        = WrfModelBackend;
  using GeometryBackend     = WrfGeometryBackend;
  using ObservationBackend  = WrfObservationBackend;
};

} // namespace metada::traits
