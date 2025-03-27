#pragma once

#include "backend_traits.hpp"
#include "mock_state_backend.hpp"
#include "mock_config_backend.hpp"
#include "mock_logger_backend.hpp"
#include "mock_model_backend.hpp"
#include "mock_geometry_backend.hpp"
#include "mock_observation_backend.hpp"

namespace metada::traits {

template<>
struct BackendTraits<MockBackendTag> {
  using ConfigBackend       = MockConfigBackend;
  using LoggerBackend       = MockLoggerBackend;
  using StateBackend        = MockStateBackend;
  using ModelBackend        = MockModelBackend;
  using GeometryBackend     = MockGeometryBackend;
  using ObservationBackend  = MockObservationBackend;
};

} // namespace metada::traits
