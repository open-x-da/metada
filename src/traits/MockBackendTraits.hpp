#pragma once

#include "BackendTraits.hpp"
#include "../backends/gmock/MockConfig.hpp"
#include "../backends/gmock/MockLogger.hpp"
#include "../backends/gmock/MockState.hpp"
//#include "mock_model_backend.hpp"
//#include "mock_geometry_backend.hpp"
//#include "mock_observation_backend.hpp"

namespace metada::traits {

struct MockBackendTag {};

template<>
struct BackendTraits<MockBackendTag> {
  using ConfigBackend       = backends::gmock::MockConfig;
  using LoggerBackend       = backends::gmock::MockLogger<ConfigBackend>;
  using StateBackend        = backends::gmock::MockState<ConfigBackend>;
  //using ModelBackend        = MockModelBackend;
  //using GeometryBackend     = MockGeometryBackend;
  //using ObservationBackend  = MockObservationBackend;
};

} // namespace metada::traits
