#pragma once

#include "BackendTraits.hpp"
#include "MockConfig.hpp"
// #include "MockLogger.hpp"
//#include "mock_logger_backend.hpp"
//#include "mock_model_backend.hpp"
//#include "mock_geometry_backend.hpp"
//#include "mock_observation_backend.hpp"

namespace metada::traits {

struct MockBackendTag {};

template<>
struct BackendTraits<MockBackendTag> {
  using ConfigBackend       = MockConfig;
  using LoggerBackend       = MockLogger;
  //using StateBackend        = MockStateBackend;
  //using ModelBackend        = MockModelBackend;
  //using GeometryBackend     = MockGeometryBackend;
  //using ObservationBackend  = MockObservationBackend;
};

} // namespace metada::traits
