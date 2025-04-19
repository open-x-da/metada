#pragma once

#include "BackendTraits.hpp"
#include "../backends/gmock/MockConfig.hpp"
#include "../backends/gmock/MockLogger.hpp"
#include "../backends/gmock/MockGeometry.hpp"
#include "../backends/gmock/MockGeometryIterator.hpp"
#include "../backends/gmock/MockState.hpp"
#include "../backends/gmock/MockModel.hpp"

namespace metada::traits {

struct MockBackendTag {};

template<>
struct BackendTraits<MockBackendTag> {
  using ConfigBackend       = backends::gmock::MockConfig;
  using LoggerBackend       = backends::gmock::MockLogger<ConfigBackend>;
  using GeometryBackend     = backends::gmock::MockGeometry<ConfigBackend>;
  using GeometryIteratorBackend = backends::gmock::MockGeometryIterator;
  using StateBackend        = backends::gmock::MockState<ConfigBackend>;
  using ModelBackend        = backends::gmock::MockModel<ConfigBackend, StateBackend>;
  //using ObservationBackend  = MockObservationBackend;
};

} // namespace metada::traits
