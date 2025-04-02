#pragma once

#include "BackendTraits.hpp"
#include "../backends/common/utils/config/json/JsonConfig.hpp"
#include "../backends/common/utils/logger/glog/GoogleLogger.hpp"

namespace metada::traits {

struct L63BackendTag {};

template<>
struct BackendTraits<L63BackendTag> {
  using ConfigBackend = backends::config::JsonConfig;
  using LoggerBackend = backends::logger::GoogleLogger<ConfigBackend>;
};

} // namespace metada::traits
