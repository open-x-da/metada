#pragma once

#include "BackendTraits.hpp"
#include "GoogleLogger.hpp"
#include "JsonConfig.hpp"

namespace metada::traits {

struct L63BackendTag {};

template<>
struct BackendTraits<L63BackendTag> {
  using ConfigBackend = backends::config::JsonConfig;
  using LoggerBackend = backends::logger::GoogleLogger;
};

} // namespace metada::traits
