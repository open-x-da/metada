#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_
#define METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_

#ifdef USE_GLOG
#include "google/GoogleLogger.h"
#else
#include "default/DefaultLogger.h"
#endif

namespace metada {
namespace framework {
namespace tools {
namespace logger {

#ifdef USE_GLOG
using LoggerBackend = backends::tools::logger::google_logger::GoogleLogger;
#else
using LoggerBackend = backends::tools::logger::default_logger::DefaultLogger;
#endif

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_