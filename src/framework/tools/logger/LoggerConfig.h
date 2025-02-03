#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERCONFIG_H_
#define METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERCONFIG_H_

#ifdef USE_GLOG
#include "google/GoogleLogger.h"
namespace metada {
namespace framework {
namespace tools {
namespace logger {
using LoggerBackend = backends::tools::logger::google_logger::GoogleLogger;
}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada
#else
#include "default/DefaultLogger.h"
namespace metada {
namespace framework {
namespace tools {
namespace logger {
using LoggerBackend = backends::tools::logger::default_logger::DefaultLogger;
}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada
#endif

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERCONFIG_H_