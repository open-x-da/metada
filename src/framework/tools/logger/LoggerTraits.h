#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_
#define METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_

namespace metada {
namespace framework {
namespace tools {
namespace logger {

template <typename Backend>
struct LoggerTraits {
  using LoggerBackend = Backend;
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_