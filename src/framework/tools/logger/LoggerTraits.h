#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_
#define METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_

namespace metada {

template <typename Backend>
struct LoggerTraits {
  using LoggerBackend = Backend;
};

}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_LOGGERTRAITS_H_