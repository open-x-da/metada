#ifndef METADA_LOGGER_H
#define METADA_LOGGER_H

#include <string>

namespace metada {
namespace logging {

template <typename Backend>
class Logger {
 public:
  Logger() = default;
  ~Logger() = default;

  void Info(const std::string& message) { Backend::LogInfo(message); }

  void Warning(const std::string& message) { Backend::LogWarning(message); }

  void Error(const std::string& message) { Backend::LogError(message); }

  void Debug(const std::string& message) { Backend::LogDebug(message); }
};

}  // namespace logging
}  // namespace metada

#endif  // METADA_LOGGER_H
