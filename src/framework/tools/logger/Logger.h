#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_LOGGER_H_
#define METADA_FRAMEWORK_TOOLS_LOGGER_LOGGER_H_

#include <string>

namespace metada {
namespace framework {
namespace tools {
namespace logger {

template <typename Backend>
class Logger {
 public:
  Logger() : backend_() {}
  ~Logger() = default;

  void Info(const std::string& message) { backend_.Info(message); }

  void Warning(const std::string& message) { backend_.Warning(message); }

  void Error(const std::string& message) { backend_.Error(message); }

  void Debug(const std::string& message) { backend_.Debug(message); }

 private:
  Backend backend_;
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_LOGGER_H_
