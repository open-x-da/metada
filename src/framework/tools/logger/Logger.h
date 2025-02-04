#ifndef METADA_FRAMEWORK_TOOLS_LOGGER_LOGGER_H_
#define METADA_FRAMEWORK_TOOLS_LOGGER_LOGGER_H_

#include <string>

namespace metada {

/**
 * @brief Generic logger class that delegates to a backend implementation
 * @tparam Backend The concrete logger backend type
 */
template <typename Backend>
class Logger {
 public:
  Logger() = default;
  ~Logger() = default;

  /**
   * @brief Log an informational message
   * @param message The message to log
   */
  void Info(const std::string& message) { backend_.Info(message); }

  /**
   * @brief Log a warning message
   * @param message The message to log
   */
  void Warning(const std::string& message) { backend_.Warning(message); }

  /**
   * @brief Log an error message
   * @param message The message to log
   */
  void Error(const std::string& message) { backend_.Error(message); }

  /**
   * @brief Log a debug message
   * @param message The message to log
   */
  void Debug(const std::string& message) { backend_.Debug(message); }

 protected:
  Backend backend_;  ///< The underlying logger backend instance
};

}  // namespace metada

#endif  // METADA_FRAMEWORK_TOOLS_LOGGER_LOGGER_H_
