#ifndef METADA_TESTS_FRAMEWORK_COMMON_UTILS_LOGGER_MOCK_LOGGER_HPP_
#define METADA_TESTS_FRAMEWORK_COMMON_UTILS_LOGGER_MOCK_LOGGER_HPP_

#include <gmock/gmock.h>

#include "utils/logger/ILogger.hpp"

namespace metada {
namespace framework {
namespace common {
namespace utils {
namespace logger {
namespace tests {

/**
 * @brief Mock logger backend for testing
 *
 * Implements ILogger interface using Google Mock to provide mock logging
 * methods. Used to verify Logger's interaction with its backend implementation.
 */
class MockLogger : public ILogger {
 public:
  MOCK_METHOD(void, Info, (const std::string& message), (override));
  MOCK_METHOD(void, Warning, (const std::string& message), (override));
  MOCK_METHOD(void, Error, (const std::string& message), (override));
  MOCK_METHOD(void, Debug, (const std::string& message), (override));

  /**
   * @brief Mock initialization method
   * @param app_name Application name (unused in mock)
   */
  static void Init([[maybe_unused]] const std::string& app_name) {}

  /**
   * @brief Mock shutdown method
   */
  static void Shutdown() {}
};

}  // namespace tests
}  // namespace logger
}  // namespace utils
}  // namespace common
}  // namespace framework
}  // namespace metada

#endif  // METADA_TESTS_FRAMEWORK_COMMON_UTILS_LOGGER_MOCK_LOGGER_HPP_