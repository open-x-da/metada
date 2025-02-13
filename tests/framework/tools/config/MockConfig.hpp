#ifndef METADA_TESTS_FRAMEWORK_TOOLS_CONFIG_MOCK_CONFIG_HPP_
#define METADA_TESTS_FRAMEWORK_TOOLS_CONFIG_MOCK_CONFIG_HPP_

#include <gmock/gmock.h>

#include "IConfig.hpp"

namespace metada {
namespace framework {
namespace tools {
namespace config {
namespace tests {

/**
 * @brief Mock configuration backend for testing
 *
 * Implements IConfig interface using Google Mock to provide mock configuration
 * methods. Used to verify Config's interaction with its backend implementation.
 */
class MockConfig : public IConfig {
 public:
  MOCK_METHOD(bool, LoadFromFile, (const std::string&), (override));
  MOCK_METHOD(bool, LoadFromString, (const std::string&), (override));
  MOCK_METHOD(ConfigValue, Get, (const std::string&), (const, override));
  MOCK_METHOD(void, Set, (const std::string&, const ConfigValue&), (override));
  MOCK_METHOD(bool, HasKey, (const std::string&), (const, override));
  MOCK_METHOD(bool, SaveToFile, (const std::string&), (const, override));
  MOCK_METHOD(std::string, ToString, (), (const, override));
  MOCK_METHOD(void, Clear, (), (override));
};

}  // namespace tests
}  // namespace config
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_TESTS_FRAMEWORK_TOOLS_CONFIG_MOCK_CONFIG_HPP_