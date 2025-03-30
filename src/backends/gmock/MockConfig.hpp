#pragma once

#include <gmock/gmock.h>

#include "utils/config/ConfigValue.hpp"
#include "utils/config/IConfig.hpp"

namespace metada::backends::gmock {

using framework::ConfigValue;
using framework::IConfig;

/**
 * @brief Mock configuration backend for testing
 *
 * Implements IConfig interface using Google Mock to provide mock configuration
 * methods. Used to verify Config's interaction with its backend implementation.
 */
class MockConfig : public IConfig {
 public:
  /**
   * @brief Default constructor
   */
  MockConfig() = default;

  /**
   * @brief Move constructor - explicitly defined for Google Mock compatibility
   */
  MockConfig(MockConfig&&) noexcept {}

  /**
   * @brief Move assignment - explicitly defined for Google Mock compatibility
   */
  MockConfig& operator=(MockConfig&&) noexcept { return *this; }

  MOCK_METHOD(bool, LoadFromFile, (const std::string&), (override));
  MOCK_METHOD(bool, LoadFromString, (const std::string&), (override));
  MOCK_METHOD(ConfigValue, Get, (const std::string&), (const, override));
  MOCK_METHOD(void, Set, (const std::string&, const ConfigValue&), (override));
  MOCK_METHOD(bool, HasKey, (const std::string&), (const, override));
  MOCK_METHOD(bool, SaveToFile, (const std::string&), (const, override));
  MOCK_METHOD(std::string, ToString, (), (const, override));
  MOCK_METHOD(void, Clear, (), (override));
};

}  // namespace metada::backends::gmock
