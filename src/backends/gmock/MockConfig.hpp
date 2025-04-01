#pragma once

#include <gmock/gmock.h>
#include <filesystem>

#include "common/utils/config/ConfigValue.hpp"

namespace metada::backends::gmock {

using framework::ConfigValue;

/**
 * @brief Mock configuration backend for testing
 *
 * Implements the configuration backend contract using Google Mock to provide mock
 * configuration methods. Used to verify Config's interaction with its backend
 * implementation.
 */
class MockConfig {
 public:
  /**
   * @brief Default constructor
   */
  MockConfig() {
    // Set up default behavior for LoadFromFile to check if file exists
    ON_CALL(*this, LoadFromFile(::testing::_))
        .WillByDefault(::testing::Invoke([](const std::string& filename) {
          return std::filesystem::exists(filename);
        }));
  }

  /**
   * @brief Move constructor - explicitly defined for Google Mock compatibility
   */
  MockConfig(MockConfig&&) noexcept {}

  /**
   * @brief Move assignment - explicitly defined for Google Mock compatibility
   */
  MockConfig& operator=(MockConfig&&) noexcept { return *this; }

  MOCK_METHOD(bool, LoadFromFile, (const std::string&));
  MOCK_METHOD(bool, LoadFromString, (const std::string&));
  MOCK_METHOD(ConfigValue, Get, (const std::string&), (const));
  MOCK_METHOD(void, Set, (const std::string&, const ConfigValue&));
  MOCK_METHOD(bool, HasKey, (const std::string&), (const));
  MOCK_METHOD(bool, SaveToFile, (const std::string&), (const));
  MOCK_METHOD(std::string, ToString, (), (const));
  MOCK_METHOD(void, Clear, ());
};

}  // namespace metada::backends::gmock
