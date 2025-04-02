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
   * @brief Disable default constructor
   */
  MockConfig() = delete;

  /**
   * @brief Default destructor
   */
  ~MockConfig() = default;

  /**
   * @brief Disable copy constructor
   */
  MockConfig(const MockConfig&) = delete;

  /**
   * @brief Disable copy assignment operator
   */
  MockConfig& operator=(const MockConfig&) = delete;

  /**
   * @brief Move constructor - explicitly defined for Google Mock compatibility
   */
  MockConfig(MockConfig&&) noexcept {}

  /**
   * @brief Move assignment - explicitly defined for Google Mock compatibility
   */
  MockConfig& operator=(MockConfig&&) noexcept { return *this; }

  /**
   * @brief Constructor that loads configuration from a file
   * @param filename Path to the configuration file
   * @throws std::runtime_error If loading fails
   */
  explicit MockConfig(const std::string& filename) { LoadFromFile(filename); }

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
