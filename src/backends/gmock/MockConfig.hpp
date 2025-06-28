#pragma once

#include <gmock/gmock.h>

#include "ConfigValue.hpp"

namespace metada::backends::gmock {

using framework::ConfigValue;

/**
 * @brief Mock configuration backend for testing
 *
 * @details This class implements the configuration backend contract using
 * Google Mock to provide mock configuration methods. It's used in unit tests to
 * verify the Config class's interaction with its backend implementation without
 * requiring actual file I/O or configuration parsing.
 *
 * The mock provides stubs for all required configuration backend methods
 * including loading, saving, getting, and setting configuration values.
 */
class MockConfig {
 public:
  /**
   * @brief Default constructor is disabled
   *
   * @details Configuration backends should always be initialized with a source.
   */
  MockConfig() = delete;

  /**
   * @brief Default destructor
   */
  ~MockConfig() = default;

  /**
   * @brief Copy constructor is disabled
   *
   * @details Configuration backends are not intended to be copied.
   */
  MockConfig(const MockConfig&) = delete;

  /**
   * @brief Copy assignment operator is disabled
   *
   * @details Configuration backends are not intended to be copied.
   */
  MockConfig& operator=(const MockConfig&) = delete;

  /**
   * @brief Move constructor - explicitly defined for Google Mock compatibility
   *
   * @details Google Mock requires move operations to be defined for proper test
   * setup.
   */
  MockConfig(MockConfig&&) noexcept {}

  /**
   * @brief Move assignment - explicitly defined for Google Mock compatibility
   *
   * @details Google Mock requires move operations to be defined for proper test
   * setup.
   *
   * @return Reference to this MockConfig instance
   */
  MockConfig& operator=(MockConfig&&) noexcept { return *this; }

  /**
   * @brief Constructor that loads configuration from a file
   *
   * @details Initializes the mock configuration by calling LoadFromFile with
   * the provided filename. In actual tests, the behavior of LoadFromFile would
   * be specified using EXPECT_CALL.
   *
   * @param filename Path to the configuration file
   * @throws std::runtime_error If loading fails (when configured in the test)
   */
  explicit MockConfig(const std::string& filename) { LoadFromFile(filename); }

  /**
   * @brief Constructor that loads configuration from a map
   *
   * @details Initializes the mock configuration by calling LoadFromMap with
   * the provided map. In actual tests, the behavior of LoadFromMap would
   * be specified using EXPECT_CALL.
   *
   * @param map The configuration map
   * @throws std::runtime_error If loading fails (when configured in the test)
   */
  explicit MockConfig(const framework::ConfigMap& map) { LoadFromMap(map); }

  MOCK_METHOD(bool, LoadFromFile, (const std::string&));
  MOCK_METHOD(bool, LoadFromString, (const std::string&));
  MOCK_METHOD(ConfigValue, Get, (const std::string&), (const));
  MOCK_METHOD(void, Set, (const std::string&, const ConfigValue&));
  MOCK_METHOD(bool, HasKey, (const std::string&), (const));
  MOCK_METHOD(bool, SaveToFile, (const std::string&), (const));
  MOCK_METHOD(std::string, ToString, (), (const));
  MOCK_METHOD(void, Clear, ());

  /**
   * @brief Create a new MockConfig object representing a subsection
   *
   * @details This method creates a new MockConfig object that represents a
   * subsection of the current configuration. In tests, the behavior of this
   * method can be specified using EXPECT_CALL.
   *
   * @param key Dot-separated path to the subsection
   * @return A new MockConfig object representing the subsection
   */
  MOCK_METHOD(MockConfig, CreateSubsection, (const std::string&), (const));

  /**
   * @brief Mock method to load configuration from a map
   *
   * @details This method is used to specify the behavior of LoadFromMap in
   * tests.
   *
   * @param map The configuration map
   */
  MOCK_METHOD(void, LoadFromMap, (const framework::ConfigMap&), ());
};

}  // namespace metada::backends::gmock
