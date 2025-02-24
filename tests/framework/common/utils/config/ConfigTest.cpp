/**
 * @file ConfigTest.cpp
 * @brief Unit tests for the Config class template
 *
 * This test suite verifies the Config class template correctly delegates
 * configuration operations to its backend implementation. It uses Google Mock
 * to create a mock config backend and validate the interactions between Config
 * and backend.
 *
 * Key test areas:
 * - Configuration loading from files and strings
 * - Getting and setting values of different types
 * - Array value handling
 * - Error cases and default values
 * - Configuration persistence
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MockConfig.hpp"
#include "utils/config/Config.hpp"

namespace metada::framework::common::utils::config::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::Throw;

/**
 * @brief Test fixture for Config class tests
 */
class ConfigTest : public ::testing::Test {
 protected:
  /** @brief Config instance with mock backend */
  Config<MockConfig> config;
};

/**
 * @brief Verify LoadFromFile() method delegation
 *
 * Tests that Config::LoadFromFile() properly delegates to backend's
 * LoadFromFile() method
 */
TEST_F(ConfigTest, LoadFromFile) {
  EXPECT_CALL(config.backend(), LoadFromFile("test.yaml"))
      .WillOnce(Return(true));
  EXPECT_TRUE(config.LoadFromFile("test.yaml"));
}

/**
 * @brief Verify LoadFromString() method delegation
 *
 * Tests that Config::LoadFromString() properly delegates to backend's
 * LoadFromString() method
 */
TEST_F(ConfigTest, LoadFromString) {
  EXPECT_CALL(config.backend(), LoadFromString("content"))
      .WillOnce(Return(true));
  EXPECT_TRUE(config.LoadFromString("content"));
}

/**
 * @brief Test Get() method with string values
 */
TEST_F(ConfigTest, GetString) {
  ConfigValue value = std::string("test");
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  EXPECT_EQ(std::get<std::string>(config.Get("key")), "test");
}

/**
 * @brief Test Get() method with integer values
 */
TEST_F(ConfigTest, GetInt) {
  ConfigValue value = 42;
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  EXPECT_EQ(std::get<int>(config.Get("key")), 42);
}

/**
 * @brief Test Get() method with double values
 */
TEST_F(ConfigTest, GetDouble) {
  ConfigValue value = 3.14;
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  EXPECT_DOUBLE_EQ(std::get<double>(config.Get("key")), 3.14);
}

/**
 * @brief Test Get() method with boolean values
 */
TEST_F(ConfigTest, GetBool) {
  ConfigValue value = true;
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  EXPECT_TRUE(std::get<bool>(config.Get("key")));
}

/**
 * @brief Test Get() method with string arrays
 */
TEST_F(ConfigTest, GetStringArray) {
  ConfigValue value = std::vector<std::string>{"a", "b", "c"};
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  auto result = std::get<std::vector<std::string>>(config.Get("key"));
  EXPECT_THAT(result, ::testing::ElementsAre("a", "b", "c"));
}

/**
 * @brief Test Get() method with integer arrays
 */
TEST_F(ConfigTest, GetIntArray) {
  ConfigValue value = std::vector<int>{1, 2, 3};
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  auto result = std::get<std::vector<int>>(config.Get("key"));
  EXPECT_THAT(result, ::testing::ElementsAre(1, 2, 3));
}

/**
 * @brief Test Get() method with double arrays
 */
TEST_F(ConfigTest, GetDoubleArray) {
  ConfigValue value = std::vector<double>{1.1, 2.2, 3.3};
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  auto result = std::get<std::vector<double>>(config.Get("key"));
  EXPECT_THAT(result, ::testing::ElementsAre(1.1, 2.2, 3.3));
}

/**
 * @brief Test Get() method with boolean arrays
 */
TEST_F(ConfigTest, GetBoolArray) {
  ConfigValue value = std::vector<bool>{true, false, true};
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  auto result = std::get<std::vector<bool>>(config.Get("key"));
  std::vector<bool> expected{true, false, true};
  EXPECT_EQ(result, expected);
}

/**
 * @brief Test Set() method with string values
 */
TEST_F(ConfigTest, SetString) {
  ConfigValue value = std::string("test");
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

/**
 * @brief Test Set() method with integer values
 */
TEST_F(ConfigTest, SetInt) {
  ConfigValue value = 42;
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

/**
 * @brief Test Set() method with double values
 */
TEST_F(ConfigTest, SetDouble) {
  ConfigValue value = 3.14;
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

/**
 * @brief Test Set() method with boolean values
 */
TEST_F(ConfigTest, SetBool) {
  ConfigValue value = true;
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

/**
 * @brief Test Set() method with string arrays
 */
TEST_F(ConfigTest, SetStringArray) {
  ConfigValue value = std::vector<std::string>{"a", "b", "c"};
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

/**
 * @brief Test Set() method with integer arrays
 */
TEST_F(ConfigTest, SetIntArray) {
  ConfigValue value = std::vector<int>{1, 2, 3};
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

/**
 * @brief Test Set() method with double arrays
 */
TEST_F(ConfigTest, SetDoubleArray) {
  ConfigValue value = std::vector<double>{1.1, 2.2, 3.3};
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

/**
 * @brief Test Set() method with boolean arrays
 */
TEST_F(ConfigTest, SetBoolArray) {
  ConfigValue value = std::vector<bool>{true, false, true};
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

/**
 * @brief Test HasKey() method for existing and non-existing keys
 */
TEST_F(ConfigTest, HasKey) {
  EXPECT_CALL(config.backend(), HasKey("existing")).WillOnce(Return(true));
  EXPECT_TRUE(config.HasKey("existing"));

  EXPECT_CALL(config.backend(), HasKey("nonexistent")).WillOnce(Return(false));
  EXPECT_FALSE(config.HasKey("nonexistent"));
}

/**
 * @brief Test SaveToFile() method delegation
 */
TEST_F(ConfigTest, SaveToFile) {
  EXPECT_CALL(config.backend(), SaveToFile("test.yaml")).WillOnce(Return(true));
  EXPECT_TRUE(config.SaveToFile("test.yaml"));
}

/**
 * @brief Test ToString() method delegation
 */
TEST_F(ConfigTest, ToString) {
  EXPECT_CALL(config.backend(), ToString()).WillOnce(Return("config content"));
  EXPECT_EQ(config.ToString(), "config content");
}

/**
 * @brief Test Clear() method delegation
 */
TEST_F(ConfigTest, Clear) {
  EXPECT_CALL(config.backend(), Clear());
  config.Clear();
}

/**
 * @brief Test LoadFromFile() failure case
 */
TEST_F(ConfigTest, LoadFromFileFailure) {
  EXPECT_CALL(config.backend(), LoadFromFile("nonexistent.yaml"))
      .WillOnce(Return(false));
  EXPECT_FALSE(config.LoadFromFile("nonexistent.yaml"));
}

/**
 * @brief Test LoadFromString() failure case
 */
TEST_F(ConfigTest, LoadFromStringFailure) {
  EXPECT_CALL(config.backend(), LoadFromString("invalid content"))
      .WillOnce(Return(false));
  EXPECT_FALSE(config.LoadFromString("invalid content"));
}

/**
 * @brief Test SaveToFile() failure case
 */
TEST_F(ConfigTest, SaveToFileFailure) {
  EXPECT_CALL(config.backend(), SaveToFile("/invalid/path/test.yaml"))
      .WillOnce(Return(false));
  EXPECT_FALSE(config.SaveToFile("/invalid/path/test.yaml"));
}

/**
 * @brief Test Get() with invalid key
 */
TEST_F(ConfigTest, GetInvalidKey) {
  EXPECT_CALL(config.backend(), Get("nonexistent"))
      .WillOnce(::testing::Throw(std::runtime_error("Key not found")));
  EXPECT_THROW(config.GetUnsafe("nonexistent"), std::runtime_error);
}

/**
 * @brief Test Get() with default value
 */
TEST_F(ConfigTest, GetDefaultValue) {
  ConfigValue defaultValue = std::string("default");
  EXPECT_CALL(config.backend(), Get("nonexistent"))
      .WillOnce(::testing::Throw(std::runtime_error("Key not found")));
  auto result = config.Get("nonexistent", defaultValue);
  EXPECT_EQ(std::get<std::string>(result), "default");
}

/**
 * @brief Test Get() safe access returning default on error
 */
TEST_F(ConfigTest, GetSafeReturnsDefaultOnError) {
  EXPECT_CALL(config.backend(), Get("nonexistent"))
      .WillOnce(::testing::Throw(std::runtime_error("Key not found")));
  ConfigValue result = config.Get("nonexistent");
  EXPECT_EQ(result,
            ConfigValue());  // Should return default-constructed ConfigValue
}

}  // namespace metada::framework::common::utils::config::tests