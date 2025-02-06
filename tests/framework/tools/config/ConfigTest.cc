#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Config.h"
#include "IConfig.h"

using ::testing::_;
using ::testing::Return;

namespace metada {
namespace framework {
namespace tools {
namespace config {

// Mock configuration backend for testing
class MockConfigBackend : public IConfig {
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

// Test fixture
class ConfigTest : public ::testing::Test {
 protected:
  Config<MockConfigBackend> config;
};

// Test LoadFromFile
TEST_F(ConfigTest, LoadFromFile) {
  EXPECT_CALL(config.backend(), LoadFromFile("test.yaml"))
      .WillOnce(Return(true));
  EXPECT_TRUE(config.LoadFromFile("test.yaml"));
}

// Test LoadFromString
TEST_F(ConfigTest, LoadFromString) {
  EXPECT_CALL(config.backend(), LoadFromString("content"))
      .WillOnce(Return(true));
  EXPECT_TRUE(config.LoadFromString("content"));
}

// Test Get with different types
TEST_F(ConfigTest, GetString) {
  ConfigValue value = std::string("test");
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  EXPECT_EQ(std::get<std::string>(config.Get("key")), "test");
}

TEST_F(ConfigTest, GetInt) {
  ConfigValue value = 42;
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  EXPECT_EQ(std::get<int>(config.Get("key")), 42);
}

TEST_F(ConfigTest, GetDouble) {
  ConfigValue value = 3.14;
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  EXPECT_DOUBLE_EQ(std::get<double>(config.Get("key")), 3.14);
}

TEST_F(ConfigTest, GetBool) {
  ConfigValue value = true;
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  EXPECT_TRUE(std::get<bool>(config.Get("key")));
}

// Test Get with arrays
TEST_F(ConfigTest, GetStringArray) {
  ConfigValue value = std::vector<std::string>{"a", "b", "c"};
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  auto result = std::get<std::vector<std::string>>(config.Get("key"));
  EXPECT_THAT(result, ::testing::ElementsAre("a", "b", "c"));
}

TEST_F(ConfigTest, GetIntArray) {
  ConfigValue value = std::vector<int>{1, 2, 3};
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  auto result = std::get<std::vector<int>>(config.Get("key"));
  EXPECT_THAT(result, ::testing::ElementsAre(1, 2, 3));
}

TEST_F(ConfigTest, GetDoubleArray) {
  ConfigValue value = std::vector<double>{1.1, 2.2, 3.3};
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  auto result = std::get<std::vector<double>>(config.Get("key"));
  EXPECT_THAT(result, ::testing::ElementsAre(1.1, 2.2, 3.3));
}

TEST_F(ConfigTest, GetBoolArray) {
  ConfigValue value = std::vector<bool>{true, false, true};
  EXPECT_CALL(config.backend(), Get("key")).WillOnce(Return(value));
  auto result = std::get<std::vector<bool>>(config.Get("key"));
  std::vector<bool> expected{true, false, true};
  EXPECT_EQ(result, expected);
}

// Test Set with different types
TEST_F(ConfigTest, SetString) {
  ConfigValue value = std::string("test");
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

TEST_F(ConfigTest, SetInt) {
  ConfigValue value = 42;
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

TEST_F(ConfigTest, SetDouble) {
  ConfigValue value = 3.14;
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

TEST_F(ConfigTest, SetBool) {
  ConfigValue value = true;
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

// Test Set with arrays
TEST_F(ConfigTest, SetStringArray) {
  ConfigValue value = std::vector<std::string>{"a", "b", "c"};
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

TEST_F(ConfigTest, SetIntArray) {
  ConfigValue value = std::vector<int>{1, 2, 3};
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

TEST_F(ConfigTest, SetDoubleArray) {
  ConfigValue value = std::vector<double>{1.1, 2.2, 3.3};
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

TEST_F(ConfigTest, SetBoolArray) {
  ConfigValue value = std::vector<bool>{true, false, true};
  EXPECT_CALL(config.backend(), Set("key", value));
  config.Set("key", value);
}

// Test HasKey
TEST_F(ConfigTest, HasKey) {
  EXPECT_CALL(config.backend(), HasKey("existing")).WillOnce(Return(true));
  EXPECT_TRUE(config.HasKey("existing"));

  EXPECT_CALL(config.backend(), HasKey("nonexistent")).WillOnce(Return(false));
  EXPECT_FALSE(config.HasKey("nonexistent"));
}

// Test SaveToFile
TEST_F(ConfigTest, SaveToFile) {
  EXPECT_CALL(config.backend(), SaveToFile("test.yaml")).WillOnce(Return(true));
  EXPECT_TRUE(config.SaveToFile("test.yaml"));
}

// Test ToString
TEST_F(ConfigTest, ToString) {
  EXPECT_CALL(config.backend(), ToString()).WillOnce(Return("config content"));
  EXPECT_EQ(config.ToString(), "config content");
}

// Test Clear
TEST_F(ConfigTest, Clear) {
  EXPECT_CALL(config.backend(), Clear());
  config.Clear();
}

// Test error cases
TEST_F(ConfigTest, LoadFromFileFailure) {
  EXPECT_CALL(config.backend(), LoadFromFile("nonexistent.yaml"))
      .WillOnce(Return(false));
  EXPECT_FALSE(config.LoadFromFile("nonexistent.yaml"));
}

TEST_F(ConfigTest, LoadFromStringFailure) {
  EXPECT_CALL(config.backend(), LoadFromString("invalid content"))
      .WillOnce(Return(false));
  EXPECT_FALSE(config.LoadFromString("invalid content"));
}

TEST_F(ConfigTest, SaveToFileFailure) {
  EXPECT_CALL(config.backend(), SaveToFile("/invalid/path/test.yaml"))
      .WillOnce(Return(false));
  EXPECT_FALSE(config.SaveToFile("/invalid/path/test.yaml"));
}

TEST_F(ConfigTest, GetInvalidKey) {
  EXPECT_CALL(config.backend(), Get("nonexistent"))
      .WillOnce(::testing::Throw(std::runtime_error("Key not found")));
  EXPECT_THROW(config.GetUnsafe("nonexistent"), std::runtime_error);
}

TEST_F(ConfigTest, GetDefaultValue) {
  ConfigValue defaultValue = std::string("default");
  EXPECT_CALL(config.backend(), Get("nonexistent"))
      .WillOnce(::testing::Throw(std::runtime_error("Key not found")));
  auto result = config.Get("nonexistent", defaultValue);
  EXPECT_EQ(std::get<std::string>(result), "default");
}

TEST_F(ConfigTest, GetSafeReturnsDefaultOnError) {
  EXPECT_CALL(config.backend(), Get("nonexistent"))
      .WillOnce(::testing::Throw(std::runtime_error("Key not found")));
  ConfigValue result = config.Get("nonexistent");
  EXPECT_EQ(result,
            ConfigValue());  // Should return default-constructed ConfigValue
}

}  // namespace config
}  // namespace tools
}  // namespace framework
}  // namespace metada