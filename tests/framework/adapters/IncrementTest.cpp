#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Increment.hpp"
#include "MockConfig.hpp"
#include "MockIncrement.hpp"
#include "MockState.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Config;
using framework::Increment;

class IncrementTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Setup common test data
    dimensions_ = {10, 20};
    test_data_ = new double[10];  // Example data
    for (int i = 0; i < 10; i++) test_data_[i] = i;
  }

  void TearDown() override { delete[] test_data_; }

  std::vector<size_t> dimensions_;
  double* test_data_;
  Config<MockConfig> config;
};

TEST_F(IncrementTest, InitializeDelegatestoBackend) {
  MockIncrement mock_backend;
  Increment<MockIncrement> increment;

  EXPECT_CALL(increment.backend(), initialize()).Times(1);

  increment.initialize();
}

TEST_F(IncrementTest, ZeroDelegatestoBackend) {
  MockIncrement mock_backend;
  Increment<MockIncrement> increment;

  EXPECT_CALL(increment.backend(), zero()).Times(1);

  increment.zero();
}

TEST_F(IncrementTest, ScaleDelegatestoBackend) {
  MockIncrement mock_backend;
  Increment<MockIncrement> increment;

  EXPECT_CALL(increment.backend(), scale(2.0)).Times(1);

  increment.scale(2.0);
}

TEST_F(IncrementTest, LinearAlgebraOperations) {
  MockIncrement mock_backend;
  Increment<MockIncrement> increment;
  MockIncrement other_backend;
  Increment<MockIncrement> other;

  EXPECT_CALL(increment.backend(), axpy(2.0, _)).Times(1);
  EXPECT_CALL(increment.backend(), dot(_)).WillOnce(Return(1.5));
  EXPECT_CALL(increment.backend(), norm()).WillOnce(Return(2.0));

  increment.axpy(2.0, other);
  EXPECT_DOUBLE_EQ(increment.dot(other), 1.5);
  EXPECT_DOUBLE_EQ(increment.norm(), 2.0);
}

TEST_F(IncrementTest, GetDataReturnsTypedPointer) {
  MockIncrement mock_backend;
  Increment<MockIncrement> increment;

  EXPECT_CALL(increment.backend(), getData()).WillOnce(Return(test_data_));

  double* data = &increment.getData<double>();
  EXPECT_EQ(data, test_data_);
}

TEST_F(IncrementTest, MetadataOperations) {
  MockIncrement mock_backend;
  Increment<MockIncrement> increment;

  EXPECT_CALL(increment.backend(), setMetadata("key1", "value1")).Times(1);
  EXPECT_CALL(increment.backend(), getMetadata("key1"))
      .WillOnce(Return("value1"));

  increment.setMetadata("key1", "value1");
  EXPECT_EQ(increment.getMetadata("key1"), "value1");
}

TEST_F(IncrementTest, IncrementInformation) {
  MockIncrement mock_backend;
  Increment<MockIncrement> increment;

  EXPECT_CALL(increment.backend(), getDimensions())
      .WillOnce(ReturnRef(dimensions_));
  EXPECT_CALL(increment.backend(), isInitialized()).WillOnce(Return(true));

  EXPECT_EQ(increment.getDimensions(), dimensions_);
  EXPECT_TRUE(increment.isInitialized());
}

TEST_F(IncrementTest, BackendAccessors) {
  Increment<MockIncrement> increment;

  // Test that we can get non-const access to backend
  MockIncrement& backend = increment.backend();
  EXPECT_CALL(backend, initialize()).Times(1);
  backend.initialize();

  // Test that we can get const access to backend
  const Increment<MockIncrement>& const_increment = increment;
  const MockIncrement& const_backend = const_increment.backend();
  EXPECT_CALL(const_backend, isInitialized()).WillOnce(Return(true));
  const_backend.isInitialized();
}

}  // namespace metada::tests