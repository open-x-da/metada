#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "MockObservation.hpp"
#include "Observation.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Const;
using ::testing::Return;
using ::testing::ReturnRef;

using metada::framework::Observation;

class ObservationTest : public ::testing::Test {
 protected:
  // Test data
  std::vector<double> location_{45.0, -120.0, 100.0};
  double timestamp_{1234567890.0};
  double temperature_{25.5};
  double uncertainty_{0.5};

  // Setup mock expectations for common operations
  void SetUp() override {
    // Setup mock to return our test data
    EXPECT_CALL(mockBackend_, isValid()).WillRepeatedly(Return(true));
    EXPECT_CALL(mockBackend_, getLocation())
        .WillRepeatedly(ReturnRef(location_));
    EXPECT_CALL(mockBackend_, getTimestamp())
        .WillRepeatedly(Return(timestamp_));
    EXPECT_CALL(mockBackend_, getQualityFlag()).WillRepeatedly(Return(1));
    EXPECT_CALL(mockBackend_, getConfidence()).WillRepeatedly(Return(0.95));
    EXPECT_CALL(mockBackend_, getObsType())
        .WillRepeatedly(Return("temperature"));
    EXPECT_CALL(mockBackend_, getSource()).WillRepeatedly(Return("satellite"));
    EXPECT_CALL(mockBackend_, getInstrument()).WillRepeatedly(Return("AMSU-A"));

    // Setup data access
    EXPECT_CALL(mockBackend_, getData())
        .WillRepeatedly(Return(static_cast<void*>(&temperature_)));
    ON_CALL(testing::Const(mockBackend_), getData())
        .WillByDefault(Return(static_cast<const void*>(&temperature_)));
    EXPECT_CALL(mockBackend_, getUncertainty())
        .WillRepeatedly(Return(static_cast<void*>(&uncertainty_)));
    ON_CALL(testing::Const(mockBackend_), getUncertainty())
        .WillByDefault(Return(static_cast<const void*>(&uncertainty_)));
    EXPECT_CALL(mockBackend_, getDataSize()).WillRepeatedly(Return(1));
  }

  MockObservation mockBackend_;
  Observation<MockObservation> observation_;
};

TEST_F(ObservationTest, InitializeCallsBackend) {
  EXPECT_CALL(mockBackend_, initialize()).Times(1);
  observation_.initialize();
}

TEST_F(ObservationTest, GetDataReturnsCorrectValue) {
  ASSERT_EQ(observation_.getData<double>(), temperature_);
}

TEST_F(ObservationTest, GetUncertaintyReturnsCorrectValue) {
  ASSERT_EQ(observation_.getUncertainty<double>(), uncertainty_);
}

TEST_F(ObservationTest, LocationIsCorrectlyReturned) {
  const auto& loc = observation_.getLocation();
  ASSERT_EQ(loc.size(), 3);
  ASSERT_EQ(loc[0], 45.0);
  ASSERT_EQ(loc[1], -120.0);
  ASSERT_EQ(loc[2], 100.0);
}

TEST_F(ObservationTest, SetLocationCallsBackend) {
  EXPECT_CALL(mockBackend_, setLocation(50.0, -110.0, 200.0)).Times(1);
  observation_.setLocation(50.0, -110.0, 200.0);
}

TEST_F(ObservationTest, FluentInterfaceWorks) {
  EXPECT_CALL(mockBackend_, setLocation(_, _, _)).Times(1);
  EXPECT_CALL(mockBackend_, setTime(_)).Times(1);
  EXPECT_CALL(mockBackend_, setQualityFlag(_)).Times(1);

  // Test chaining of methods
  observation_.setLocation(50.0, -110.0, 200.0)
      .setTime(1234567891.0)
      .setQualityFlag(2);
}

TEST_F(ObservationTest, AttributesWorkCorrectly) {
  EXPECT_CALL(mockBackend_, setAttribute("key1", "value1")).Times(1);
  EXPECT_CALL(mockBackend_, getAttribute("key1")).WillOnce(Return("value1"));
  EXPECT_CALL(mockBackend_, hasAttribute("key1")).WillOnce(Return(true));

  observation_.setAttribute("key1", "value1");
  ASSERT_EQ(observation_.getAttribute("key1"), "value1");
  ASSERT_TRUE(observation_.hasAttribute("key1"));
}

TEST_F(ObservationTest, ThrowsOnInvalidAccess) {
  EXPECT_CALL(mockBackend_, isValid()).WillOnce(Return(false));
  EXPECT_THROW(observation_.getData<double>(), std::runtime_error);
}

}  // namespace metada::tests