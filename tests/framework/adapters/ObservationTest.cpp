#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "MockConfig.hpp"
#include "MockObservation.hpp"
#include "Observation.hpp"
#include "utils/config/Config.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Const;
using ::testing::Return;
using ::testing::ReturnRef;

using metada::framework::Config;
using metada::framework::Observation;

class ObservationTest : public ::testing::Test {
 protected:
  // Test data
  std::vector<std::string> variableNames_{"temperature", "humidity"};
  std::vector<size_t> dimensions_{10, 5};
  std::vector<std::vector<double>> locations_{{45.0, -120.0, 100.0},
                                              {46.0, -121.0, 200.0}};
  std::vector<double> times_{1234567890.0, 1234567891.0};
  std::vector<int> qualityFlags_{1, 2};
  std::vector<double> confidenceValues_{0.95, 0.85};
  double temperature_{25.5};
  double uncertainty_{0.5};

  // Mock config and observation
  Config<MockConfig> mockConfig_;
  MockObservation mockBackend_{mockConfig_};
  Observation<MockObservation> observation_{mockBackend_};

  // Setup mock expectations for common operations
  void SetUp() override {
    // Setup mock to return our test data
    EXPECT_CALL(mockBackend_, isValid()).WillRepeatedly(Return(true));
    EXPECT_CALL(mockBackend_, getVariableNames())
        .WillRepeatedly(ReturnRef(variableNames_));
    EXPECT_CALL(mockBackend_, getDimensions())
        .WillRepeatedly(ReturnRef(dimensions_));
    EXPECT_CALL(mockBackend_, getLocations())
        .WillRepeatedly(ReturnRef(locations_));
    EXPECT_CALL(mockBackend_, getTimes()).WillRepeatedly(ReturnRef(times_));
    EXPECT_CALL(mockBackend_, getQualityFlags())
        .WillRepeatedly(ReturnRef(qualityFlags_));
    EXPECT_CALL(mockBackend_, getConfidenceValues())
        .WillRepeatedly(ReturnRef(confidenceValues_));

    // Setup data access
    EXPECT_CALL(mockBackend_, getData())
        .WillRepeatedly(Return(static_cast<void*>(&temperature_)));
    ON_CALL(testing::Const(mockBackend_), getData())
        .WillByDefault(Return(static_cast<const void*>(&temperature_)));
    EXPECT_CALL(mockBackend_, getUncertainty())
        .WillRepeatedly(Return(static_cast<void*>(&uncertainty_)));
    ON_CALL(testing::Const(mockBackend_), getUncertainty())
        .WillByDefault(Return(static_cast<const void*>(&uncertainty_)));
    EXPECT_CALL(mockBackend_, getSize()).WillRepeatedly(Return(2));
  }
};

TEST_F(ObservationTest, InitializeCallsBackend) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize()).Times(1);

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Call initialize
  testObs.initialize();
}

TEST_F(ObservationTest, InitializeWithConfigCallsBackend) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize(testing::_)).Times(1);

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Call initialize with config
  testObs.initialize(testConfig.backend());
}

TEST_F(ObservationTest, GetDataReturnsCorrectValue) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize()).Times(1);  // Called in constructor
  EXPECT_CALL(testMock, isValid()).WillRepeatedly(Return(true));
  EXPECT_CALL(testMock, getData())
      .WillRepeatedly(Return(static_cast<void*>(&temperature_)));

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Test the method
  ASSERT_EQ(testObs.getData<double>(), temperature_);
}

TEST_F(ObservationTest, GetUncertaintyReturnsCorrectValue) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize()).Times(1);  // Called in constructor
  EXPECT_CALL(testMock, getUncertainty())
      .WillRepeatedly(Return(static_cast<void*>(&uncertainty_)));

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Test the method
  ASSERT_EQ(testObs.getUncertainty<double>(), uncertainty_);
}

TEST_F(ObservationTest, LocationsAreCorrectlyReturned) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize()).Times(1);  // Called in constructor
  EXPECT_CALL(testMock, getLocations()).WillRepeatedly(ReturnRef(locations_));

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Test the method
  const auto& locs = testObs.getLocations();
  ASSERT_EQ(locs.size(), 2);
  ASSERT_EQ(locs[0][0], 45.0);
  ASSERT_EQ(locs[0][1], -120.0);
  ASSERT_EQ(locs[0][2], 100.0);
  ASSERT_EQ(locs[1][0], 46.0);
  ASSERT_EQ(locs[1][1], -121.0);
  ASSERT_EQ(locs[1][2], 200.0);
}

TEST_F(ObservationTest, SetLocationsCallsBackend) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize()).Times(1);  // Called in constructor
  EXPECT_CALL(testMock, setLocations(_)).Times(1);

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Test the method
  std::vector<std::vector<double>> newLocations{{50.0, -110.0, 200.0},
                                                {51.0, -111.0, 300.0}};
  testObs.setLocations(newLocations);
}

TEST_F(ObservationTest, FluentInterfaceWorks) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize()).Times(1);  // Called in constructor
  EXPECT_CALL(testMock, setLocations(_)).Times(1);
  EXPECT_CALL(testMock, setTimes(_)).Times(1);
  EXPECT_CALL(testMock, setQualityFlags(_)).Times(1);

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Test chaining of methods
  std::vector<std::vector<double>> newLocations{{50.0, -110.0, 200.0},
                                                {51.0, -111.0, 300.0}};
  std::vector<double> newTimes{1234567892.0, 1234567893.0};
  std::vector<int> newFlags{3, 4};

  testObs.setLocations(newLocations)
      .setTimes(newTimes)
      .setQualityFlags(newFlags);
}

TEST_F(ObservationTest, MetadataOperationsWork) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize()).Times(1);  // Called in constructor

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Test the methods
  testObs.setMetadata("key1", "value1");
  EXPECT_CALL(testMock, getMetadata("key1")).WillOnce(Return("value1"));
  EXPECT_CALL(testMock, hasMetadata("key1")).WillOnce(Return(true));

  ASSERT_EQ(testObs.getMetadata("key1"), "value1");
  ASSERT_TRUE(testObs.hasMetadata("key1"));
}

TEST_F(ObservationTest, ThrowsOnInvalidAccess) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize()).Times(1);  // Called in constructor
  EXPECT_CALL(testMock, isValid()).WillOnce(Return(false));

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Test the method
  EXPECT_THROW(testObs.getData<double>(), std::runtime_error);
}

TEST_F(ObservationTest, CopyConstructionWorks) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize()).Times(1);  // Called in constructor
  EXPECT_CALL(testMock, isValid()).WillRepeatedly(Return(true));
  EXPECT_CALL(testMock, getData())
      .WillRepeatedly(Return(static_cast<void*>(&temperature_)));

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Create a copied observation
  EXPECT_CALL(testMock, copyFrom(_)).Times(1);
  Observation<MockObservation> copiedObs(testObs);

  // Verify the copied observation works correctly
  EXPECT_CALL(testMock, getData())
      .WillRepeatedly(Return(static_cast<void*>(&temperature_)));
  ASSERT_EQ(copiedObs.getData<double>(), temperature_);
}

TEST_F(ObservationTest, MoveConstructionWorks) {
  // Create a new mock for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock(testConfig);
  EXPECT_CALL(testMock, initialize()).Times(1);  // Called in constructor
  EXPECT_CALL(testMock, isValid()).WillRepeatedly(Return(true));
  EXPECT_CALL(testMock, getData())
      .WillRepeatedly(Return(static_cast<void*>(&temperature_)));

  // Create a new observation with this mock
  Observation<MockObservation> testObs(testMock);

  // Create a moved observation
  Observation<MockObservation> movedObs(std::move(testObs));

  // Verify the moved observation works correctly
  ASSERT_EQ(movedObs.getData<double>(), temperature_);
}

TEST_F(ObservationTest, ArithmeticOperationsWork) {
  // Create mocks for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock1(testConfig);
  MockObservation testMock2(testConfig);

  // Setup expectations
  EXPECT_CALL(testMock1, initialize()).Times(1);
  EXPECT_CALL(testMock2, initialize()).Times(1);
  EXPECT_CALL(testMock1, add(_)).Times(1);
  EXPECT_CALL(testMock1, subtract(_)).Times(1);
  EXPECT_CALL(testMock1, multiply(2.0)).Times(1);

  // Create observations
  Observation<MockObservation> obs1(testMock1);
  Observation<MockObservation> obs2(testMock2);

  // Test arithmetic operations
  obs1.add(obs2);
  obs1.subtract(obs2);
  obs1.multiply(2.0);
}

TEST_F(ObservationTest, OperatorOverloadsWork) {
  // Create mocks for this specific test
  Config<MockConfig> testConfig;
  MockObservation testMock1(testConfig);
  MockObservation testMock2(testConfig);

  // Setup expectations
  EXPECT_CALL(testMock1, initialize()).Times(1);
  EXPECT_CALL(testMock2, initialize()).Times(1);

  // For operator+
  EXPECT_CALL(testMock1, copyFrom(_)).Times(3);  // Once for each operator test
  EXPECT_CALL(testMock1, add(_)).Times(1);
  EXPECT_CALL(testMock1, subtract(_)).Times(1);
  EXPECT_CALL(testMock1, multiply(2.0)).Times(1);

  // Create observations
  Observation<MockObservation> obs1(testMock1);
  Observation<MockObservation> obs2(testMock2);

  // Test operator overloads
  Observation<MockObservation> result1 = obs1 + obs2;
  Observation<MockObservation> result2 = obs1 - obs2;
  Observation<MockObservation> result3 = obs1 * 2.0;
}

}  // namespace metada::tests