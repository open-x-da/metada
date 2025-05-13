/**
 * @file ObsIOTest.cpp
 * @brief Unit tests for the ObsIO adapter class
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This test suite verifies the functionality of the ObsIO adapter class,
 * which provides a unified interface for reading and writing observation data
 * in different file formats. The tests cover:
 *
 * - Construction with parameters
 * - Move semantics
 * - Reading observations
 * - Writing observations
 *
 * The test suite uses Google Test/Mock framework for mocking and assertions.
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

#include "Config.hpp"
#include "DateTime.hpp"
#include "MockBackendTraits.hpp"
#include "ObsIO.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Config;
using framework::ObsIO;
using framework::ObsRecord;

/**
 * @brief Test fixture for ObsIO tests
 *
 * This fixture sets up the necessary objects for testing the ObsIO adapter,
 * including initialization parameters and sample observation records.
 */
class ObsIOTest : public ::testing::Test {
 protected:
  std::string config_file_;

  // Sample observation records for testing
  std::vector<ObsRecord> records_;

  /**
   * @brief Set up the test environment
   *
   * Creates sample data needed for testing.
   */
  void SetUp() override {
    // Get the directory where the test file is located
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    config_file_ = (test_dir / "test_config.yaml").string();

    // Create sample observation records
    ObsRecord record1;
    record1.type = "temperature";
    record1.value = 25.5;
    record1.station_id = "STATION_001";
    record1.longitude = -75.0;
    record1.latitude = 40.0;
    record1.elevation = 100.0;
    record1.datetime = DateTime();
    record1.report_type = "1";
    record1.input_report_type = "0";
    record1.instrument_type = "0";
    record1.qc_marker = 0;

    ObsRecord record2;
    record2.type = "pressure";
    record2.value = 1013.2;
    record2.station_id = "STATION_001";
    record2.longitude = -75.0;
    record2.latitude = 40.0;
    record2.elevation = 100.0;
    record2.datetime = DateTime();
    record2.report_type = "1";
    record2.input_report_type = "0";
    record2.instrument_type = "0";
    record2.qc_marker = 0;

    records_.push_back(record1);
    records_.push_back(record2);
  }

  /**
   * @brief Clean up the test environment
   */
  void TearDown() override { records_.clear(); }
};

/**
 * @brief Test that ObsIO correctly follows NonCopyable pattern
 *
 * Verifies that ObsIO adheres to the NonCopyable design pattern
 * by ensuring it cannot be copied but can be moved. This is important
 * for resource management and ownership semantics.
 */
TEST_F(ObsIOTest, NonCopyableButMovable) {
  using ObsIOType = ObsIO<traits::MockBackendTag>;

  // Static assertions to verify type traits
  static_assert(!std::is_default_constructible_v<ObsIOType>,
                "ObsIO should not be default constructible");
  static_assert(!std::is_copy_constructible_v<ObsIOType>,
                "ObsIO should not be copy constructible");
  static_assert(!std::is_copy_assignable_v<ObsIOType>,
                "ObsIO should not be copy assignable");
  static_assert(std::is_move_constructible_v<ObsIOType>,
                "ObsIO should be move constructible");
  static_assert(std::is_move_assignable_v<ObsIOType>,
                "ObsIO should be move assignable");
}

/**
 * @brief Test construction with parameters
 *
 * Verifies that an ObsIO object can be properly constructed with initialization
 * parameters.
 */
TEST_F(ObsIOTest, ConstructionFromConfig) {
  // Construct the ObsIO object
  auto config = Config<traits::MockBackendTag>(config_file_);
  ObsIO<traits::MockBackendTag> obsIO(std::move(config));

  // Verify the backend is accessible
  auto& backend = obsIO.backend();
  EXPECT_NE(&backend, nullptr);
}

/**
 * @brief Test move semantics
 *
 * Verifies that ObsIO correctly implements move construction and move
 * assignment, transferring ownership of resources and ensuring the moved-from
 * object is valid.
 */
TEST_F(ObsIOTest, MoveSemantics) {
  // Create an ObsIO object
  auto config = Config<traits::MockBackendTag>(config_file_);
  ObsIO<traits::MockBackendTag> obsIO1(std::move(config));

  // Test move constructor
  auto config2 = Config<traits::MockBackendTag>(config_file_);
  ObsIO<traits::MockBackendTag> obsIO2(std::move(config2));

  // Test move assignment
  auto config3 = Config<traits::MockBackendTag>(config_file_);
  ObsIO<traits::MockBackendTag> obsIO3(std::move(config3));
  obsIO3 = std::move(obsIO2);

  // Verify the final object's backend is accessible
  auto& backend = obsIO3.backend();
  EXPECT_NE(&backend, nullptr);
}

/**
 * @brief Test read method
 *
 * Verifies that the read method correctly delegates to the backend
 * and handles errors appropriately.
 */
TEST_F(ObsIOTest, ReadObservations) {
  // Create an ObsIO object
  auto config = Config<traits::MockBackendTag>(config_file_);
  ObsIO<traits::MockBackendTag> obsIO(std::move(config));
  auto& backend = obsIO.backend();

  // Setup expectation
  EXPECT_CALL(backend, read()).WillOnce(Return(records_));

  // Call the method and verify
  std::vector<ObsRecord> result = obsIO.read();
  EXPECT_EQ(result.size(), records_.size());
  EXPECT_EQ(result[0].type, records_[0].type);
  EXPECT_EQ(result[0].value, records_[0].value);
  EXPECT_EQ(result[1].type, records_[1].type);
  EXPECT_EQ(result[1].value, records_[1].value);
}

/**
 * @brief Test read method with backend errors
 *
 * Verifies that the read method properly propagates exceptions from the
 * backend.
 */
TEST_F(ObsIOTest, ReadObservationsError) {
  // Create an ObsIO object
  auto config = Config<traits::MockBackendTag>(config_file_);
  ObsIO<traits::MockBackendTag> obsIO(std::move(config));
  auto& backend = obsIO.backend();

  // Setup expectation to throw an exception
  EXPECT_CALL(backend, read())
      .WillOnce(::testing::Throw(std::runtime_error("Test error")));

  // Call the method and verify exception is propagated
  EXPECT_THROW(obsIO.read(), std::runtime_error);
}

/**
 * @brief Test write method
 *
 * Verifies that the write method correctly delegates to the backend
 * and handles errors appropriately.
 */
TEST_F(ObsIOTest, WriteObservations) {
  // Create an ObsIO object
  auto config = Config<traits::MockBackendTag>(config_file_);
  ObsIO<traits::MockBackendTag> obsIO(std::move(config));
  auto& backend = obsIO.backend();

  // Setup expectation
  EXPECT_CALL(backend, write(_)).Times(1);

  // Call the method
  obsIO.write(records_);
}

/**
 * @brief Test write method with backend errors
 *
 * Verifies that the write method properly propagates exceptions from the
 * backend.
 */
TEST_F(ObsIOTest, WriteObservationsError) {
  // Create an ObsIO object
  auto config = Config<traits::MockBackendTag>(config_file_);
  ObsIO<traits::MockBackendTag> obsIO(std::move(config));
  auto& backend = obsIO.backend();

  // Setup expectation to throw an exception
  EXPECT_CALL(backend, write(_))
      .WillOnce(::testing::Throw(std::runtime_error("Test error")));

  // Call the method and verify exception is propagated
  EXPECT_THROW(obsIO.write(records_), std::runtime_error);
}

}  // namespace metada::tests