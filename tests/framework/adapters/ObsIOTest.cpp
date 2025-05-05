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
 * - Format information retrieval
 * - Capability checks
 *
 * The test suite uses Google Test/Mock framework for mocking and assertions.
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

#include "DateTime.hpp"
#include "MockBackendTraits.hpp"
#include "ObsIO.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::ObservationRecord;
using framework::ObsIO;

/**
 * @brief Test fixture for ObsIO tests
 *
 * This fixture sets up the necessary objects for testing the ObsIO adapter,
 * including initialization parameters and sample observation records.
 */
class ObsIOTest : public ::testing::Test {
 protected:
  // Initialization parameters for ObsIO
  std::string params_;

  // Sample observation records for testing
  std::vector<ObservationRecord> records_;

  /**
   * @brief Set up the test environment
   *
   * Creates sample data needed for testing.
   */
  void SetUp() override {
    // Set up initialization parameters
    params_ = "config=path/to/config.yaml";

    // Create sample observation records
    ObservationRecord record1;
    record1.type = "temperature";
    record1.value = 25.5;
    record1.location = "STATION_001";
    record1.datetime = DateTime();
    record1.qc_marker = 0;

    ObservationRecord record2;
    record2.type = "pressure";
    record2.value = 1013.2;
    record2.location = "STATION_001";
    record2.datetime = DateTime();
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
 * @brief Test construction with parameters
 *
 * Verifies that an ObsIO object can be properly constructed with initialization
 * parameters.
 */
TEST_F(ObsIOTest, Construction) {
  // Construct the ObsIO object
  ObsIO<traits::MockBackendTag> obsIO(params_);

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
  ObsIO<traits::MockBackendTag> obsIO1(params_);

  // Test move constructor
  ObsIO<traits::MockBackendTag> obsIO2(std::move(obsIO1));

  // Test move assignment
  ObsIO<traits::MockBackendTag> obsIO3(params_);
  obsIO3 = std::move(obsIO2);

  // Verify the final object's backend is accessible
  auto& backend = obsIO3.backend();
  EXPECT_NE(&backend, nullptr);
}

/**
 * @brief Test canRead method
 *
 * Verifies that the canRead method correctly delegates to the backend.
 */
TEST_F(ObsIOTest, CanRead) {
  // Create an ObsIO object
  ObsIO<traits::MockBackendTag> obsIO(params_);
  auto& backend = obsIO.backend();

  // Setup expectation
  const std::string filename = "observations.bufr";
  EXPECT_CALL(backend, canRead(filename)).WillOnce(Return(true));

  // Call the method and verify
  EXPECT_TRUE(obsIO.canRead(filename));
}

/**
 * @brief Test canWrite method
 *
 * Verifies that the canWrite method correctly delegates to the backend.
 */
TEST_F(ObsIOTest, CanWrite) {
  // Create an ObsIO object
  ObsIO<traits::MockBackendTag> obsIO(params_);
  auto& backend = obsIO.backend();

  // Setup expectation
  EXPECT_CALL(backend, canWrite()).WillOnce(Return(true));

  // Call the method and verify
  EXPECT_TRUE(obsIO.canWrite());
}

/**
 * @brief Test getFormatName method
 *
 * Verifies that the getFormatName method correctly delegates to the backend.
 */
TEST_F(ObsIOTest, GetFormatName) {
  // Create an ObsIO object
  ObsIO<traits::MockBackendTag> obsIO(params_);
  auto& backend = obsIO.backend();

  // Setup expectation
  const std::string formatName = "BUFR";
  EXPECT_CALL(backend, getFormatName()).WillOnce(Return(formatName));

  // Call the method and verify
  EXPECT_EQ(obsIO.getFormatName(), formatName);
}

/**
 * @brief Test getFileExtensions method
 *
 * Verifies that the getFileExtensions method correctly delegates to the
 * backend.
 */
TEST_F(ObsIOTest, GetFileExtensions) {
  // Create an ObsIO object
  ObsIO<traits::MockBackendTag> obsIO(params_);
  auto& backend = obsIO.backend();

  // Setup expectation
  const std::vector<std::string> extensions = {".bufr", ".BUFR", ".bfr"};
  EXPECT_CALL(backend, getFileExtensions()).WillOnce(Return(extensions));

  // Call the method and verify
  EXPECT_EQ(obsIO.getFileExtensions(), extensions);
}

/**
 * @brief Test readObservations method
 *
 * Verifies that the readObservations method correctly delegates to the backend
 * and handles errors appropriately.
 */
TEST_F(ObsIOTest, ReadObservations) {
  // Create an ObsIO object
  ObsIO<traits::MockBackendTag> obsIO(params_);
  auto& backend = obsIO.backend();

  // Setup expectations
  const std::string filename = "observations.bufr";
  EXPECT_CALL(backend, canRead(filename)).WillOnce(Return(true));
  EXPECT_CALL(backend, read(filename)).WillOnce(Return(records_));

  // Call the method and verify
  std::vector<ObservationRecord> result = obsIO.read(filename);
  EXPECT_EQ(result.size(), records_.size());
  EXPECT_EQ(result[0].type, records_[0].type);
  EXPECT_EQ(result[0].value, records_[0].value);
  EXPECT_EQ(result[1].type, records_[1].type);
  EXPECT_EQ(result[1].value, records_[1].value);
}

/**
 * @brief Test readObservations method with unsupported file
 *
 * Verifies that the readObservations method throws an exception when trying to
 * read an unsupported file format.
 */
TEST_F(ObsIOTest, ReadObservationsUnsupportedFormat) {
  // Create an ObsIO object
  ObsIO<traits::MockBackendTag> obsIO(params_);
  auto& backend = obsIO.backend();

  // Setup expectation
  const std::string filename = "observations.unsupported";
  EXPECT_CALL(backend, canRead(filename)).WillOnce(Return(false));

  // Call the method and verify exception
  EXPECT_THROW(obsIO.read(filename), std::runtime_error);
}

/**
 * @brief Test writeObservations method
 *
 * Verifies that the writeObservations method correctly delegates to the backend
 * and handles errors appropriately.
 */
TEST_F(ObsIOTest, WriteObservations) {
  // Create an ObsIO object
  ObsIO<traits::MockBackendTag> obsIO(params_);
  auto& backend = obsIO.backend();

  // Setup expectations
  const std::string filename = "output.bufr";
  EXPECT_CALL(backend, canWrite()).WillOnce(Return(true));
  EXPECT_CALL(backend, write(filename, _)).Times(1);

  // Call the method
  obsIO.write(filename, records_);
}

/**
 * @brief Test writeObservations method with unsupported backend
 *
 * Verifies that the writeObservations method throws an exception when trying to
 * write with a backend that doesn't support writing.
 */
TEST_F(ObsIOTest, WriteObservationsUnsupported) {
  // Create an ObsIO object
  ObsIO<traits::MockBackendTag> obsIO(params_);
  auto& backend = obsIO.backend();

  // Setup expectation
  EXPECT_CALL(backend, canWrite()).WillOnce(Return(false));

  // Call the method and verify exception
  EXPECT_THROW(obsIO.write("output.bufr", records_), std::runtime_error);
}

}  // namespace metada::tests