/**
 * @file GeometryTest.cpp
 * @brief Unit tests for the Geometry adapter and IGeometry interface
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This test suite verifies the proper functioning of the Geometry adapter
 * and the implementation of the IGeometry interface by MockGeometry.
 * The tests ensure that the adapter correctly delegates operations to the
 * backend implementation.
 *
 * Core functionality:
 * - Initialization and construction
 * - Dimension and resolution management
 * - Domain bounds and grid spacing
 * - Coordinate transformations
 *
 * Advanced features:
 * - Grid iteration
 * - Compatibility checking
 * - Error handling
 * - Copy/move semantics
 * - Configuration loading
 *
 * The test suite uses Google Test/Mock framework for mocking and assertions.
 *
 * @see IGeometry
 * @see Geometry
 * @see MockGeometry
 */

// Standard C++ includes
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Framework includes
#include "Geometry.hpp"
#include "GeometryIterator.hpp"
#include "GeometryPointIterator.hpp"
#include "IGeometry.hpp"
#include "IGeometryIterator.hpp"
#include "MockGeometry.hpp"

// Legacy includes
#include "AppTraits.hpp"
#include "ApplicationContext.hpp"
#include "MockConfig.hpp"
#include "MockLogger.hpp"
#include "MockModel.hpp"
#include "MockObsOperator.hpp"
#include "MockObservation.hpp"
#include "MockState.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::DoAll;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::SetArgReferee;

using framework::Config;
using framework::Geometry;
using framework::GeometryIterator;
using framework::IGeometry;
using framework::Logger;
using framework::runs::ApplicationContext;

// Define test traits
using Traits = AppTraits<MockLogger, MockConfig, MockGeometry, MockState,
                         MockModel, MockObservation, MockObsOperator>;

/**
 * @brief Test fixture for Geometry adapter tests
 */
class GeometryTest : public ::testing::Test {
 protected:
  /**
   * @brief Application context for testing
   */
  std::unique_ptr<ApplicationContext<Traits>> app_context_;

  /**
   * @brief Test dimensions and resolutions
   */
  const size_t default_dimensions_ = 3;
  const std::vector<size_t> default_resolution_ = {10, 20, 30};
  const std::vector<double> default_min_bounds_ = {-1.0, -2.0, -3.0};
  const std::vector<double> default_max_bounds_ = {1.0, 2.0, 3.0};

  /**
   * @brief Test backends and adapters
   */
  std::unique_ptr<Traits::GeometryType> mock_backend1_;
  std::unique_ptr<Traits::GeometryType> mock_backend2_;
  std::unique_ptr<Geometry<Traits::GeometryType>> geometry1_;
  std::unique_ptr<Geometry<Traits::GeometryType>> geometry2_;
  std::unique_ptr<Traits::GeometryType> direct_mock_geometry_;

  /**
   * @brief Initialize test setup
   */
  void SetUp() override {
    // Create application context
    app_context_ = std::make_unique<ApplicationContext<Traits>>("GeometryTest");

    // Configure mock objects
    auto& config = app_context_->getConfig();

    // Set mock configuration for geometry
    ON_CALL(config.backend(), Get("geometry.dimensions"))
        .WillByDefault(Return(static_cast<int>(default_dimensions_)));

    // Resolution for each dimension
    for (size_t i = 0; i < default_dimensions_; ++i) {
      std::string key = "geometry.resolution." + std::to_string(i);
      ON_CALL(config.backend(), Get(key))
          .WillByDefault(Return(static_cast<int>(default_resolution_[i])));

      // Min and max bounds
      std::string min_key = "geometry.min_bounds." + std::to_string(i);
      std::string max_key = "geometry.max_bounds." + std::to_string(i);

      ON_CALL(config.backend(), Get(min_key))
          .WillByDefault(Return(static_cast<double>(default_min_bounds_[i])));

      ON_CALL(config.backend(), Get(max_key))
          .WillByDefault(Return(static_cast<double>(default_max_bounds_[i])));
    }

    // Create mock backends for the Geometry adapter
    mock_backend1_ = std::make_unique<Traits::GeometryType>();
    mock_backend1_->setTestDimensions(default_dimensions_);
    mock_backend1_->setTestResolution(default_resolution_);
    mock_backend1_->setTestBounds(default_min_bounds_, default_max_bounds_);

    // Set up expectations for the first mock backend
    EXPECT_CALL(*mock_backend1_, getDimensions())
        .WillRepeatedly(Return(default_dimensions_));
    EXPECT_CALL(*mock_backend1_, getResolution())
        .WillRepeatedly(Return(default_resolution_));
    EXPECT_CALL(*mock_backend1_, getTotalPoints())
        .WillRepeatedly(Return(default_resolution_[0] * default_resolution_[1] *
                               default_resolution_[2]));

    for (size_t i = 0; i < default_dimensions_; ++i) {
      mock_backend1_->setTestDomainBounds(i, default_min_bounds_[i],
                                          default_max_bounds_[i]);
      EXPECT_CALL(*mock_backend1_, getDomainBounds(i))
          .WillRepeatedly(Return(std::array<double, 2>{
              default_min_bounds_[i], default_max_bounds_[i]}));
    }

    // Create another mock backend for a simple 1D geometry
    mock_backend2_ = std::make_unique<Traits::GeometryType>();
    mock_backend2_->setTestDimensions(1);
    mock_backend2_->setTestResolution({10});
    mock_backend2_->setTestBounds({0.0}, {1.0});

    EXPECT_CALL(*mock_backend2_, getDimensions()).WillRepeatedly(Return(1));
    EXPECT_CALL(*mock_backend2_, getResolution())
        .WillRepeatedly(Return(std::vector<size_t>{10}));
    EXPECT_CALL(*mock_backend2_, getTotalPoints()).WillRepeatedly(Return(10));

    // Create the geometry objects
    geometry1_ =
        std::make_unique<Geometry<Traits::GeometryType>>(mock_backend1_.get());
    geometry2_ =
        std::make_unique<Geometry<Traits::GeometryType>>(mock_backend2_.get());

    // Create a separate mock geometry for direct interface testing
    direct_mock_geometry_ = std::make_unique<Traits::GeometryType>();
    direct_mock_geometry_->setTestDimensions(default_dimensions_);
    direct_mock_geometry_->setTestResolution(default_resolution_);
    direct_mock_geometry_->setTestBounds(default_min_bounds_,
                                         default_max_bounds_);
  }

  /**
   * @brief Clean up test resources
   */
  void TearDown() override {
    direct_mock_geometry_.reset();
    geometry1_.reset();
    geometry2_.reset();
    mock_backend1_.reset();
    mock_backend2_.reset();
    app_context_.reset();
  }

  /**
   * @brief Helper to get config reference
   */
  Config<MockConfig>& getConfig() { return app_context_->getConfig(); }
};

/**
 * @brief Test Geometry adapter for configuration and basic properties
 */
TEST_F(GeometryTest, BasicProperties) {
  // Set up expectations for direct mock
  EXPECT_CALL(*direct_mock_geometry_, getDimensions())
      .WillRepeatedly(Return(default_dimensions_));
  EXPECT_CALL(*direct_mock_geometry_, getResolution())
      .WillRepeatedly(Return(default_resolution_));
  EXPECT_CALL(*direct_mock_geometry_, getTotalPoints())
      .WillRepeatedly(Return(default_resolution_[0] * default_resolution_[1] *
                             default_resolution_[2]));

  // Test adapter method: getDimensions
  EXPECT_EQ(geometry1_->getDimensions(), default_dimensions_);
  EXPECT_EQ(direct_mock_geometry_->getDimensions(), default_dimensions_);

  // Test adapter method: getResolution
  EXPECT_EQ(geometry1_->getResolution(), default_resolution_);
  EXPECT_EQ(direct_mock_geometry_->getResolution(), default_resolution_);

  // Test adapter method: getTotalPoints
  size_t expected_total =
      default_resolution_[0] * default_resolution_[1] * default_resolution_[2];
  EXPECT_EQ(geometry1_->getTotalPoints(), expected_total);
}

/**
 * @brief Test the Geometry adapter properly forwards calls to the backend
 */
TEST_F(GeometryTest, AdapterDelegation) {
  // Create a mock geometry backend for testing delegation
  auto mockBackend = std::make_unique<Traits::GeometryType>();

  // Set up expectations for methods that should be called on the backend
  EXPECT_CALL(*mockBackend, getDimensions()).WillOnce(Return(3));

  std::vector<size_t> testResolution{10, 10, 10};
  EXPECT_CALL(*mockBackend, getResolution()).WillOnce(Return(testResolution));

  EXPECT_CALL(*mockBackend, getTotalPoints()).WillOnce(Return(1000));

  std::vector<double> testSpacing{0.1, 0.1, 0.1};
  EXPECT_CALL(*mockBackend, getGridSpacing()).WillOnce(Return(testSpacing));

  // Create a Geometry adapter with the mock backend
  Geometry<Traits::GeometryType> geometry(mockBackend.get());

  // Verify adapter delegates to the backend
  EXPECT_EQ(geometry.getDimensions(), 3);

  // Use EXPECT_TRUE with container equality for vector comparisons
  auto resolution = geometry.getResolution();
  EXPECT_TRUE(resolution == testResolution);

  EXPECT_EQ(geometry.getTotalPoints(), 1000);

  // Use EXPECT_TRUE with container equality for vector comparisons
  auto spacing = geometry.getGridSpacing();
  EXPECT_TRUE(spacing == testSpacing);
}

}  // namespace metada::tests