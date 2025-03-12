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

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "AppTraits.hpp"
#include "ApplicationContext.hpp"
#include "Geometry.hpp"
#include "IGeometry.hpp"
#include "MockConfig.hpp"
#include "MockGeometry.hpp"
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
    mock_backend1_ = std::make_unique<MockGeometry>();
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

    // Create geometry adapters with mock backends
    /*
    geometry1_ =
        std::make_unique<Geometry<Traits::GeometryType>>(*mock_backend1_);
    geometry2_ =
        std::make_unique<Geometry<Traits::GeometryType>>(*mock_backend2_);

    // Create a separate mock geometry for direct interface testing
    direct_mock_geometry_ = std::make_unique<Traits::GeometryType>();
    direct_mock_geometry_->setTestDimensions(default_dimensions_);
    direct_mock_geometry_->setTestResolution(default_resolution_);
    direct_mock_geometry_->setTestBounds(default_min_bounds_,
                                         default_max_bounds_);
    */
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

  /**
   * @brief Helper to create a geometry with default test values
   * @return A Geometry adapter with MockGeometry backend
   */
  /*
  std::unique_ptr<Geometry<Traits::GeometryType>> createGeometry() {
    auto mock_backend = std::make_unique<Traits::GeometryType>();
    mock_backend->setTestDimensions(default_dimensions_);
    mock_backend->setTestResolution(default_resolution_);
    mock_backend->setTestBounds(default_min_bounds_, default_max_bounds_);

    // Set up minimal expectations
    EXPECT_CALL(*mock_backend, getDimensions())
        .WillRepeatedly(Return(default_dimensions_));
    EXPECT_CALL(*mock_backend, getResolution())
        .WillRepeatedly(Return(default_resolution_));

    for (size_t i = 0; i < default_dimensions_; ++i) {
      EXPECT_CALL(*mock_backend, getDomainBounds(i))
          .WillRepeatedly(Return(std::array<double, 2>{
              default_min_bounds_[i], default_max_bounds_[i]}));
    }

    return std::make_unique<Geometry<Traits::GeometryType>>(*mock_backend);
  }*/
};

/**
 * @brief Test Geometry adapter for configuration and basic properties
 */
/*
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
*/

/**
 * @brief Test Geometry adapter for domain bounds
 */
/*
TEST_F(GeometryTest, DomainBounds) {
  // Set up mock expectations for domain bounds on direct mock
  for (size_t i = 0; i < default_dimensions_; ++i) {
    direct_mock_geometry_->setTestDomainBounds(i, default_min_bounds_[i],
                                               default_max_bounds_[i]);
  }

  // Test getDomainBounds adapter method
  for (size_t i = 0; i < default_dimensions_; ++i) {
    auto bounds = geometry1_->getDomainBounds(i);
    EXPECT_DOUBLE_EQ(bounds[0], default_min_bounds_[i]);
    EXPECT_DOUBLE_EQ(bounds[1], default_max_bounds_[i]);

    auto mock_bounds = direct_mock_geometry_->getDomainBounds(i);
    EXPECT_DOUBLE_EQ(mock_bounds[0], default_min_bounds_[i]);
    EXPECT_DOUBLE_EQ(mock_bounds[1], default_max_bounds_[i]);
  }

  // Test adapter method: modifying domain bounds
  double new_min = -5.0;
  double new_max = 5.0;

  // Set up expectation for direct mock to handle setDomainBounds
  EXPECT_CALL(*direct_mock_geometry_, setDomainBounds(0, new_min, new_max))
      .Times(1);
  direct_mock_geometry_->setDomainBounds(0, new_min, new_max);

  // Set up expectation for the backend mock that will be used via the adapter
  EXPECT_CALL(*mock_backend1_, setDomainBounds(0, new_min, new_max)).Times(1);

  // Test the adapter's delegation to the backend
  geometry1_->setDomainBounds(0, new_min, new_max);

  // Set up the expected return value for the updated bounds
  EXPECT_CALL(*mock_backend1_, getDomainBounds(0))
      .WillOnce(Return(std::array<double, 2>{new_min, new_max}));

  auto updated_bounds = geometry1_->getDomainBounds(0);
  EXPECT_DOUBLE_EQ(updated_bounds[0], new_min);
  EXPECT_DOUBLE_EQ(updated_bounds[1], new_max);
}
*/

/**
 * @brief Test Geometry adapter for coordinate transformations
 */
/*
TEST_F(GeometryTest, CoordinateTransformations) {
  // Create a 2D mock backend for simpler testing
  auto mock_backend = std::make_unique<Traits::GeometryType>();
  mock_backend->setTestDimensions(2);
  mock_backend->setTestResolution({5, 5});
  mock_backend->setTestBounds({0.0, 0.0}, {1.0, 1.0});
  mock_backend->setTestGridSpacing({0.25, 0.25});

  // Set up expectations for dimension and resolution
  EXPECT_CALL(*mock_backend, getDimensions()).WillRepeatedly(Return(2));
  EXPECT_CALL(*mock_backend, getResolution())
      .WillRepeatedly(Return(std::vector<size_t>{5, 5}));

  // Wrap the mock in the Geometry adapter
  auto geometry = Geometry<Traits::GeometryType>(*mock_backend);

  // Set up direct mock geometry with the same specs
  direct_mock_geometry_->setTestDimensions(2);
  direct_mock_geometry_->setTestResolution({5, 5});
  direct_mock_geometry_->setTestBounds({0.0, 0.0}, {1.0, 1.0});
  direct_mock_geometry_->setTestGridSpacing({0.25, 0.25});

  // Test adapter for coordinate transformations
  std::vector<size_t> indices = {2, 3};
  std::vector<double> expected_coords = {0.5, 0.75};

  // Setup mock expectations for direct mock
  EXPECT_CALL(*direct_mock_geometry_, getCoordinates(indices))
      .WillOnce(Return(expected_coords));

  // Setup expectation for the backend mock
  EXPECT_CALL(*mock_backend, getCoordinates(indices))
      .WillOnce(Return(expected_coords));

  // Test the adapter's delegation to the backend
  auto coords = geometry.getCoordinates(indices);
  EXPECT_DOUBLE_EQ(coords[0], expected_coords[0]);
  EXPECT_DOUBLE_EQ(coords[1], expected_coords[1]);

  // And verify the direct mock behaves the same way
  auto mock_coords = direct_mock_geometry_->getCoordinates(indices);
  EXPECT_DOUBLE_EQ(mock_coords[0], expected_coords[0]);
  EXPECT_DOUBLE_EQ(mock_coords[1], expected_coords[1]);

  // Test the inverse mapping (getIndices)
  EXPECT_CALL(*direct_mock_geometry_, getIndices(expected_coords))
      .WillOnce(Return(indices));

  EXPECT_CALL(*mock_backend, getIndices(expected_coords))
      .WillOnce(Return(indices));

  auto reverse_indices = geometry.getIndices(expected_coords);
  EXPECT_EQ(reverse_indices, indices);

  auto mock_indices = direct_mock_geometry_->getIndices(expected_coords);
  EXPECT_EQ(mock_indices, indices);
}
*/

/**
 * @brief Test Geometry adapter for point containment
 */
/*
TEST_F(GeometryTest, PointContainment) {
  // Create a 2D mock backend for testing
  auto mock_backend = std::make_unique<Traits::GeometryType>();
  mock_backend->setTestDimensions(2);
  mock_backend->setTestResolution({5, 5});
  mock_backend->setTestBounds({-1.0, -2.0}, {1.0, 2.0});

  // Set up expectations for dimension and resolution
  EXPECT_CALL(*mock_backend, getDimensions()).WillRepeatedly(Return(2));
  EXPECT_CALL(*mock_backend, getResolution())
      .WillRepeatedly(Return(std::vector<size_t>{5, 5}));

  // Wrap the mock in the Geometry adapter
  auto geometry = Geometry<Traits::GeometryType>(*mock_backend);

  // Set up direct mock with same specs
  direct_mock_geometry_->setTestDimensions(2);
  direct_mock_geometry_->setTestResolution({5, 5});
  direct_mock_geometry_->setTestBounds({-1.0, -2.0}, {1.0, 2.0});

  // Test points inside domain
  std::vector<double> point_inside = {0.0, 0.0};
  EXPECT_CALL(*direct_mock_geometry_, containsPoint(point_inside))
      .WillOnce(Return(true));

  EXPECT_CALL(*mock_backend, containsPoint(point_inside))
      .WillOnce(Return(true));

  EXPECT_TRUE(geometry.containsPoint(point_inside));
  EXPECT_TRUE(direct_mock_geometry_->containsPoint(point_inside));

  // Test points outside domain
  std::vector<double> point_outside = {2.0, 0.0};
  EXPECT_CALL(*direct_mock_geometry_, containsPoint(point_outside))
      .WillOnce(Return(false));

  EXPECT_CALL(*mock_backend, containsPoint(point_outside))
      .WillOnce(Return(false));

  EXPECT_FALSE(geometry.containsPoint(point_outside));
  EXPECT_FALSE(direct_mock_geometry_->containsPoint(point_outside));
}
*/

/**
 * @brief Test Geometry adapter for grid iteration
 */
/*
TEST_F(GeometryTest, GridIteration) {
  // Create a small 2D mock backend for simpler testing
  auto mock_backend = std::make_unique<Traits::GeometryType>();
  mock_backend->setTestDimensions(2);
  mock_backend->setTestResolution({3, 2});
  mock_backend->setTestBounds({0.0, 0.0}, {1.0, 1.0});

  // Setup required calls for the backend
  EXPECT_CALL(*mock_backend, getDimensions()).WillRepeatedly(Return(2));
  EXPECT_CALL(*mock_backend, getResolution())
      .WillRepeatedly(Return(std::vector<size_t>{3, 2}));

  // Create the adapter with our backend
  auto geometry = Geometry<Traits::GeometryType>(*mock_backend);

  // Expected coordinates from iterating through all grid points
  std::vector<std::vector<double>> expected_coords = {
      {0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {0.5, 1.0}, {1.0, 1.0}};

  // Set up the mock backend to return the expected coordinates in order
  for (size_t i = 0; i < expected_coords.size(); ++i) {
    std::vector<size_t> position;
    if (i < 3) {
      position = {i, 0};  // First row
    } else {
      position = {i - 3, 1};  // Second row
    }

    EXPECT_CALL(*mock_backend, getCoordinates(position))
        .WillOnce(Return(expected_coords[i]));
  }

  // Test iteration through the adapter
  size_t count = 0;
  for (const auto& coords : geometry) {
    if (count < expected_coords.size()) {
      EXPECT_DOUBLE_EQ(coords[0], expected_coords[count][0]);
      EXPECT_DOUBLE_EQ(coords[1], expected_coords[count][1]);
    }
    count++;
  }
  EXPECT_EQ(count, 6);  // Verify we visited all points
}
*/

/**
 * @brief Test Geometry adapter for clone and compatibility
 */
/*
TEST_F(GeometryTest, CloneAndCompatibility) {
  // Setup expectations for the backend
  EXPECT_CALL(*mock_backend1_, getDimensions())
      .WillRepeatedly(Return(default_dimensions_));
  EXPECT_CALL(*mock_backend1_, getResolution())
      .WillRepeatedly(Return(default_resolution_));

  // Test the clone method
  auto original = createGeometry();
  auto cloned = original->clone();

  // Verify the clone has the same properties
  EXPECT_EQ(cloned->getDimensions(), original->getDimensions());
  EXPECT_EQ(cloned->getResolution(), original->getResolution());

  // Test compatibility checking
  // Same properties should be compatible
  auto compatible = createGeometry();

  // Set up isCompatible to return true for same properties
  EXPECT_CALL(*mock_backend1_, isCompatible(testing::_)).WillOnce(Return(true));

  EXPECT_TRUE(original->isCompatible(*compatible));

  // Different properties should not be compatible
  auto incompatible_backend = std::make_unique<Traits::GeometryType>();
  incompatible_backend->setTestDimensions(2);
  incompatible_backend->setTestResolution({5, 5});

  EXPECT_CALL(*incompatible_backend, getDimensions()).WillRepeatedly(Return(2));
  EXPECT_CALL(*incompatible_backend, getResolution())
      .WillRepeatedly(Return(std::vector<size_t>{5, 5}));

  auto incompatible =
      std::make_unique<Geometry<Traits::GeometryType>>(*incompatible_backend);

  // Set up isCompatible to return false for different properties
  EXPECT_CALL(*mock_backend1_, isCompatible(testing::_))
      .WillOnce(Return(false));

  EXPECT_FALSE(original->isCompatible(*incompatible));
}
*/

/**
 * @brief Test Geometry adapter error handling
 */
/*
TEST_F(GeometryTest, ErrorHandling) {
  // Set up the mock backend to throw errors in specific scenarios
  EXPECT_CALL(*mock_backend1_, getDomainBounds(default_dimensions_))
      .WillOnce(
          testing::Throw(std::out_of_range("Dimension index out of range")));

  EXPECT_CALL(*mock_backend1_, setDomainBounds(0, 1.0, 0.0))
      .WillOnce(testing::Throw(
          std::invalid_argument("Min value must be less than max value")));

  // Test adapter error handling for out-of-bounds dimensions
  EXPECT_THROW(
      {
        geometry1_->getDomainBounds(default_dimensions_);  // Out of range
      },
      std::out_of_range);

  // Test adapter error handling for invalid bounds
  EXPECT_THROW(
      {
        geometry1_->setDomainBounds(0, 1.0, 0.0);  // Min > Max
      },
      std::invalid_argument);

  // Create a separate mock that will throw on getCoordinates
  auto error_mock = std::make_unique<Traits::GeometryType>();
  error_mock->setTestDimensions(2);
  error_mock->setTestResolution({5, 5});

  EXPECT_CALL(*error_mock, getDimensions()).WillRepeatedly(Return(2));
  EXPECT_CALL(*error_mock, getResolution())
      .WillRepeatedly(Return(std::vector<size_t>{5, 5}));

  EXPECT_CALL(*error_mock, getCoordinates(std::vector<size_t>{0, 5}))
      .WillOnce(testing::Throw(std::out_of_range("Index out of range")));

  auto error_geometry = Geometry<Traits::GeometryType>(*error_mock);

  // Test adapter error handling for invalid coordinates
  EXPECT_THROW(
      {
        error_geometry.getCoordinates({0, 5});  // Second index out of range
      },
      std::out_of_range);
}
*/

/**
 * @brief Test the Geometry adapter properly forwards calls to the backend
 */
/*
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
  Geometry<Traits::GeometryType> geometry(*mockBackend);

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
*/

}  // namespace metada::tests