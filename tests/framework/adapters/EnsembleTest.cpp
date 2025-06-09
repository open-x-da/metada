#include <gmock/gmock.h>

#include <filesystem>
#include <fstream>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "MockBackendTraits.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "State.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::MockBackendTag;
using StateType = fwk::State<BackendTag>;
using Ensemble = fwk::Ensemble<BackendTag>;
using ConfigType = fwk::Config<BackendTag>;
using GeometryType = fwk::Geometry<BackendTag>;

namespace metada::tests {

using ::testing::InSequence;
using ::testing::Invoke;
using ::testing::NiceMock;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;

class EnsembleTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create a temporary config file
    temp_config_path_ =
        std::filesystem::temp_directory_path() / "test_config.json";
    std::ofstream config_file(temp_config_path_);
    config_file << R"({"ensemble_size": 64})";
    config_file.close();

    // Initialize config with the temporary file
    config_ = std::make_unique<ConfigType>(temp_config_path_.string());

    // Set up default mock behavior
    ON_CALL(config_->backend(), Get("ensemble_size"))
        .WillByDefault(Return(fwk::ConfigValue(64)));

    // Create mock geometry
    geometry_ = std::make_unique<GeometryType>(*config_);
  }

  void TearDown() override {
    // Clean up the temporary config file
    std::filesystem::remove(temp_config_path_);
    config_.reset();
    geometry_.reset();
  }

  std::filesystem::path temp_config_path_;
  std::unique_ptr<ConfigType> config_;
  std::unique_ptr<GeometryType> geometry_;
};

TEST_F(EnsembleTest, InitializationCreatesCorrectNumberOfMembers) {
  Ensemble ensemble(*config_, *geometry_);
  EXPECT_EQ(ensemble.Size(), 64);
}

TEST_F(EnsembleTest, MemberAccessIsValid) {
  Ensemble ensemble(*config_, *geometry_);
  EXPECT_NO_THROW(ensemble.GetMember(0));
  EXPECT_NO_THROW(ensemble.GetMember(63));
}

TEST_F(EnsembleTest, ComputeMean) {
  Ensemble ensemble(*config_, *geometry_);
  ensemble.ComputeMean();
  StateType& mean = ensemble.Mean();
  EXPECT_TRUE(mean.isInitialized());
}

TEST_F(EnsembleTest, OutOfBoundsAccessThrows) {
  Ensemble ensemble(*config_, *geometry_);
  EXPECT_THROW(ensemble.GetMember(64), std::out_of_range);
  EXPECT_THROW(ensemble.GetPerturbation(64), std::out_of_range);
}

TEST_F(EnsembleTest, ComputePerturbationsCreatesValidPerturbations) {
  Ensemble ensemble(*config_, *geometry_);
  ensemble.ComputePerturbations();
  EXPECT_EQ(ensemble.Size(), 64);
  for (size_t i = 0; i < ensemble.Size(); ++i) {
    EXPECT_TRUE(ensemble.GetPerturbation(i).isInitialized());
  }
}

}  // namespace metada::tests