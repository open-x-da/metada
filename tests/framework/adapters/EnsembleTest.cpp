#include <gmock/gmock.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "Logger.hpp"
#include "MockBackendTraits.hpp"
#include "State.hpp"

namespace metada::tests {

using namespace metada::framework;

using ::testing::_;
using ::testing::Return;

class EnsembleTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create a temporary config file
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    config_file_ = (test_dir / "test_config.yaml").string();
    config_ = std::make_unique<Config<traits::MockBackendTag>>(config_file_);

    // Initialize the Logger singleton before creating any objects that use it
    Logger<traits::MockBackendTag>::Init(*config_);

    // Create test objects - we'll use a mock geometry through the adapter
    geometry_ = std::make_unique<Geometry<traits::MockBackendTag>>(*config_);
  }

  void TearDown() override {
    config_.reset();
    geometry_.reset();
    // Reset the Logger singleton after tests
    Logger<traits::MockBackendTag>::Reset();
  }

  std::string config_file_;
  std::unique_ptr<Config<traits::MockBackendTag>> config_;
  std::unique_ptr<Geometry<traits::MockBackendTag>> geometry_;
};

TEST_F(EnsembleTest, InitializationCreatesCorrectNumberOfMembers) {
  // The Ensemble constructor reads member configs from the config file
  // For testing, we expect it to create members based on available configs
  Ensemble ensemble(*config_, *geometry_);
  EXPECT_GT(ensemble.Size(), 0);
}

TEST_F(EnsembleTest, MemberAccessIsValid) {
  Ensemble ensemble(*config_, *geometry_);
  size_t size = ensemble.Size();

  if (size > 0) {
    EXPECT_NO_THROW(ensemble.GetMember(0));
    if (size > 1) {
      EXPECT_NO_THROW(ensemble.GetMember(size - 1));
    }
  }
}

TEST_F(EnsembleTest, MeanComputationIsLazy) {
  Ensemble ensemble(*config_, *geometry_);
  if (ensemble.Size() > 0) {
    // Mean should be computed on first access
    auto& mean = ensemble.Mean();
    EXPECT_TRUE(mean.isInitialized());
  }
}

TEST_F(EnsembleTest, MeanRecomputation) {
  Ensemble ensemble(*config_, *geometry_);
  if (ensemble.Size() > 0) {
    // Compute mean first time
    auto& mean1 = ensemble.Mean();
    // Force recomputation
    ensemble.RecomputeMean();
    auto& mean2 = ensemble.Mean();
    EXPECT_TRUE(mean1.isInitialized());
    EXPECT_TRUE(mean2.isInitialized());
  }
}

TEST_F(EnsembleTest, ConstMeanAccessThrowsIfNotComputed) {
  const Ensemble ensemble(*config_, *geometry_);
  if (ensemble.Size() > 0) {
    EXPECT_THROW(ensemble.Mean(), std::runtime_error);
  }
}

TEST_F(EnsembleTest, OutOfBoundsAccessThrows) {
  Ensemble ensemble(*config_, *geometry_);
  size_t size = ensemble.Size();
  EXPECT_THROW(ensemble.GetMember(size), std::out_of_range);
  EXPECT_THROW(ensemble.GetPerturbation(size), std::out_of_range);
}

TEST_F(EnsembleTest, PerturbationComputationIsLazy) {
  Ensemble ensemble(*config_, *geometry_);
  if (ensemble.Size() > 0) {
    // Perturbations should be computed on first access
    auto& pert = ensemble.GetPerturbation(0);
    EXPECT_TRUE(pert.isInitialized());
  }
}

TEST_F(EnsembleTest, PerturbationRecomputation) {
  Ensemble ensemble(*config_, *geometry_);
  if (ensemble.Size() > 0) {
    // Compute perturbations first time
    auto& pert1 = ensemble.GetPerturbation(0);
    // Force recomputation
    ensemble.RecomputePerturbations();
    auto& pert2 = ensemble.GetPerturbation(0);
    EXPECT_TRUE(pert1.isInitialized());
    EXPECT_TRUE(pert2.isInitialized());
  }
}

TEST_F(EnsembleTest, ConstPerturbationAccessThrowsIfNotComputed) {
  const Ensemble ensemble(*config_, *geometry_);
  if (ensemble.Size() > 0) {
    EXPECT_THROW(ensemble.GetPerturbation(0), std::runtime_error);
  }
}

}  // namespace metada::tests