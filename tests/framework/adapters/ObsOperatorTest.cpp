#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

#include "AppTraits.hpp"
#include "ApplicationContext.hpp"
#include "MockConfig.hpp"
#include "MockIncrement.hpp"
#include "MockLogger.hpp"
#include "MockObsOperator.hpp"
#include "MockObservation.hpp"
#include "MockState.hpp"
#include "ObsOperator.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Increment;
using framework::Observation;
using framework::ObsOperator;
using framework::State;
using framework::runs::ApplicationContext;

using Traits = AppTraits<MockLogger, MockConfig, MockState, MockIncrement,
                         MockObservation, MockObsOperator>;

class ObsOperatorTest : public ::testing::Test {
 protected:
  /** @brief Application context instance */
  std::unique_ptr<ApplicationContext<Traits>> context_;

  // Test data
  std::vector<std::string> state_vars_;
  std::vector<std::string> obs_vars_;

  // Adapters
  std::unique_ptr<State<Traits::StateType>> stateAdapter_;
  std::unique_ptr<Observation<Traits::ObservationType>> obsAdapter_;
  std::unique_ptr<Increment<Traits::IncrementType>> incrementAdapter_;

  // System under test
  std::unique_ptr<ObsOperator<Traits::StateType, Traits::IncrementType,
                              Traits::ObservationType, Traits::ObsOperatorType>>
      obsOperator_;

  void SetUp() override {
    context_ = std::make_unique<ApplicationContext<Traits>>("ObsOperatorTest");

    // Setup test data
    state_vars_ = {"temperature", "pressure"};
    obs_vars_ = {"radiance"};
  }

  void TearDown() override {
    context_.reset();
    state_vars_.clear();
    obs_vars_.clear();
  }
};

}  // namespace metada::tests
