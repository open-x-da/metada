#ifndef METADA_TESTS_FRAMEWORK_REPR_MOCK_MODEL_HPP_
#define METADA_TESTS_FRAMEWORK_REPR_MOCK_MODEL_HPP_

#include <gmock/gmock.h>

#include "IModel.hpp"

namespace metada {
namespace framework {
namespace repr {
namespace tests {

class MockModel : public IModel {
 public:
  MOCK_METHOD(void, initialize, (), (override));
  MOCK_METHOD(void, finalize, (), (override));

  MOCK_METHOD(void, step, (const IState& state1, IState& state2),
              (const, override));
  MOCK_METHOD(void, stepTL,
              (const IIncrement& increment1, IIncrement& increment2),
              (const, override));
  MOCK_METHOD(void, stepAD,
              (const IIncrement& increment2, IIncrement& increment1),
              (const, override));

  MOCK_METHOD(void, setParameter, (const std::string& name, double value),
              (override));
  MOCK_METHOD(double, getParameter, (const std::string& name),
              (const, override));

  MOCK_METHOD(const std::vector<std::string>&, getRequiredStateVariables, (),
              (const, override));
  MOCK_METHOD(bool, isInitialized, (), (const, override));
};

}  // namespace tests
}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_TESTS_FRAMEWORK_REPR_MOCK_MODEL_HPP_