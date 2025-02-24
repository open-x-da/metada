#pragma once

#include <gmock/gmock.h>

#include <memory>

#include "IModel.hpp"
#include "MockState.hpp"
#include "State.hpp"

namespace metada::framework::interfaces::tests {

class MockModel : public IModel {
 public:
  MOCK_METHOD(void, initialize, (), (override));
  MOCK_METHOD(void, finalize, (), (override));
  MOCK_METHOD(void, step, (double dt), (override));

  MOCK_METHOD(void, setState, (const IState&), (override));
  MOCK_METHOD(const IState&, getState, (), (const, override));

  MOCK_METHOD(void, setParameter, (const std::string& name, double value),
              (override));
  MOCK_METHOD(double, getParameter, (const std::string& name),
              (const, override));

  MOCK_METHOD(std::string, getName, (), (const, override));
  MOCK_METHOD(std::string, getVersion, (), (const, override));
  MOCK_METHOD(double, getCurrentTime, (), (const, override));
  MOCK_METHOD(bool, isInitialized, (), (const, override));
};

}  // namespace metada::framework::interfaces::tests
