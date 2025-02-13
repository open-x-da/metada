#ifndef METADA_TESTS_FRAMEWORK_REPR_MOCK_STATE_HPP_
#define METADA_TESTS_FRAMEWORK_REPR_MOCK_STATE_HPP_

#include <gmock/gmock.h>

#include "IConfig.hpp"
#include "IState.hpp"

using metada::framework::tools::config::IConfig;

namespace metada {
namespace framework {
namespace repr {
namespace tests {

class MockState : public IState {
 public:
  // Core state operations
  MOCK_METHOD(void, initialize, (const IConfig& config), (override));
  MOCK_METHOD(void, reset, (), (override));
  MOCK_METHOD(void, validate, (), (const, override));
  MOCK_METHOD(bool, isInitialized, (), (const, override));

  MOCK_METHOD(void*, getData, (), (override));
  MOCK_METHOD(const void*, getData, (), (const, override));

  MOCK_METHOD(void, setMetadata,
              (const std::string& key, const std::string& value), (override));
  MOCK_METHOD(std::string, getMetadata, (const std::string& key),
              (const, override));

  MOCK_METHOD(const std::vector<std::string>&, getVariableNames, (),
              (const, override));
  MOCK_METHOD(const std::vector<size_t>&, getDimensions, (), (const, override));

  // Copy operations
  MOCK_METHOD(void, copyFrom, (const IState& other), (override));
  MOCK_METHOD(void, moveFrom, (IState && other), (override));
  MOCK_METHOD(bool, equals, (const IState& other), (const, override));
};

}  // namespace tests
}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_TESTS_FRAMEWORK_REPR_MOCK_STATE_HPP_