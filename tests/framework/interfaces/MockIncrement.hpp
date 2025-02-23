#ifndef METADA_TESTS_FRAMEWORK_REPR_MOCK_INCREMENT_HPP_
#define METADA_TESTS_FRAMEWORK_REPR_MOCK_INCREMENT_HPP_

#include <gmock/gmock.h>

#include "IIncrement.hpp"

namespace metada {
namespace framework {
namespace repr {
namespace tests {

class MockIncrement : public IIncrement {
 public:
  MOCK_METHOD(void, initialize, (), (override));
  MOCK_METHOD(void, zero, (), (override));
  MOCK_METHOD(void, scale, (double alpha), (override));

  MOCK_METHOD(void, axpy, (double alpha, const IIncrement& other), (override));
  MOCK_METHOD(double, dot, (const IIncrement& other), (const, override));
  MOCK_METHOD(double, norm, (), (const, override));

  MOCK_METHOD(void, addToState, (IState & state), (const, override));
  MOCK_METHOD(void, differenceFromStates,
              (const IState& state1, const IState& state2), (override));

  MOCK_METHOD(void*, getData, (), (override));
  MOCK_METHOD(const void*, getData, (), (const, override));

  MOCK_METHOD(void, setMetadata,
              (const std::string& key, const std::string& value), (override));
  MOCK_METHOD(std::string, getMetadata, (const std::string& key),
              (const, override));

  MOCK_METHOD(const std::vector<size_t>&, getDimensions, (), (const, override));
  MOCK_METHOD(bool, isInitialized, (), (const, override));
};

}  // namespace tests
}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_TESTS_FRAMEWORK_REPR_MOCK_INCREMENT_HPP_