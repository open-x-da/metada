#pragma once

#include <gmock/gmock.h>

#include "IIncrement.hpp"

namespace metada {

namespace framework {
class IState;
}  // namespace framework

namespace tests {

using metada::framework::IIncrement;
using metada::framework::IState;

class MockIncrement : public IIncrement {
 public:
  MOCK_METHOD(void, initialize, (), (override));
  MOCK_METHOD(void, zero, (), (override));
  MOCK_METHOD(void, scale, (double alpha), (override));

  MOCK_METHOD(void, axpy, (double alpha, const IIncrement& other), (override));
  MOCK_METHOD(double, dot, (const IIncrement& other), (const, override));
  MOCK_METHOD(double, norm, (), (const, override));

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
}  // namespace metada
