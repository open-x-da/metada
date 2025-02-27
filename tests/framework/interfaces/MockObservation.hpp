#pragma once

#include <gmock/gmock.h>

#include "IObservation.hpp"

namespace metada::tests {

using metada::framework::IObservation;

/**
 * @brief Mock implementation of IObservation for testing
 */
class MockObservation : public IObservation {
 public:
  // Lifecycle management
  MOCK_METHOD(void, initialize, (), (override));
  MOCK_METHOD(void, validate, (), (const, override));
  MOCK_METHOD(bool, isValid, (), (const, override));

  // Data access
  MOCK_METHOD(void*, getData, (), (override));
  MOCK_METHOD(const void*, getData, (), (const, override));
  MOCK_METHOD(void*, getUncertainty, (), (override));
  MOCK_METHOD(const void*, getUncertainty, (), (const, override));
  MOCK_METHOD(size_t, getDataSize, (), (const, override));

  // Spatiotemporal metadata
  MOCK_METHOD(void, setLocation, (double, double, double), (override));
  MOCK_METHOD(void, setTime, (double), (override));
  MOCK_METHOD(const std::vector<double>&, getLocation, (), (const, override));
  MOCK_METHOD(double, getTimestamp, (), (const, override));

  // Quality control
  MOCK_METHOD(void, setQualityFlag, (int), (override));
  MOCK_METHOD(int, getQualityFlag, (), (const, override));
  MOCK_METHOD(void, setConfidence, (double), (override));
  MOCK_METHOD(double, getConfidence, (), (const, override));

  // Extensible attributes
  MOCK_METHOD(void, setAttribute, (const std::string&, const std::string&),
              (override));
  MOCK_METHOD(std::string, getAttribute, (const std::string&),
              (const, override));
  MOCK_METHOD(bool, hasAttribute, (const std::string&), (const, override));
  MOCK_METHOD((std::map<std::string, std::string>), getAllAttributes, (),
              (const, override));

  // Observation metadata
  MOCK_METHOD(void, setObsType, (const std::string&), (override));
  MOCK_METHOD(std::string, getObsType, (), (const, override));
  MOCK_METHOD(void, setSource, (const std::string&), (override));
  MOCK_METHOD(std::string, getSource, (), (const, override));
  MOCK_METHOD(void, setInstrument, (const std::string&), (override));
  MOCK_METHOD(std::string, getInstrument, (), (const, override));
};

}  // namespace metada::tests