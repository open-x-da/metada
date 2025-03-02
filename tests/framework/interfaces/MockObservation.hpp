#pragma once

#include <gmock/gmock.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "IObservation.hpp"
#include "MockConfig.hpp"
#include "utils/config/Config.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using metada::framework::Config;
using metada::framework::IConfig;
using metada::framework::IObservation;

/**
 * @brief Mock implementation of IObservation for testing
 */
class MockObservation : public IObservation {
 private:
  const Config<MockConfig>& config_;
  bool initialized_{false};
  std::vector<std::string> variableNames_;
  std::vector<size_t> dimensions_;
  std::vector<std::vector<double>> locations_;
  std::vector<double> times_;
  std::vector<int> qualityFlags_;
  std::vector<double> confidenceValues_;
  std::map<std::string, std::string> metadata_;

 public:
  /**
   * @brief Constructor that initializes observation from config
   */
  explicit MockObservation(const Config<MockConfig>& config) : config_(config) {
    // Set up default behaviors for mocked methods
    ON_CALL(*this, isInitialized()).WillByDefault([this]() {
      return initialized_;
    });
    ON_CALL(*this, getVariableNames()).WillByDefault(ReturnRef(variableNames_));
    ON_CALL(*this, getDimensions()).WillByDefault(ReturnRef(dimensions_));
    ON_CALL(*this, getLocations()).WillByDefault(ReturnRef(locations_));
    ON_CALL(*this, getTimes()).WillByDefault(ReturnRef(times_));
    ON_CALL(*this, getQualityFlags()).WillByDefault(ReturnRef(qualityFlags_));
    ON_CALL(*this, getConfidenceValues())
        .WillByDefault(ReturnRef(confidenceValues_));

    ON_CALL(*this, setMetadata(testing::_, testing::_))
        .WillByDefault(
            [this](const std::string& key, const std::string& value) {
              metadata_[key] = value;
            });
    ON_CALL(*this, getMetadata(testing::_))
        .WillByDefault([this](const std::string& key) {
          auto it = metadata_.find(key);
          if (it != metadata_.end()) {
            return it->second;
          }
          throw std::runtime_error("Metadata key not found: " + key);
        });
    ON_CALL(*this, hasMetadata(testing::_))
        .WillByDefault([this](const std::string& key) {
          return metadata_.find(key) != metadata_.end();
        });

    ON_CALL(*this, hasVariable(testing::_))
        .WillByDefault([this](const std::string& name) {
          return std::find(variableNames_.begin(), variableNames_.end(),
                           name) != variableNames_.end();
        });

    ON_CALL(*this, setLocations(testing::_))
        .WillByDefault(
            [this](const std::vector<std::vector<double>>& locations) {
              locations_ = locations;
            });
    ON_CALL(*this, setTimes(testing::_))
        .WillByDefault([this](const std::vector<double>& timestamps) {
          times_ = timestamps;
        });
    ON_CALL(*this, setQualityFlags(testing::_))
        .WillByDefault(
            [this](const std::vector<int>& flags) { qualityFlags_ = flags; });
    ON_CALL(*this, setConfidenceValues(testing::_))
        .WillByDefault([this](const std::vector<double>& values) {
          confidenceValues_ = values;
        });

    initialize(config.backend());
  }

  /**
   * @brief Copy constructor
   */
  MockObservation(const MockObservation& other)
      : config_(other.config_),
        initialized_(other.initialized_),
        variableNames_(other.variableNames_),
        dimensions_(other.dimensions_),
        locations_(other.locations_),
        times_(other.times_),
        qualityFlags_(other.qualityFlags_),
        confidenceValues_(other.confidenceValues_),
        metadata_(other.metadata_) {}

  /**
   * @brief Move constructor
   */
  MockObservation(MockObservation&& other) noexcept
      : config_(other.config_),
        initialized_(other.initialized_),
        variableNames_(std::move(other.variableNames_)),
        dimensions_(std::move(other.dimensions_)),
        locations_(std::move(other.locations_)),
        times_(std::move(other.times_)),
        qualityFlags_(std::move(other.qualityFlags_)),
        confidenceValues_(std::move(other.confidenceValues_)),
        metadata_(std::move(other.metadata_)) {
    other.initialized_ = false;
  }

  /**
   * @brief Copy assignment operator
   */
  MockObservation& operator=(const MockObservation& other) {
    if (this != &other) {
      // We can't reassign config_ since it's a reference
      initialized_ = other.initialized_;
      variableNames_ = other.variableNames_;
      dimensions_ = other.dimensions_;
      locations_ = other.locations_;
      times_ = other.times_;
      qualityFlags_ = other.qualityFlags_;
      confidenceValues_ = other.confidenceValues_;
      metadata_ = other.metadata_;
    }
    return *this;
  }

  /**
   * @brief Move assignment operator
   */
  MockObservation& operator=(MockObservation&& other) noexcept {
    if (this != &other) {
      // We can't reassign config_ since it's a reference
      initialized_ = other.initialized_;
      variableNames_ = std::move(other.variableNames_);
      dimensions_ = std::move(other.dimensions_);
      locations_ = std::move(other.locations_);
      times_ = std::move(other.times_);
      qualityFlags_ = std::move(other.qualityFlags_);
      confidenceValues_ = std::move(other.confidenceValues_);
      metadata_ = std::move(other.metadata_);
      other.initialized_ = false;
    }
    return *this;
  }

  // Lifecycle management
  MOCK_METHOD(void, initialize, (const IConfig& config), (override));
  MOCK_METHOD(void, reset, (), (override));
  MOCK_METHOD(void, validate, (), (const, override));
  MOCK_METHOD(bool, isValid, (), (const, override));
  MOCK_METHOD(bool, isInitialized, (), (const, override));

  // Copy/Move operations
  MOCK_METHOD(void, copyFrom, (const IObservation& other), (override));
  MOCK_METHOD(void, moveFrom, (IObservation && other), (override));
  MOCK_METHOD(bool, equals, (const IObservation& other), (const, override));

  // Data access
  MOCK_METHOD(void*, getData, (), (override));
  MOCK_METHOD(const void*, getData, (), (const, override));
  MOCK_METHOD(void*, getUncertainty, (), (override));
  MOCK_METHOD(const void*, getUncertainty, (), (const, override));
  MOCK_METHOD(size_t, getSize, (), (const, override));

  // Metadata operations
  MOCK_METHOD(void, setMetadata,
              (const std::string& key, const std::string& value), (override));
  MOCK_METHOD(std::string, getMetadata, (const std::string& key),
              (const, override));
  MOCK_METHOD(bool, hasMetadata, (const std::string& key), (const, override));

  // Observation information
  MOCK_METHOD(const std::vector<std::string>&, getVariableNames, (),
              (const, override));
  MOCK_METHOD(bool, hasVariable, (const std::string& name), (const, override));
  MOCK_METHOD(const std::vector<size_t>&, getDimensions, (), (const, override));

  // Spatiotemporal metadata
  MOCK_METHOD(void, setLocations,
              (const std::vector<std::vector<double>>& locations), (override));
  MOCK_METHOD(void, setTimes, (const std::vector<double>& timestamps),
              (override));
  MOCK_METHOD(const std::vector<std::vector<double>>&, getLocations, (),
              (const, override));
  MOCK_METHOD(const std::vector<double>&, getTimes, (), (const, override));

  // Arithmetic operations
  MOCK_METHOD(void, add, (const IObservation& other), (override));
  MOCK_METHOD(void, subtract, (const IObservation& other), (override));
  MOCK_METHOD(void, multiply, (double scalar), (override));

  // Quality control
  MOCK_METHOD(void, setQualityFlags, (const std::vector<int>& flags),
              (override));
  MOCK_METHOD(const std::vector<int>&, getQualityFlags, (), (const, override));
  MOCK_METHOD(void, setConfidenceValues, (const std::vector<double>& values),
              (override));
  MOCK_METHOD(const std::vector<double>&, getConfidenceValues, (),
              (const, override));

  // Get the config
  const Config<MockConfig>& config() const { return config_; }
};

}  // namespace metada::tests