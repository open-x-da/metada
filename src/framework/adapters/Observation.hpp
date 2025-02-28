#pragma once

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "IObservation.hpp"

namespace metada::framework {

/**
 * @brief Adapter class for observation implementations in data assimilation
 * systems
 *
 * This template class provides a type-safe interface for handling observational
 * data using a backend that implements the IObservation interface.
 *
 * Features:
 * - Type-safe data access through templates
 * - Delegation to backend implementation
 * - Comprehensive error handling
 * - Fluent interface for chaining operations
 *
 * @tparam ObservationBackend The backend implementation type
 */
template <typename ObservationBackend>
class Observation {
 private:
  ObservationBackend& backend_;  ///< Backend implementation instance

 public:
  /** @brief Default constructor */
  Observation() = default;

  /** @brief Constructor with initialization */
  Observation(ObservationBackend& backend) : backend_(backend) { initialize(); }

  /** @brief Constructor with initialization */
  Observation(ObservationBackend&& backend) : backend_(std::move(backend)) {
    initialize();
  }

  /** @brief Access to backend instance */
  ObservationBackend& backend() { return backend_; }
  const ObservationBackend& backend() const { return backend_; }

  // Lifecycle management
  Observation& initialize() {
    backend_.initialize();
    return *this;
  }

  void validate() const { backend_.validate(); }

  bool isValid() const { return backend_.isValid(); }

  // Type-safe data access
  template <typename T>
  T& getData() {
    if (!isValid()) {
      throw std::runtime_error("Accessing data from invalid observation");
    }
    return *static_cast<T*>(backend_.getData());
  }

  template <typename T>
  const T& getData() const {
    if (!isValid()) {
      throw std::runtime_error("Accessing data from invalid observation");
    }
    return *static_cast<const T*>(backend_.getData());
  }

  template <typename T>
  T& getUncertainty() {
    return *static_cast<T*>(backend_.getUncertainty());
  }

  template <typename T>
  const T& getUncertainty() const {
    return *static_cast<const T*>(backend_.getUncertainty());
  }

  size_t getDataSize() const { return backend_.getDataSize(); }

  // Spatiotemporal metadata with fluent interface
  Observation& setLocation(double lat, double lon, double elevation) {
    backend_.setLocation(lat, lon, elevation);
    return *this;
  }

  Observation& setTime(double timestamp) {
    backend_.setTime(timestamp);
    return *this;
  }

  const std::vector<double>& getLocation() const {
    return backend_.getLocation();
  }

  double getTimestamp() const { return backend_.getTimestamp(); }

  // Quality control with fluent interface
  Observation& setQualityFlag(int flag) {
    backend_.setQualityFlag(flag);
    return *this;
  }

  int getQualityFlag() const { return backend_.getQualityFlag(); }

  Observation& setConfidence(double value) {
    backend_.setConfidence(value);
    return *this;
  }

  double getConfidence() const { return backend_.getConfidence(); }

  // Extensible attributes with fluent interface
  Observation& setAttribute(const std::string& key, const std::string& value) {
    backend_.setAttribute(key, value);
    return *this;
  }

  std::string getAttribute(const std::string& key) const {
    return backend_.getAttribute(key);
  }

  bool hasAttribute(const std::string& key) const {
    return backend_.hasAttribute(key);
  }

  std::map<std::string, std::string> getAllAttributes() const {
    return backend_.getAllAttributes();
  }

  // Observation metadata with fluent interface
  Observation& setObsType(const std::string& type) {
    backend_.setObsType(type);
    return *this;
  }

  std::string getObsType() const { return backend_.getObsType(); }

  Observation& setSource(const std::string& source) {
    backend_.setSource(source);
    return *this;
  }

  std::string getSource() const { return backend_.getSource(); }

  Observation& setInstrument(const std::string& instrument) {
    backend_.setInstrument(instrument);
    return *this;
  }

  std::string getInstrument() const { return backend_.getInstrument(); }
};

}  // namespace metada::framework