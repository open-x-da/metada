#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "IObservation.hpp"

namespace metada::framework::repr {

/**
 * @brief Main observation class template providing a generic interface to
 * observation implementations
 *
 * This class template provides a static interface for handling observational
 * data using a backend specified by the ObservationBackend template parameter.
 * The backend must implement the IObservation interface.
 *
 * @tparam ObservationBackend The observation backend type that implements
 * IObservation
 */
template <typename ObservationBackend>
class Observation {
 private:
  ObservationBackend backend_;  ///< Instance of the observation backend

 public:
  /** @brief Default constructor */
  Observation() = default;

  /** @brief Get direct access to the backend instance */
  ObservationBackend& backend() { return backend_; }

  /** @brief Get const access to the backend instance */
  const ObservationBackend& backend() const { return backend_; }

  // Core observation operations
  void initialize() { backend_.initialize(); }
  void validate() const { backend_.validate(); }
  bool isValid() const { return backend_.isValid(); }

  // Data access
  template <typename T>
  T& getValue() {
    return *static_cast<T*>(backend_.getValue());
  }

  template <typename T>
  const T& getValue() const {
    return *static_cast<const T*>(backend_.getValue());
  }

  template <typename T>
  T& getError() {
    return *static_cast<T*>(backend_.getError());
  }

  template <typename T>
  const T& getError() const {
    return *static_cast<const T*>(backend_.getError());
  }

  // Metadata management
  void setLocation(double lat, double lon, double height) {
    backend_.setLocation(lat, lon, height);
  }

  void setTime(double timestamp) { backend_.setTime(timestamp); }
  void setQualityFlag(int flag) { backend_.setQualityFlag(flag); }

  // Observation attributes
  void setAttribute(const std::string& key, const std::string& value) {
    backend_.setAttribute(key, value);
  }

  std::string getAttribute(const std::string& key) const {
    return backend_.getAttribute(key);
  }

  // Location access
  const std::vector<double>& getLocation() const {
    return backend_.getLocation();
  }
  double getTimestamp() const { return backend_.getTimestamp(); }
  int getQualityFlag() const { return backend_.getQualityFlag(); }
};

}  // namespace metada::framework::repr