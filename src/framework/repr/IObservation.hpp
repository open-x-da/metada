#ifndef METADA_FRAMEWORK_REPR_IOBSERVATION_HPP_
#define METADA_FRAMEWORK_REPR_IOBSERVATION_HPP_

#include <map>
#include <string>
#include <vector>

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Abstract interface for observation implementations
 *
 * This interface defines the contract that all observation implementations must
 * follow. It provides a unified API for handling observational data.
 *
 * Key features:
 * - Multiple observation types and sources
 * - Observation errors/uncertainties
 * - Quality control flags
 * - Spatial and temporal metadata
 */
class IObservation {
 public:
  virtual ~IObservation() = default;

  // Core observation operations
  virtual void initialize() = 0;
  virtual void validate() const = 0;
  virtual bool isValid() const = 0;

  // Data access
  virtual void* getValue() = 0;
  virtual const void* getValue() const = 0;
  virtual void* getError() = 0;
  virtual const void* getError() const = 0;

  // Metadata management
  virtual void setLocation(double lat, double lon, double height) = 0;
  virtual void setTime(double timestamp) = 0;
  virtual void setQualityFlag(int flag) = 0;

  // Observation attributes
  virtual void setAttribute(const std::string& key,
                            const std::string& value) = 0;
  virtual std::string getAttribute(const std::string& key) const = 0;

  // Location access
  virtual const std::vector<double>& getLocation() const = 0;
  virtual double getTimestamp() const = 0;
  virtual int getQualityFlag() const = 0;
};

}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_FRAMEWORK_REPR_IOBSERVATION_HPP_