#pragma once

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace metada::framework {

/**
 * @brief Abstract interface for observation implementations in data
 * assimilation systems
 *
 * This interface defines the contract for handling observational data,
 * providing a unified API for DA systems to interact with observation
 * databases.
 *
 * Key features:
 * - Type-safe value and error access
 * - Comprehensive metadata for spatiotemporal context
 * - Quality control and uncertainty quantification
 * - Flexible attribute system for extensibility
 */
class IObservation {
 public:
  virtual ~IObservation() = default;

  // Lifecycle management
  virtual void initialize() = 0;
  virtual void validate() const = 0;
  virtual bool isValid() const = 0;

  // Data access - with explicit type information
  virtual void* getData() = 0;
  virtual const void* getData() const = 0;
  virtual void* getUncertainty() = 0;
  virtual const void* getUncertainty() const = 0;
  virtual size_t getDataSize() const = 0;

  // Spatiotemporal metadata
  virtual void setLocation(double lat, double lon, double elevation) = 0;
  virtual void setTime(double timestamp) = 0;
  virtual const std::vector<double>& getLocation() const = 0;
  virtual double getTimestamp() const = 0;

  // Quality control
  virtual void setQualityFlag(int flag) = 0;
  virtual int getQualityFlag() const = 0;
  virtual void setConfidence(double value) = 0;
  virtual double getConfidence() const = 0;

  // Extensible attributes
  virtual void setAttribute(const std::string& key,
                            const std::string& value) = 0;
  virtual std::string getAttribute(const std::string& key) const = 0;
  virtual bool hasAttribute(const std::string& key) const = 0;
  virtual std::map<std::string, std::string> getAllAttributes() const = 0;

  // Observation metadata
  virtual void setObsType(const std::string& type) = 0;
  virtual std::string getObsType() const = 0;
  virtual void setSource(const std::string& source) = 0;
  virtual std::string getSource() const = 0;
  virtual void setInstrument(const std::string& instrument) = 0;
  virtual std::string getInstrument() const = 0;
};

}  // namespace metada::framework