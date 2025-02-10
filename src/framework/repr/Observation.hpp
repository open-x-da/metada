#ifndef METADA_FRAMEWORK_REPR_OBSERVATION_HPP_
#define METADA_FRAMEWORK_REPR_OBSERVATION_HPP_

#include <memory>
#include <string>
#include <vector>
#include <map>

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Base class for observation representations
 * 
 * Provides a generic interface for handling observational data, supporting:
 * - Multiple observation types and sources
 * - Observation errors/uncertainties
 * - Quality control flags
 * - Spatial and temporal metadata
 * 
 * @tparam T Type of observation values (typically double)
 */
template<typename T>
class Observation {
public:
    virtual ~Observation() = default;

    // Core observation operations
    virtual void initialize() = 0;
    virtual void validate() const = 0;
    virtual bool isValid() const = 0;

    // Data access
    virtual T& getValue() = 0;
    virtual const T& getValue() const = 0;
    virtual T& getError() = 0;
    virtual const T& getError() const = 0;

    // Metadata management
    virtual void setLocation(double lat, double lon, double height) = 0;
    virtual void setTime(double timestamp) = 0;
    virtual void setQualityFlag(int flag) = 0;
    
    // Observation attributes
    virtual void setAttribute(const std::string& key, const std::string& value) = 0;
    virtual std::string getAttribute(const std::string& key) const = 0;

protected:
    T value_;                    // Observation value
    T error_;                    // Observation error/uncertainty
    std::vector<double> location_;  // Spatial coordinates
    double timestamp_{0.0};      // Observation time
    int quality_flag_{0};        // Quality control flag
    std::map<std::string, std::string> attributes_;  // Additional metadata
};

}}} // namespace metada::framework::repr

#endif // METADA_FRAMEWORK_REPR_OBSERVATION_HPP_ 