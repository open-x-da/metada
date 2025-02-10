#ifndef METADA_FRAMEWORK_REPR_STATE_HPP_
#define METADA_FRAMEWORK_REPR_STATE_HPP_

#include <vector>
#include <string>

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Base class for representing model state variables
 * 
 * Provides a generic interface for handling state variables in scientific models.
 * This includes support for:
 * - Multiple variable types (scalar, vector, matrix)
 * - Named dimensions
 * - Metadata handling
 * - Basic mathematical operations
 */
template<typename T>
class State {
public:
    virtual ~State() = default;

    // Core state operations
    virtual void initialize() = 0;
    virtual void reset() = 0;
    virtual void validate() const = 0;

    // Data access
    virtual T& getData() = 0;
    virtual const T& getData() const = 0;

    // Metadata
    virtual void setMetadata(const std::string& key, const std::string& value) = 0;
    virtual std::string getMetadata(const std::string& key) const = 0;

protected:
    std::vector<std::string> variable_names_;
    std::vector<size_t> dimensions_;
};

}}} // namespace metada::framework::repr

#endif // METADA_FRAMEWORK_REPR_STATE_HPP_ 