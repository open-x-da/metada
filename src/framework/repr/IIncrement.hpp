#ifndef METADA_FRAMEWORK_REPR_IINCREMENT_HPP_
#define METADA_FRAMEWORK_REPR_IINCREMENT_HPP_

#include <string>
#include <vector>
#include "State.hpp"

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Abstract interface for increment implementations
 *
 * This interface defines the contract that all increment implementations must follow.
 * It provides a unified API for handling state increments/perturbations in scientific models.
 *
 * Key features:
 * - Addition/subtraction with states
 * - Scaling operations
 * - Norm calculations
 * - Linear algebra operations
 */
class IIncrement {
public:
    virtual ~IIncrement() = default;

    // Core increment operations
    virtual void initialize() = 0;
    virtual void zero() = 0;
    virtual void scale(double alpha) = 0;

    // Linear algebra operations
    virtual void axpy(double alpha, const IIncrement& other) = 0;  // this += alpha * other
    virtual double dot(const IIncrement& other) const = 0;         // inner product
    virtual double norm() const = 0;                               // L2 norm

    // State operations
    virtual void addToState(IState& state) const = 0;             // state += increment
    virtual void differenceFromStates(const IState& state1,       // increment = state1 - state2
                                    const IState& state2) = 0;

    // Data access
    virtual void* getData() = 0;
    virtual const void* getData() const = 0;

    // Metadata
    virtual void setMetadata(const std::string& key, const std::string& value) = 0;
    virtual std::string getMetadata(const std::string& key) const = 0;

    // Increment information
    virtual const std::vector<size_t>& getDimensions() const = 0;
    virtual bool isInitialized() const = 0;
};

}}} // namespace metada::framework::repr

#endif // METADA_FRAMEWORK_REPR_IINCREMENT_HPP_ 