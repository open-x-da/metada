#ifndef METADA_FRAMEWORK_REPR_INCREMENT_HPP_
#define METADA_FRAMEWORK_REPR_INCREMENT_HPP_

#include "State.hpp"
#include <memory>
#include <string>
#include <vector>

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Base class for state increment representations
 * 
 * Provides a generic interface for handling state increments/perturbations,
 * supporting:
 * - Addition/subtraction with states
 * - Scaling operations
 * - Norm calculations
 * - Linear algebra operations
 * 
 * @tparam T Type of increment values (typically double)
 */
template<typename T>
class Increment {
public:
    virtual ~Increment() = default;

    // Core increment operations
    virtual void initialize() = 0;
    virtual void zero() = 0;
    virtual void scale(double alpha) = 0;

    // Linear algebra operations
    virtual void axpy(double alpha, const Increment<T>& other) = 0;  // this += alpha * other
    virtual double dot(const Increment<T>& other) const = 0;         // inner product
    virtual double norm() const = 0;                                 // L2 norm

    // State operations
    virtual void addToState(State<T>& state) const = 0;             // state += increment
    virtual void differenceFromStates(const State<T>& state1,       // increment = state1 - state2
                                    const State<T>& state2) = 0;

    // Data access
    virtual T& getData() = 0;
    virtual const T& getData() const = 0;

    // Metadata
    virtual void setMetadata(const std::string& key, const std::string& value) = 0;
    virtual std::string getMetadata(const std::string& key) const = 0;

protected:
    T data_;                     // Increment data
    std::vector<size_t> dimensions_;  // Dimensions of the increment
    bool initialized_{false};     // Initialization flag
};

}}} // namespace metada::framework::repr

#endif // METADA_FRAMEWORK_REPR_INCREMENT_HPP_ 