#ifndef METADA_FRAMEWORK_REPR_STATE_HPP_
#define METADA_FRAMEWORK_REPR_STATE_HPP_

#include "IState.hpp"
#include <vector>
#include <string>

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Main state class template providing a generic interface to state implementations
 * 
 * @tparam StateBackend The state backend type that implements IState
 */
template<typename StateBackend>
class State {
private:
    StateBackend& backend_;  ///< Instance of the state backend

public:
    /** @brief Constructor that takes a backend reference */
    explicit State(StateBackend& backend) : backend_(backend) {}

    /** @brief Copy constructor */
    explicit State(const State& other) : backend_(other.backend_) {}

    /** @brief Copy assignment operator */
    State& operator=(const State& other) {
        if (this != &other) {
            backend_.copyFrom(other.backend_);
        }
        return *this;
    }

    /** @brief Equality operator */
    bool operator==(const State& other) const {
        return backend_.equals(other.backend_);
    }

    /** @brief Inequality operator */
    bool operator!=(const State& other) const {
        return !(*this == other);
    }

    /** @brief Get direct access to the backend instance */
    StateBackend& backend() { return backend_; }

    /** @brief Get const access to the backend instance */
    const StateBackend& backend() const { return backend_; }

    // Data access
    template<typename T>
    T& getData() { return *static_cast<T*>(backend_.getData()); }
    
    template<typename T>
    const T& getData() const { return *static_cast<const T*>(backend_.getData()); }

    // Metadata operations
    void setMetadata(const std::string& key, const std::string& value) {
        backend_.setMetadata(key, value);
    }

    std::string getMetadata(const std::string& key) const {
        return backend_.getMetadata(key);
    }

    // State information
    const std::vector<std::string>& getVariableNames() const {
        return backend_.getVariableNames();
    }

    const std::vector<size_t>& getDimensions() const {
        return backend_.getDimensions();
    }
};

}}} // namespace metada::framework::repr

#endif // METADA_FRAMEWORK_REPR_STATE_HPP_ 