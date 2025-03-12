#pragma once

#include <memory>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include "state_c_api.h"

namespace metada {
namespace backends {
namespace interfaces {

/**
 * @brief C++ wrapper for the Lorenz63 Fortran State class
 */
class Lorenz63State {
public:
    /**
     * @brief Construct a new Lorenz63State with given initial values
     * 
     * @param x Initial x value
     * @param y Initial y value
     * @param z Initial z value
     */
    Lorenz63State(float x = 0.0f, float y = 0.0f, float z = 0.0f) 
        : ptr_(state_create(x, y, z), state_deleter) {
        if (!ptr_) {
            throw std::runtime_error("Failed to create Lorenz63 state");
        }
    }

    /**
     * @brief Get the x component
     * 
     * @return float x value
     */
    float getX() const {
        return state_get_x(ptr_.get());
    }

    /**
     * @brief Get the y component
     * 
     * @return float y value
     */
    float getY() const {
        return state_get_y(ptr_.get());
    }

    /**
     * @brief Get the z component
     * 
     * @return float z value
     */
    float getZ() const {
        return state_get_z(ptr_.get());
    }

    /**
     * @brief Set the x component
     * 
     * @param x New x value
     */
    void setX(float x) {
        state_set_x(ptr_.get(), x);
    }

    /**
     * @brief Set the y component
     * 
     * @param y New y value
     */
    void setY(float y) {
        state_set_y(ptr_.get(), y);
    }

    /**
     * @brief Set the z component
     * 
     * @param z New z value
     */
    void setZ(float z) {
        state_set_z(ptr_.get(), z);
    }

    /**
     * @brief Get all components as a vector
     * 
     * @return std::vector<float> {x, y, z}
     */
    std::vector<float> getComponents() const {
        return {getX(), getY(), getZ()};
    }

    /**
     * @brief Calculate distance between this state and another
     * 
     * @param other Another Lorenz63State
     * @return float Euclidean distance
     */
    float distance(const Lorenz63State& other) const {
        return state_distance(ptr_.get(), other.ptr_.get());
    }

    /**
     * @brief Get the underlying Fortran pointer (for internal use)
     * 
     * @return void* Raw pointer to Fortran state
     */
    void* getPtr() const {
        return ptr_.get();
    }

    /**
     * @brief Friend declaration for the stream output operator
     */
    friend std::ostream& operator<<(std::ostream& os, const Lorenz63State& state);

private:
    // Custom deleter for Fortran state
    static void state_deleter(void* ptr) {
        if (ptr) {
            state_destroy(ptr);
        }
    }

    // Smart pointer to manage Fortran state
    std::shared_ptr<void> ptr_;
};

/**
 * @brief Stream output operator for Lorenz63State
 * 
 * @param os Output stream
 * @param state Lorenz63State to print
 * @return std::ostream& Reference to the output stream
 */
inline std::ostream& operator<<(std::ostream& os, const Lorenz63State& state) {
    os << "Lorenz63State(" 
       << std::fixed << std::setprecision(6) << state.getX() << ", " 
       << std::fixed << std::setprecision(6) << state.getY() << ", " 
       << std::fixed << std::setprecision(6) << state.getZ() << ")";
    return os;
}

} // namespace interfaces
} // namespace backends
} // namespace metada 