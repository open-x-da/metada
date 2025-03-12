#pragma once

#include <memory>
#include <iostream>
#include <iomanip>
#include "geometry_c_api.h"
#include "../../state/lorenz63/State.hpp"

namespace metada::backends::lorenz63 {

/**
 * @brief C++ wrapper for the Lorenz63 Fortran Geometry class
 * This class defines the boundaries of the phase space for the Lorenz63 system
 */
class Geometry {
public:
    /**
     * @brief Construct a new Geometry with given boundaries
     * 
     * @param x_min Minimum x value
     * @param x_max Maximum x value
     * @param y_min Minimum y value
     * @param y_max Maximum y value
     * @param z_min Minimum z value
     * @param z_max Maximum z value
     */
    Geometry(float x_min = -30.0f, float x_max = 30.0f, 
                     float y_min = -30.0f, float y_max = 30.0f, 
                     float z_min = 0.0f, float z_max = 60.0f)
        : ptr_(geometry_create(x_min, x_max, y_min, y_max, z_min, z_max), geometry_deleter) {
        if (!ptr_) {
            throw std::runtime_error("Failed to create Lorenz63 geometry");
        }
    }

    /**
     * @brief Check if a state is within the geometry boundaries
     * 
     * @param state The state to check
     * @return bool True if the state is within boundaries, false otherwise
     */
    bool containsPoint(const State& state) const {
        return geometry_contains_point(ptr_.get(), state.getPtr()) != 0;
    }

    /**
     * @brief Get the range of x-axis
     * 
     * @param min_val Output minimum x value
     * @param max_val Output maximum x value
     */
    void getXRange(float& min_val, float& max_val) const {
        geometry_get_x_range(ptr_.get(), &min_val, &max_val);
    }

    /**
     * @brief Get the range of y-axis
     * 
     * @param min_val Output minimum y value
     * @param max_val Output maximum y value
     */
    void getYRange(float& min_val, float& max_val) const {
        geometry_get_y_range(ptr_.get(), &min_val, &max_val);
    }

    /**
     * @brief Get the range of z-axis
     * 
     * @param min_val Output minimum z value
     * @param max_val Output maximum z value
     */
    void getZRange(float& min_val, float& max_val) const {
        geometry_get_z_range(ptr_.get(), &min_val, &max_val);
    }

    /**
     * @brief Get the underlying Fortran pointer (for internal use)
     * 
     * @return void* Raw pointer to Fortran geometry
     */
    void* getPtr() const {
        return ptr_.get();
    }

    /**
     * @brief Friend declaration for the stream output operator
     */
    friend std::ostream& operator<<(std::ostream& os, const Geometry& geometry);

private:
    // Custom deleter for Fortran geometry
    static void geometry_deleter(void* ptr) {
        if (ptr) {
            geometry_destroy(ptr);
        }
    }

    // Smart pointer to manage Fortran geometry
    std::shared_ptr<void> ptr_;
};

/**
 * @brief Stream output operator for Geometry
 * 
 * @param os Output stream
 * @param geometry Geometry to print
 * @return std::ostream& Reference to the output stream
 */
inline std::ostream& operator<<(std::ostream& os, const Geometry& geometry) {
    float x_min, x_max, y_min, y_max, z_min, z_max;
    geometry.getXRange(x_min, x_max);
    geometry.getYRange(y_min, y_max);
    geometry.getZRange(z_min, z_max);
    
    os << "Geometry(X:[" << x_min << ", " << x_max << "], "
       << "Y:[" << y_min << ", " << y_max << "], "
       << "Z:[" << z_min << ", " << z_max << "])";
    return os;
}

} // namespace metada::backends::lorenz63 