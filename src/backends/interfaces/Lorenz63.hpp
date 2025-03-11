#ifndef LORENZ63_HPP
#define LORENZ63_HPP

#include <memory>
#include <vector>
#include <stdexcept>
#include "lorenz63_fortran.h"

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
 * @brief C++ wrapper for the Lorenz63 Fortran Model class
 */
class Lorenz63Model {
public:
    /**
     * @brief Construct a new Lorenz63Model with given parameters
     * 
     * @param sigma Prandtl number (default 10.0)
     * @param rho Rayleigh number (default 28.0)
     * @param beta Geometric factor (default 8/3)
     * @param dt Time step (default 0.01)
     */
    Lorenz63Model(float sigma = 10.0f, float rho = 28.0f, float beta = 8.0f/3.0f, float dt = 0.01f)
        : ptr_(lorenz63_model_create(sigma, rho, beta, dt), model_deleter) {
        if (!ptr_) {
            throw std::runtime_error("Failed to create Lorenz63 model");
        }
    }

    /**
     * @brief Get the parameters of the model
     * 
     * @param sigma Output Prandtl number
     * @param rho Output Rayleigh number
     * @param beta Output Geometric factor
     * @param dt Output Time step
     */
    void getParameters(float& sigma, float& rho, float& beta, float& dt) const {
        lorenz63_model_get_params(ptr_.get(), &sigma, &rho, &beta, &dt);
    }

    /**
     * @brief Set the parameters of the model
     * 
     * @param sigma Prandtl number
     * @param rho Rayleigh number
     * @param beta Geometric factor
     * @param dt Time step
     */
    void setParameters(float sigma, float rho, float beta, float dt) {
        lorenz63_model_set_params(ptr_.get(), sigma, rho, beta, dt);
    }

    /**
     * @brief Get the underlying Fortran pointer (for internal use)
     * 
     * @return void* Raw pointer to Fortran model
     */
    void* getPtr() const {
        return ptr_.get();
    }

private:
    // Custom deleter for Fortran model
    static void model_deleter(void* ptr) {
        if (ptr) {
            lorenz63_model_destroy(ptr);
        }
    }

    // Smart pointer to manage Fortran model
    std::shared_ptr<void> ptr_;
};

/**
 * @brief C++ wrapper for the Lorenz63 Fortran RK4 Integrator class
 */
class Lorenz63Integrator {
public:
    /**
     * @brief Construct a new Lorenz63Integrator with given model
     * 
     * @param model Reference to Lorenz63Model
     */
    Lorenz63Integrator(const Lorenz63Model& model)
        : ptr_(rk4_integrator_create(model.getPtr()), integrator_deleter) {
        if (!ptr_) {
            throw std::runtime_error("Failed to create RK4 integrator");
        }
    }

    /**
     * @brief Perform one step of integration
     * 
     * @param state Lorenz63State to update
     */
    void step(Lorenz63State& state) {
        rk4_integrator_step(ptr_.get(), state.getPtr());
    }

    /**
     * @brief Run the integration for multiple steps
     * 
     * @param state Lorenz63State to update
     * @param numSteps Number of steps to run
     */
    void run(Lorenz63State& state, int numSteps) {
        for (int i = 0; i < numSteps; ++i) {
            step(state);
        }
    }

private:
    // Custom deleter for Fortran integrator
    static void integrator_deleter(void* ptr) {
        if (ptr) {
            rk4_integrator_destroy(ptr);
        }
    }

    // Smart pointer to manage Fortran integrator
    std::shared_ptr<void> ptr_;
};

} // namespace interfaces
} // namespace backends
} // namespace metada

#endif // LORENZ63_HPP 