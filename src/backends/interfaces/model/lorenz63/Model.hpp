#pragma once

#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>

#include "IModel.hpp"
#include "IState.hpp"
#include "State.hpp"
#include "model_c_api.h"

namespace metada::backends::lorenz63 {

// Forward declaration
class Integrator;

// Helper struct to manage void pointers with shared_ptr
struct VoidPointerWrapper {
  void* ptr;
  explicit VoidPointerWrapper(void* p) : ptr(p) {}
};

/**
 * @brief C++ wrapper for the Lorenz63 Fortran Model class that implements
 * IModel
 */
class Model : public framework::IModel {
 public:
  /**
   * @brief Construct a new Lorenz63Model with configuration
   *
   * @param config Configuration object containing model parameters
   */
  template <typename T>
  explicit Model(const framework::Config<T>& config)
      : framework::IModel(config), initialized_(false), ptr_(nullptr) {
    initializeFromConfig(config);
  }

  /**
   * @brief Initialize the model with the given configuration
   *
   * @param config Configuration object with model parameters
   * @throws std::runtime_error if initialization fails
   */
  void initialize(const framework::IConfig& config) override {
    try {
      // Extract parameters from config
      float sigma = std::get<float>(config.Get("model.parameters.sigma"));
      float rho = std::get<float>(config.Get("model.parameters.rho"));
      float beta = std::get<float>(config.Get("model.parameters.beta"));
      float dt = std::get<float>(config.Get("model.integrator.dt"));

      // Create the model
      void* model_ptr = lorenz63_model_create(sigma, rho, beta, dt);
      if (!model_ptr) {
        throw std::runtime_error("Failed to create Lorenz63 model");
      }

      // Use the wrapper to manage the void pointer
      ptr_ = std::shared_ptr<VoidPointerWrapper>(
          new VoidPointerWrapper(model_ptr), [](VoidPointerWrapper* w) {
            if (w && w->ptr) lorenz63_model_destroy(w->ptr);
            delete w;
          });

      // Store parameters for later access
      parameters_["sigma"] = std::to_string(sigma);
      parameters_["rho"] = std::to_string(rho);
      parameters_["beta"] = std::to_string(beta);
      parameters_["dt"] = std::to_string(dt);
      parameters_["num_steps"] = std::to_string(static_cast<int>(
          std::get<int>(config.Get("model.integrator.num_steps"))));

      initialized_ = true;
    } catch (const std::exception& e) {
      throw std::runtime_error(
          std::string("Failed to initialize Lorenz63 model: ") + e.what());
    }
  }

  /**
   * @brief Initialize the model with the given configuration (templated
   * version)
   *
   * @param config Configuration object with model parameters
   * @throws std::runtime_error if initialization fails
   */
  template <typename T>
  void initializeFromConfig(const framework::Config<T>& config) {
    // Just delegate to the non-template version
    initialize(config);
  }

  /**
   * @brief Reset the model to its initial state
   */
  void reset() override {
    if (!initialized_) {
      throw std::runtime_error("Cannot reset uninitialized model");
    }

    // Reset parameters to their initial values
    float sigma = std::stof(parameters_["sigma"]);
    float rho = std::stof(parameters_["rho"]);
    float beta = std::stof(parameters_["beta"]);
    float dt = std::stof(parameters_["dt"]);

    lorenz63_model_set_params(getPtr(), sigma, rho, beta, dt);
  }

  /**
   * @brief Finalize the model and release resources
   */
  void finalize() override {
    ptr_.reset();
    initialized_ = false;
  }

  /**
   * @brief Check if the model is initialized
   *
   * @return true if the model is initialized, false otherwise
   */
  bool isInitialized() const override { return initialized_; }

  /**
   * @brief Get a parameter from the model
   *
   * @param name The name of the parameter
   * @return The value of the parameter as a string
   * @throws std::out_of_range if the parameter doesn't exist
   */
  std::string getParameter(const std::string& name) const override {
    auto it = parameters_.find(name);
    if (it == parameters_.end()) {
      throw std::out_of_range("Parameter not found: " + name);
    }
    return it->second;
  }

  /**
   * @brief Set a parameter in the model
   *
   * @param name The name of the parameter
   * @param value The value of the parameter
   * @throws std::invalid_argument if the parameter is invalid
   */
  void setParameter(const std::string& name,
                    const std::string& value) override {
    if (!initialized_) {
      throw std::runtime_error("Cannot set parameter on uninitialized model");
    }

    try {
      float floatValue = std::stof(value);

      // Update internal parameter store
      parameters_[name] = value;

      // Update model parameters
      float sigma, rho, beta, dt;
      getParameters(sigma, rho, beta, dt);

      if (name == "sigma")
        sigma = floatValue;
      else if (name == "rho")
        rho = floatValue;
      else if (name == "beta")
        beta = floatValue;
      else if (name == "dt")
        dt = floatValue;
      else if (name != "num_steps") {
        throw std::invalid_argument("Unknown parameter: " + name);
      }

      lorenz63_model_set_params(getPtr(), sigma, rho, beta, dt);
    } catch (const std::invalid_argument&) {
      throw std::invalid_argument("Invalid value for parameter: " + name);
    }
  }

  /**
   * @brief Run the model from initial state to final state
   *
   * @param initialState The initial state of the model
   * @param finalState The final state after model execution
   * @param startTime The start time for the model run (unused)
   * @param endTime The end time for the model run (unused)
   * @throws std::runtime_error if the model run fails
   */
  void run(const framework::IState& initialState, framework::IState& finalState,
           double startTime, double endTime) override;

  /**
   * @brief Get the parameters of the model
   *
   * @param sigma Output Prandtl number
   * @param rho Output Rayleigh number
   * @param beta Output Geometric factor
   * @param dt Output Time step
   */
  void getParameters(float& sigma, float& rho, float& beta, float& dt) const {
    if (!initialized_) {
      throw std::runtime_error("Cannot get parameters of uninitialized model");
    }
    lorenz63_model_get_params(getPtr(), &sigma, &rho, &beta, &dt);
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
    if (!initialized_) {
      throw std::runtime_error("Cannot set parameters of uninitialized model");
    }
    lorenz63_model_set_params(getPtr(), sigma, rho, beta, dt);

    // Update parameter store
    parameters_["sigma"] = std::to_string(sigma);
    parameters_["rho"] = std::to_string(rho);
    parameters_["beta"] = std::to_string(beta);
    parameters_["dt"] = std::to_string(dt);
  }

  /**
   * @brief Get the underlying Fortran pointer (for internal use)
   *
   * @return void* Raw pointer to Fortran model
   */
  void* getPtr() const {
    if (!initialized_) {
      throw std::runtime_error("Cannot get pointer to uninitialized model");
    }
    return ptr_->ptr;
  }

 private:
  // Smart pointer to manage Fortran model using the wrapper
  std::shared_ptr<VoidPointerWrapper> ptr_;

  // Tracking initialization state
  bool initialized_;

  // Store for model parameters
  std::unordered_map<std::string, std::string> parameters_;
};

/**
 * @brief C++ wrapper for the Lorenz63 Fortran RK4 Integrator class
 */
class Integrator {
 public:
  /**
   * @brief Construct a new Integrator with given model
   *
   * @param model Reference to Model
   */
  Integrator(const Model& model)
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
  void step(State& state) { rk4_integrator_step(ptr_.get(), state.getPtr()); }

  /**
   * @brief Run the integration for multiple steps
   *
   * @param state Lorenz63State to update
   * @param numSteps Number of steps to run
   */
  void run(State& state, int numSteps) {
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

// Implementation of Model::run
inline void Model::run(const framework::IState& initialState,
                       [[maybe_unused]] framework::IState& finalState,
                       [[maybe_unused]] double startTime,
                       [[maybe_unused]] double endTime) {
  if (!initialized_) {
    throw std::runtime_error("Cannot run uninitialized model");
  }

  try {
    // Cast to concrete Lorenz63 state types
    const State& lorenzInitialState = dynamic_cast<const State&>(initialState);

    // Copy initial state to final state
    std::unique_ptr<State> lorenzFinalState(
        dynamic_cast<State*>(lorenzInitialState.clone().release()));

    // Create integrator and run simulation
    Integrator integrator(*this);
    int numSteps = std::stoi(parameters_["num_steps"]);
    integrator.run(*lorenzFinalState, numSteps);

  } catch (const std::bad_cast&) {
    throw std::runtime_error(
        "Incompatible state type provided to Lorenz63 model");
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("Model run failed: ") + e.what());
  }
}

}  // namespace metada::backends::lorenz63