#pragma once

#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "IState.hpp"
#include "state_c_api.h"

namespace metada::backends::lorenz63 {

/**
 * @brief C++ wrapper for the Lorenz63 Fortran State class that implements
 * IState interface
 */
class State : public framework::IState {
 public:
  /**
   * @brief Default constructor required by IState
   */
  State() : ptr_(state_create(0.0f, 0.0f, 0.0f), state_deleter) {
    if (!ptr_) {
      throw std::runtime_error("Failed to create Lorenz63 state");
    }
    initializeVariableNames();
  }

  /**
   * @brief Construct a new Lorenz63State with given initial values from config
   *
   * @param config Configuration object with initial values
   */
  template <typename T>
  State(const metada::framework::Config<T>& config)
      : ptr_(state_create(
                 std::get<float>(config.Get("model.initial_conditions.x")),
                 std::get<float>(config.Get("model.initial_conditions.y")),
                 std::get<float>(config.Get("model.initial_conditions.z"))),
             state_deleter) {
    if (!ptr_) {
      throw std::runtime_error("Failed to create Lorenz63 state");
    }
    initializeVariableNames();
  }

  /**
   * @brief Get the x component
   *
   * @return float x value
   */
  float getX() const { return state_get_x(ptr_.get()); }

  /**
   * @brief Get the y component
   *
   * @return float y value
   */
  float getY() const { return state_get_y(ptr_.get()); }

  /**
   * @brief Get the z component
   *
   * @return float z value
   */
  float getZ() const { return state_get_z(ptr_.get()); }

  /**
   * @brief Set the x component
   *
   * @param x New x value
   */
  void setX(float x) { state_set_x(ptr_.get(), x); }

  /**
   * @brief Set the y component
   *
   * @param y New y value
   */
  void setY(float y) { state_set_y(ptr_.get(), y); }

  /**
   * @brief Set the z component
   *
   * @param z New z value
   */
  void setZ(float z) { state_set_z(ptr_.get(), z); }

  /**
   * @brief Get all components as a vector
   *
   * @return std::vector<float> {x, y, z}
   */
  std::vector<float> getComponents() const { return {getX(), getY(), getZ()}; }

  /**
   * @brief Calculate distance between this state and another
   *
   * @param other Another Lorenz63State
   * @return float Euclidean distance
   */
  float distance(const State& other) const {
    return state_distance(ptr_.get(), other.ptr_.get());
  }

  /**
   * @brief Get the underlying Fortran pointer (for internal use)
   *
   * @return void* Raw pointer to Fortran state
   */
  void* getPtr() const { return ptr_.get(); }

  // IState interface implementation

  /**
   * @brief Initialize state from configuration
   * @param config Configuration object containing initialization parameters
   * @throws std::runtime_error If initialization fails
   */
  void initialize(
      [[maybe_unused]] const metada::framework::IConfig& config) override {
    // Extract x, y, z from config and initialize
    // This is a placeholder implementation; actual config parsing would depend
    // on IConfig's API For now, just reset to default values
    reset();
  }

  /**
   * @brief Reset the state to initial values (0,0,0)
   * @throws std::runtime_error If reset operation fails
   */
  void reset() override {
    setX(0.0f);
    setY(0.0f);
    setZ(0.0f);
  }

  /**
   * @brief Set all values to zero
   */
  void zero() override {
    setX(0.0f);
    setY(0.0f);
    setZ(0.0f);
  }

  /**
   * @brief Validate the state consistency
   * @throws std::runtime_error If validation fails
   */
  void validate() const override {
    // For Lorenz63, all real values are valid
    // Could implement checks for NaN or Inf if needed
    if (std::isnan(getX()) || std::isnan(getY()) || std::isnan(getZ()) ||
        std::isinf(getX()) || std::isinf(getY()) || std::isinf(getZ())) {
      throw std::runtime_error("Invalid state: contains NaN or Inf values");
    }
  }

  /**
   * @brief Copy state data from another instance
   * @param other Source state to copy from
   * @throws std::runtime_error If copy operation fails
   */
  void copyFrom(const IState& other) override {
    try {
      const State& otherState = dynamic_cast<const State&>(other);
      setX(otherState.getX());
      setY(otherState.getY());
      setZ(otherState.getZ());

      // Copy metadata
      metadata_ = otherState.metadata_;
    } catch (const std::bad_cast&) {
      throw std::runtime_error("Cannot copy from incompatible state type");
    }
  }

  /**
   * @brief Move state data from another instance
   * @param other Source state to move from
   * @throws std::runtime_error If move operation fails
   */
  void moveFrom(IState&& other) override {
    try {
      State&& otherState = dynamic_cast<State&&>(other);
      setX(otherState.getX());
      setY(otherState.getY());
      setZ(otherState.getZ());

      // Move metadata
      metadata_ = std::move(otherState.metadata_);
    } catch (const std::bad_cast&) {
      throw std::runtime_error("Cannot move from incompatible state type");
    }
  }

  /**
   * @brief Compare equality with another state
   * @param other State to compare with
   * @return true if states are equal, false otherwise
   */
  bool equals(const IState& other) const override {
    try {
      const State& otherState = dynamic_cast<const State&>(other);
      const float EPSILON = 1e-6f;
      return std::abs(getX() - otherState.getX()) < EPSILON &&
             std::abs(getY() - otherState.getY()) < EPSILON &&
             std::abs(getZ() - otherState.getZ()) < EPSILON;
    } catch (const std::bad_cast&) {
      return false;
    }
  }

  /**
   * @brief Get raw pointer to underlying data
   * @return Void pointer to data
   */
  void* getData() override { return ptr_.get(); }

  /**
   * @brief Get const raw pointer to underlying data
   * @return Const void pointer to data
   */
  const void* getData() const override { return ptr_.get(); }

  /**
   * @brief Set metadata value for given key
   * @param key Metadata key to set
   * @param value Metadata value to associate with key
   */
  void setMetadata(const std::string& key, const std::string& value) override {
    metadata_[key] = value;
  }

  /**
   * @brief Get metadata value for given key
   * @param key Metadata key to retrieve
   * @return Metadata value associated with key
   * @throws std::out_of_range If key does not exist
   */
  std::string getMetadata(const std::string& key) const override {
    auto it = metadata_.find(key);
    if (it == metadata_.end()) {
      throw std::out_of_range("Metadata key not found: " + key);
    }
    return it->second;
  }

  /**
   * @brief Get names of state variables
   * @return Const reference to vector containing variable names
   */
  const std::vector<std::string>& getVariableNames() const override {
    return variableNames_;
  }

  /**
   * @brief Check if the state contains a specific variable
   * @param name Name of the variable to check
   * @return true if the variable exists in the state, false otherwise
   */
  bool hasVariable(const std::string& name) const override {
    for (const auto& varName : variableNames_) {
      if (varName == name) {
        return true;
      }
    }
    return false;
  }

  /**
   * @brief Get dimensions of state space
   * @return Const reference to vector containing dimension sizes
   */
  const std::vector<size_t>& getDimensions() const override {
    return dimensions_;
  }

  /**
   * @brief Add another state to this one
   * @param other State to add
   * @throws std::runtime_error If states are incompatible
   */
  void add(const IState& other) override {
    try {
      const State& otherState = dynamic_cast<const State&>(other);
      setX(getX() + otherState.getX());
      setY(getY() + otherState.getY());
      setZ(getZ() + otherState.getZ());
    } catch (const std::bad_cast&) {
      throw std::runtime_error("Cannot add incompatible state type");
    }
  }

  /**
   * @brief Subtract another state from this one
   * @param other State to subtract
   * @throws std::runtime_error If states are incompatible
   */
  void subtract(const IState& other) override {
    try {
      const State& otherState = dynamic_cast<const State&>(other);
      setX(getX() - otherState.getX());
      setY(getY() - otherState.getY());
      setZ(getZ() - otherState.getZ());
    } catch (const std::bad_cast&) {
      throw std::runtime_error("Cannot subtract incompatible state type");
    }
  }

  /**
   * @brief Multiply this state by a scalar
   * @param scalar Value to multiply by
   */
  void multiply(double scalar) override {
    setX(static_cast<float>(getX() * scalar));
    setY(static_cast<float>(getY() * scalar));
    setZ(static_cast<float>(getZ() * scalar));
  }

  /**
   * @brief Calculate the dot product of this state with another state
   * @param other State to calculate dot product with
   * @return Resulting dot product value
   * @throws std::runtime_error If dot product operation fails
   */
  double dot(const IState& other) const override {
    try {
      const State& otherState = dynamic_cast<const State&>(other);
      return getX() * otherState.getX() + getY() * otherState.getY() +
             getZ() * otherState.getZ();
    } catch (const std::bad_cast&) {
      throw std::runtime_error(
          "Cannot calculate dot product with incompatible state type");
    }
  }

  /**
   * @brief Calculate the norm of this state
   * @return Resulting norm value
   */
  double norm() const override {
    return std::sqrt(getX() * getX() + getY() * getY() + getZ() * getZ());
  }

  /**
   * @brief Friend declaration for the stream output operator
   */
  friend std::ostream& operator<<(std::ostream& os, const State& state);

 private:
  /**
   * @brief Initialize variable names for the Lorenz63 model
   */
  void initializeVariableNames() {
    variableNames_ = {"x", "y", "z"};
    dimensions_ = {3};  // Single dimension of size 3
  }

  // Custom deleter for Fortran state
  static void state_deleter(void* ptr) {
    if (ptr) {
      state_destroy(ptr);
    }
  }

  // Smart pointer to manage Fortran state
  std::shared_ptr<void> ptr_;

  // Store variable names and dimensions
  std::vector<std::string> variableNames_;
  std::vector<size_t> dimensions_;

  // Metadata storage
  std::unordered_map<std::string, std::string> metadata_;
};

/**
 * @brief Stream output operator for State
 *
 * @param os Output stream
 * @param state State to print
 * @return std::ostream& Reference to the output stream
 */
inline std::ostream& operator<<(std::ostream& os, const State& state) {
  os << "State(" << std::fixed << std::setprecision(6) << state.getX() << ", "
     << std::fixed << std::setprecision(6) << state.getY() << ", " << std::fixed
     << std::setprecision(6) << state.getZ() << ")";
  return os;
}

}  // namespace metada::backends::lorenz63