/**
 * @file IState.hpp
 * @brief Interface defining the contract for state implementations
 * @ingroup repr
 *
 * This header provides the abstract interface that all state implementations
 * must follow. It defines a unified API for handling state variables in
 * scientific models, including:
 * - Initialization and validation
 * - Data access and manipulation
 * - Metadata management
 * - State information queries
 */

#ifndef METADA_SRC_FRAMEWORK_REPR_ISTATE_HPP_
#define METADA_SRC_FRAMEWORK_REPR_ISTATE_HPP_

#include <string>
#include <vector>

#include "IConfig.hpp"

using metada::framework::tools::config::IConfig;

namespace metada {
namespace framework {
namespace repr {

/**
 * @brief Abstract interface for state implementations
 *
 * This interface defines the contract that all state implementations must
 * follow. It provides a unified API for handling state variables in scientific
 * models.
 *
 * Key features:
 * - State initialization from configuration
 * - Core state operations (reset, validate)
 * - Copy/move semantics
 * - Raw data access
 * - Metadata management
 * - State information queries
 */
class IState {
 public:
  /** @brief Virtual destructor for proper cleanup */
  virtual ~IState() = default;

  /**
   * @brief Initialize state from configuration
   * @param config Configuration object containing initialization parameters
   * @throws std::runtime_error If initialization fails
   */
  virtual void initialize(const IConfig& config) = 0;

  /**
   * @brief Reset the state to initial values
   * @throws std::runtime_error If reset operation fails
   */
  virtual void reset() = 0;

  /**
   * @brief Validate the state consistency
   * @throws std::runtime_error If validation fails
   */
  virtual void validate() const = 0;

  /**
   * @brief Copy state data from another instance
   * @param other Source state to copy from
   * @throws std::runtime_error If copy operation fails
   */
  virtual void copyFrom(const IState& other) = 0;

  /**
   * @brief Move state data from another instance
   * @param other Source state to move from
   * @throws std::runtime_error If move operation fails
   */
  virtual void moveFrom(IState&& other) = 0;

  /**
   * @brief Compare equality with another state
   * @param other State to compare with
   * @return true if states are equal, false otherwise
   */
  virtual bool equals(const IState& other) const = 0;

  /**
   * @brief Get raw pointer to underlying data
   * @return Void pointer to data
   * @throws std::runtime_error If data access fails
   */
  virtual void* getData() = 0;

  /**
   * @brief Get const raw pointer to underlying data
   * @return Const void pointer to data
   * @throws std::runtime_error If data access fails
   */
  virtual const void* getData() const = 0;

  /**
   * @brief Set metadata value for given key
   * @param key Metadata key to set
   * @param value Metadata value to associate with key
   * @throws std::runtime_error If metadata operation fails
   */
  virtual void setMetadata(const std::string& key,
                           const std::string& value) = 0;

  /**
   * @brief Get metadata value for given key
   * @param key Metadata key to retrieve
   * @return Metadata value associated with key
   * @throws std::out_of_range If key does not exist
   */
  virtual std::string getMetadata(const std::string& key) const = 0;

  /**
   * @brief Get names of state variables
   * @return Const reference to vector containing variable names
   */
  virtual const std::vector<std::string>& getVariableNames() const = 0;

  /**
   * @brief Get dimensions of state space
   * @return Const reference to vector containing dimension sizes
   */
  virtual const std::vector<size_t>& getDimensions() const = 0;
};

}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_SRC_FRAMEWORK_REPR_ISTATE_HPP_