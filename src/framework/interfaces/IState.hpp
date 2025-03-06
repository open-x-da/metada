/**
 * @file IState.hpp
 * @brief Interface defining the contract for state implementations
 * @ingroup repr
 *
 * @details
 * This header provides the abstract interface that all state implementations
 * must follow. It defines a unified API for handling state variables in
 * scientific models.
 *
 * The interface includes:
 * - Initialization and validation of state
 * - Data access and manipulation through type-safe methods
 * - Metadata management for storing state properties
 * - State information queries for dimensions and variables
 * - Copy/move semantics for state transfer
 *
 * Implementations must provide thread-safe access to state data and proper
 * exception handling for all operations.
 */

#pragma once

#include <string>
#include <vector>

namespace metada::framework {

class IConfig;  // Forward declaration

/**
 * @brief Abstract interface for state implementations
 *
 * @details
 * This interface defines the contract that all state implementations must
 * follow. It provides a unified API for handling state variables in scientific
 * models.
 *
 * The interface is designed to support:
 * - Flexible state initialization from configuration
 * - Core state operations (reset, validate) with proper error handling
 * - Copy/move semantics for efficient state transfer
 * - Type-safe raw data access with const correctness
 * - Extensible metadata management through key-value pairs
 * - Comprehensive state information queries
 *
 * Implementations should ensure:
 * - Thread safety for all operations
 * - Proper exception handling and error reporting
 * - Efficient memory management
 * - Type safety for data access
 * - Consistent state validation
 */
class IState {
 protected:
  /**
   * @brief Default constructor is protected to enforce initialization with a
   * config
   * @details Only derived classes can use this constructor. States should
   * always be initialized with a configuration object.
   */
  IState() = default;

 public:
  /**
   * @brief Virtual destructor for proper cleanup
   * @details Ensures proper cleanup of derived class resources
   */
  virtual ~IState() = default;

  /**
   * @brief Default copy constructor
   * @details Allows copying of IState-derived objects
   */
  IState(const IState&) = default;

  /**
   * @brief Default copy assignment operator
   * @details Allows copy assignment of IState-derived objects
   */
  IState& operator=(const IState&) = default;

  /**
   * @brief Default move constructor
   * @details Allows moving of IState-derived objects
   */
  IState(IState&&) = default;

  /**
   * @brief Default move assignment operator
   * @details Allows move assignment of IState-derived objects
   */
  IState& operator=(IState&&) = default;

  /**
   * @brief Initialize state from configuration
   * @param config Configuration object containing initialization parameters
   * @throws std::runtime_error If initialization fails
   * @details Implementations should validate configuration parameters and
   * initialize internal state accordingly
   */
  virtual void initialize(const IConfig& config) = 0;

  /**
   * @brief Reset the state to initial values
   * @throws std::runtime_error If reset operation fails
   * @details Restores state to initial values defined during initialization
   */
  virtual void reset() = 0;

  /**
   * @brief Set all values to zero
   * @return Reference to this state
   */
  virtual void zero() = 0;

  /**
   * @brief Validate the state consistency
   * @throws std::runtime_error If validation fails
   * @details Verifies internal state consistency and validity of all values
   */
  virtual void validate() const = 0;

  /**
   * @brief Copy state data from another instance
   * @param other Source state to copy from
   * @throws std::runtime_error If copy operation fails
   * @details Performs deep copy of state data and metadata
   */
  virtual void copyFrom(const IState& other) = 0;

  /**
   * @brief Move state data from another instance
   * @param other Source state to move from
   * @throws std::runtime_error If move operation fails
   * @details Transfers ownership of state data efficiently
   */
  virtual void moveFrom(IState&& other) = 0;

  /**
   * @brief Compare equality with another state
   * @param other State to compare with
   * @return true if states are equal, false otherwise
   * @details Performs deep comparison of state data and metadata
   */
  virtual bool equals(const IState& other) const = 0;

  /**
   * @brief Get raw pointer to underlying data
   * @return Void pointer to data
   * @throws std::runtime_error If data access fails
   * @details Provides direct access to underlying state data
   */
  virtual void* getData() = 0;

  /**
   * @brief Get const raw pointer to underlying data
   * @return Const void pointer to data
   * @throws std::runtime_error If data access fails
   * @details Provides read-only access to underlying state data
   */
  virtual const void* getData() const = 0;

  /**
   * @brief Set metadata value for given key
   * @param key Metadata key to set
   * @param value Metadata value to associate with key
   * @throws std::runtime_error If metadata operation fails
   * @details Associates metadata value with specified key
   */
  virtual void setMetadata(const std::string& key,
                           const std::string& value) = 0;

  /**
   * @brief Get metadata value for given key
   * @param key Metadata key to retrieve
   * @return Metadata value associated with key
   * @throws std::out_of_range If key does not exist
   * @details Retrieves metadata value for specified key
   */
  virtual std::string getMetadata(const std::string& key) const = 0;

  /**
   * @brief Get names of state variables
   * @return Const reference to vector containing variable names
   * @details Returns ordered list of state variable identifiers
   */
  virtual const std::vector<std::string>& getVariableNames() const = 0;

  /**
   * @brief Check if the state contains a specific variable
   * @param name Name of the variable to check
   * @return true if the variable exists in the state, false otherwise
   * @details Verifies if a named variable is present in the state
   */
  virtual bool hasVariable(const std::string& name) const = 0;

  /**
   * @brief Get dimensions of state space
   * @return Const reference to vector containing dimension sizes
   * @details Returns ordered list of dimension sizes for state space
   */
  virtual const std::vector<size_t>& getDimensions() const = 0;

  /**
   * @brief Add another state to this one
   * @param other State to add
   * @throws std::runtime_error If states are incompatible
   */
  virtual void add(const IState& other) = 0;

  /**
   * @brief Subtract another state from this one
   * @param other State to subtract
   * @throws std::runtime_error If states are incompatible
   */
  virtual void subtract(const IState& other) = 0;

  /**
   * @brief Multiply this state by a scalar
   * @param scalar Value to multiply by
   * @throws std::runtime_error If multiplication fails
   */
  virtual void multiply(double scalar) = 0;
};

}  // namespace metada::framework