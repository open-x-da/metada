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
 public:
  /**
   * @brief Default constructor
   * @details Constructs an empty state
   */
  IState() = default;

  /**
   * @brief Virtual destructor for proper cleanup
   * @details Ensures proper cleanup of derived class resources
   */
  virtual ~IState() = default;

  /**
   * @brief Copy constructor
   * @details Allows copying of IState-derived objects
   */
  IState(const IState&) = delete;

  /**
   * @brief Copy assignment operator
   * @details Allows copy assignment of IState-derived objects
   */
  IState& operator=(const IState&) = delete;

  /**
   * @brief Move constructor
   * @details Allows moving of IState-derived objects
   */
  IState(IState&&) = delete;

  /**
   * @brief Move assignment operator
   * @details Allows move assignment of IState-derived objects
   */
  IState& operator=(IState&&) = delete;

  /**
   * @brief Initialize state from configuration
   * @param config Configuration object containing initialization parameters
   * @throws std::runtime_error If initialization fails
   * @details Implementations should validate configuration parameters and
   * initialize internal state accordingly
   */
  virtual void initialize() = 0;

  /**
   * @brief Set all values to zero
   * @return Reference to this state
   */
  virtual void zero() = 0;

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
  virtual const std::vector<size_t>& getDimensions(
      const std::string& name) const = 0;

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

  /**
   * @brief Calculate the dot product of this state with another state
   * @param other State to calculate dot product with
   * @return Resulting dot product value
   * @throws std::runtime_error If dot product operation fails
   */
  virtual double dot(const IState& other) const = 0;

  /**
   * @brief Calculate the norm of this state
   * @return Resulting norm value
   * @throws std::runtime_error If norm operation fails
   */
  virtual double norm() const = 0;
};

}  // namespace metada::framework