#pragma once

namespace metada::framework {

/**
 * @brief Interface for observation implementations in data assimilation systems
 *
 * @details
 * This interface defines the contract for handling observational data in data
 * assimilation systems. Observations represent measurements of the real system
 * and are a fundamental component in data assimilation algorithms.
 *
 * The interface is designed to support:
 * - Observation initialization and quality control
 * - Core mathematical operations (add, subtract, multiply)
 * - Equality comparison
 *
 * Implementations should ensure:
 * - Thread safety for all operations
 * - Proper exception handling and error reporting
 * - Efficient memory management
 * - Type safety for data access
 * - Consistent behavior across operations
 */
class IObservation {
 public:
  /**
   * @brief Default constructor
   * @details Constructs an empty observation
   */
  IObservation() = default;

  /**
   * @brief Virtual destructor for proper cleanup
   * @details Ensures proper cleanup of derived class resources
   */
  virtual ~IObservation() = default;

  /**
   * @brief Copy constructor
   * @details Allows copying of IObservation-derived objects
   */
  IObservation(const IObservation&) = delete;

  /**
   * @brief Copy assignment operator
   * @details Allows copy assignment of IObservation-derived objects
   */
  IObservation& operator=(const IObservation&) = delete;

  /**
   * @brief Move constructor
   * @details Allows moving of IObservation-derived objects
   */
  IObservation(IObservation&&) = delete;

  /**
   * @brief Move assignment operator
   * @details Allows move assignment of IObservation-derived objects
   */
  IObservation& operator=(IObservation&&) = delete;

  /**
   * @brief Initialize the observation
   * @throws std::runtime_error If initialization fails
   */
  virtual void initialize() = 0;

  /**
   * @brief Apply quality control to the observation
   */
  virtual void applyQC() = 0;

  /**
   * @brief Check if this observation equals another
   * @param other Observation to compare with
   * @return true if equal, false otherwise
   */
  virtual bool equals(const IObservation& other) const = 0;

  /**
   * @brief Add another observation to this one
   * @param other Observation to add
   */
  virtual void add(const IObservation& other) = 0;

  /**
   * @brief Subtract another observation from this one
   * @param other Observation to subtract
   */
  virtual void subtract(const IObservation& other) = 0;

  /**
   * @brief Multiply this observation by a scalar
   * @param scalar Scalar value
   */
  virtual void multiply(double scalar) = 0;
};

}  // namespace metada::framework