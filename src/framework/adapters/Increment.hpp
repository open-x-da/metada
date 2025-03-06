#pragma once

#include <memory>
#include <string>
#include <vector>

namespace metada::framework {

/**
 * @brief Increment class representing perturbations or differences between
 * entities
 *
 * @details
 * This class template provides a representation of an increment as a difference
 * between entities in data assimilation systems. The increment is a core
 * concept in data assimilation, representing perturbations, differences, or
 * innovation vectors.
 *
 * Common use cases include:
 * - Analysis increments (x_a - x_b)
 * - Observation innovations (y - H(x))
 * - Ensemble perturbations
 * - Iterative minimization steps
 *
 * The EntityType template parameter must support the following operations:
 * - += operator for adding another entity
 * - *= operator for scaling by a scalar value
 * - -= operator for computing differences
 * - zero() method for resetting values
 * - norm() method for computing the L2 norm
 * - dot() method for computing inner products
 * - getData<T>() method for accessing the underlying data
 * - setMetadata/getMetadata for managing metadata
 * - getDimensions() for size/shape information
 * - isInitialized() for checking validity
 *
 * @tparam EntityType The type of entity this increment operates on
 */
template <typename EntityType>
class Increment {
 private:
  EntityType entity_;        ///< The underlying entity
  bool initialized_{false};  ///< Initialization flag

  // Constructor that takes an entity
  explicit Increment(const EntityType& entity)
      : entity_(entity), initialized_(true) {}

 public:
  /** @brief Default constructor */
  Increment() = default;

  /**
   * @brief Constructor from two entities (computes the difference)
   *
   * @param entity1 The first entity
   * @param entity2 The second entity
   * @details Creates an increment as the difference between entity1 and entity2
   * (entity1 - entity2)
   */
  Increment(const EntityType& entity1, const EntityType& entity2)
      : entity_(entity1), initialized_(true) {
    entity_ -= entity2;
  }

  /**
   * @brief Get access to the underlying entity
   *
   * @return Reference to the underlying entity
   */
  EntityType& entity() { return entity_; }

  /**
   * @brief Get const access to the underlying entity
   *
   * @return Const reference to the underlying entity
   */
  const EntityType& entity() const { return entity_; }

  /**
   * @brief Set all elements to zero
   *
   * @return Reference to this increment
   */
  Increment& zero() {
    entity_.zero();
    return *this;
  }

  /**
   * @brief Scale the increment by a factor
   *
   * @param alpha The scaling factor
   * @return Reference to this increment
   */
  Increment& scale(double alpha) {
    entity_ *= alpha;
    return *this;
  }

  /**
   * @brief Add another increment scaled by a factor (this += alpha * other)
   *
   * This is the AXPY operation (A times X Plus Y) common in linear algebra and
   * extensively used in iterative minimization algorithms.
   *
   * @param alpha The scaling factor
   * @param other The increment to add
   * @return Reference to this increment
   */
  Increment& axpy(double alpha, const Increment& other) {
    EntityType scaled = other.entity_;
    scaled *= alpha;
    entity_ += scaled;
    return *this;
  }

  /**
   * @brief Compute dot product with another increment
   *
   * The inner product is essential for algorithms like conjugate gradient and
   * for computing norms and projections.
   *
   * @param other The other increment
   * @return The dot product value
   */
  double dot(const Increment& other) const {
    return entity_.dot(other.entity_);
  }

  /**
   * @brief Compute the L2 norm of the increment
   *
   * The norm is used to measure the size of the increment and for convergence
   * criteria in iterative algorithms.
   *
   * @return The L2 norm value
   */
  double norm() const { return entity_.norm(); }

  /**
   * @brief Apply this increment to an entity
   *
   * Adds this increment to the provided entity. This is used, for example,
   * to update the background state with an analysis increment.
   *
   * @param target The entity to apply this increment to
   * @return Reference to the modified entity
   */
  template <typename TargetType>
  TargetType& applyTo(TargetType& target) const {
    target += entity_;
    return target;
  }

  /**
   * @brief Create an increment as a copy of an entity
   */
  static Increment createFromEntity(const EntityType& entity) {
    return Increment(entity);
  }

  /**
   * @brief Create an increment as the difference between entities
   */
  static Increment createFromDifference(const EntityType& first,
                                        const EntityType& second) {
    Increment result = createFromEntity(first);
    result.entity_ -= second;
    return result;
  }

  /**
   * @brief Get typed access to the underlying data
   *
   * @tparam T The type to cast the data to
   * @return Reference to the data with the requested type
   */
  template <typename T>
  T& getData() {
    return entity_.template getData<T>();
  }

  /**
   * @brief Get const typed access to the underlying data
   *
   * @tparam T The type to cast the data to
   * @return Const reference to the data with the requested type
   */
  template <typename T>
  const T& getData() const {
    return entity_.template getData<T>();
  }

  /**
   * @brief Set metadata value
   *
   * @param key The metadata key
   * @param value The metadata value
   */
  void setMetadata(const std::string& key, const std::string& value) {
    entity_.setMetadata(key, value);
  }

  /**
   * @brief Get metadata value
   *
   * @param key The metadata key
   * @return The metadata value
   */
  std::string getMetadata(const std::string& key) const {
    return entity_.getMetadata(key);
  }

  /**
   * @brief Get the dimensions of the increment
   *
   * @return Vector of dimensions
   */
  const std::vector<size_t>& getDimensions() const {
    return entity_.getDimensions();
  }

  /**
   * @brief Check if the increment is initialized
   *
   * @return true if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_ && entity_.isInitialized(); }

  /**
   * @brief Addition operator
   *
   * Creates a new increment by adding another increment to this one.
   * Used in combining increments from different sources.
   *
   * @param other The increment to add
   * @return A new increment containing the sum
   */
  Increment operator+(const Increment& other) const {
    Increment result(*this);
    result.entity_ += other.entity_;
    return result;
  }

  /**
   * @brief Multiplication operator
   *
   * Creates a new increment by scaling this increment by a scalar value.
   *
   * @param scalar The scalar to multiply by
   * @return A new increment containing the product
   */
  Increment operator*(double scalar) const {
    Increment result(*this);
    result.entity_ *= scalar;
    return result;
  }

  /**
   * @brief Addition assignment operator
   *
   * Adds another increment to this one in-place.
   *
   * @param other The increment to add
   * @return Reference to this increment
   */
  Increment& operator+=(const Increment& other) {
    entity_ += other.entity_;
    return *this;
  }

  /**
   * @brief Multiplication assignment operator
   *
   * Scales this increment by a scalar value in-place.
   *
   * @param scalar The scalar to multiply by
   * @return Reference to this increment
   */
  Increment& operator*=(double scalar) {
    entity_ *= scalar;
    return *this;
  }

  /**
   * @brief Non-member multiplication operator (scalar * Increment)
   *
   * Allows scalar multiplication with the scalar on the left side.
   *
   * @param scalar The scalar to multiply by
   * @param increment The increment to multiply
   * @return A new increment containing the product
   */
  friend Increment operator*(double scalar, const Increment& increment) {
    return increment * scalar;
  }
};
}  // namespace metada::framework
