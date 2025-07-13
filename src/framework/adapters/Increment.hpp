#pragma once

#include <string>

namespace metada::framework {

/**
 * @class Increment
 * @brief Represents perturbations or differences between entities in data
 * assimilation systems
 *
 * @details
 * The Increment class template provides a representation of differences between
 * entities, which is a core concept in data assimilation. Common use cases
 * include:
 * - Analysis increments (x_a - x_b)
 * - Observation innovations (y - H(x))
 * - Ensemble perturbations
 * - Iterative minimization steps
 *
 * @tparam EntityType The type of entity this increment operates on. Must
 * support:
 * - Arithmetic operators (+=, -=, *=)
 * - zero() method for resetting values
 * - norm() method for computing L2 norm
 * - dot() method for inner products
 * - getData<T>() method for data access
 * - Metadata methods (setMetadata/getMetadata)
 * - size() for total size information
 * - isInitialized() for validity checks
 */
template <typename EntityType>
class Increment {
 private:
  EntityType entity_;        ///< The underlying entity
  bool initialized_{false};  ///< Initialization flag

  /**
   * @brief Private constructor taking an entity
   * @param entity The entity to initialize with
   */
  explicit Increment(const EntityType& entity)
      : entity_(entity.clone()), initialized_(true) {}

 public:
  /** @brief Default constructor */
  Increment() = default;

  /**
   * @brief Constructs increment as difference between two entities
   * @param entity1 First entity
   * @param entity2 Second entity
   * @details Creates increment as (entity1 - entity2)
   */
  Increment(const EntityType& entity1, const EntityType& entity2)
      : entity_(entity1.clone()), initialized_(true) {
    entity_ -= entity2;
  }

  // Core accessors
  /**
   * @brief Get access to the underlying entity
   * @return Reference to the entity
   */
  EntityType& entity() { return entity_; }

  /**
   * @brief Get const access to the underlying entity
   * @return Const reference to the entity
   */
  const EntityType& entity() const { return entity_; }

  /**
   * @brief Gets typed access to underlying data
   * @tparam T Type to cast data to
   * @return Reference to data as type T
   */
  template <typename T>
  T& getData() {
    return entity_.template getData<T>();
  }

  /**
   * @brief Gets const typed access to underlying data
   * @tparam T Type to cast data to
   * @return Const reference to data as type T
   */
  template <typename T>
  const T& getData() const {
    return entity_.template getData<T>();
  }

  /**
   * @brief Gets total size of increment
   * @return Total number of elements
   */
  size_t size() const { return entity_.size(); }

  /**
   * @brief Checks if increment is initialized
   * @return true if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_ && entity_.isInitialized(); }

  // Metadata operations
  /**
   * @brief Sets metadata value
   * @param key Metadata key
   * @param value Metadata value
   */
  void setMetadata(const std::string& key, const std::string& value) {
    entity_.setMetadata(key, value);
  }

  /**
   * @brief Gets metadata value
   * @param key Metadata key
   * @return Metadata value
   */
  std::string getMetadata(const std::string& key) const {
    return entity_.getMetadata(key);
  }

  // Core mathematical operations
  /**
   * @brief Sets all elements to zero
   * @return Reference to this increment
   */
  Increment& zero() {
    entity_.zero();
    return *this;
  }

  /**
   * @brief Scales increment by a factor
   * @param alpha Scaling factor
   * @return Reference to this increment
   */
  Increment& scale(double alpha) {
    entity_ *= alpha;
    return *this;
  }

  /**
   * @brief Performs AXPY operation: this += alpha * other
   * @param alpha Scaling factor
   * @param other Increment to add
   * @return Reference to this increment
   */
  Increment& axpy(double alpha, const Increment& other) {
    auto scaled = other.entity_.clone();
    scaled *= alpha;
    entity_ += scaled;
    return *this;
  }

  /**
   * @brief Fills the increment with random values in [-0.5, 0.5]
   */
  void randomize() {
    auto& data = getData<std::vector<double>>();
    for (auto& v : data) {
      v = (double(rand()) / RAND_MAX - 0.5);
    }
  }

  /**
   * @brief Computes dot product with another increment
   * @param other Increment to compute dot product with
   * @return Dot product value
   */
  double dot(const Increment& other) const {
    const auto& data1 = getData<std::vector<double>>();
    const auto& data2 = other.getData<std::vector<double>>();
    double result = 0.0;
    for (size_t i = 0; i < data1.size(); ++i) {
      result += data1[i] * data2[i];
    }
    return result;
  }

  /**
   * @brief Computes L2 norm of increment
   * @return L2 norm value
   */
  double norm() const { return entity_.norm(); }

  // Factory methods
  /**
   * @brief Creates increment as copy of entity
   * @param entity Entity to copy
   * @return New increment
   */
  static Increment createFromEntity(const EntityType& entity) {
    return Increment(entity);
  }

  /**
   * @brief Creates increment as difference between entities
   * @param first First entity
   * @param second Second entity
   * @return New increment representing (first - second)
   */
  static Increment createFromDifference(const EntityType& first,
                                        const EntityType& second) {
    Increment result = createFromEntity(first);
    result.entity_ -= second;
    return result;
  }

  /**
   * @brief Applies increment to target entity
   * @tparam TargetType Type of target entity
   * @param target Entity to apply increment to
   * @return Reference to modified target
   */
  template <typename TargetType>
  TargetType& applyTo(TargetType& target) const {
    target += entity_;
    return target;
  }

  // Operator overloads
  /**
   * @brief Adds two increments
   * @param other Increment to add
   * @return New increment containing sum
   */
  Increment operator+(const Increment& other) const {
    Increment result(*this);
    result.entity_ += other.entity_;
    return result;
  }

  /**
   * @brief Multiplies increment by scalar
   * @param scalar Scalar multiplier
   * @return New increment containing product
   */
  Increment operator*(double scalar) const {
    Increment result(*this);
    result.entity_ *= scalar;
    return result;
  }

  /**
   * @brief Adds increment in-place
   * @param other Increment to add
   * @return Reference to this increment
   */
  Increment& operator+=(const Increment& other) {
    entity_ += other.entity_;
    return *this;
  }

  /**
   * @brief Multiplies by scalar in-place
   * @param scalar Scalar multiplier
   * @return Reference to this increment
   */
  Increment& operator*=(double scalar) {
    entity_ *= scalar;
    return *this;
  }

  /**
   * @brief Non-member scalar multiplication
   * @param scalar Scalar multiplier
   * @param increment Increment to multiply
   * @return New increment containing product
   */
  friend Increment operator*(double scalar, const Increment& increment) {
    return increment * scalar;
  }
};
}  // namespace metada::framework
