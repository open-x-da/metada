/**
 * @file MockBackgroundErrorCovariance.hpp
 * @brief Mock background error covariance implementation for testing
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This file provides a mock implementation of background error covariance
 * operations for testing purposes. The mock implementation uses Google Mock
 * to provide test-friendly interfaces for setting expectations and verifying
 * interactions in unit and integration tests.
 */

#pragma once

#include <gmock/gmock.h>
#include <memory>

// Forward declarations
class MockConfig;
class MockGeometry;
template<typename ConfigBackend, typename GeometryBackend>
class MockState;

namespace metada::backends::gmock {

/**
 * @brief Mock background error covariance backend implementation
 * 
 * @details This mock class provides a test-friendly interface for background
 * error covariance operations. It can be configured to return predefined
 * values or simulate specific behaviors for testing scenarios.
 */
class MockBackgroundErrorCovariance {
 public:
  /** @brief Default constructor */
  MockBackgroundErrorCovariance() = default;

  /** @brief Constructor with config backend */
  template<typename ConfigBackend>
  explicit MockBackgroundErrorCovariance([[maybe_unused]] const ConfigBackend& config_backend) {
    // Mock constructor - no actual initialization needed for testing
  }

  /** @brief Virtual destructor */
  virtual ~MockBackgroundErrorCovariance() = default;

  /** @brief Mock method for computing quadratic form */
  template<typename StateType>
  double computeQuadraticFormDiagonal([[maybe_unused]] const StateType& state) const {
    return 1.0; // Mock return value
  }

  /** @brief Mock method for applying inverse diagonal */
  template<typename StateType>
  void applyInverseDiagonal([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for applying diagonal */
  template<typename StateType>
  void applyDiagonal([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for applying square root diagonal */
  template<typename StateType>
  void applySquareRootDiagonal([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for computing ensemble quadratic form */
  template<typename StateType>
  double computeQuadraticFormEnsemble([[maybe_unused]] const StateType& state) const {
    return 1.0; // Mock return value
  }

  /** @brief Mock method for applying inverse ensemble */
  template<typename StateType>
  void applyInverseEnsemble([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for applying ensemble */
  template<typename StateType>
  void applyEnsemble([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for applying square root ensemble */
  template<typename StateType>
  void applySquareRootEnsemble([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for computing parametric quadratic form */
  template<typename StateType>
  double computeQuadraticFormParametric([[maybe_unused]] const StateType& state) const {
    return 1.0; // Mock return value
  }

  /** @brief Mock method for applying inverse parametric */
  template<typename StateType>
  void applyInverseParametric([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for applying parametric */
  template<typename StateType>
  void applyParametric([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for applying square root parametric */
  template<typename StateType>
  void applySquareRootParametric([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for computing hybrid quadratic form */
  template<typename StateType>
  double computeQuadraticFormHybrid([[maybe_unused]] const StateType& state) const {
    return 1.0; // Mock return value
  }

  /** @brief Mock method for applying inverse hybrid */
  template<typename StateType>
  void applyInverseHybrid([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for applying hybrid */
  template<typename StateType>
  void applyHybrid([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for applying square root hybrid */
  template<typename StateType>
  void applySquareRootHybrid([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for computing full quadratic form */
  template<typename StateType>
  double computeQuadraticFormFull([[maybe_unused]] const StateType& state) const {
    return 1.0; // Mock return value
  }

  /** @brief Mock method for applying inverse full */
  template<typename StateType>
  void applyInverseFull([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for applying full */
  template<typename StateType>
  void applyFull([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }

  /** @brief Mock method for applying square root full */
  template<typename StateType>
  void applySquareRootFull([[maybe_unused]] const StateType& state, [[maybe_unused]] StateType& result) const {
    // Mock implementation
  }
};

} // namespace metada::backends::gmock 