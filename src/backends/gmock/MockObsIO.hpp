/**
 * @file MockObsIO.hpp
 * @brief Mock implementation of ObsIO interface for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the ObsIO interface using Google
 * Mock. It allows testing code that depends on ObsIO by providing mock
 * implementations of all interface methods that can be configured with
 * expectations and behaviors.
 *
 * The mock implementation supports:
 * - Setting expectations on method calls
 * - Configuring return values and behaviors
 * - Verifying interaction patterns
 * - Testing error conditions
 *
 * @see ObsIOConcepts
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <string>
#include <vector>

#include "ObsRecord.hpp"

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of ObsIO for testing
 *
 * @details
 * Provides mock methods for all ObsIO interface operations:
 *
 * - read() - Read observations from the configured data source
 * - write() - Write observations to the configured data destination
 *
 * @note All mock methods use Google Mock's MOCK_METHOD macro to enable
 * setting expectations and verifying calls.
 */
template <typename ConfigBackend>
class MockObsIO {
 public:
  /**
   * @brief Default constructor is deleted
   *
   * MockObsIO objects must be constructed with initialization parameters.
   */
  MockObsIO() = delete;

  /**
   * @brief Destructor
   */
  ~MockObsIO() = default;

  /**
   * @brief Copy constructor is deleted
   */
  MockObsIO(const MockObsIO&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  MockObsIO& operator=(const MockObsIO&) = delete;

  /**
   * @brief Move constructor
   *
   * @param other The other MockObsIO object to move from
   */
  MockObsIO(MockObsIO&& other) noexcept : config_(std::move(other.config_)) {}

  /**
   * @brief Move assignment operator
   *
   * @param other The other MockObsIO object to move from
   * @return Reference to this MockObsIO object
   */
  MockObsIO& operator=([[maybe_unused]] MockObsIO&& other) noexcept {
    if (this != &other) {
      config_ = std::move(other.config_);
    }
    return *this;
  }

  /**
   * @brief Constructor with initialization parameters
   *
   * @param params Initialization parameters as a string
   */
  explicit MockObsIO(ConfigBackend&& config) : config_(std::move(config)) {}

  /**
   * @brief Read observations from the configured data source
   *
   * @return Vector of observation records
   */
  MOCK_METHOD(std::vector<framework::ObsRecord>, read, ());

  /**
   * @brief Write observations to the configured data destination
   *
   * @param records Vector of observation records to write
   */
  MOCK_METHOD(void, write, (const std::vector<framework::ObsRecord>& records));

 private:
  ConfigBackend config_; /**< Initialization parameters */
};

}  // namespace metada::backends::gmock