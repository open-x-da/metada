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

#include "DateTime.hpp"
#include "ObsIOConcepts.hpp"

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of ObsIO for testing
 *
 * @details
 * Provides mock methods for all ObsIO interface operations, organized
 * into the following categories:
 *
 * @par File Operations
 * - readObservations() - Read observations from a file
 * - writeObservations() - Write observations to a file
 * - canRead() - Check if a file can be read
 * - canWrite() - Check if observations can be written
 *
 * @par Format Information
 * - getFormatName() - Get the name of the format
 * - getFileExtensions() - Get the file extensions supported by this backend
 *
 * @note All mock methods use Google Mock's MOCK_METHOD macro to enable
 * setting expectations and verifying calls.
 */
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
  MockObsIO(MockObsIO&& other) noexcept : params_(std::move(other.params_)) {}

  /**
   * @brief Move assignment operator
   *
   * @param other The other MockObsIO object to move from
   * @return Reference to this MockObsIO object
   */
  MockObsIO& operator=(MockObsIO&& other) noexcept {
    if (this != &other) {
      params_ = std::move(other.params_);
    }
    return *this;
  }

  /**
   * @brief Constructor with initialization parameters
   *
   * @param params Initialization parameters as a string
   */
  explicit MockObsIO(const std::string& params) : params_(params) {}

  /**
   * @brief Read observations from a file
   *
   * @param filename Path to the file to read
   * @return Vector of observation records
   */
  MOCK_METHOD(std::vector<framework::ObservationRecord>, read,
              (const std::string& filename));

  /**
   * @brief Check if a file can be read by this backend
   *
   * @param filename Path to the file to check
   * @return True if the file can be read, false otherwise
   */
  MOCK_METHOD(bool, canRead, (const std::string& filename), (const));

  /**
   * @brief Write observations to a file
   *
   * @param filename Path to the file to write
   * @param records Vector of observation records to write
   */
  MOCK_METHOD(void, write,
              (const std::string& filename,
               const std::vector<framework::ObservationRecord>& records));

  /**
   * @brief Check if this backend can write observations
   *
   * @return True if observations can be written, false otherwise
   */
  MOCK_METHOD(bool, canWrite, (), (const));

  /**
   * @brief Get the name of the format supported by this backend
   *
   * @return The format name
   */
  MOCK_METHOD(std::string, getFormatName, (), (const));

  /**
   * @brief Get the file extensions supported by this backend
   *
   * @return Vector of supported file extensions
   */
  MOCK_METHOD(std::vector<std::string>, getFileExtensions, (), (const));

 private:
  std::string params_; /**< Initialization parameters */
};

}  // namespace metada::backends::gmock