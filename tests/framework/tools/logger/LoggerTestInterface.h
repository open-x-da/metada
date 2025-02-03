#ifndef METADA_TESTS_FRAMEWORK_TOOLS_LOGGER_LOGGERTESTINTERFACE_H_
#define METADA_TESTS_FRAMEWORK_TOOLS_LOGGER_LOGGERTESTINTERFACE_H_

#include "Logger.h"

namespace metada {
namespace framework {
namespace tools {
namespace logger {

/**
 * @brief Test interface for Logger class to enable testing with mock backends
 *
 * This interface exists to provide access to the protected backend_ member of
 * the Logger class during testing. The Logger class intentionally keeps its
 * backend private/protected to maintain encapsulation, but for testing we need
 * to:
 *
 * 1. Access the mock backend to set expectations on method calls
 * 2. Verify that Logger properly delegates to its backend implementation
 *
 * This test interface provides a backend() method to expose the backend
 * instance while keeping the production Logger class properly encapsulated.
 *
 * @see Logger
 * @see LoggerTest for usage examples with mock backends
 */
template <typename Backend>
class LoggerTestInterface : public Logger<Backend> {
 public:
  /**
   * @brief Get reference to the underlying backend implementation
   * @return Reference to the backend instance
   */
  Backend& backend() { return this->backend_; }
};

}  // namespace logger
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_TESTS_FRAMEWORK_TOOLS_LOGGER_LOGGERTESTINTERFACE_H_