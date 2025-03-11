#include "utils/logger/LogStream.hpp"

#include "utils/logger/ILogger.hpp"

namespace metada::framework {

void LogStream::Flush() {
  const std::string message = stream_.str();

  // Use the unified LogMessage method instead of individual severity methods
  logger_.LogMessage(level_, message);

  // Clear the stream for potential reuse
  stream_.str("");
  stream_.clear();
}

}  // namespace metada::framework