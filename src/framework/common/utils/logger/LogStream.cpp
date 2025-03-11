#include "utils/logger/LogStream.hpp"

#include "utils/logger/ILogger.hpp"

namespace metada::framework {

void LogStream::Flush() {
  const std::string message = stream_.str();

  switch (level_) {
    case LogLevel::Debug:
      logger_.Debug(message);
      break;
    case LogLevel::Info:
      logger_.Info(message);
      break;
    case LogLevel::Warning:
      logger_.Warning(message);
      break;
    case LogLevel::Error:
      logger_.Error(message);
      break;
  }

  // Clear the stream for potential reuse
  stream_.str("");
  stream_.clear();
}

}  // namespace metada::framework