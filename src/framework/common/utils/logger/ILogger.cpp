#include "utils/logger/ILogger.hpp"

#include "utils/logger/LogStream.hpp"

namespace metada::framework {

LogStream ILogger::InfoStream() {
  return LogStream(*this, LogLevel::Info);
}

LogStream ILogger::WarningStream() {
  return LogStream(*this, LogLevel::Warning);
}

LogStream ILogger::ErrorStream() {
  return LogStream(*this, LogLevel::Error);
}

LogStream ILogger::DebugStream() {
  return LogStream(*this, LogLevel::Debug);
}

}  // namespace metada::framework