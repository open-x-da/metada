#include "DateTime.hpp"

#include <format>
#include <regex>
#include <sstream>

namespace metada {

DateTime::DateTime() noexcept : m_timePoint(std::chrono::system_clock::now()) {
  updateComponents();
}

DateTime::DateTime(int year, int month, int day, int hour, int minute,
                   int second) noexcept {
  // Create date components
  auto y = std::chrono::year{year};
  auto m = std::chrono::month{static_cast<unsigned>(month)};
  auto d = std::chrono::day{static_cast<unsigned>(day)};
  m_date = std::chrono::year_month_day{y, m, d};

  // Create time components
  m_time = std::chrono::hh_mm_ss<std::chrono::seconds>{
      std::chrono::hours{hour} + std::chrono::minutes{minute} +
      std::chrono::seconds{second}};

  // Convert to timepoint
  auto dp = std::chrono::sys_days(m_date);
  m_timePoint = dp + m_time.hours() + m_time.minutes() + m_time.seconds();
}

DateTime::DateTime(
    const std::chrono::system_clock::time_point& timepoint) noexcept
    : m_timePoint(timepoint) {
  updateComponents();
}

DateTime::DateTime(const std::string& iso8601_string) {
  // Regular expression for ISO8601 format: YYYY-MM-DDThh:mm:ss[.sss]Z
  // This will also accept formats like YYYY-MM-DDThh:mm:ss[.sss]+00:00
  static const std::regex iso8601_regex(
      R"(^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})(?:\.\d+)?(?:Z|[+-]\d{2}:\d{2})$)");

  std::smatch matches;
  if (!std::regex_match(iso8601_string, matches, iso8601_regex)) {
    throw std::invalid_argument("Invalid ISO8601 format: " + iso8601_string);
  }

  // Extract components
  int year = std::stoi(matches[1]);
  int month = std::stoi(matches[2]);
  int day = std::stoi(matches[3]);
  int hour = std::stoi(matches[4]);
  int minute = std::stoi(matches[5]);
  int second = std::stoi(matches[6]);

  // Validate components
  if (month < 1 || month > 12 || day < 1 || day > 31 || hour < 0 || hour > 23 ||
      minute < 0 || minute > 59 || second < 0 || second > 59) {
    throw std::invalid_argument(
        "Invalid date/time components in ISO8601 string: " + iso8601_string);
  }

  // Create date components
  auto y = std::chrono::year{year};
  auto m = std::chrono::month{static_cast<unsigned>(month)};
  auto d = std::chrono::day{static_cast<unsigned>(day)};
  m_date = std::chrono::year_month_day{y, m, d};

  // Create time components
  m_time = std::chrono::hh_mm_ss<std::chrono::seconds>{
      std::chrono::hours{hour} + std::chrono::minutes{minute} +
      std::chrono::seconds{second}};

  // Convert to timepoint
  auto dp = std::chrono::sys_days(m_date);
  m_timePoint = dp + m_time.hours() + m_time.minutes() + m_time.seconds();

  // Note: This implementation ignores the timezone information and treats
  // the time as if it were in the local timezone. For a more complete
  // implementation, timezone adjustment would be needed.
}

int DateTime::year() const noexcept {
  return static_cast<int>(m_date.year());
}

int DateTime::month() const noexcept {
  return static_cast<int>(static_cast<unsigned>(m_date.month()));
}

int DateTime::day() const noexcept {
  return static_cast<int>(static_cast<unsigned>(m_date.day()));
}

int DateTime::hour() const noexcept {
  return static_cast<int>(m_time.hours().count());
}

int DateTime::minute() const noexcept {
  return static_cast<int>(m_time.minutes().count());
}

int DateTime::second() const noexcept {
  return static_cast<int>(m_time.seconds().count());
}

std::chrono::system_clock::time_point DateTime::timePoint() const noexcept {
  return m_timePoint;
}

std::string DateTime::format(std::string_view fmt) const {
  try {
    // Create a vector of arguments for std::vformat
    std::vector<std::string> args = {
        std::to_string(year()),   std::to_string(month()),
        std::to_string(day()),    std::to_string(hour()),
        std::to_string(minute()), std::to_string(second())};

    // Use custom formatting since we can't validate fmt at compile time
    std::string result = std::string(fmt);
    // Replace placeholders with values
    size_t pos = 0;
    for (const auto& arg : args) {
      pos = result.find("%", pos);
      if (pos != std::string::npos) {
        result.replace(pos, 2, arg);
        pos += arg.length();
      }
    }
    return result;
  } catch (const std::exception&) {
    // Fallback to a default format if the provided format is invalid
    return std::format("{:04d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}", year(),
                       month(), day(), hour(), minute(), second());
  }
}

std::string DateTime::iso8601() const {
  // Format: YYYY-MM-DDThh:mm:ssZ
  return std::format("{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}Z", year(),
                     month(), day(), hour(), minute(), second());
}

DateTime DateTime::now() noexcept {
  return DateTime{std::chrono::system_clock::now()};
}

void DateTime::updateComponents() noexcept {
  // Convert to date
  auto days = std::chrono::time_point_cast<std::chrono::days>(m_timePoint);
  m_date = std::chrono::year_month_day{
      std::chrono::sys_days{days.time_since_epoch()}};

  // Convert to time
  auto timeOfDay = m_timePoint - days;
  auto seconds = std::chrono::duration_cast<std::chrono::seconds>(timeOfDay);
  m_time = std::chrono::hh_mm_ss<std::chrono::seconds>{seconds};
}

}  // namespace metada