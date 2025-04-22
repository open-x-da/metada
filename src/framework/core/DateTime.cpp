#include "DateTime.hpp"

#include <format>

namespace metada {
namespace core {

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

}  // namespace core
}  // namespace metada