#include "Duration.hpp"

#include <format>
#include <regex>
#include <sstream>

namespace metada::core {

Duration::Duration() noexcept : m_duration(std::chrono::seconds(0)) {}

Duration::Duration(int64_t seconds) noexcept
    : m_duration(std::chrono::seconds(seconds)) {}

Duration::Duration(const std::string& duration_string) {
  // Initialize to zero
  m_duration = std::chrono::seconds(0);

  if (duration_string.empty()) {
    return;
  }

  // Regular expressions for different duration parts
  static const std::regex days_regex(R"((\d+)\s*d)");
  static const std::regex hours_regex(R"((\d+)\s*h)");
  static const std::regex minutes_regex(
      R"((\d+)\s*m(?!s))");  // m not followed by s
  static const std::regex seconds_regex(R"((\d+)\s*s)");
  static const std::regex integer_regex(
      R"(^(\d+)$)");  // Just a number, assume seconds

  std::smatch match;

  // Parse days
  std::string str = duration_string;
  if (std::regex_search(str, match, days_regex)) {
    int64_t days = std::stoll(match[1]);
    m_duration += std::chrono::seconds(days * 24 * 60 * 60);
  }

  // Parse hours
  if (std::regex_search(str, match, hours_regex)) {
    int64_t hours = std::stoll(match[1]);
    m_duration += std::chrono::seconds(hours * 60 * 60);
  }

  // Parse minutes
  if (std::regex_search(str, match, minutes_regex)) {
    int64_t minutes = std::stoll(match[1]);
    m_duration += std::chrono::seconds(minutes * 60);
  }

  // Parse seconds
  if (std::regex_search(str, match, seconds_regex)) {
    int64_t seconds = std::stoll(match[1]);
    m_duration += std::chrono::seconds(seconds);
  }

  // If it's just a number, assume seconds
  if (std::regex_match(str, match, integer_regex)) {
    int64_t seconds = std::stoll(match[1]);
    m_duration = std::chrono::seconds(seconds);
    return;
  }

  // If no valid parts found and it's not just a number, it's an error
  if (m_duration.count() == 0 && !str.empty()) {
    throw std::invalid_argument("Invalid duration format: " + duration_string);
  }
}

int64_t Duration::totalSeconds() const noexcept {
  return m_duration.count();
}

int64_t Duration::days() const noexcept {
  return m_duration.count() / (24 * 60 * 60);
}

int Duration::hours() const noexcept {
  return (m_duration.count() % (24 * 60 * 60)) / (60 * 60);
}

int Duration::minutes() const noexcept {
  return (m_duration.count() % (60 * 60)) / 60;
}

int Duration::seconds() const noexcept {
  return m_duration.count() % 60;
}

std::chrono::seconds Duration::asChrono() const noexcept {
  return m_duration;
}

std::string Duration::toString() const {
  std::stringstream ss;
  bool hasPrevious = false;

  if (days() > 0) {
    ss << days() << "d";
    hasPrevious = true;
  }

  if (hours() > 0 || hasPrevious) {
    if (hasPrevious) ss << " ";
    ss << hours() << "h";
    hasPrevious = true;
  }

  if (minutes() > 0 || hasPrevious) {
    if (hasPrevious) ss << " ";
    ss << minutes() << "m";
    hasPrevious = true;
  }

  // Always include seconds (even when they are zero)
  if (hasPrevious) ss << " ";
  ss << seconds() << "s";

  return ss.str();
}

Duration Duration::fromDays(int64_t days) noexcept {
  return Duration(days * 24 * 60 * 60);
}

Duration Duration::fromHours(int64_t hours) noexcept {
  return Duration(hours * 60 * 60);
}

Duration Duration::fromMinutes(int64_t minutes) noexcept {
  return Duration(minutes * 60);
}

Duration Duration::fromSeconds(int64_t seconds) noexcept {
  return Duration(seconds);
}

Duration Duration::operator+(const Duration& other) const noexcept {
  return Duration(m_duration + other.m_duration);
}

Duration Duration::operator-(const Duration& other) const noexcept {
  return Duration(m_duration - other.m_duration);
}

Duration Duration::operator*(int64_t factor) const noexcept {
  return Duration(m_duration * factor);
}

Duration Duration::operator/(int64_t divisor) const {
  if (divisor == 0) {
    throw std::invalid_argument("Division by zero in Duration");
  }
  return Duration(m_duration / divisor);
}

Duration& Duration::operator+=(const Duration& other) noexcept {
  m_duration += other.m_duration;
  return *this;
}

Duration& Duration::operator-=(const Duration& other) noexcept {
  m_duration -= other.m_duration;
  return *this;
}

Duration& Duration::operator*=(int64_t factor) noexcept {
  m_duration *= factor;
  return *this;
}

Duration& Duration::operator/=(int64_t divisor) {
  if (divisor == 0) {
    throw std::invalid_argument("Division by zero in Duration");
  }
  m_duration /= divisor;
  return *this;
}

Duration Duration::operator-() const noexcept {
  return Duration(-m_duration.count());
}

std::ostream& operator<<(std::ostream& os, const Duration& duration) {
  return os << duration.toString();
}

}  // namespace metada::core