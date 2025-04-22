#pragma once

#include <chrono>
#include <compare>
#include <concepts>
#include <format>
#include <stdexcept>
#include <string>
#include <string_view>

namespace metada::core {

/**
 * @brief A DateTime class utilizing C++20 features for date and time handling
 */
class DateTime {
 public:
  // Default constructor: current time
  DateTime() noexcept;

  // Construct from year, month, day, hour, minute, second
  DateTime(int year, int month, int day, int hour = 0, int minute = 0,
           int second = 0) noexcept;

  // Construct from timepoint
  explicit DateTime(
      const std::chrono::system_clock::time_point& timepoint) noexcept;

  // Construct from ISO8601 string (e.g., "2023-04-15T12:30:45Z")
  explicit DateTime(const std::string& iso8601_string);

  // Get components
  [[nodiscard]] int year() const noexcept;
  [[nodiscard]] int month() const noexcept;
  [[nodiscard]] int day() const noexcept;
  [[nodiscard]] int hour() const noexcept;
  [[nodiscard]] int minute() const noexcept;
  [[nodiscard]] int second() const noexcept;

  // Get as time_point
  [[nodiscard]] std::chrono::system_clock::time_point timePoint()
      const noexcept;

  // Format the date/time using std::format
  [[nodiscard]] std::string format(
      std::string_view fmt = "%Y-%m-%d %H:%M:%S") const;

  // Format as ISO8601 string (YYYY-MM-DDThh:mm:ssZ)
  [[nodiscard]] std::string iso8601() const;

  // Static methods
  [[nodiscard]] static DateTime now() noexcept;

  // Custom comparison operators
  bool operator==(const DateTime& other) const noexcept {
    return m_timePoint == other.m_timePoint;
  }

  bool operator!=(const DateTime& other) const noexcept {
    return m_timePoint != other.m_timePoint;
  }

  bool operator<(const DateTime& other) const noexcept {
    return m_timePoint < other.m_timePoint;
  }

  bool operator<=(const DateTime& other) const noexcept {
    return m_timePoint <= other.m_timePoint;
  }

  bool operator>(const DateTime& other) const noexcept {
    return m_timePoint > other.m_timePoint;
  }

  bool operator>=(const DateTime& other) const noexcept {
    return m_timePoint >= other.m_timePoint;
  }

  // Operator overloads for immutable duration arithmetic
  template <typename Rep, typename Period>
  [[nodiscard]] DateTime operator+(
      const std::chrono::duration<Rep, Period>& duration) const noexcept {
    DateTime result(*this);
    result.m_timePoint += duration;
    result.updateComponents();
    return result;
  }

  template <typename Rep, typename Period>
  [[nodiscard]] DateTime operator-(
      const std::chrono::duration<Rep, Period>& duration) const noexcept {
    DateTime result(*this);
    result.m_timePoint -= duration;
    result.updateComponents();
    return result;
  }

  // Compound assignment operators for mutable operations
  template <typename Rep, typename Period>
  DateTime& operator+=(
      const std::chrono::duration<Rep, Period>& duration) noexcept {
    m_timePoint += duration;
    updateComponents();
    return *this;
  }

  template <typename Rep, typename Period>
  DateTime& operator-=(
      const std::chrono::duration<Rep, Period>& duration) noexcept {
    m_timePoint -= duration;
    updateComponents();
    return *this;
  }

  // Compute difference between two DateTimes
  [[nodiscard]] std::chrono::seconds operator-(
      const DateTime& other) const noexcept {
    return std::chrono::duration_cast<std::chrono::seconds>(m_timePoint -
                                                            other.m_timePoint);
  }

 private:
  std::chrono::system_clock::time_point m_timePoint;
  std::chrono::year_month_day m_date;
  std::chrono::hh_mm_ss<std::chrono::seconds> m_time;

  // Update m_date and m_time from m_timePoint
  void updateComponents() noexcept;

  // Add/subtract durations (private implementation methods)
  template <typename Rep, typename Period>
  DateTime& add(const std::chrono::duration<Rep, Period>& duration) noexcept {
    m_timePoint += duration;
    updateComponents();
    return *this;
  }

  template <typename Rep, typename Period>
  DateTime& subtract(
      const std::chrono::duration<Rep, Period>& duration) noexcept {
    m_timePoint -= duration;
    updateComponents();
    return *this;
  }
};

}  // namespace metada::core
