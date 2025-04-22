#pragma once

#include <chrono>
#include <compare>
#include <concepts>
#include <format>
#include <string>
#include <string_view>

namespace metada {
namespace core {

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

  // Add/subtract durations
  template <typename Rep, typename Period>
  DateTime& add(const std::chrono::duration<Rep, Period>& duration) noexcept;

  template <typename Rep, typename Period>
  DateTime& subtract(
      const std::chrono::duration<Rep, Period>& duration) noexcept;

  // Format the date/time using std::format
  [[nodiscard]] std::string format(
      std::string_view fmt = "%Y-%m-%d %H:%M:%S") const;

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

 private:
  std::chrono::system_clock::time_point m_timePoint;
  std::chrono::year_month_day m_date;
  std::chrono::hh_mm_ss<std::chrono::seconds> m_time;

  // Update m_date and m_time from m_timePoint
  void updateComponents() noexcept;
};

//
// Template implementations
//

template <typename Rep, typename Period>
DateTime& DateTime::add(
    const std::chrono::duration<Rep, Period>& duration) noexcept {
  m_timePoint += duration;
  updateComponents();
  return *this;
}

template <typename Rep, typename Period>
DateTime& DateTime::subtract(
    const std::chrono::duration<Rep, Period>& duration) noexcept {
  m_timePoint -= duration;
  updateComponents();
  return *this;
}

}  // namespace core
}  // namespace metada