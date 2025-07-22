#pragma once

#include <chrono>
#include <compare>
#include <stdexcept>
#include <string>

namespace metada {

class DateTime;  // Forward declaration

/**
 * @brief A Duration class wrapping C++20 chrono duration functionality
 */
class Duration {
 public:
  // Default constructor: zero duration
  Duration() noexcept;

  // Construct from seconds
  explicit Duration(int64_t seconds) noexcept;

  // Construct from chrono duration
  template <typename Rep, typename Period>
  explicit Duration(const std::chrono::duration<Rep, Period>& duration) noexcept
      : m_duration(std::chrono::duration_cast<std::chrono::seconds>(duration)) {
  }

  // Construct from string (format: "1d 2h 3m 4s" or variations)
  explicit Duration(const std::string& duration_string);

  // Get components
  [[nodiscard]] int64_t totalSeconds() const noexcept;
  [[nodiscard]] int64_t days() const noexcept;
  [[nodiscard]] int hours() const noexcept;
  [[nodiscard]] int minutes() const noexcept;
  [[nodiscard]] int seconds() const noexcept;

  // Get as chrono duration
  [[nodiscard]] std::chrono::seconds asChrono() const noexcept;

  // Format as string
  [[nodiscard]] std::string toString() const;

  // Static factory methods
  [[nodiscard]] static Duration fromDays(int64_t days) noexcept;
  [[nodiscard]] static Duration fromHours(int64_t hours) noexcept;
  [[nodiscard]] static Duration fromHoursF(
      double hours) noexcept;  // For float DHR values
  [[nodiscard]] static Duration fromMinutes(int64_t minutes) noexcept;
  [[nodiscard]] static Duration fromSeconds(int64_t seconds) noexcept;

  // Arithmetic operators
  Duration operator+(const Duration& other) const noexcept;
  Duration operator-(const Duration& other) const noexcept;
  Duration operator*(int64_t factor) const noexcept;
  Duration operator/(int64_t divisor) const;

  // Compound assignment operators
  Duration& operator+=(const Duration& other) noexcept;
  Duration& operator-=(const Duration& other) noexcept;
  Duration& operator*=(int64_t factor) noexcept;
  Duration& operator/=(int64_t divisor);

  // Unary operators
  Duration operator-() const noexcept;

  // Comparison operators
  bool operator==(const Duration& other) const noexcept = default;
  auto operator<=>(const Duration& other) const noexcept = default;

 private:
  std::chrono::seconds m_duration;
};

// Operator overload for stream insertion
std::ostream& operator<<(std::ostream& os, const Duration& duration);

}  // namespace metada