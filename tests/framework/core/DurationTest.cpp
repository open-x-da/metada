#include <gtest/gtest.h>

#include <sstream>

#include "DateTime.hpp"
#include "Duration.hpp"

using namespace metada::core;

TEST(DurationTest, DefaultConstructor) {
  Duration d;
  EXPECT_EQ(d.totalSeconds(), 0);
  EXPECT_EQ(d.days(), 0);
  EXPECT_EQ(d.hours(), 0);
  EXPECT_EQ(d.minutes(), 0);
  EXPECT_EQ(d.seconds(), 0);
}

TEST(DurationTest, SecondsConstructor) {
  Duration d(3665);  // 1h 1m 5s
  EXPECT_EQ(d.totalSeconds(), 3665);
  EXPECT_EQ(d.days(), 0);
  EXPECT_EQ(d.hours(), 1);
  EXPECT_EQ(d.minutes(), 1);
  EXPECT_EQ(d.seconds(), 5);
}

TEST(DurationTest, ChronoDurationConstructor) {
  auto chrono_hours = std::chrono::hours(2);
  Duration d(chrono_hours);
  EXPECT_EQ(d.totalSeconds(), 7200);
  EXPECT_EQ(d.hours(), 2);
}

TEST(DurationTest, StringConstructor) {
  // Test basic string parsing
  Duration d1("1d 2h 3m 4s");
  EXPECT_EQ(d1.days(), 1);
  EXPECT_EQ(d1.hours(), 2);
  EXPECT_EQ(d1.minutes(), 3);
  EXPECT_EQ(d1.seconds(), 4);
  EXPECT_EQ(d1.totalSeconds(), 93784);  // 1d 2h 3m 4s in seconds

  // Test with spaces
  Duration d2("1d   2h   3m   4s");
  EXPECT_EQ(d2.totalSeconds(), 93784);

  // Test with different order
  Duration d3("4s 3m 2h 1d");
  EXPECT_EQ(d3.totalSeconds(), 93784);

  // Test with missing parts
  Duration d4("2h 30m");
  EXPECT_EQ(d4.hours(), 2);
  EXPECT_EQ(d4.minutes(), 30);
  EXPECT_EQ(d4.totalSeconds(), 9000);

  // Test with just a number (seconds)
  Duration d5("3600");
  EXPECT_EQ(d5.hours(), 1);
  EXPECT_EQ(d5.totalSeconds(), 3600);

  // Test empty string
  Duration d6("");
  EXPECT_EQ(d6.totalSeconds(), 0);

  // Test invalid format
  EXPECT_THROW(Duration("invalid"), std::invalid_argument);
}

TEST(DurationTest, StaticFactoryMethods) {
  Duration d1 = Duration::fromDays(1);
  EXPECT_EQ(d1.totalSeconds(), 86400);
  EXPECT_EQ(d1.days(), 1);

  Duration d2 = Duration::fromHours(2);
  EXPECT_EQ(d2.totalSeconds(), 7200);
  EXPECT_EQ(d2.hours(), 2);

  Duration d3 = Duration::fromMinutes(30);
  EXPECT_EQ(d3.totalSeconds(), 1800);
  EXPECT_EQ(d3.minutes(), 30);

  Duration d4 = Duration::fromSeconds(45);
  EXPECT_EQ(d4.totalSeconds(), 45);
  EXPECT_EQ(d4.seconds(), 45);
}

TEST(DurationTest, ArithmeticOperations) {
  Duration d1 = Duration::fromHours(1);
  Duration d2 = Duration::fromMinutes(30);

  // Addition
  Duration sum = d1 + d2;
  EXPECT_EQ(sum.hours(), 1);
  EXPECT_EQ(sum.minutes(), 30);
  EXPECT_EQ(sum.totalSeconds(), 5400);

  // Subtraction
  Duration diff = d1 - d2;
  EXPECT_EQ(diff.minutes(), 30);
  EXPECT_EQ(diff.totalSeconds(), 1800);

  // Multiplication
  Duration mul = d2 * 2;
  EXPECT_EQ(mul.hours(), 1);
  EXPECT_EQ(mul.totalSeconds(), 3600);

  // Division
  Duration div = d1 / 2;
  EXPECT_EQ(div.minutes(), 30);
  EXPECT_EQ(div.totalSeconds(), 1800);

  // Division by zero
  EXPECT_THROW(d1 / 0, std::invalid_argument);
}

TEST(DurationTest, CompoundAssignment) {
  Duration d1 = Duration::fromHours(1);
  Duration d2 = Duration::fromMinutes(30);

  // +=
  d1 += d2;
  EXPECT_EQ(d1.hours(), 1);
  EXPECT_EQ(d1.minutes(), 30);
  EXPECT_EQ(d1.totalSeconds(), 5400);

  // -=
  d1 -= d2;
  EXPECT_EQ(d1.totalSeconds(), 3600);
  EXPECT_EQ(d1.hours(), 1);
  EXPECT_EQ(d1.minutes(), 0);

  // *=
  d2 *= 2;
  EXPECT_EQ(d2.hours(), 1);
  EXPECT_EQ(d2.totalSeconds(), 3600);

  // /=
  d2 /= 2;
  EXPECT_EQ(d2.minutes(), 30);
  EXPECT_EQ(d2.totalSeconds(), 1800);

  // /= with zero
  EXPECT_THROW(d2 /= 0, std::invalid_argument);
}

TEST(DurationTest, UnaryOperator) {
  Duration d1 = Duration::fromHours(1);
  Duration d2 = -d1;

  EXPECT_EQ(d2.totalSeconds(), -3600);
  EXPECT_EQ(d1.totalSeconds(), 3600);  // Original unchanged
}

TEST(DurationTest, Comparison) {
  Duration d1 = Duration::fromHours(1);
  Duration d2 = Duration::fromMinutes(60);
  Duration d3 = Duration::fromMinutes(30);

  EXPECT_EQ(d1, d2);
  EXPECT_NE(d1, d3);
  EXPECT_GT(d1, d3);
  EXPECT_LT(d3, d1);
  EXPECT_GE(d1, d2);
  EXPECT_LE(d2, d1);
}

TEST(DurationTest, ToString) {
  Duration d1("1d 2h 3m 4s");
  EXPECT_EQ(d1.toString(), "1d 2h 3m 4s");

  Duration d2 = Duration::fromHours(1);
  EXPECT_EQ(d2.toString(), "1h 0m 0s");

  Duration d3 = Duration::fromSeconds(0);
  EXPECT_EQ(d3.toString(), "0s");

  Duration d4 = Duration::fromSeconds(65);
  EXPECT_EQ(d4.toString(), "1m 5s");
}

TEST(DurationTest, StreamOperator) {
  Duration d("1d 2h 3m 4s");
  std::stringstream ss;
  ss << d;
  EXPECT_EQ(ss.str(), d.toString());
}

TEST(DurationTest, DateTimeIntegration) {
  // Test DateTime - DateTime = Duration
  DateTime dt1(2023, 4, 15, 12, 0, 0);
  DateTime dt2(2023, 4, 15, 13, 30, 45);
  Duration diff = dt2 - dt1;
  EXPECT_EQ(diff.hours(), 1);
  EXPECT_EQ(diff.minutes(), 30);
  EXPECT_EQ(diff.seconds(), 45);

  // Test DateTime + Duration = DateTime
  DateTime dt3 = dt1 + diff;
  EXPECT_EQ(dt3, dt2);

  // Test DateTime - Duration = DateTime
  DateTime dt4 = dt2 - diff;
  EXPECT_EQ(dt4, dt1);

  // Test DateTime += Duration
  DateTime dt5 = dt1;
  dt5 += diff;
  EXPECT_EQ(dt5, dt2);

  // Test DateTime -= Duration
  DateTime dt6 = dt2;
  dt6 -= diff;
  EXPECT_EQ(dt6, dt1);

  // Ensure Duration works with larger time spans
  DateTime start(2023, 1, 1);
  DateTime end(2023, 12, 31);
  Duration yearDiff = end - start;
  EXPECT_EQ(yearDiff.days(), 364);  // Not leap year

  // Parse from string and operate
  Duration parsedDuration("2d 5h");
  DateTime result = dt1 + parsedDuration;
  EXPECT_EQ(result.day(), 17);
  EXPECT_EQ(result.hour(), 17);
}