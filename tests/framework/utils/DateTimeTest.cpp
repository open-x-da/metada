#include <gtest/gtest.h>

#include <chrono>
#include <sstream>
#include <thread>

#include "DateTime.hpp"

using namespace metada::utils;

TEST(DateTimeTest, DefaultConstructor) {
  DateTime dt;
  // Current time should be close to now
  auto now = std::chrono::system_clock::now();
  auto diff =
      std::chrono::duration_cast<std::chrono::seconds>(now - dt.timePoint())
          .count();
  EXPECT_LE(std::abs(diff), 1);  // Allow 1 second tolerance
}

TEST(DateTimeTest, ComponentConstructor) {
  DateTime dt(2023, 4, 15, 12, 30, 45);
  EXPECT_EQ(dt.year(), 2023);
  EXPECT_EQ(dt.month(), 4);
  EXPECT_EQ(dt.day(), 15);
  EXPECT_EQ(dt.hour(), 12);
  EXPECT_EQ(dt.minute(), 30);
  EXPECT_EQ(dt.second(), 45);
}

TEST(DateTimeTest, PartialComponentConstructor) {
  // Test with just year, month, day (time should default to 00:00:00)
  DateTime date_only(2023, 4, 15);
  EXPECT_EQ(date_only.year(), 2023);
  EXPECT_EQ(date_only.month(), 4);
  EXPECT_EQ(date_only.day(), 15);
  EXPECT_EQ(date_only.hour(), 0);
  EXPECT_EQ(date_only.minute(), 0);
  EXPECT_EQ(date_only.second(), 0);

  // Test with year, month, day, hour (minutes and seconds should default to
  // 00:00)
  DateTime with_hour(2023, 4, 15, 12);
  EXPECT_EQ(with_hour.year(), 2023);
  EXPECT_EQ(with_hour.month(), 4);
  EXPECT_EQ(with_hour.day(), 15);
  EXPECT_EQ(with_hour.hour(), 12);
  EXPECT_EQ(with_hour.minute(), 0);
  EXPECT_EQ(with_hour.second(), 0);

  // Test with year, month, day, hour, minute (seconds should default to 00)
  DateTime with_minute(2023, 4, 15, 12, 30);
  EXPECT_EQ(with_minute.year(), 2023);
  EXPECT_EQ(with_minute.month(), 4);
  EXPECT_EQ(with_minute.day(), 15);
  EXPECT_EQ(with_minute.hour(), 12);
  EXPECT_EQ(with_minute.minute(), 30);
  EXPECT_EQ(with_minute.second(), 0);

  // Test ISO8601 output with incomplete time
  EXPECT_EQ(date_only.iso8601(), "2023-04-15T00:00:00Z");
  EXPECT_EQ(with_hour.iso8601(), "2023-04-15T12:00:00Z");
  EXPECT_EQ(with_minute.iso8601(), "2023-04-15T12:30:00Z");
}

TEST(DateTimeTest, TimePointConstructor) {
  // Create a specific time point
  std::chrono::system_clock::time_point tp =
      std::chrono::system_clock::from_time_t(0) + std::chrono::hours(12);

  DateTime dt(tp);
  EXPECT_EQ(dt.timePoint(), tp);
}

TEST(DateTimeTest, ISO8601Constructor) {
  // Test valid ISO8601 string
  DateTime dt1("2023-04-15T12:30:45Z");
  EXPECT_EQ(dt1.year(), 2023);
  EXPECT_EQ(dt1.month(), 4);
  EXPECT_EQ(dt1.day(), 15);
  EXPECT_EQ(dt1.hour(), 12);
  EXPECT_EQ(dt1.minute(), 30);
  EXPECT_EQ(dt1.second(), 45);

  // Test with timezone offset
  DateTime dt2("2023-04-15T12:30:45+00:00");
  EXPECT_EQ(dt2.year(), 2023);
  EXPECT_EQ(dt2.month(), 4);
  EXPECT_EQ(dt2.day(), 15);
  EXPECT_EQ(dt2.hour(), 12);
  EXPECT_EQ(dt2.minute(), 30);
  EXPECT_EQ(dt2.second(), 45);

  // Test with milliseconds
  DateTime dt3("2023-04-15T12:30:45.123Z");
  EXPECT_EQ(dt3.year(), 2023);
  EXPECT_EQ(dt3.month(), 4);
  EXPECT_EQ(dt3.day(), 15);
  EXPECT_EQ(dt3.hour(), 12);
  EXPECT_EQ(dt3.minute(), 30);
  EXPECT_EQ(dt3.second(), 45);

  // Test invalid format exceptions
  EXPECT_THROW(DateTime("2023-04-15"), std::invalid_argument);
  EXPECT_THROW(DateTime("2023-04-15T12:30"), std::invalid_argument);
  EXPECT_THROW(DateTime("not-a-date"), std::invalid_argument);

  // Test invalid date components
  EXPECT_THROW(DateTime("2023-13-15T12:30:45Z"),
               std::invalid_argument);  // Invalid month
  EXPECT_THROW(DateTime("2023-04-32T12:30:45Z"),
               std::invalid_argument);  // Invalid day
  EXPECT_THROW(DateTime("2023-04-15T24:30:45Z"),
               std::invalid_argument);  // Invalid hour
}

TEST(DateTimeTest, OperatorPlus) {
  DateTime dt(2023, 4, 15, 12, 30, 45);

  // Test adding days
  DateTime nextDay = dt + std::chrono::days(1);
  EXPECT_EQ(nextDay.year(), 2023);
  EXPECT_EQ(nextDay.month(), 4);
  EXPECT_EQ(nextDay.day(), 16);
  EXPECT_EQ(nextDay.hour(), 12);
  EXPECT_EQ(nextDay.minute(), 30);
  EXPECT_EQ(nextDay.second(), 45);

  // Test adding hours
  DateTime laterHour = dt + std::chrono::hours(3);
  EXPECT_EQ(laterHour.hour(), 15);

  // Test adding minutes
  DateTime laterMinute = dt + std::chrono::minutes(45);
  EXPECT_EQ(laterMinute.hour(), 13);
  EXPECT_EQ(laterMinute.minute(), 15);

  // Test that original DateTime remains unchanged
  EXPECT_EQ(dt.day(), 15);
  EXPECT_EQ(dt.hour(), 12);
  EXPECT_EQ(dt.minute(), 30);
}

TEST(DateTimeTest, OperatorMinus) {
  DateTime dt(2023, 4, 15, 12, 30, 45);

  // Test subtracting days
  DateTime previousDay = dt - std::chrono::days(1);
  EXPECT_EQ(previousDay.year(), 2023);
  EXPECT_EQ(previousDay.month(), 4);
  EXPECT_EQ(previousDay.day(), 14);
  EXPECT_EQ(previousDay.hour(), 12);
  EXPECT_EQ(previousDay.minute(), 30);
  EXPECT_EQ(previousDay.second(), 45);

  // Test subtracting hours
  DateTime earlierHour = dt - std::chrono::hours(3);
  EXPECT_EQ(earlierHour.hour(), 9);

  // Test subtracting minutes
  DateTime earlierMinute = dt - std::chrono::minutes(45);
  EXPECT_EQ(earlierMinute.hour(), 11);
  EXPECT_EQ(earlierMinute.minute(), 45);

  // Test that original DateTime remains unchanged
  EXPECT_EQ(dt.day(), 15);
  EXPECT_EQ(dt.hour(), 12);
  EXPECT_EQ(dt.minute(), 30);
}

TEST(DateTimeTest, OperatorDateDifference) {
  DateTime dt1(2023, 4, 15, 12, 30, 45);
  DateTime dt2(2023, 4, 15, 14, 45, 15);
  DateTime dt3(2023, 4, 16, 12, 30, 45);

  // Test difference in hours
  auto diff1 = dt2 - dt1;
  EXPECT_EQ(diff1.totalSeconds(), 8070);  // 2h 14m 30s = 8070 seconds

  // Test difference in days
  auto diff2 = dt3 - dt1;
  EXPECT_EQ(diff2.totalSeconds(), 86400);  // 24h = 86400 seconds

  // Test negative difference
  auto diff3 = dt1 - dt2;
  EXPECT_EQ(diff3.totalSeconds(), -8070);  // -2h 14m 30s = -8070 seconds
}

TEST(DateTimeTest, CompoundAssignmentOperators) {
  // Test +=
  DateTime dt1(2023, 4, 15, 12, 30, 45);
  dt1 += std::chrono::hours(1);
  EXPECT_EQ(dt1.hour(), 13);
  EXPECT_EQ(dt1.minute(), 30);

  dt1 += std::chrono::minutes(45);
  EXPECT_EQ(dt1.hour(), 14);
  EXPECT_EQ(dt1.minute(), 15);

  // Test -=
  DateTime dt2(2023, 4, 15, 12, 30, 45);
  dt2 -= std::chrono::hours(3);
  EXPECT_EQ(dt2.hour(), 9);
  EXPECT_EQ(dt2.minute(), 30);

  dt2 -= std::chrono::minutes(45);
  EXPECT_EQ(dt2.hour(), 8);
  EXPECT_EQ(dt2.minute(), 45);

  // Test sequential operators
  DateTime dt3(2023, 4, 15, 12, 0, 0);
  dt3 += std::chrono::hours(3);
  dt3 += std::chrono::minutes(30);
  EXPECT_EQ(dt3.hour(), 15);
  EXPECT_EQ(dt3.minute(), 30);

  DateTime dt4(2023, 4, 15, 12, 0, 0);
  dt4 -= std::chrono::hours(3);
  dt4 -= std::chrono::minutes(30);
  EXPECT_EQ(dt4.hour(), 8);
  EXPECT_EQ(dt4.minute(), 30);
}

TEST(DateTimeTest, ChainedOperatorExpressions) {
  // Base date and various durations
  DateTime base(2023, 4, 15, 12, 0, 0);
  auto oneDay = std::chrono::days(1);
  auto threeHours = std::chrono::hours(3);
  auto ninetyMinutes = std::chrono::minutes(90);

  // Test d = a + b (single addition chain)
  DateTime result1 = base + oneDay;
  EXPECT_EQ(result1.year(), 2023);
  EXPECT_EQ(result1.month(), 4);
  EXPECT_EQ(result1.day(), 16);
  EXPECT_EQ(result1.hour(), 12);

  // Test d = a + b + c (multiple addition chain)
  DateTime result2 = base + oneDay + threeHours;
  EXPECT_EQ(result2.day(), 16);
  EXPECT_EQ(result2.hour(), 15);

  // Test d = a - b (single subtraction chain)
  DateTime result3 = base - threeHours;
  EXPECT_EQ(result3.day(), 15);
  EXPECT_EQ(result3.hour(), 9);

  // Test d = a - b - c (multiple subtraction chain)
  DateTime result4 = base - threeHours - ninetyMinutes;
  EXPECT_EQ(result4.hour(), 7);
  EXPECT_EQ(result4.minute(), 30);

  // Test d = a + b - c (mixed operators)
  DateTime result5 = base + oneDay - threeHours;
  EXPECT_EQ(result5.day(), 16);
  EXPECT_EQ(result5.hour(), 9);

  // Test d = a - b + c (mixed operators, different order)
  DateTime result6 = base - threeHours + ninetyMinutes;
  EXPECT_EQ(result6.day(), 15);
  EXPECT_EQ(result6.hour(), 10);
  EXPECT_EQ(result6.minute(), 30);

  // Test DateTime difference calculation with expressions
  auto dateA = base + oneDay;
  auto dateB = base - threeHours;
  auto duration = dateA - dateB;
  EXPECT_EQ(duration.totalSeconds(),
            86400 + 10800);  // 1 day + 3 hours in seconds

  // Verify base date hasn't changed (immutability check)
  EXPECT_EQ(base.day(), 15);
  EXPECT_EQ(base.hour(), 12);
  EXPECT_EQ(base.minute(), 0);
}

TEST(DateTimeTest, FormatDefault) {
  DateTime dt(2023, 4, 15, 12, 30, 45);
  std::string formatted = dt.format();
  EXPECT_TRUE(formatted.find("2023") != std::string::npos);
  EXPECT_TRUE(formatted.find("4") != std::string::npos);
  EXPECT_TRUE(formatted.find("15") != std::string::npos);
  EXPECT_TRUE(formatted.find("12") != std::string::npos);
  EXPECT_TRUE(formatted.find("30") != std::string::npos);
  EXPECT_TRUE(formatted.find("45") != std::string::npos);
}

TEST(DateTimeTest, ISO8601Format) {
  // Test basic ISO8601 formatting
  DateTime dt1(2023, 4, 15, 12, 30, 45);
  EXPECT_EQ(dt1.iso8601(), "2023-04-15T12:30:45Z");

  // Test single-digit values (should be zero-padded)
  DateTime dt2(2023, 1, 5, 9, 5, 5);
  EXPECT_EQ(dt2.iso8601(), "2023-01-05T09:05:05Z");

  // Test round-trip (parse ISO8601 string and convert back)
  DateTime dt3("2023-12-31T23:59:59Z");
  EXPECT_EQ(dt3.iso8601(), "2023-12-31T23:59:59Z");
}

TEST(DateTimeTest, FormatCustom) {
  DateTime dt(2023, 4, 15, 12, 30, 45);
  std::string formatted = dt.format("%Y/%m/%d");
  EXPECT_TRUE(formatted.find("2023") != std::string::npos);
  EXPECT_TRUE(formatted.find("4") != std::string::npos);
  EXPECT_TRUE(formatted.find("15") != std::string::npos);
}

TEST(DateTimeTest, StaticNow) {
  DateTime dt = DateTime::now();
  auto now = std::chrono::system_clock::now();
  auto diff =
      std::chrono::duration_cast<std::chrono::seconds>(now - dt.timePoint())
          .count();
  EXPECT_LE(std::abs(diff), 1);  // Allow 1 second tolerance
}

TEST(DateTimeTest, Comparison) {
  DateTime dt1(2023, 4, 15, 12, 30, 45);
  DateTime dt2(2023, 4, 15, 12, 30, 45);
  DateTime dt3(2023, 4, 15, 12, 30, 46);

  EXPECT_EQ(dt1, dt2);
  EXPECT_NE(dt1, dt3);
  EXPECT_LT(dt1, dt3);
  EXPECT_LE(dt1, dt2);
  EXPECT_GT(dt3, dt1);
  EXPECT_GE(dt2, dt1);
}

TEST(DateTimeTest, EdgeCases) {
  // Test leap year
  DateTime leapDay(2020, 2, 29);
  EXPECT_EQ(leapDay.year(), 2020);
  EXPECT_EQ(leapDay.month(), 2);
  EXPECT_EQ(leapDay.day(), 29);

  // Test wrapping around midnight
  DateTime evening(2023, 4, 15, 23, 59, 59);
  evening += std::chrono::seconds(1);
  EXPECT_EQ(evening.day(), 16);
  EXPECT_EQ(evening.hour(), 0);
  EXPECT_EQ(evening.minute(), 0);
  EXPECT_EQ(evening.second(), 0);

  // Test month boundary
  DateTime endOfMonth(2023, 4, 30, 12, 0, 0);
  endOfMonth += std::chrono::days(1);
  EXPECT_EQ(endOfMonth.month(), 5);
  EXPECT_EQ(endOfMonth.day(), 1);
}

TEST(DateTimeTest, StreamInsertionOperator) {
  DateTime dt(2023, 4, 15, 12, 30, 45);

  // Test stream insertion operator (<<)
  std::stringstream ss;
  ss << dt;

  // Should match the default format output
  std::string expected = dt.format();
  EXPECT_EQ(ss.str(), expected);

  // Test in a more complex stream expression
  std::stringstream ss2;
  ss2 << "DateTime: " << dt << " is valid";
  EXPECT_EQ(ss2.str(), "DateTime: " + expected + " is valid");
}