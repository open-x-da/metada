#include <gtest/gtest.h>

#include <chrono>
#include <thread>

#include "DateTime.hpp"

using namespace metada::core;

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

TEST(DateTimeTest, TimePointConstructor) {
  // Create a specific time point
  std::chrono::system_clock::time_point tp =
      std::chrono::system_clock::from_time_t(0) + std::chrono::hours(12);

  DateTime dt(tp);
  EXPECT_EQ(dt.timePoint(), tp);
}

TEST(DateTimeTest, GetComponents) {
  DateTime dt(2023, 4, 15, 12, 30, 45);
  EXPECT_EQ(dt.year(), 2023);
  EXPECT_EQ(dt.month(), 4);
  EXPECT_EQ(dt.day(), 15);
  EXPECT_EQ(dt.hour(), 12);
  EXPECT_EQ(dt.minute(), 30);
  EXPECT_EQ(dt.second(), 45);
}

TEST(DateTimeTest, AddDuration) {
  DateTime dt(2023, 4, 15, 12, 30, 45);

  dt.add(std::chrono::days(1));
  EXPECT_EQ(dt.day(), 16);

  dt.add(std::chrono::hours(1));
  EXPECT_EQ(dt.hour(), 13);

  dt.add(std::chrono::minutes(30));
  EXPECT_EQ(dt.minute(), 0);
  EXPECT_EQ(dt.hour(), 14);
}

TEST(DateTimeTest, SubtractDuration) {
  DateTime dt(2023, 4, 15, 12, 30, 45);

  dt.subtract(std::chrono::hours(1));
  EXPECT_EQ(dt.hour(), 11);

  dt.subtract(std::chrono::minutes(31));
  EXPECT_EQ(dt.hour(), 10);
  EXPECT_EQ(dt.minute(), 59);
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
  evening.add(std::chrono::seconds(1));
  EXPECT_EQ(evening.day(), 16);
  EXPECT_EQ(evening.hour(), 0);
  EXPECT_EQ(evening.minute(), 0);
  EXPECT_EQ(evening.second(), 0);

  // Test month boundary
  DateTime endOfMonth(2023, 4, 30, 12, 0, 0);
  endOfMonth.add(std::chrono::days(1));
  EXPECT_EQ(endOfMonth.month(), 5);
  EXPECT_EQ(endOfMonth.day(), 1);
}