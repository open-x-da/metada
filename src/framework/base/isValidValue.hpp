#pragma once
#include <cmath>
#include <limits>

namespace metada {

constexpr double R8BFMS = 10.0E10;  // BUFR missing value for real*8

// High-performance check for BUFR missing value or NaN
inline bool isValidValue(double value) noexcept {
  constexpr double epsilon = std::numeric_limits<double>::epsilon() * 100.0;
  // Use fast math: avoid function call overhead, use std::abs and std::isnan
  return (std::abs(value - R8BFMS) > std::abs(R8BFMS) * epsilon) &&
         !std::isnan(value);
}

}  // namespace metada
