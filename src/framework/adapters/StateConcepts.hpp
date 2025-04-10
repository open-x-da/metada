#pragma once

#include <concepts>
#include <memory>
#include <string>
#include <vector>

namespace metada::framework {

// Individual concepts for state backend requirements
template <typename T>
concept HasClone = requires(const T& t) {
  { t.clone() } -> std::convertible_to<std::unique_ptr<T>>;
};

template <typename T>
concept HasGetData = requires(T& t, const T& ct) {
  { t.getData() } -> std::same_as<void*>;
  { ct.getData() } -> std::same_as<const void*>;
};

template <typename T>
concept HasGetVariableNames = requires(const T& t) {
  { t.getVariableNames() } -> std::same_as<const std::vector<std::string>&>;
};

template <typename T>
concept HasGetDimensions = requires(const T& t) {
  {
    t.getDimensions(std::string{})
  } -> std::same_as<const std::vector<size_t>&>;
};

/**
 * @brief Combined concept for all state backend requirements
 *
 * @details This concept enforces the contract that any state backend
 * must implement. It provides compile-time validation of the required methods
 * and their signatures, ensuring that backends can be properly used with the
 * State class.
 *
 * The concept requires specific method signatures for:
 * - Cloning the state
 * - Accessing state data
 * - Getting variable names
 * - Getting dimensions of state variables
 */
template <typename T>
concept StateBackendType = HasClone<T> && HasGetData<T> &&
                           HasGetVariableNames<T> && HasGetDimensions<T>;

}  // namespace metada::framework