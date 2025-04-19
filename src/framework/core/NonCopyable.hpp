#pragma once

namespace metada::framework {

/**
 * @brief Utility class to make derived classes non-copyable
 *
 * Classes can inherit from this to prevent copying while still allowing
 * move semantics. This is cleaner than explicitly deleting copy operations
 * in each class.
 *
 * Example usage:
 * @code
 * class MyResource : public NonCopyable {
 *   // This class is now non-copyable but still movable
 * };
 * @endcode
 *
 * This is similar to boost::noncopyable but doesn't require Boost.
 */
class NonCopyable {
 protected:
  // Allow construction and destruction
  NonCopyable() = default;
  ~NonCopyable() = default;

  // Prevent copying
  NonCopyable(const NonCopyable&) = delete;
  NonCopyable& operator=(const NonCopyable&) = delete;

  // Allow moving
  NonCopyable(NonCopyable&&) = default;
  NonCopyable& operator=(NonCopyable&&) = default;
};

}  // namespace metada::framework