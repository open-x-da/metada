#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

namespace metada::framework {

namespace tools::config {
template <typename U>
class Config;
}
using tools::config::Config;
/**
 * @brief Template class for ensemble representation
 *
 * Provides implementation for ensemble operations using a specific State type.
 * Designed for use in ensemble-based data assimilation methods like LETKF.
 *
 * @tparam StateBackend The type of State used in the ensemble
 */
template <typename T>
class Ensemble {
 private:
  std::vector<T> members_;
  T mean_;
  size_t size_;
  std::vector<T> perturbations_;

 public:
  // Constructor with config
  template <typename U>
  explicit Ensemble(const Config<U>& config, size_t size)
      : members_(), mean_(config), size_(size), perturbations_() {
    // Initialize members
    members_.reserve(size);
    for (size_t i = 0; i < size; ++i) {
      members_.emplace_back(config);
    }
  }

  // Member access
  T& getMember(size_t index) {
    if (index >= size_) {
      throw std::out_of_range("Ensemble member index out of range");
    }
    return members_[index];
  }

  const T& getMember(size_t index) const {
    if (index >= size_) {
      throw std::out_of_range("Ensemble member index out of range");
    }
    return members_[index];
  }

  size_t getSize() const { return size_; }

  // Statistical operations
  void computeMean() {
    // Reset mean state
    mean_.reset();

    // Sum all members using State's += operator
    for (size_t i = 0; i < size_; ++i) {
      mean_ += members_[i];
    }

    // Divide by ensemble size
    mean_ *= (1.0 / static_cast<double>(size_));
  }

  T& getMean() { return mean_; }

  const T& getMean() const { return mean_; }

  void computePerturbations() {
    // Ensure mean is computed
    computeMean();

    // Resize perturbations if needed - this copies mean_ to each slot
    perturbations_.resize(size_, mean_);

    // Compute perturbations for each member in-place
    for (size_t i = 0; i < size_; ++i) {
      // Perturbation = -(Mean - Member) = Member - Mean
      perturbations_[i] -= members_[i];  // Now contains (Mean - Member)
      perturbations_[i] *= -1.0;         // Convert to (Member - Mean)
    }
  }

  T& getPerturbation(size_t index) {
    if (index >= size_) {
      throw std::out_of_range("Perturbation index out of range");
    }
    return perturbations_[index];
  }

  const T& getPerturbation(size_t index) const {
    if (index >= size_) {
      throw std::out_of_range("Perturbation index out of range");
    }
    return perturbations_[index];
  }
};

}  // namespace metada::framework