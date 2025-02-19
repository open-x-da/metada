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
    double* mean_data = &mean_.template getData<double>();
    const auto dim = mean_.getDimensions()[0];

    // Sum all members
    for (const auto& member : members_) {
      const double* member_data = &member.template getData<double>();
      for (size_t i = 0; i < dim; ++i) {
        mean_data[i] += member_data[i];
      }
    }

    // Divide by ensemble size
    for (size_t i = 0; i < dim; ++i) {
      mean_data[i] /= static_cast<double>(size_);
    }
  }

  const T& getMean() const { return mean_; }

  void computePerturbations() {
    // Ensure mean is computed
    computeMean();

    // Resize perturbations if needed
    perturbations_.resize(size_, mean_);
    const auto dim = mean_.getDimensions()[0];

    // Compute perturbations for each member
    for (size_t i = 0; i < size_; ++i) {
      auto& pert_data = perturbations_[i].template getData<double>();
      const auto& member_data = members_[i].template getData<double>();
      const auto& mean_data = mean_.template getData<double>();

      for (size_t j = 0; j < dim; ++j) {
        pert_data[j] = member_data[j] - mean_data[j];
      }
    }
  }

  const T& getPerturbation(size_t index) const {
    if (index >= size_) {
      throw std::out_of_range("Perturbation index out of range");
    }
    return perturbations_[index];
  }
};

}  // namespace metada::framework