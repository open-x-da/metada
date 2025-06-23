#pragma once

#include <memory>
#include <stdexcept>
#include <vector>

#include "SimpleGeometry.hpp"
#include "SimpleState.hpp"

namespace metada::backends::simple {

class SimpleEnsemble {
 public:
  // Prevent default construction and copy operations
  SimpleEnsemble() = delete;
  SimpleEnsemble(const SimpleEnsemble&) = delete;
  SimpleEnsemble& operator=(const SimpleEnsemble&) = delete;

  // Move constructor/assignment
  SimpleEnsemble(SimpleEnsemble&& other) noexcept
      : geometry_(other.geometry_),
        members_(std::move(other.members_)),
        mean_(std::move(other.mean_)),
        size_(other.size_),
        perturbations_(std::move(other.perturbations_)) {}

  SimpleEnsemble& operator=(SimpleEnsemble&& other) noexcept {
    if (this != &other) {
      // geometry_ is const reference, so we can't move it
      // members_, mean_, and perturbations_ can be moved
      members_ = std::move(other.members_);
      mean_ = std::move(other.mean_);
      size_ = other.size_;
      perturbations_ = std::move(other.perturbations_);
    }
    return *this;
  }

  ~SimpleEnsemble() = default;

  // Construct from config and geometry
  template <typename ConfigBackend>
  SimpleEnsemble(const ConfigBackend& config, const SimpleGeometry& geometry)
      : geometry_(geometry) {
    // Get state array from config
    const auto& state_vec = config.Get("state").asArray();
    size_ = state_vec.size();
    members_.reserve(size_);

    // Initialize each member from its state config
    for (size_t i = 0; i < size_; ++i) {
      members_.emplace_back(state_vec[i], geometry_);
    }
  }

  // Member access
  SimpleState& GetMember(size_t index) {
    if (index >= size_) {
      throw std::out_of_range("Ensemble member index out of range");
    }
    return members_[index];
  }

  const SimpleState& GetMember(size_t index) const {
    if (index >= size_) {
      throw std::out_of_range("Ensemble member index out of range");
    }
    return members_[index];
  }

  size_t Size() const { return size_; }

  // Ensemble operations
  void ComputeMean() {
    mean_ = members_[0].clone();  // Clone first member
    mean_->zero();
    for (size_t i = 0; i < size_; ++i) {
      mean_->add(members_[i]);
    }
    // Scale by 1/size
    for (auto& value : *mean_) {
      value.second *= (1.0 / static_cast<double>(size_));
    }
  }

  SimpleState& Mean() { return *mean_; }
  const SimpleState& Mean() const { return *mean_; }

  void ComputePerturbations() {
    ComputeMean();
    perturbations_.clear();
    perturbations_.reserve(size_);
    for (size_t i = 0; i < size_; ++i) {
      perturbations_.push_back(members_[i].clone());
      perturbations_.back()->add(*mean_);  // member - mean
    }
  }

  SimpleState& GetPerturbation(size_t index) {
    if (index >= size_) {
      throw std::out_of_range("Perturbation index out of range");
    }
    return *perturbations_[index];
  }

  const SimpleState& GetPerturbation(size_t index) const {
    if (index >= size_) {
      throw std::out_of_range("Perturbation index out of range");
    }
    return *perturbations_[index];
  }

 private:
  const SimpleGeometry& geometry_;
  std::vector<SimpleState> members_;
  std::unique_ptr<SimpleState> mean_;
  size_t size_;
  std::vector<std::unique_ptr<SimpleState>> perturbations_;
};

}  // namespace metada::backends::simple