#pragma once

#include <memory>
#include <stdexcept>
#include <vector>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "Logger.hpp"
#include "NonCopyable.hpp"
#include "State.hpp"
#include "StateConcepts.hpp"

namespace metada::framework {

/**
 * @brief Forward declaration of Config class
 */
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

/**
 * @brief Forward declaration of Geometry class
 */
template <typename BackendTag>
  requires GeometryBackendType<BackendTag>
class Geometry;

/**
 * @brief Adapter class for ensemble of states in data assimilation systems
 *
 * @details This class provides a type-safe interface for managing an ensemble
 * of State<BackendTag> objects. It supports ensemble operations such as mean
 * and perturbation computation, and is designed for use in ensemble-based data
 * assimilation methods like LETKF.
 *
 * @tparam BackendTag The backend tag type that must satisfy StateBackendType
 */
template <typename BackendTag>
  requires StateBackendType<BackendTag>
class Ensemble : public NonCopyable {
 public:
  using StateType = State<BackendTag>;

  /**
   * @brief Construct an ensemble of states with the given config and geometry
   * @param config Configuration object for state initialization
   * @param geometry Geometry object for state initialization
   */
  explicit Ensemble(const Config<BackendTag>& config,
                    const Geometry<BackendTag>& geometry)
      : config_(config),
        geometry_(geometry),
        members_(),
        size_(0),
        perturbations_() {
    logger_.Info() << "Ensemble starting construction";

    const auto member_configs = config.GetSubsectionsFromVector("members");
    size_ = member_configs.size();
    members_.reserve(size_);
    for (const auto& member_config : member_configs) {
      members_.emplace_back(member_config.GetSubsection("state"), geometry);
    }

    logger_.Info() << "Ensemble constructed with " << size_ << " members";
  }

  /**
   * @brief Get mutable access to an ensemble member
   * @param index Index of the member
   * @return Reference to the member state
   * @throws std::out_of_range if index is invalid
   */
  StateType& GetMember(size_t index) {
    if (index >= size_) {
      throw std::out_of_range("Ensemble member index out of range");
    }
    return members_[index];
  }

  /**
   * @brief Get const access to an ensemble member
   * @param index Index of the member
   * @return Const reference to the member state
   * @throws std::out_of_range if index is invalid
   */
  const StateType& GetMember(size_t index) const {
    if (index >= size_) {
      throw std::out_of_range("Ensemble member index out of range");
    }
    return members_[index];
  }

  /**
   * @brief Get the number of ensemble members
   * @return Ensemble size
   */
  size_t Size() const { return size_; }

  /**
   * @brief Force recomputation of the mean state
   */
  void RecomputeMean() {
    if (!mean_) {
      mean_ = std::make_unique<StateType>(members_[0].clone());
    }
    mean_->zero();
    for (size_t i = 0; i < size_; ++i) {
      *mean_ += members_[i];
    }
    *mean_ *= (1.0 / static_cast<double>(size_));
  }

  /**
   * @brief Get the mean state of the ensemble, computing it if necessary
   * @return Reference to the mean state
   */
  StateType& Mean() {
    if (!mean_) {
      RecomputeMean();
    }
    return *mean_;
  }

  /**
   * @brief Get const access to the mean state
   * @return Const reference to the mean state
   * @throws std::runtime_error if mean hasn't been computed
   */
  const StateType& Mean() const {
    if (!mean_) {
      throw std::runtime_error("Mean has not been computed");
    }
    return *mean_;
  }

  /**
   * @brief Force recomputation of all perturbations
   */
  void RecomputePerturbations() {
    RecomputeMean();
    perturbations_.clear();
    perturbations_.reserve(size_);
    for (size_t i = 0; i < size_; ++i) {
      perturbations_.push_back(
          std::make_unique<StateType>(members_[i].clone()));
      *perturbations_.back() -= *mean_;
    }
  }

  /**
   * @brief Get a perturbation, computing it if necessary
   * @param index Index of the perturbation
   * @return Reference to the perturbation state
   * @throws std::out_of_range if index is invalid
   */
  StateType& GetPerturbation(size_t index) {
    if (index >= size_) {
      throw std::out_of_range("Perturbation index out of range");
    }
    if (perturbations_.empty()) {
      RecomputePerturbations();
    }
    return *perturbations_[index];
  }

  /**
   * @brief Get const access to a perturbation
   * @param index Index of the perturbation
   * @return Const reference to the perturbation state
   * @throws std::out_of_range if index is invalid
   * @throws std::runtime_error if perturbations haven't been computed
   */
  const StateType& GetPerturbation(size_t index) const {
    if (index >= size_) {
      throw std::out_of_range("Perturbation index out of range");
    }
    if (perturbations_.empty()) {
      throw std::runtime_error("Perturbations have not been computed");
    }
    return *perturbations_[index];
  }

 private:
  const Config<BackendTag>& config_;
  const Geometry<BackendTag>& geometry_;
  std::vector<StateType> members_;
  std::unique_ptr<StateType> mean_;
  size_t size_;
  std::vector<std::unique_ptr<StateType>> perturbations_;
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework