#pragma once

#include <stdexcept>
#include <vector>

#include "IEnsemble.hpp"
#include "State.hpp"

namespace metada {
namespace framework {

/**
 * @brief Template class for ensemble representation
 *
 * Provides implementation for ensemble operations using a specific State type.
 * Designed for use in ensemble-based data assimilation methods like LETKF.
 *
 * @tparam StateBackend The type of State used in the ensemble
 */
template <typename StateBackend>
class Ensemble {
 private:
  std::vector<State<StateBackend>> members_;
  State<StateBackend> mean_;
  size_t ensemble_size_;
  std::vector<State<StateBackend>> perturbations_;

  // Helper method for matrix inversion (needs implementation)
  std::vector<std::vector<double>> invertMatrix(
      const std::vector<std::vector<double>>& matrix) {
    // Implement matrix inversion (e.g., using LU decomposition)
    // Return inverted matrix
    return std::vector<std::vector<double>>();
  }

 public:
  // Constructor with config
  template <typename T>
  explicit Ensemble(const tools::config::Config<T>& config,
                    size_t ensemble_size)
      : members_(),
        mean_(config),
        ensemble_size_(ensemble_size),
        perturbations_() {
    // Initialize members
    members_.reserve(ensemble_size);
    for (size_t i = 0; i < ensemble_size; ++i) {
      members_.emplace_back(config);
    }
  }

  // Member access
  State<StateBackend>& getMember(size_t index) {
    if (index >= ensemble_size_) {
      throw std::out_of_range("Ensemble member index out of range");
    }
    return members_[index];
  }

  const State<StateBackend>& getMember(size_t index) const {
    if (index >= ensemble_size_) {
      throw std::out_of_range("Ensemble member index out of range");
    }
    return members_[index];
  }

  size_t getSize() const { return ensemble_size_; }

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
      mean_data[i] /= static_cast<double>(ensemble_size_);
    }
  }

  const State<StateBackend>& getMean() const { return mean_; }

  void computePerturbations() {
    // Ensure mean is computed
    computeMean();

    // Resize perturbations if needed
    perturbations_.resize(ensemble_size_, mean_);
    const auto dim = mean_.getDimensions()[0];

    // Compute perturbations for each member
    for (size_t i = 0; i < ensemble_size_; ++i) {
      auto& pert_data = perturbations_[i].template getData<double>();
      const auto& member_data = members_[i].template getData<double>();
      const auto& mean_data = mean_.template getData<double>();

      for (size_t j = 0; j < dim; ++j) {
        pert_data[j] = member_data[j] - mean_data[j];
      }
    }
  }

  const State<StateBackend>& getPerturbation(size_t index) const {
    if (index >= ensemble_size_) {
      throw std::out_of_range("Perturbation index out of range");
    }
    return perturbations_[index];
  }

  // LETKF operations
  void inflate(double factor) {
    if (factor <= 0.0) {
      throw std::invalid_argument("Inflation factor must be positive");
    }

    // Compute mean and perturbations
    computeMean();
    computePerturbations();

    // Apply inflation to perturbations and update members
    for (size_t i = 0; i < ensemble_size_; ++i) {
      auto& member_data = members_[i].template getData<double>();
      const auto& pert_data = perturbations_[i].template getData<double>();
      const auto& mean_data = mean_.template getData<double>();

      for (size_t j = 0; j < mean_.getDimensions()[0]; ++j) {
        member_data[j] = mean_data[j] + factor * pert_data[j];
      }
    }
  }

  void transform(const std::vector<std::vector<double>>& transform_matrix) {
    if (transform_matrix.size() != ensemble_size_ ||
        transform_matrix[0].size() != ensemble_size_) {
      throw std::invalid_argument("Invalid transform matrix dimensions");
    }

    // Create temporary storage for transformed members
    std::vector<State<StateBackend>> transformed_members = members_;

    // Apply transformation
    for (size_t i = 0; i < ensemble_size_; ++i) {
      auto& new_member_data = transformed_members[i].template getData<double>();
      std::fill(new_member_data, new_member_data + mean_.getDimensions()[0],
                0.0);

      for (size_t j = 0; j < ensemble_size_; ++j) {
        const auto& member_data = members_[j].template getData<double>();
        for (size_t k = 0; k < mean_.getDimensions()[0]; ++k) {
          new_member_data[k] += transform_matrix[i][j] * member_data[k];
        }
      }
    }

    // Update members
    members_ = std::move(transformed_members);
  }

  void localizeCovariance(const std::vector<double>& localization_weights) {
    if (localization_weights.size() != mean_.getDimensions()[0]) {
      throw std::invalid_argument("Invalid localization weights dimension");
    }

    // Apply localization to perturbations
    computePerturbations();
    for (auto& pert : perturbations_) {
      auto& pert_data = pert.template getData<double>();
      for (size_t i = 0; i < mean_.getDimensions()[0]; ++i) {
        pert_data[i] *= std::sqrt(localization_weights[i]);
      }
    }

    // Update members using localized perturbations
    for (size_t i = 0; i < ensemble_size_; ++i) {
      auto& member_data = members_[i].template getData<double>();
      const auto& pert_data = perturbations_[i].template getData<double>();
      const auto& mean_data = mean_.template getData<double>();

      for (size_t j = 0; j < mean_.getDimensions()[0]; ++j) {
        member_data[j] = mean_data[j] + pert_data[j];
      }
    }
  }

  // Compute ensemble covariance matrix in observation space
  std::vector<std::vector<double>> computeObservationCovariance(
      const std::vector<double>& obs_perturbations) const {
    if (obs_perturbations.size() != ensemble_size_) {
      throw std::invalid_argument("Invalid observation perturbation size");
    }

    std::vector<std::vector<double>> cov(ensemble_size_,
                                         std::vector<double>(ensemble_size_));

    // Compute (HX')^T(HX')
    for (size_t i = 0; i < ensemble_size_; ++i) {
      for (size_t j = 0; j < ensemble_size_; ++j) {
        cov[i][j] = obs_perturbations[i] * obs_perturbations[j] /
                    (ensemble_size_ - 1.0);
      }
    }
    return cov;
  }

  // Compute analysis weights for LETKF
  std::vector<std::vector<double>> computeAnalysisWeights(
      const std::vector<double>& innovations,
      const std::vector<double>& obs_errors) {
    if (innovations.size() != ensemble_size_ ||
        obs_errors.size() != ensemble_size_) {
      throw std::invalid_argument(
          "Invalid innovation or observation error size");
    }

    // Compute (HPbH^T + R)^(-1)
    auto cov = computeObservationCovariance(innovations);
    for (size_t i = 0; i < ensemble_size_; ++i) {
      cov[i][i] += obs_errors[i];
    }

    // Compute analysis weights using matrix inversion
    return invertMatrix(cov);  // Need to implement matrix inversion
  }

  // Update ensemble using LETKF weights
  void updateWithWeights(const std::vector<std::vector<double>>& weights) {
    if (weights.size() != ensemble_size_ ||
        weights[0].size() != ensemble_size_) {
      throw std::invalid_argument("Invalid weights matrix dimensions");
    }

    // Apply weights to update ensemble members
    std::vector<State<StateBackend>> updated_members = members_;
    for (size_t i = 0; i < ensemble_size_; ++i) {
      auto& new_member = updated_members[i];
      new_member.reset();

      for (size_t j = 0; j < ensemble_size_; ++j) {
        const auto& member = members_[j];
        const auto& member_data = member.template getData<double>();
        auto& new_data = new_member.template getData<double>();

        for (size_t k = 0; k < new_member.getDimensions()[0]; ++k) {
          new_data[k] += weights[i][j] * member_data[k];
        }
      }
    }

    members_ = std::move(updated_members);
  }
};

}  // namespace framework
}  // namespace metada