#pragma once

#include <memory>
#include <vector>

#include "IState.hpp"

namespace metada {
namespace framework {
/**
 * @brief Interface for ensemble representations
 *
 * Defines core functionality needed for ensemble-based data assimilation:
 * - Ensemble member management
 * - Statistical operations (mean, perturbations)
 * - Ensemble transformations
 */
class IEnsemble {
 public:
  virtual ~IEnsemble() = default;

  // Ensemble management
  virtual void initialize(const tools::config::IConfig& config) = 0;
  virtual size_t getSize() const = 0;
  virtual void resize(size_t new_size) = 0;

  // Member access
  virtual IState& getMember(size_t index) = 0;
  virtual const IState& getMember(size_t index) const = 0;

  // Statistical operations
  virtual void computeMean() = 0;
  virtual const IState& getMean() const = 0;
  virtual void computePerturbations() = 0;
  virtual const IState& getPerturbation(size_t index) const = 0;

  // LETKF specific operations
  virtual void inflate(double factor) = 0;
  virtual void transform(
      const std::vector<std::vector<double>>& transform_matrix) = 0;
  virtual void localizeCovariance(
      const std::vector<double>& localization_weights) = 0;

  // Validation
  virtual void validate() const = 0;
  virtual bool isValid() const = 0;
};

}  // namespace framework
}  // namespace metada