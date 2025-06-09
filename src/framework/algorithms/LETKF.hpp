#pragma once
#include "Ensemble.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"

namespace metada::framework {

/**
 * @brief Local Ensemble Transform Kalman Filter (LETKF) implementation.
 *
 * This class implements the LETKF algorithm for ensemble-based data
 * assimilation. It updates an ensemble of states using observations and an
 * observation operator.
 *
 * @tparam BackendTag The backend tag type that must satisfy the required
 * concepts
 */
template <typename BackendTag>
class LETKF {
 public:
  /**
   * @brief Construct a LETKF object.
   * @param ensemble Reference to the ensemble to be updated.
   * @param obs Reference to the observation object.
   * @param obs_op Reference to the observation operator.
   * @param inflation Covariance inflation factor.
   */
  LETKF(Ensemble<BackendTag>& ensemble, const Observation<BackendTag>& obs,
        const ObsOperator<BackendTag>& obs_op, double inflation);

  /**
   * @brief Perform the LETKF analysis step, updating the ensemble.
   */
  void analyse();

 private:
  Ensemble<BackendTag>& ensemble_;
  const Observation<BackendTag>& obs_;
  const ObsOperator<BackendTag>& obs_op_;
  double inflation_;
};

}  // namespace metada::framework