#include "LETKF.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

#include "MockBackendTraits.hpp"

namespace metada::framework {

/*
 * LETKF Algorithm: Step-by-Step Documentation
 *
 * The Local Ensemble Transform Kalman Filter (LETKF) is an efficient
 * ensemble-based data assimilation algorithm. This documentation outlines the
 * serial LETKF algorithm as implemented in this file, including the
 * mathematical formulas.
 *
 * Notation:
 *   - X^b: Background ensemble matrix (state_dim × ens_size)
 *   - Y^b: Background ensemble in observation space (obs_dim × ens_size)
 *   - y^o: Observation vector (obs_dim)
 *   - R: Observation error covariance (obs_dim × obs_dim)
 *   - ρ: Inflation factor
 *
 * Step-by-Step Implementation:
 *
 * 1. Extract Ensemble Matrix
 *    - X^b = [x_1^b, x_2^b, ..., x_k^b] ∈ ℝ^{state_dim × ens_size}
 *
 * 2. Compute Ensemble Mean and Anomalies
 *    - \bar{x}^b = (1/k) ∑_{i=1}^k x_i^b
 *    - X' = X^b - \bar{x}^b (column-wise)
 *    - Apply inflation: X' ← ρ X'
 *
 * 3. Propagate Ensemble to Observation Space
 *    - For each member: y_i^b = H(x_i^b)
 *    - Y^b = [y_1^b, ..., y_k^b] ∈ ℝ^{obs_dim × ens_size}
 *    - \bar{y}^b = (1/k) ∑_{i=1}^k y_i^b
 *    - Y' = Y^b - \bar{y}^b (column-wise)
 *
 * 4. Compute Analysis Weights in Ensemble Space
 *    - P_e = [ (k-1)I + (Y')^T R^{-1} Y' ]^{-1}
 *    - w^a = P_e (Y')^T R^{-1} (y^o - \bar{y}^b)
 *    - W^a = sqrtm[ (k-1) P_e ] (matrix square root, e.g., via Cholesky)
 *
 * 5. Update Ensemble
 *    - For each member j:
 *        x_j^a = \bar{x}^b + X' W^a_{:,j} + X' w^a
 *      or, in matrix form:
 *        X^a = \bar{x}^b + X' W^a + X' w^a
 *
 * 6. Write Back to Ensemble
 *    - Update each ensemble member with the new analysis state.
 *
 * References:
 *   - Hunt et al. (2007), "Efficient Data Assimilation for Spatiotemporal
 * Chaos: A Local Ensemble Transform Kalman Filter"
 *   -
 * https://en.wikipedia.org/wiki/Ensemble_Kalman_filter#Local_ensemble_transform_Kalman_filter_(LETKF)
 */

/**
 * @brief LETKF constructor
 * @param ensemble Reference to the ensemble to be updated
 * @param obs Reference to the observation object
 * @param obs_op Reference to the observation operator
 * @param inflation Covariance inflation factor
 */
template <typename BackendTag>
LETKF<BackendTag>::LETKF(Ensemble<BackendTag>& ensemble,
                         Observation<BackendTag>& obs,
                         const ObsOperator<BackendTag>& obs_op,
                         double inflation)
    : ensemble_(ensemble), obs_(obs), obs_op_(obs_op), inflation_(inflation) {}

/**
 * @brief Perform the LETKF analysis step, updating the ensemble.
 */
template <typename BackendTag>
void LETKF<BackendTag>::analyse() {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  const int ens_size = ensemble_.Size();
  const int state_dim =
      ensemble_.GetMember(0).template getData<std::vector<double>>().size();
  const int obs_dim = obs_.template getData<std::vector<double>>().size();

  // 1. Build Xb (state_dim x ens_size)
  MatrixXd Xb(state_dim, ens_size);
  for (int i = 0; i < ens_size; ++i) {
    const auto& member_data =
        ensemble_.GetMember(i).template getData<std::vector<double>>();
    Xb.col(i) = Eigen::Map<const VectorXd>(member_data.data(), state_dim);
  }

  // 2. Compute mean and anomalies
  VectorXd xb_mean = Xb.rowwise().mean();
  MatrixXd Xb_pert = Xb.colwise() - xb_mean;
  Xb_pert *= inflation_;  // Apply inflation

  // 3. Propagate to obs space
  MatrixXd Yb(obs_dim, ens_size);
  for (int i = 0; i < ens_size; ++i) {
    obs_op_.apply(ensemble_.GetMember(i), obs_);
    const auto& member_data = obs_.template getData<std::vector<double>>();
    Yb.col(i) = Eigen::Map<const VectorXd>(member_data.data(), obs_dim);
  }

  VectorXd yb_mean = Yb.rowwise().mean();
  MatrixXd Yb_pert = Yb.colwise() - yb_mean;

  // 4. Get observation vector and R
  const auto& obs_data = obs_.template getData<std::vector<double>>();
  VectorXd yo = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);
  const auto& R_data = obs_.getCovariance();
  MatrixXd R = Eigen::Map<const MatrixXd>(R_data.data(), obs_dim, obs_dim);

  // 5. Compute innovation
  VectorXd d = yo - yb_mean;

  // 6. Compute analysis weights
  MatrixXd Pa = (Yb_pert.transpose() * R.inverse() * Yb_pert +
                 (ens_size - 1) * MatrixXd::Identity(ens_size, ens_size))
                    .inverse();
  MatrixXd Wa = Pa * Yb_pert.transpose() * R.inverse() * d;

  // 7. Update ensemble
  MatrixXd Xa = Xb + Xb_pert * Wa;
  for (int i = 0; i < ens_size; ++i) {
    auto& member = ensemble_.GetMember(i);
    auto& data = member.template getData<std::vector<double>>();
    Eigen::Map<VectorXd>(data.data(), state_dim) = Xa.col(i);
  }
}

// Explicit template instantiations
template class LETKF<traits::MockBackendTag>;

}  // namespace metada::framework