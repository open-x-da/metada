#include "LETKF.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

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
                         const Observation<BackendTag>& obs,
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

  const int ens_size = ensemble_.size();
  // You may need to add state_size() and size() methods to Ensemble/Observation
  const int state_dim = ensemble_.member(0).asVector().size();
  const int obs_dim = obs_.asVector().size();

  // 1. Build Xb (state_dim x ens_size)
  MatrixXd Xb(state_dim, ens_size);
  for (int i = 0; i < ens_size; ++i)
    Xb.col(i) =
        ensemble_.member(i).asVector();  // asVector() returns Eigen::VectorXd

  // 2. Compute mean and anomalies
  VectorXd xb_mean = Xb.rowwise().mean();
  MatrixXd Xb_pert = Xb.colwise() - xb_mean;
  Xb_pert *= inflation_;  // Apply inflation

  // 3. Propagate to obs space
  MatrixXd Yb(obs_dim, ens_size);
  for (int i = 0; i < ens_size; ++i)
    Yb.col(i) = obs_op_.apply(ensemble_.member(i));  // returns Eigen::VectorXd

  VectorXd yb_mean = Yb.rowwise().mean();
  MatrixXd Yb_pert = Yb.colwise() - yb_mean;

  // 4. Get observation vector and R
  VectorXd yo = obs_.asVector();   // asVector() returns Eigen::VectorXd
  MatrixXd R = obs_.covariance();  // returns Eigen::MatrixXd

  // 5. Compute analysis weights in ensemble space
  MatrixXd Rinv = R.inverse();
  MatrixXd Pe = ((ens_size - 1) * MatrixXd::Identity(ens_size, ens_size) +
                 Yb_pert.transpose() * Rinv * Yb_pert)
                    .inverse();
  VectorXd wa = Pe * (Yb_pert.transpose() * Rinv * (yo - yb_mean));
  MatrixXd Wa = Pe.llt().matrixL();  // Lower-triangular Cholesky (sqrt)

  // 6. Update ensemble
  MatrixXd Xa(state_dim, ens_size);
  for (int i = 0; i < ens_size; ++i)
    Xa.col(i) = xb_mean + Xb_pert * Wa.col(i) + Xb_pert * wa;

  // 7. Write back to ensemble
  for (int i = 0; i < ens_size; ++i)
    ensemble_.member(i).fromVector(
        Xa.col(i));  // fromVector() sets state from Eigen::VectorXd
}

// Explicit template instantiation for common backend(s) can be added here if
// needed

}  // namespace metada::framework