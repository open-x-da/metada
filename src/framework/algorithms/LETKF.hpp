#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

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
 * The LETKF algorithm follows these steps:
 * 1. Build forecast quantities:
 *    \[
 *    \bar{\mathbf x}^f = \frac{1}{N}\sum_{i=1}^{N}\mathbf x_i^f
 *    \]
 *    \[
 *    \mathbf X'^{\,f} = \bigl[\,\mathbf x_1^f-\bar{\mathbf
 * x}^f\,\;\cdots\;\mathbf x_N^f-\bar{\mathbf x}^f\bigr]
 *    \]
 *    \[
 *    \bar{\mathbf y}^f = \frac{1}{N}\sum_{i=1}^{N}\mathbf y_i^f
 *    \]
 *    \[
 *    \mathbf Y'^{\,f} = \bigl[\,\mathbf y_1^f-\bar{\mathbf
 * y}^f\,\;\cdots\;\mathbf y_N^f-\bar{\mathbf y}^f\bigr]
 *    \]
 *
 * 2. Compute innovation:
 *    \[
 *    \mathbf d = \mathbf y^{\text{obs}}-\bar{\mathbf y}^f
 *    \]
 *
 * 3. Compute ensemble-space gain matrices:
 *    \[
 *    \mathbf A^{-1} = \alpha (N-1)\mathbf I + \mathbf Y'^{f\,\! \top}\,\mathbf
 * R^{-1}\,\mathbf Y'^{\,f}
 *    \]
 *    \[
 *    \mathbf{P} = \mathbf{A}^{-1}
 *    \]
 *    \[
 *    \bar{\mathbf{w}}^a =
 * \mathbf{P}\,(\mathbf{Y}^{f})'^{\!\top}\mathbf{R}^{-1}\mathbf{d}
 *    \]
 *    \[
 *    \mathbf{W}^a = \sqrt{N-1}\;\mathbf{P}^{1/2}\,\mathbf{Q}
 *    \]
 *
 * 4. Analysis update in state space:
 *    \[
 *    \bar{\mathbf x}^a = \bar{\mathbf x}^f + \mathbf X'^{\,f}\,\bar{\mathbf
 * w}^a
 *    \]
 *    \[
 *    \mathbf X'^{\,a} = \mathbf X'^{\,f}\,\mathbf W^a
 *    \]
 *    \[
 *    \mathbf X^{\,a} = \bar{\mathbf x}^a\mathbf 1^{\top} + \mathbf X'^{\,a}
 *    \]
 *
 * For more details, see Hunt et al. (2007) "Efficient Data Assimilation for
 * Spatiotemporal Chaos: A Local Ensemble Transform Kalman Filter"
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
  LETKF(Ensemble<BackendTag>& ensemble, Observation<BackendTag>& obs,
        const ObsOperator<BackendTag>& obs_op, double inflation)
      : ensemble_(ensemble),
        obs_(obs),
        obs_op_(obs_op),
        inflation_(inflation) {}

  /**
   * @brief Perform the LETKF analysis step, updating the ensemble.
   */
  void analyse() {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int ens_size = ensemble_.Size();
    const int state_dim =
        ensemble_.GetMember(0).template getData<std::vector<double>>().size();
    const int obs_dim = obs_.template getData<std::vector<double>>().size();

    // 1. Build forecast quantities
    // Build Xb (state_dim x ens_size)
    MatrixXd Xb(state_dim, ens_size);
    for (int i = 0; i < ens_size; ++i) {
      const auto& member_data =
          ensemble_.GetMember(i).template getData<std::vector<double>>();
      Xb.col(i) = Eigen::Map<const VectorXd>(member_data.data(), state_dim);
    }

    // Compute mean and anomalies
    VectorXd xb_mean = Xb.rowwise().mean();
    MatrixXd Xb_pert = Xb.colwise() - xb_mean;
    Xb_pert *= inflation_;  // Apply inflation
    std::cout << "Xb_pert: " << Xb_pert << std::endl;

    // Propagate to obs space
    MatrixXd Yb(obs_dim, ens_size);
    for (int i = 0; i < ens_size; ++i) {
      const auto& obs_data =
          obs_op_.template apply(ensemble_.GetMember(i), obs_);
      Yb.col(i) = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);
    }
    std::cout << "Yb: " << Yb << std::endl;

    // Compute mean and anomalies in observation space
    VectorXd yb_mean = Yb.rowwise().mean();
    MatrixXd Yb_pert = Yb.colwise() - yb_mean;
    std::cout << "Yb_pert: " << Yb_pert << std::endl;

    // 2. Compute innovation
    const auto& obs_data = obs_.template getData<std::vector<double>>();
    VectorXd yo = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);
    VectorXd d = yo - yb_mean;
    std::cout << "yo: " << yo << std::endl;
    std::cout << "yb_mean: " << yb_mean << std::endl;
    std::cout << "d: " << d << std::endl;

    // 3. Compute ensemble-space gain matrices
    const auto& R_data = obs_.getCovariance();
    MatrixXd R = Eigen::Map<const MatrixXd>(R_data.data(), obs_dim, obs_dim);

    // Compute analysis error in ensemble space
    MatrixXd Pa = (Yb_pert.transpose() * R.inverse() * Yb_pert +
                   (ens_size - 1) * MatrixXd::Identity(ens_size, ens_size))
                      .inverse();
    std::cout << "Pa: " << Pa << std::endl;

    // Compute weights for the mean
    MatrixXd wa = Pa * Yb_pert.transpose() * R.inverse() * d;
    std::cout << "wa: " << wa << std::endl;

    // Compute square-root for anomalies (deterministic LETKF with Q = Identity)
    MatrixXd Pa_sqrt = Pa.llt().matrixL();
    MatrixXd Wa = std::sqrt(ens_size - 1) * Pa_sqrt;
    std::cout << "Wa: " << Wa << std::endl;

    // 4. Analysis update in state space
    // Update mean
    VectorXd xa_mean = xb_mean + Xb_pert * wa;
    std::cout << "xa_mean: " << xa_mean << std::endl;

    // Update anomalies
    MatrixXd Xa_pert = Xb_pert * Wa;
    std::cout << "Xa_pert: " << Xa_pert << std::endl;

    // Full analysis ensemble
    MatrixXd Xa = xa_mean.replicate(1, ens_size) + Xa_pert;
    std::cout << "Xa: " << Xa << std::endl;

    // Update ensemble members
    for (int i = 0; i < ens_size; ++i) {
      auto& member = ensemble_.GetMember(i);
      auto& data = member.template getData<std::vector<double>>();
      Eigen::Map<VectorXd>(data.data(), state_dim) = Xa.col(i);
    }
  }

 private:
  Ensemble<BackendTag>& ensemble_;
  Observation<BackendTag>& obs_;
  const ObsOperator<BackendTag>& obs_op_;
  double inflation_;
};

}  // namespace metada::framework