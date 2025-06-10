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
      obs_op_.template apply(ensemble_.GetMember(i), obs_);
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

 private:
  Ensemble<BackendTag>& ensemble_;
  Observation<BackendTag>& obs_;
  const ObsOperator<BackendTag>& obs_op_;
  double inflation_;
};

}  // namespace metada::framework