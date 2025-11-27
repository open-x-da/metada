#pragma once
#include <Eigen/Dense>
#include <cmath>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "Logger.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"

namespace metada::framework {

/**
 * @brief Ensemble Transform Kalman Filter (ETKF) implementation.
 *
 * This class implements the ETKF algorithm for ensemble-based data
 * assimilation. It updates an ensemble of states using observations and an
 * observation operator.
 *
 * The ETKF algorithm follows these steps:
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
class ETKF {
 public:
  /**
   * @brief Construct a ETKF object.
   * @param ensemble Reference to the ensemble to be updated.
   * @param obs Reference to the observation object.
   * @param obs_op Reference to the observation operator.
   * @param config Configuration object containing inflation factor.
   */
  ETKF(Ensemble<BackendTag>& ensemble, Observation<BackendTag>& obs,
       const ObsOperator<BackendTag>& obs_op, const Config<BackendTag>& config)
      : ensemble_(ensemble),
        obs_(obs),
        obs_op_(obs_op),
        inflation_(config.Get("inflation").asFloat()),
        output_base_file_(config.Get("output_base_file").asString()),
        format_(config.Get("format").asString()) {
    logger_.Info() << "ETKF constructed";
  }

  /**
   * @brief Perform the ETKF analysis step, updating the ensemble.
   */
  void Analyse() {
    logger_.Info() << "ETKF analysis started";
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int ens_size = ensemble_.Size();
    const int state_dim = ensemble_.GetMember(0).size();
    const int obs_dim = obs_.size();

    // 1. Build forecast quantities using Ensemble methods
    // Compute mean and perturbations using Ensemble class
    ensemble_.RecomputeMean();
    ensemble_.RecomputePerturbations();

    const auto& mean_state = ensemble_.Mean();
    const auto mean_data = mean_state.template getDataPtr<double>();
    VectorXd xb_mean = Eigen::Map<const VectorXd>(mean_data, state_dim);

    // Build Xb_pert (state_dim x ens_size) from perturbations
    MatrixXd Xb_pert(state_dim, ens_size);
    for (int i = 0; i < ens_size; ++i) {
      const auto& pert_data =
          ensemble_.GetPerturbation(i).template getDataPtr<double>();
      Xb_pert.col(i) = Eigen::Map<const VectorXd>(pert_data, state_dim);
    }
    Xb_pert *= inflation_;  // Apply inflation

    // Propagate to obs space
    MatrixXd Yb(obs_dim, ens_size);
    for (int i = 0; i < ens_size; ++i) {
      const auto& obs_data = obs_op_.apply(ensemble_.GetMember(i), obs_);
      Yb.col(i) = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);
    }

    // Compute mean and anomalies in observation space
    VectorXd yb_mean = Yb.rowwise().mean();
    MatrixXd Yb_pert = Yb.colwise() - yb_mean;

    // 2. Compute innovation
    const auto obs_data = obs_.getObservationValues();
    VectorXd yo = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);
    VectorXd d = yo - yb_mean;

    // 3. Compute ensemble-space gain matrices
    const auto& R_data = obs_.getCovariance();
    // R_data is now a vector of variances, create diagonal matrix
    VectorXd R_diag = Eigen::Map<const VectorXd>(R_data.data(), obs_dim);
    MatrixXd R = R_diag.asDiagonal();

    // Compute analysis error in ensemble space
    MatrixXd Pa = (Yb_pert.transpose() * R.inverse() * Yb_pert +
                   (ens_size - 1) * MatrixXd::Identity(ens_size, ens_size))
                      .inverse();

    // Compute weights for the mean
    MatrixXd wa = Pa * Yb_pert.transpose() * R.inverse() * d;

    // Compute square-root for anomalies (deterministic ETKF with Q = Identity)
    MatrixXd Pa_sqrt = Pa.llt().matrixL();
    MatrixXd Wa = std::sqrt(ens_size - 1) * Pa_sqrt;

    // 4. Analysis update in state space
    // Update mean
    VectorXd xa_mean = xb_mean + Xb_pert * wa;

    // Update anomalies
    MatrixXd Xa_pert = Xb_pert * Wa;

    // Full analysis ensemble
    MatrixXd Xa = xa_mean.replicate(1, ens_size) + Xa_pert;

    // Update ensemble members
    for (int i = 0; i < ens_size; ++i) {
      auto& member = ensemble_.GetMember(i);
      auto data = member.template getDataPtr<double>();
      Eigen::Map<VectorXd>(data, state_dim) = Xa.col(i);
    }

    logger_.Info() << "ETKF analysis completed";
  }

  /**
   * @brief Save the analyzed ensemble to files
   * @param base_filename Base filename for saving (without extension)
   */
  void saveEnsemble() const {
    logger_.Info() << "ETKF saving ensemble";
    const int ens_size = ensemble_.Size();

    // Save analysis mean using the ensemble's mean state
    const auto& mean_state = ensemble_.Mean();
    mean_state.saveToFile(output_base_file_ + "_mean." + format_);

    // Save individual ensemble members
    for (int i = 0; i < ens_size; ++i) {
      const auto& member = ensemble_.GetMember(i);
      member.saveToFile(output_base_file_ + "_member_" + std::to_string(i) +
                        "." + format_);
    }
    logger_.Info() << "ETKF ensemble saved";
  }

 private:
  Ensemble<BackendTag>& ensemble_;
  Observation<BackendTag>& obs_;
  const ObsOperator<BackendTag>& obs_op_;
  double inflation_;
  std::string output_base_file_;
  std::string format_ = "txt";  // Default format
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework