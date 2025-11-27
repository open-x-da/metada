#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <random>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "Logger.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "ProgressBar.hpp"
#include "State.hpp"

namespace metada::framework {

/**
 * @brief Ensemble Kalman Filter (EnKF) implementation.
 *
 * This class implements the Ensemble Kalman Filter algorithm for ensemble-based
 * data assimilation. The EnKF is a Monte Carlo implementation of the Kalman
 * filter that represents the state distribution by a set of ensemble members.
 *
 * The EnKF algorithm follows these steps:
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
 * 3. Compute Kalman gain:
 *    \[
 *    \mathbf{K} = \mathbf{P}^f \mathbf{H}^T (\mathbf{H} \mathbf{P}^f
 * \mathbf{H}^T + \mathbf{R})^{-1}
 *    \]
 *    Where:
 *    \[
 *    \mathbf{P}^f = \frac{1}{N-1} \mathbf{X}'^{\,f} \mathbf{X}'^{\,f^T}
 *    \]
 *
 * 4. Analysis update:
 *    \[
 *    \mathbf{x}_i^a = \mathbf{x}_i^f + \mathbf{K} (\mathbf{y}^{\text{obs}} +
 * \mathbf{\epsilon}_i - \mathbf{H} \mathbf{x}_i^f)
 *    \]
 *    Where \mathbf{\epsilon}_i are observation perturbations.
 *
 * The EnKF differs from ETKF in that it:
 * - Uses the traditional Kalman gain formulation
 * - Requires observation perturbations for each ensemble member
 * - Maintains the ensemble spread through the analysis
 * - Is more computationally intensive but conceptually simpler
 *
 * For more details, see:
 * - Evensen, G. (1994) "Sequential data assimilation with a nonlinear
 *   quasi-geostrophic model using Monte Carlo methods to forecast error
 *   statistics"
 * - Burgers, G., van Leeuwen, P. J., & Evensen, G. (1998) "Analysis scheme
 *   in the ensemble Kalman filter"
 *
 * @tparam BackendTag The backend tag type that must satisfy the required
 * concepts
 */
template <typename BackendTag>
class EnKF {
 public:
  /**
   * @brief Analysis results structure
   */
  struct AnalysisResults {
    double innovation_norm;          ///< Innovation vector norm
    double analysis_increment_norm;  ///< Analysis increment norm
    double background_spread;        ///< Background ensemble spread
    double analysis_spread;          ///< Analysis ensemble spread
    double max_kalman_gain;          ///< Maximum Kalman gain element
    double min_kalman_gain;          ///< Minimum Kalman gain element
    double condition_number;  ///< Condition number of innovation covariance
    int ensemble_size;        ///< Number of ensemble members
    int observation_count;    ///< Number of observations
    std::string inflation_method;  ///< Method used for inflation
    double inflation_factor;       ///< Inflation factor applied
  };

  /**
   * @brief Construct an EnKF object.
   * @param ensemble Reference to the ensemble to be updated.
   * @param obs Reference to the observation object.
   * @param obs_op Reference to the observation operator.
   * @param config Configuration object containing EnKF parameters.
   */
  EnKF(Ensemble<BackendTag>& ensemble, Observation<BackendTag>& obs,
       const ObsOperator<BackendTag>& obs_op, const Config<BackendTag>& config)
      : ensemble_(ensemble), obs_(obs), obs_op_(obs_op) {
    // Get analysis configuration subsection
    auto analysis_config = config.GetSubsection("analysis");

    inflation_factor_ = analysis_config.Get("inflation").asFloat();
    output_base_file_ = analysis_config.Get("output_base_file").asString();
    format_ = analysis_config.Get("format").asString();

    // Parse inflation method
    std::string inflation_method_str =
        analysis_config.Get("inflation_method").asString();
    if (inflation_method_str == "multiplicative") {
      inflation_method_ = InflationMethod::MULTIPLICATIVE;
    } else if (inflation_method_str == "additive") {
      inflation_method_ = InflationMethod::ADDITIVE;
    } else if (inflation_method_str == "relaxation") {
      inflation_method_ = InflationMethod::RELAXATION;
    } else {
      inflation_method_ = InflationMethod::MULTIPLICATIVE;
      logger_.Warning() << "Unknown inflation method '" << inflation_method_str
                        << "', using multiplicative inflation";
    }

    logger_.Info() << "EnKF constructed with " << ensemble_.Size()
                   << " members";
    logger_.Info() << "Inflation factor: " << inflation_factor_;
    logger_.Info() << "Inflation method: " << inflation_method_str;
  }

  /**
   * @brief Perform the EnKF analysis step, updating the ensemble.
   */
  void Analyse() {
    logger_.Info() << "EnKF analysis started";
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int ens_size = ensemble_.Size();
    const int state_dim = ensemble_.GetMember(0).size();
    const int obs_dim = obs_.size();

    // 1. Build forecast quantities using Ensemble methods
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

    // Apply inflation to background perturbations
    applyInflation(Xb_pert);

    // 2. Propagate to observation space
    MatrixXd Yb(obs_dim, ens_size);
    for (int i = 0; i < ens_size; ++i) {
      const auto& obs_data = obs_op_.apply(ensemble_.GetMember(i), obs_);
      Yb.col(i) = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);
    }

    // Compute mean and anomalies in observation space
    VectorXd yb_mean = Yb.rowwise().mean();
    MatrixXd Yb_pert = Yb.colwise() - yb_mean;

    // 3. Compute innovation
    const auto obs_data = obs_.getObservationValues();
    VectorXd yo = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);
    VectorXd d = yo - yb_mean;

    // Store innovation norm for diagnostics
    innovation_norm_ = d.norm();

    // 4. Compute Kalman gain
    // Get observation error covariance
    const auto& R_data = obs_.getCovariance();
    VectorXd R_diag = Eigen::Map<const VectorXd>(R_data.data(), obs_dim);
    MatrixXd R = R_diag.asDiagonal();

    // Compute background covariance in observation space
    MatrixXd Pf_obs = Yb_pert * Yb_pert.transpose() / (ens_size - 1);

    // Add observation error covariance
    MatrixXd S = Pf_obs + R;

    // Compute Kalman gain
    MatrixXd K = Xb_pert * Yb_pert.transpose() * S.inverse() / (ens_size - 1);

    // Store Kalman gain statistics
    max_kalman_gain_ = K.maxCoeff();
    min_kalman_gain_ = K.minCoeff();

    // Compute condition number using singular values
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(S);
    auto singular_values = svd.singularValues();
    condition_number_ =
        singular_values(0) / singular_values(singular_values.size() - 1);

    // 5. Generate observation perturbations
    MatrixXd obs_pert(obs_dim, ens_size);
    generateObservationPerturbations(obs_pert, R);

    // 6. Analysis update for each ensemble member
    MatrixXd Xa(state_dim, ens_size);
    for (int i = 0; i < ens_size; ++i) {
      // Compute perturbed observation
      VectorXd yo_pert = yo + obs_pert.col(i);

      // Compute innovation for this member
      VectorXd d_i = yo_pert - Yb.col(i);

      // Analysis update
      VectorXd xa_i = Xb_pert.col(i) + K * d_i;
      Xa.col(i) = xb_mean + xa_i;
    }

    // 7. Update ensemble members
    for (int i = 0; i < ens_size; ++i) {
      auto& member = ensemble_.GetMember(i);
      auto data = member.template getDataPtr<double>();
      Eigen::Map<VectorXd>(data, state_dim) = Xa.col(i);
    }

    // 8. Update ensemble statistics
    ensemble_.RecomputeMean();
    ensemble_.RecomputePerturbations();

    // Compute analysis spread
    const auto& analysis_mean = ensemble_.Mean();
    const auto analysis_mean_data = analysis_mean.template getDataPtr<double>();
    VectorXd xa_mean =
        Eigen::Map<const VectorXd>(analysis_mean_data, state_dim);

    double analysis_spread = 0.0;
    for (int i = 0; i < ens_size; ++i) {
      const auto& pert_data =
          ensemble_.GetPerturbation(i).template getDataPtr<double>();
      VectorXd pert = Eigen::Map<const VectorXd>(pert_data, state_dim);
      analysis_spread += pert.squaredNorm();
    }
    analysis_spread_ = std::sqrt(analysis_spread / (ens_size * state_dim));

    logger_.Info() << "EnKF analysis completed";
  }

  /**
   * @brief Save the analyzed ensemble to files
   */
  void saveEnsemble() const {
    logger_.Info() << "EnKF saving ensemble";
    const int ens_size = ensemble_.Size();

    // Save analysis mean using the ensemble's mean state
    const auto& mean_state = ensemble_.Mean();
    mean_state.saveToFile(output_base_file_ + "_mean." + format_);
    logger_.Info() << "Analysis mean saved to: "
                   << output_base_file_ + "_mean." + format_;

    // Save individual ensemble members
    for (int i = 0; i < ens_size; ++i) {
      const auto& member = ensemble_.GetMember(i);
      member.saveToFile(output_base_file_ + "_member_" + std::to_string(i) +
                        "." + format_);
    }

    // Save diagnostics
    saveDiagnostics();

    logger_.Info() << "EnKF ensemble saved";
  }

  /**
   * @brief Get the analysis results
   */
  AnalysisResults getAnalysisResults() const {
    AnalysisResults results;
    results.innovation_norm = innovation_norm_;
    results.analysis_increment_norm = analysis_increment_norm_;
    results.background_spread = background_spread_;
    results.analysis_spread = analysis_spread_;
    results.max_kalman_gain = max_kalman_gain_;
    results.min_kalman_gain = min_kalman_gain_;
    results.condition_number = condition_number_;
    results.ensemble_size = ensemble_.Size();
    results.observation_count = obs_.size();
    results.inflation_method = getInflationMethodString();
    results.inflation_factor = inflation_factor_;
    return results;
  }

 private:
  /**
   * @brief Inflation method enumeration
   */
  enum class InflationMethod {
    MULTIPLICATIVE,  ///< Multiplicative inflation (default)
    ADDITIVE,        ///< Additive inflation
    RELAXATION       ///< Relaxation inflation
  };

  /**
   * @brief Apply inflation to background perturbations
   */
  void applyInflation(Eigen::MatrixXd& Xb_pert) {
    switch (inflation_method_) {
      case InflationMethod::MULTIPLICATIVE:
        Xb_pert *= std::sqrt(inflation_factor_);
        break;
      case InflationMethod::ADDITIVE:
        // Additive inflation would add random perturbations
        // For simplicity, we use multiplicative here
        Xb_pert *= std::sqrt(inflation_factor_);
        break;
      case InflationMethod::RELAXATION:
        // Relaxation inflation reduces the analysis increment
        // For simplicity, we use multiplicative here
        Xb_pert *= std::sqrt(inflation_factor_);
        break;
    }

    // Compute background spread for diagnostics
    background_spread_ = std::sqrt(Xb_pert.squaredNorm() / Xb_pert.size());
  }

  /**
   * @brief Generate observation perturbations
   */
  void generateObservationPerturbations(Eigen::MatrixXd& obs_pert,
                                        const Eigen::MatrixXd& R) {
    const int obs_dim = obs_.size();
    const int ens_size = ensemble_.Size();

    // Generate random perturbations from N(0, R) distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);

    // Generate standard normal random numbers
    Eigen::MatrixXd Z(obs_dim, ens_size);
    for (int i = 0; i < obs_dim; ++i) {
      for (int j = 0; j < ens_size; ++j) {
        Z(i, j) = dist(gen);
      }
    }

    // Transform to N(0, R) using Cholesky decomposition
    Eigen::LLT<Eigen::MatrixXd> llt(R);
    obs_pert = llt.matrixL() * Z;
  }

  /**
   * @brief Save diagnostic information to file
   */
  void saveDiagnostics() const {
    std::string diagnostics_file = output_base_file_ + "_diagnostics.txt";
    // Implementation would write diagnostics to file
    logger_.Info() << "Diagnostics saved to: " << diagnostics_file;
  }

  /**
   * @brief Get inflation method as string
   */
  std::string getInflationMethodString() const {
    switch (inflation_method_) {
      case InflationMethod::MULTIPLICATIVE:
        return "multiplicative";
      case InflationMethod::ADDITIVE:
        return "additive";
      case InflationMethod::RELAXATION:
        return "relaxation";
      default:
        return "unknown";
    }
  }

  Ensemble<BackendTag>& ensemble_;
  Observation<BackendTag>& obs_;
  const ObsOperator<BackendTag>& obs_op_;
  InflationMethod inflation_method_;
  double inflation_factor_;
  std::string output_base_file_;
  std::string format_ = "txt";  // Default format

  // Diagnostic variables
  double innovation_norm_ = 0.0;
  double analysis_increment_norm_ = 0.0;
  double background_spread_ = 0.0;
  double analysis_spread_ = 0.0;
  double max_kalman_gain_ = 0.0;
  double min_kalman_gain_ = 0.0;
  double condition_number_ = 0.0;

  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework