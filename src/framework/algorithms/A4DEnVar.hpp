#pragma once
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "Logger.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "ProgressBar.hpp"
#include "State.hpp"

namespace metada::framework {

/**
 * @brief Analytical Four-Dimensional Ensemble-Variational (4DEnVar)
 * implementation.
 *
 * This class implements the analytical four-dimensional ensemble-variational
 * data assimilation algorithm as described in Liang et al. (2021). The A4DEnVar
 * combines the strengths of 4D-Var (temporal consistency) and ensemble methods
 * (flow-dependent error covariance) with an analytical solution that avoids
 * iterative minimization.
 *
 * The A4DEnVar algorithm follows these steps:
 * 1. Build ensemble perturbations across time windows:
 *    \[
 *    \mathbf{X}'_i(t) = \mathbf{x}_i(t) - \bar{\mathbf{x}}(t)
 *    \]
 *    where \(\mathbf{x}_i(t)\) is the i-th ensemble member at time t.
 *
 * 2. Compute ensemble-based background error covariance:
 *    \[
 *    \mathbf{P}^b = \frac{1}{N-1} \sum_{i=1}^{N} \mathbf{X}'_i \mathbf{X}'_i^T
 *    \]
 *
 * 3. Formulate the 4DEnVar cost function:
 *    \[
 *    J(\mathbf{x}_0) = \frac{1}{2} (\mathbf{x}_0 - \mathbf{x}_0^b)^T
 * \mathbf{P}^{b-1} (\mathbf{x}_0 - \mathbf{x}_0^b)
 *    + \frac{1}{2} \sum_{k=1}^{K} (\mathbf{y}_k - \mathbf{H}_k \mathbf{x}_k)^T
 * \mathbf{R}_k^{-1} (\mathbf{y}_k - \mathbf{H}_k \mathbf{x}_k)
 *    \]
 *
 * 4. Analytical solution for the analysis increment:
 *    \[
 *    \delta\mathbf{x}_0^a = \mathbf{P}^b \mathbf{H}^T (\mathbf{H} \mathbf{P}^b
 * \mathbf{H}^T + \mathbf{R})^{-1} \mathbf{d}
 *    \]
 *    where \(\mathbf{d}\) is the innovation vector and \(\mathbf{H}\) is the
 *    tangent linear observation operator.
 *
 * 5. Update the ensemble:
 *    \[
 *    \mathbf{x}_i^a = \mathbf{x}_i^b + \delta\mathbf{x}_0^a
 *    \]
 *
 * The A4DEnVar provides several advantages:
 * - Analytical solution avoids iterative minimization
 * - Flow-dependent covariance from ensemble
 * - Temporal consistency across multiple time windows
 * - Computationally efficient compared to standard 4D-Var
 *
 * For more details, see:
 * - Liang, K., et al. (2021). "An analytical four-dimensional
 * ensemble-variational data assimilation scheme." Journal of Advances in
 * Modeling Earth Systems, 13, e2020MS002314.
 *
 * @tparam BackendTag The backend tag type that must satisfy the required
 * concepts
 */
template <typename BackendTag>
class A4DEnVar {
 public:
  /**
   * @brief Analysis results structure
   */
  struct AnalysisResults {
    double cost_function_value;      ///< Final cost function value
    double background_cost;          ///< Background cost component
    double observation_cost;         ///< Observation cost component
    double innovation_norm;          ///< Innovation vector norm
    double analysis_increment_norm;  ///< Analysis increment norm
    double background_spread;        ///< Background ensemble spread
    double analysis_spread;          ///< Analysis ensemble spread
    double max_kalman_gain;          ///< Maximum Kalman gain element
    double min_kalman_gain;          ///< Minimum Kalman gain element
    double condition_number;      ///< Condition number of innovation covariance
    int ensemble_size;            ///< Number of ensemble members
    int observation_count;        ///< Total number of observations
    int time_windows;             ///< Number of time windows
    std::string covariance_type;  ///< Type of covariance used
    double localization_radius;   ///< Localization radius used
    std::string localization_function;  ///< Localization function used
    double inflation_factor;            ///< Inflation factor applied
    std::string inflation_method;       ///< Method used for inflation
  };

  /**
   * @brief Construct an A4DEnVar object.
   * @param ensemble Reference to the ensemble to be updated.
   * @param observations Vector of observation objects for different time
   * windows.
   * @param obs_operators Vector of observation operators for different time
   * windows.
   * @param config Configuration object containing A4DEnVar parameters.
   */
  A4DEnVar(
      Ensemble<BackendTag>& ensemble,
      const std::vector<std::unique_ptr<Observation<BackendTag>>>& observations,
      const std::vector<std::unique_ptr<ObsOperator<BackendTag>>>&
          obs_operators,
      const Config<BackendTag>& config)
      : ensemble_(ensemble),
        observations_(observations),
        obs_operators_(obs_operators) {
    // Get analysis configuration subsection
    auto analysis_config = config.GetSubsection("analysis");

    inflation_factor_ = analysis_config.Get("inflation").asFloat();
    localization_radius_ = analysis_config.Get("localization_radius").asFloat();
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

    // Parse localization function
    std::string loc_function_str =
        analysis_config.Get("localization_function").asString();
    if (loc_function_str == "gaussian") {
      localization_function_ = LocalizationFunction::GAUSSIAN;
    } else if (loc_function_str == "exponential") {
      localization_function_ = LocalizationFunction::EXPONENTIAL;
    } else if (loc_function_str == "cutoff") {
      localization_function_ = LocalizationFunction::CUTOFF;
    } else if (loc_function_str == "gaspari_cohn") {
      localization_function_ = LocalizationFunction::GASPARI_COHN;
    } else {
      localization_function_ = LocalizationFunction::GAUSSIAN;
      logger_.Warning() << "Unknown localization function '" << loc_function_str
                        << "', using Gaussian localization";
    }

    // Parse covariance type
    std::string cov_type_str =
        analysis_config.Get("covariance_type").asString();
    if (cov_type_str == "ensemble") {
      covariance_type_ = CovarianceType::ENSEMBLE;
    } else if (cov_type_str == "hybrid") {
      covariance_type_ = CovarianceType::HYBRID;
    } else if (cov_type_str == "static") {
      covariance_type_ = CovarianceType::STATIC;
    } else {
      covariance_type_ = CovarianceType::ENSEMBLE;
      logger_.Warning() << "Unknown covariance type '" << cov_type_str
                        << "', using ensemble covariance";
    }

    logger_.Info() << "A4DEnVar constructed with " << ensemble_.Size()
                   << " members and " << observations_.size()
                   << " time windows";
    logger_.Info() << "Inflation factor: " << inflation_factor_;
    logger_.Info() << "Localization radius: " << localization_radius_;
    logger_.Info() << "Localization function: " << loc_function_str;
    logger_.Info() << "Covariance type: " << cov_type_str;
  }

  /**
   * @brief Perform the A4DEnVar analysis step, updating the ensemble.
   */
  void Analyse() {
    logger_.Info() << "A4DEnVar analysis started";
    logger_.Info() << "Number of observations: " << observations_.size();
    logger_.Info() << "Number of obs operators: " << obs_operators_.size();
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int ens_size = ensemble_.Size();
    const int state_dim = ensemble_.GetMember(0).size();
    const int time_windows = observations_.size();

    // 1. Build ensemble perturbations across time windows
    std::vector<MatrixXd> X_pert_windows;
    std::vector<VectorXd> x_mean_windows;

    for (int t = 0; t < time_windows; ++t) {
      // Compute ensemble mean for this time window
      ensemble_.RecomputeMean();
      const auto& mean_state = ensemble_.Mean();
      const auto mean_data = mean_state.template getDataPtr<double>();
      VectorXd x_mean = Eigen::Map<const VectorXd>(mean_data, state_dim);
      x_mean_windows.push_back(x_mean);

      // Build perturbation matrix for this time window
      MatrixXd X_pert(state_dim, ens_size);
      for (int i = 0; i < ens_size; ++i) {
        const auto& pert_data =
            ensemble_.GetPerturbation(i).template getDataPtr<double>();
        X_pert.col(i) = Eigen::Map<const VectorXd>(pert_data, state_dim);
      }
      X_pert_windows.push_back(X_pert);
    }

    // 2. Compute ensemble-based background error covariance
    MatrixXd P_b = computeEnsembleCovariance(X_pert_windows[0]);

    // Apply inflation
    applyInflation(P_b);

    // Apply localization
    MatrixXd P_b_local = applyLocalization(P_b);

    // 3. Build observation operators and innovation vectors
    std::vector<MatrixXd> H_matrices;
    std::vector<VectorXd> innovations;
    std::vector<MatrixXd> R_matrices;
    int total_obs = 0;

    logger_.Info() << "Building observation operators for " << time_windows
                   << " time windows";
    logger_.Info() << "State dimension: " << state_dim;

    for (int t = 0; t < time_windows; ++t) {
      const auto& obs = *observations_[t];
      const auto& obs_op = *obs_operators_[t];
      const int obs_dim = obs.size();
      logger_.Info() << "Time window " << t << ": obs_dim = " << obs_dim;

      // Get observation data
      const auto obs_data = obs.getObservationValues();
      VectorXd y_o = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);

      // Get observation error covariance
      const auto& R_data = obs.getCovariance();
      VectorXd R_diag = Eigen::Map<const VectorXd>(R_data.data(), obs_dim);
      MatrixXd R = R_diag.asDiagonal();

      // Compute innovation for ensemble mean
      const auto& simulated_obs = obs_op.apply(ensemble_.Mean(), obs);
      VectorXd y_b = Eigen::Map<const VectorXd>(simulated_obs.data(), obs_dim);
      VectorXd innovation = y_o - y_b;

      // Store matrices and vectors
      // Create H matrix that maps from state space to observation space
      // For A4DEnVar, we need to ensure H has the correct dimensions
      MatrixXd H = MatrixXd::Zero(obs_dim, state_dim);

      // Create a simple linear mapping - in practice, this should be the actual
      // tangent linear of the observation operator
      // For now, we'll use a simple interpolation-like mapping
      for (int i = 0; i < obs_dim; ++i) {
        // Map each observation to a corresponding state element
        // This is a simplified assumption - in practice, H should come from the
        // obs operator
        int state_idx = i % state_dim;  // Simple modulo mapping
        H(i, state_idx) = 1.0;
      }

      // Ensure H matrix is valid
      if (H.rows() != obs_dim || H.cols() != state_dim) {
        logger_.Error() << "H matrix dimension mismatch: " << H.rows() << "x"
                        << H.cols() << " expected " << obs_dim << "x"
                        << state_dim;
        throw std::runtime_error("H matrix dimension mismatch");
      }

      H_matrices.push_back(H);
      innovations.push_back(innovation);
      R_matrices.push_back(R);
      total_obs += obs_dim;
    }

    // 4. Formulate the 4DEnVar problem
    logger_.Info() << "Building 4D observation operator...";
    MatrixXd H_4d = buildFourDObservationOperator(H_matrices);
    logger_.Info() << "Building 4D innovation vector...";
    VectorXd d_4d = buildFourDInnovationVector(innovations);
    logger_.Info() << "Building 4D observation covariance...";
    MatrixXd R_4d = buildFourDObservationCovariance(R_matrices);

    // Debug: Print matrix dimensions
    logger_.Info() << "Matrix dimensions:";
    logger_.Info() << "H_4d: " << H_4d.rows() << " x " << H_4d.cols();
    logger_.Info() << "P_b_local: " << P_b_local.rows() << " x "
                   << P_b_local.cols();
    logger_.Info() << "R_4d: " << R_4d.rows() << " x " << R_4d.cols();
    logger_.Info() << "d_4d: " << d_4d.size();
    logger_.Info() << "Total observations: " << total_obs;
    logger_.Info() << "Time windows: " << time_windows;

    // 5. Compute analytical solution
    // Check matrix dimensions before multiplication
    logger_.Info() << "Computing S = H * P * H^T + R";
    logger_.Info() << "H_4d: " << H_4d.rows() << " x " << H_4d.cols();
    logger_.Info() << "P_b_local: " << P_b_local.rows() << " x "
                   << P_b_local.cols();
    logger_.Info() << "R_4d: " << R_4d.rows() << " x " << R_4d.cols();

    MatrixXd S, K;
    try {
      // Compute S = H * P * H^T + R step by step to catch dimension issues
      MatrixXd HP = H_4d * P_b_local;
      logger_.Info() << "HP: " << HP.rows() << " x " << HP.cols();
      MatrixXd HPHt = HP * H_4d.transpose();
      logger_.Info() << "HPHt: " << HPHt.rows() << " x " << HPHt.cols();
      S = HPHt + R_4d;
      logger_.Info() << "S: " << S.rows() << " x " << S.cols();

      K = P_b_local * H_4d.transpose() * S.inverse();
    } catch (const std::exception& e) {
      logger_.Error() << "Matrix multiplication error: " << e.what();
      throw;
    }

    // Store Kalman gain statistics
    max_kalman_gain_ = K.maxCoeff();
    min_kalman_gain_ = K.minCoeff();

    // Compute condition number using singular values
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(S);
    auto singular_values = svd.singularValues();
    condition_number_ =
        singular_values(0) / singular_values(singular_values.size() - 1);

    // 6. Compute analysis increment
    VectorXd delta_x = K * d_4d;

    // Store analysis increment norm
    analysis_increment_norm_ = delta_x.norm();

    // 7. Update ensemble members
    for (int i = 0; i < ens_size; ++i) {
      auto& member = ensemble_.GetMember(i);
      auto data = member.template getDataPtr<double>();
      VectorXd x_current = Eigen::Map<VectorXd>(data, state_dim);
      VectorXd x_updated = x_current + delta_x;
      Eigen::Map<VectorXd>(data, state_dim) = x_updated;
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

    // Compute cost function components
    computeCostFunctionComponents(P_b_local, d_4d, R_4d);

    logger_.Info() << "A4DEnVar analysis completed";
  }

  /**
   * @brief Save the analyzed ensemble to files
   */
  void saveEnsemble() const {
    logger_.Info() << "A4DEnVar saving ensemble";
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

    logger_.Info() << "A4DEnVar ensemble saved";
  }

  /**
   * @brief Get the analysis results
   */
  AnalysisResults getAnalysisResults() const {
    AnalysisResults results;
    results.cost_function_value = cost_function_value_;
    results.background_cost = background_cost_;
    results.observation_cost = observation_cost_;
    results.innovation_norm = innovation_norm_;
    results.analysis_increment_norm = analysis_increment_norm_;
    results.background_spread = background_spread_;
    results.analysis_spread = analysis_spread_;
    results.max_kalman_gain = max_kalman_gain_;
    results.min_kalman_gain = min_kalman_gain_;
    results.condition_number = condition_number_;
    results.ensemble_size = ensemble_.Size();
    results.observation_count = total_observation_count_;
    results.time_windows = observations_.size();
    results.covariance_type = getCovarianceTypeString();
    results.localization_radius = localization_radius_;
    results.localization_function = getLocalizationFunctionString();
    results.inflation_factor = inflation_factor_;
    results.inflation_method = getInflationMethodString();
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
   * @brief Localization function type enumeration
   */
  enum class LocalizationFunction {
    GAUSSIAN,     ///< Gaussian localization function
    EXPONENTIAL,  ///< Exponential localization function
    CUTOFF,       ///< Cutoff localization function
    GASPARI_COHN  ///< Gaspari-Cohn localization function
  };

  /**
   * @brief Covariance type enumeration
   */
  enum class CovarianceType {
    ENSEMBLE,  ///< Pure ensemble covariance
    HYBRID,    ///< Hybrid ensemble-static covariance
    STATIC     ///< Static climatological covariance
  };

  /**
   * @brief Compute ensemble-based background error covariance
   */
  Eigen::MatrixXd computeEnsembleCovariance(
      const Eigen::MatrixXd& X_pert) const {
    const int ens_size = X_pert.cols();
    return X_pert * X_pert.transpose() / (ens_size - 1);
  }

  /**
   * @brief Apply inflation to background covariance
   */
  void applyInflation(Eigen::MatrixXd& P_b) {
    switch (inflation_method_) {
      case InflationMethod::MULTIPLICATIVE:
        P_b *= inflation_factor_;
        break;
      case InflationMethod::ADDITIVE:
        // Additive inflation would add random perturbations
        // For simplicity, we use multiplicative here
        P_b *= inflation_factor_;
        break;
      case InflationMethod::RELAXATION:
        // Relaxation inflation reduces the analysis increment
        // For simplicity, we use multiplicative here
        P_b *= inflation_factor_;
        break;
    }

    // Compute background spread for diagnostics
    background_spread_ = std::sqrt(P_b.trace() / P_b.rows());
  }

  /**
   * @brief Apply localization to covariance matrix
   */
  Eigen::MatrixXd applyLocalization(const Eigen::MatrixXd& cov) const {
    using Eigen::MatrixXd;

    const int dim = cov.rows();
    MatrixXd localized_cov = cov;

    // Create localization matrix (simplified - in practice would use distance
    // information)
    MatrixXd L = MatrixXd::Ones(dim, dim);

    // Apply localization function
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        double distance = std::abs(i - j) / static_cast<double>(dim);
        L(i, j) = computeLocalizationFunction(distance);
      }
    }

    return cov.cwiseProduct(L);
  }

  /**
   * @brief Compute localization function value
   */
  double computeLocalizationFunction(double distance) const {
    double normalized_distance = distance / localization_radius_;

    switch (localization_function_) {
      case LocalizationFunction::GAUSSIAN:
        return std::exp(-0.5 * normalized_distance * normalized_distance);

      case LocalizationFunction::EXPONENTIAL:
        return std::exp(-normalized_distance);

      case LocalizationFunction::CUTOFF:
        return (normalized_distance <= 1.0) ? 1.0 : 0.0;

      case LocalizationFunction::GASPARI_COHN:
        return computeGaspariCohnFunction(normalized_distance);

      default:
        return std::exp(-0.5 * normalized_distance * normalized_distance);
    }
  }

  /**
   * @brief Compute Gaspari-Cohn localization function
   */
  double computeGaspariCohnFunction(double r) const {
    if (r >= 2.0) return 0.0;
    if (r >= 1.0) {
      double z = r - 1.0;
      return ((-0.25 * z + 0.5) * z + 0.625) * z + 0.125;
    } else {
      double z = r;
      return (((-0.25 * z + 0.5) * z + 0.625) * z - 5.0) * z + 4.0;
    }
  }

  /**
   * @brief Build 4D observation operator matrix
   */
  Eigen::MatrixXd buildFourDObservationOperator(
      const std::vector<Eigen::MatrixXd>& H_matrices) const {
    logger_.Info() << "buildFourDObservationOperator: " << H_matrices.size()
                   << " matrices";

    int total_obs = 0;
    for (size_t i = 0; i < H_matrices.size(); ++i) {
      const auto& H = H_matrices[i];
      logger_.Info() << "H[" << i << "]: " << H.rows() << " x " << H.cols();
      total_obs += H.rows();
    }

    if (H_matrices.empty()) {
      logger_.Error() << "No H matrices provided";
      return Eigen::MatrixXd::Zero(1, 1);
    }

    int state_dim = H_matrices[0].cols();
    logger_.Info() << "Total obs: " << total_obs
                   << ", state_dim: " << state_dim;

    Eigen::MatrixXd H_4d = Eigen::MatrixXd::Zero(total_obs, state_dim);
    logger_.Info() << "H_4d created: " << H_4d.rows() << " x " << H_4d.cols();

    int row_offset = 0;
    for (size_t i = 0; i < H_matrices.size(); ++i) {
      const auto& H = H_matrices[i];
      logger_.Info() << "Copying H[" << i << "] to block at row " << row_offset;
      H_4d.block(row_offset, 0, H.rows(), H.cols()) = H;
      row_offset += H.rows();
    }

    return H_4d;
  }

  /**
   * @brief Build 4D innovation vector
   */
  Eigen::VectorXd buildFourDInnovationVector(
      const std::vector<Eigen::VectorXd>& innovations) const {
    int total_obs = 0;
    for (const auto& d : innovations) {
      total_obs += d.size();
    }

    Eigen::VectorXd d_4d(total_obs);
    int offset = 0;
    for (const auto& d : innovations) {
      d_4d.segment(offset, d.size()) = d;
      offset += d.size();
    }

    return d_4d;
  }

  /**
   * @brief Build 4D observation error covariance matrix
   */
  Eigen::MatrixXd buildFourDObservationCovariance(
      const std::vector<Eigen::MatrixXd>& R_matrices) const {
    int total_obs = 0;
    for (const auto& R : R_matrices) {
      total_obs += R.rows();
    }

    Eigen::MatrixXd R_4d = Eigen::MatrixXd::Zero(total_obs, total_obs);

    int offset = 0;
    for (const auto& R : R_matrices) {
      R_4d.block(offset, offset, R.rows(), R.cols()) = R;
      offset += R.rows();
    }

    return R_4d;
  }

  /**
   * @brief Compute cost function components
   */
  void computeCostFunctionComponents(const Eigen::MatrixXd& P_b,
                                     const Eigen::VectorXd& d,
                                     const Eigen::MatrixXd& R) {
    // Background cost
    background_cost_ = 0.5 * d.transpose() * P_b.inverse() * d;

    // Observation cost
    observation_cost_ = 0.5 * d.transpose() * R.inverse() * d;

    // Total cost
    cost_function_value_ = background_cost_ + observation_cost_;

    // Innovation norm
    innovation_norm_ = d.norm();
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

  /**
   * @brief Get localization function as string
   */
  std::string getLocalizationFunctionString() const {
    switch (localization_function_) {
      case LocalizationFunction::GAUSSIAN:
        return "gaussian";
      case LocalizationFunction::EXPONENTIAL:
        return "exponential";
      case LocalizationFunction::CUTOFF:
        return "cutoff";
      case LocalizationFunction::GASPARI_COHN:
        return "gaspari_cohn";
      default:
        return "unknown";
    }
  }

  /**
   * @brief Get covariance type as string
   */
  std::string getCovarianceTypeString() const {
    switch (covariance_type_) {
      case CovarianceType::ENSEMBLE:
        return "ensemble";
      case CovarianceType::HYBRID:
        return "hybrid";
      case CovarianceType::STATIC:
        return "static";
      default:
        return "unknown";
    }
  }

  Ensemble<BackendTag>& ensemble_;
  const std::vector<std::unique_ptr<Observation<BackendTag>>>& observations_;
  const std::vector<std::unique_ptr<ObsOperator<BackendTag>>>& obs_operators_;
  InflationMethod inflation_method_;
  double inflation_factor_;
  LocalizationFunction localization_function_;
  double localization_radius_;
  CovarianceType covariance_type_;
  std::string output_base_file_;
  std::string format_ = "txt";  // Default format

  // Diagnostic variables
  double cost_function_value_ = 0.0;
  double background_cost_ = 0.0;
  double observation_cost_ = 0.0;
  double innovation_norm_ = 0.0;
  double analysis_increment_norm_ = 0.0;
  double background_spread_ = 0.0;
  double analysis_spread_ = 0.0;
  double max_kalman_gain_ = 0.0;
  double min_kalman_gain_ = 0.0;
  double condition_number_ = 0.0;
  int total_observation_count_ = 0;

  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework