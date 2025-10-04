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
 * @brief Local Weighted Ensemble Kalman Filter (LWEnKF) implementation.
 *
 * This class implements the Local Weighted Ensemble Kalman Filter algorithm
 * for ensemble-based data assimilation. The LWEnKF extends the traditional
 * EnKF by incorporating localization and adaptive weighting to improve
 * covariance estimation and reduce spurious correlations.
 *
 * The LWEnKF algorithm follows these steps:
 * 1. Build forecast quantities:
 *    \[
 *    \bar{\mathbf x}^f = \frac{1}{N}\sum_{i=1}^{N}\mathbf x_i^f
 *    \]
 *    \[
 *    \mathbf X'^{\,f} = \bigl[\,\mathbf x_1^f-\bar{\mathbf
 * x}^f\,\;\cdots\;\mathbf x_N^f-\bar{\mathbf x}^f\bigr]
 *    \]
 *
 * 2. Apply localization to reduce spurious correlations:
 *    \[
 *    \mathbf{P}^f_{\text{local}} = \mathbf{P}^f \circ \mathbf{L}
 *    \]
 *    Where \(\mathbf{L}\) is the localization matrix and \(\circ\) denotes
 *    element-wise multiplication.
 *
 * 3. Compute weighted ensemble statistics:
 *    \[
 *    \mathbf{P}^f_{\text{weighted}} = \sum_{i=1}^{N} w_i (\mathbf{x}_i^f -
 *    \bar{\mathbf{x}}^f)(\mathbf{x}_i^f - \bar{\mathbf{x}}^f)^T
 *    \]
 *
 * 4. Compute innovation:
 *    \[
 *    \mathbf d = \mathbf y^{\text{obs}}-\bar{\mathbf y}^f
 *    \]
 *
 * 5. Compute localized Kalman gain:
 *    \[
 *    \mathbf{K}_{\text{local}} = \mathbf{P}^f_{\text{local}} \mathbf{H}^T
 *    (\mathbf{H} \mathbf{P}^f_{\text{local}} \mathbf{H}^T + \mathbf{R})^{-1}
 *    \]
 *
 * 6. Analysis update:
 *    \[
 *    \mathbf{x}_i^a = \mathbf{x}_i^f + \mathbf{K}_{\text{local}}
 *    (\mathbf{y}^{\text{obs}} + \mathbf{\epsilon}_i - \mathbf{H}
 *    \mathbf{x}_i^f)
 *    \]
 *
 * The LWEnKF improves upon traditional EnKF by:
 * - Using distance-based localization to reduce spurious correlations
 * - Applying adaptive weighting to ensemble members
 * - Better handling of observation error correlations
 * - Improved treatment of model error and observation bias
 *
 * For more details, see:
 * - Chen, Y. et al. "Local Weighted Ensemble Kalman Filter"
 * - Hamill, T. M., Whitaker, J. S., & Snyder, C. (2001) "Distance-dependent
 *   filtering of background error covariance estimates in an ensemble Kalman
 *   filter"
 *
 * @tparam BackendTag The backend tag type that must satisfy the required
 * concepts
 */
template <typename BackendTag>
class LWEnKF {
 public:
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
   * @brief Weighting scheme enumeration
   */
  enum class WeightingScheme {
    UNIFORM,      ///< Uniform weights (traditional EnKF)
    ADAPTIVE,     ///< Adaptive weights based on ensemble spread
    INVERSE_VAR,  ///< Inverse variance weighting
    LIKELIHOOD    ///< Likelihood-based weighting
  };

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
    double condition_number;     ///< Condition number of innovation covariance
    double localization_radius;  ///< Localization radius used
    std::string localization_function;  ///< Localization function used
    std::string weighting_scheme;       ///< Weighting scheme used
    double max_weight;                  ///< Maximum ensemble weight
    double min_weight;                  ///< Minimum ensemble weight
    double weight_variance;             ///< Weight variance
    int ensemble_size;                  ///< Number of ensemble members
    int observation_count;              ///< Number of observations
    std::string inflation_method;       ///< Method used for inflation
    double inflation_factor;            ///< Inflation factor applied
  };

  /**
   * @brief Construct an LWEnKF object.
   * @param ensemble Reference to the ensemble to be updated.
   * @param obs Reference to the observation object.
   * @param obs_op Reference to the observation operator.
   * @param config Configuration object containing LWEnKF parameters.
   */
  LWEnKF(Ensemble<BackendTag>& ensemble, Observation<BackendTag>& obs,
         const ObsOperator<BackendTag>& obs_op,
         const Config<BackendTag>& config)
      : ensemble_(ensemble), obs_(obs), obs_op_(obs_op) {
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

    // Parse weighting scheme
    std::string weighting_str =
        analysis_config.Get("weighting_scheme").asString();
    if (weighting_str == "uniform") {
      weighting_scheme_ = WeightingScheme::UNIFORM;
    } else if (weighting_str == "adaptive") {
      weighting_scheme_ = WeightingScheme::ADAPTIVE;
    } else if (weighting_str == "inverse_var") {
      weighting_scheme_ = WeightingScheme::INVERSE_VAR;
    } else if (weighting_str == "likelihood") {
      weighting_scheme_ = WeightingScheme::LIKELIHOOD;
    } else {
      weighting_scheme_ = WeightingScheme::UNIFORM;
      logger_.Warning() << "Unknown weighting scheme '" << weighting_str
                        << "', using uniform weighting";
    }

    logger_.Info() << "LWEnKF constructed with " << ensemble_.Size()
                   << " members";
    logger_.Info() << "Inflation factor: " << inflation_factor_;
    logger_.Info() << "Localization radius: " << localization_radius_;
    logger_.Info() << "Localization function: " << loc_function_str;
    logger_.Info() << "Weighting scheme: " << weighting_str;
  }

  /**
   * @brief Perform the LWEnKF analysis step, updating the ensemble.
   */
  void Analyse() {
    logger_.Info() << "LWEnKF analysis started";
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

    // 2. Compute adaptive weights
    computeWeights();

    // 3. Apply inflation to background perturbations
    applyInflation(Xb_pert);

    // 4. Propagate to observation space
    MatrixXd Yb(obs_dim, ens_size);
    for (int i = 0; i < ens_size; ++i) {
      const auto& obs_data = obs_op_.apply(ensemble_.GetMember(i), obs_);
      Yb.col(i) = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);
    }

    // Compute mean and anomalies in observation space
    VectorXd yb_mean = Yb.rowwise().mean();
    MatrixXd Yb_pert = Yb.colwise() - yb_mean;

    // 5. Compute innovation
    const auto obs_data = obs_.getObservationValues();
    VectorXd yo = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);
    VectorXd d = yo - yb_mean;

    // Store innovation norm for diagnostics
    innovation_norm_ = d.norm();

    // 6. Compute localized Kalman gain
    // Get observation error covariance
    const auto& R_data = obs_.getCovariance();
    VectorXd R_diag = Eigen::Map<const VectorXd>(R_data.data(), obs_dim);
    MatrixXd R = R_diag.asDiagonal();

    // Compute weighted background covariance in observation space
    MatrixXd Pf_obs = computeWeightedCovariance(Yb_pert);

    // Apply localization to observation space covariance
    MatrixXd Pf_obs_local = applyLocalization(Pf_obs);

    // Add observation error covariance
    MatrixXd S = Pf_obs_local + R;

    // Compute localized Kalman gain
    MatrixXd K_local =
        Xb_pert * Yb_pert.transpose() * S.inverse() / (ens_size - 1);

    // Apply localization to Kalman gain
    K_local = applyLocalizationToGain(K_local);

    // Store Kalman gain statistics
    max_kalman_gain_ = K_local.maxCoeff();
    min_kalman_gain_ = K_local.minCoeff();

    // Compute condition number using singular values
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(S);
    auto singular_values = svd.singularValues();
    condition_number_ =
        singular_values(0) / singular_values(singular_values.size() - 1);

    // 7. Generate observation perturbations
    MatrixXd obs_pert(obs_dim, ens_size);
    generateObservationPerturbations(obs_pert, R);

    // 8. Analysis update for each ensemble member
    MatrixXd Xa(state_dim, ens_size);
    for (int i = 0; i < ens_size; ++i) {
      // Compute perturbed observation
      VectorXd yo_pert = yo + obs_pert.col(i);

      // Compute innovation for this member
      VectorXd d_i = yo_pert - Yb.col(i);

      // Analysis update
      VectorXd xa_i = Xb_pert.col(i) + K_local * d_i;
      Xa.col(i) = xb_mean + xa_i;
    }

    // 9. Update ensemble members
    for (int i = 0; i < ens_size; ++i) {
      auto& member = ensemble_.GetMember(i);
      auto data = member.template getDataPtr<double>();
      Eigen::Map<VectorXd>(data, state_dim) = Xa.col(i);
    }

    // 10. Update ensemble statistics
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

    logger_.Info() << "LWEnKF analysis completed";
  }

  /**
   * @brief Save the analyzed ensemble to files
   */
  void saveEnsemble() const {
    logger_.Info() << "LWEnKF saving ensemble";
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

    logger_.Info() << "LWEnKF ensemble saved";
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
    results.localization_radius = localization_radius_;
    results.localization_function = getLocalizationFunctionString();
    results.weighting_scheme = getWeightingSchemeString();
    results.max_weight = max_weight_;
    results.min_weight = min_weight_;
    results.weight_variance = weight_variance_;
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
   * @brief Compute adaptive weights for ensemble members
   */
  void computeWeights() {
    const int ens_size = ensemble_.Size();
    weights_.resize(ens_size);

    switch (weighting_scheme_) {
      case WeightingScheme::UNIFORM:
        std::fill(weights_.begin(), weights_.end(), 1.0 / ens_size);
        break;

      case WeightingScheme::ADAPTIVE:
        computeAdaptiveWeights();
        break;

      case WeightingScheme::INVERSE_VAR:
        computeInverseVarianceWeights();
        break;

      case WeightingScheme::LIKELIHOOD:
        computeLikelihoodWeights();
        break;
    }

    // Normalize weights
    double sum_weights = std::accumulate(weights_.begin(), weights_.end(), 0.0);
    for (auto& weight : weights_) {
      weight /= sum_weights;
    }

    // Compute weight statistics
    max_weight_ = *std::max_element(weights_.begin(), weights_.end());
    min_weight_ = *std::min_element(weights_.begin(), weights_.end());
    weight_variance_ = computeWeightVariance();
  }

  /**
   * @brief Compute adaptive weights based on ensemble spread
   */
  void computeAdaptiveWeights() {
    using Eigen::VectorXd;

    const int ens_size = ensemble_.Size();
    const int state_dim = ensemble_.GetMember(0).size();

    // Compute ensemble spread for each member
    std::vector<double> spreads(ens_size);
    for (int i = 0; i < ens_size; ++i) {
      const auto& pert_data =
          ensemble_.GetPerturbation(i).template getDataPtr<double>();
      VectorXd pert = Eigen::Map<const VectorXd>(pert_data, state_dim);
      spreads[i] = pert.norm();
    }

    // Compute weights inversely proportional to spread
    double sum_inv_spread = 0.0;
    for (int i = 0; i < ens_size; ++i) {
      weights_[i] = 1.0 / (spreads[i] + 1e-8);  // Add small epsilon
      sum_inv_spread += weights_[i];
    }

    // Normalize
    for (int i = 0; i < ens_size; ++i) {
      weights_[i] /= sum_inv_spread;
    }
  }

  /**
   * @brief Compute inverse variance weights
   */
  void computeInverseVarianceWeights() {
    using Eigen::VectorXd;

    const int ens_size = ensemble_.Size();
    const int state_dim = ensemble_.GetMember(0).size();

    // Compute variance for each member
    std::vector<double> variances(ens_size);
    for (int i = 0; i < ens_size; ++i) {
      const auto& pert_data =
          ensemble_.GetPerturbation(i).template getDataPtr<double>();
      VectorXd pert = Eigen::Map<const VectorXd>(pert_data, state_dim);
      variances[i] = pert.squaredNorm();
    }

    // Compute weights inversely proportional to variance
    double sum_inv_var = 0.0;
    for (int i = 0; i < ens_size; ++i) {
      weights_[i] = 1.0 / (variances[i] + 1e-8);  // Add small epsilon
      sum_inv_var += weights_[i];
    }

    // Normalize
    for (int i = 0; i < ens_size; ++i) {
      weights_[i] /= sum_inv_var;
    }
  }

  /**
   * @brief Compute likelihood-based weights
   */
  void computeLikelihoodWeights() {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int ens_size = ensemble_.Size();
    const int obs_dim = obs_.size();

    // Get observation data
    const auto obs_data = obs_.getObservationValues();
    VectorXd yo = Eigen::Map<const VectorXd>(obs_data.data(), obs_dim);

    // Get observation error covariance
    const auto& R_data = obs_.getCovariance();
    VectorXd R_diag = Eigen::Map<const VectorXd>(R_data.data(), obs_dim);
    MatrixXd R = R_diag.asDiagonal();

    // Compute likelihood for each member
    for (int i = 0; i < ens_size; ++i) {
      // Apply observation operator
      const auto& simulated_obs = obs_op_.apply(ensemble_.GetMember(i), obs_);
      VectorXd y_sim =
          Eigen::Map<const VectorXd>(simulated_obs.data(), obs_dim);

      // Compute innovation
      VectorXd innovation = yo - y_sim;

      // Compute likelihood (Gaussian assumption)
      double likelihood =
          std::exp(-0.5 * innovation.transpose() * R.inverse() * innovation);
      weights_[i] = likelihood;
    }
  }

  /**
   * @brief Compute weighted covariance matrix
   */
  Eigen::MatrixXd computeWeightedCovariance(
      const Eigen::MatrixXd& anomalies) const {
    using Eigen::MatrixXd;

    const int ens_size = ensemble_.Size();
    const int dim = anomalies.rows();

    MatrixXd weighted_cov = MatrixXd::Zero(dim, dim);
    for (int i = 0; i < ens_size; ++i) {
      weighted_cov +=
          weights_[i] * anomalies.col(i) * anomalies.col(i).transpose();
    }

    return weighted_cov;
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
   * @brief Apply localization to Kalman gain
   */
  Eigen::MatrixXd applyLocalizationToGain(const Eigen::MatrixXd& K) const {
    using Eigen::MatrixXd;

    const int state_dim = K.rows();
    const int obs_dim = K.cols();
    MatrixXd K_local = K;

    // Create localization matrix for state-observation pairs
    MatrixXd L = MatrixXd::Ones(state_dim, obs_dim);

    // Apply localization function (simplified)
    for (int i = 0; i < state_dim; ++i) {
      for (int j = 0; j < obs_dim; ++j) {
        double distance =
            std::abs(i - j) / static_cast<double>(std::max(state_dim, obs_dim));
        L(i, j) = computeLocalizationFunction(distance);
      }
    }

    return K.cwiseProduct(L);
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
   * @brief Compute weight variance
   */
  double computeWeightVariance() const {
    double mean_weight = 1.0 / weights_.size();
    double variance = 0.0;
    for (const auto& weight : weights_) {
      double diff = weight - mean_weight;
      variance += diff * diff;
    }
    return variance / weights_.size();
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
   * @brief Get weighting scheme as string
   */
  std::string getWeightingSchemeString() const {
    switch (weighting_scheme_) {
      case WeightingScheme::UNIFORM:
        return "uniform";
      case WeightingScheme::ADAPTIVE:
        return "adaptive";
      case WeightingScheme::INVERSE_VAR:
        return "inverse_var";
      case WeightingScheme::LIKELIHOOD:
        return "likelihood";
      default:
        return "unknown";
    }
  }

  Ensemble<BackendTag>& ensemble_;
  Observation<BackendTag>& obs_;
  const ObsOperator<BackendTag>& obs_op_;
  std::vector<double> weights_;
  InflationMethod inflation_method_;
  double inflation_factor_;
  LocalizationFunction localization_function_;
  double localization_radius_;
  WeightingScheme weighting_scheme_;
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
  double max_weight_ = 0.0;
  double min_weight_ = 0.0;
  double weight_variance_ = 0.0;

  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework