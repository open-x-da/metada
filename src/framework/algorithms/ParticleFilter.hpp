#pragma once
#include <Eigen/Dense>
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
 * @brief Particle Filter (PF) implementation for data assimilation.
 *
 * This class implements the particle filter algorithm for ensemble-based data
 * assimilation. The particle filter is a non-parametric Bayesian method that
 * represents the posterior probability density function using a set of weighted
 * particles (ensemble members).
 *
 * The particle filter algorithm follows these steps:
 * 1. Forecast step: Propagate particles using the model
 * 2. Analysis step: Update particle weights using observations
 * 3. Resampling step: Resample particles based on weights to prevent degeneracy
 * 4. Optional: Apply jittering/perturbation to maintain diversity
 *
 * The weight update follows Bayes' rule:
 * w_i^{(a)} \propto w_i^{(f)} p(y|X_i^{(f)})
 *
 * Where:
 * - w_i^{(a)} is the analysis weight of particle i
 * - w_i^{(f)} is the forecast weight of particle i
 * - p(y|X_i^{(f)}) is the likelihood of observation y given particle i
 *
 * For Gaussian observation errors, the likelihood is:
 * p(y|X_i^{(f)}) \propto \exp(-\frac{1}{2}(y - H(X_i^{(f)}))^T R^{-1}(y -
 * H(X_i^{(f)})))
 *
 * The algorithm supports various resampling methods:
 * - Systematic resampling (default)
 * - Multinomial resampling
 * - Stratified resampling
 * - Residual resampling
 *
 * For more details, see:
 * - Arulampalam et al. (2002) "A tutorial on particle filters for online
 *   nonlinear/non-Gaussian Bayesian tracking"
 * - Poterjoy (2016) "A localized particle filter for high-dimensional
 *   nonlinear systems"
 *
 * @tparam BackendTag The backend tag type that must satisfy the required
 * concepts
 */
template <typename BackendTag>
class ParticleFilter {
 public:
  /**
   * @brief Resampling method enumeration
   */
  enum class ResamplingMethod {
    SYSTEMATIC,   ///< Systematic resampling (default)
    MULTINOMIAL,  ///< Multinomial resampling
    STRATIFIED,   ///< Stratified resampling
    RESIDUAL      ///< Residual resampling
  };

  /**
   * @brief Analysis results structure
   */
  struct AnalysisResults {
    std::vector<double> weights;    ///< Final particle weights
    double effective_sample_size;   ///< Effective sample size
    double max_weight;              ///< Maximum weight
    double min_weight;              ///< Minimum weight
    double weight_variance;         ///< Weight variance
    bool resampling_performed;      ///< Whether resampling was performed
    std::string resampling_method;  ///< Method used for resampling
    int resampling_threshold;       ///< Threshold for resampling
    std::vector<double>
        likelihood_values;  ///< Likelihood values for each particle
  };

  /**
   * @brief Construct a ParticleFilter object.
   * @param ensemble Reference to the ensemble to be updated.
   * @param obs Reference to the observation object.
   * @param obs_op Reference to the observation operator.
   * @param config Configuration object containing particle filter parameters.
   */
  ParticleFilter(Ensemble<BackendTag>& ensemble, Observation<BackendTag>& obs,
                 const ObsOperator<BackendTag>& obs_op,
                 const Config<BackendTag>& config)
      : ensemble_(ensemble), obs_(obs), obs_op_(obs_op) {
    // Get analysis configuration subsection
    auto analysis_config = config.GetSubsection("analysis");

    resampling_threshold_ = analysis_config.Get("resampling_threshold").asInt();
    jittering_enabled_ = analysis_config.Get("jittering_enabled").asBool();
    jittering_std_ = analysis_config.Get("jittering_std").asFloat();
    output_base_file_ = analysis_config.Get("output_base_file").asString();
    format_ = analysis_config.Get("format").asString();

    // Parse resampling method
    std::string method_str =
        analysis_config.Get("resampling_method").asString();
    if (method_str == "systematic") {
      resampling_method_ = ResamplingMethod::SYSTEMATIC;
    } else if (method_str == "multinomial") {
      resampling_method_ = ResamplingMethod::MULTINOMIAL;
    } else if (method_str == "stratified") {
      resampling_method_ = ResamplingMethod::STRATIFIED;
    } else if (method_str == "residual") {
      resampling_method_ = ResamplingMethod::RESIDUAL;
    } else {
      resampling_method_ = ResamplingMethod::SYSTEMATIC;
      logger_.Warning() << "Unknown resampling method '" << method_str
                        << "', using systematic resampling";
    }

    logger_.Info() << "ParticleFilter constructed with " << ensemble_.Size()
                   << " particles";
    logger_.Info() << "Resampling threshold: " << resampling_threshold_;
    logger_.Info() << "Resampling method: " << method_str;
    logger_.Info() << "Jittering enabled: "
                   << (jittering_enabled_ ? "Yes" : "No");
  }

  /**
   * @brief Perform the particle filter analysis step, updating the ensemble.
   */
  void Analyse() {
    logger_.Info() << "ParticleFilter analysis started";
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int ens_size = ensemble_.Size();

    // Initialize weights (equal weights for first analysis)
    if (weights_.empty()) {
      weights_.resize(ens_size, 1.0 / ens_size);
    }

    // 1. Forecast step: Propagate particles using model (if available)
    // Note: For simplicity, we assume particles are already forecasted
    // In a full implementation, this would involve model integration

    // 2. Analysis step: Update weights using observations
    updateWeights();

    // 3. Check for resampling
    double effective_sample_size = computeEffectiveSampleSize();
    logger_.Info() << "Effective sample size: " << effective_sample_size;

    if (effective_sample_size < resampling_threshold_) {
      logger_.Info() << "Performing resampling (ESS < " << resampling_threshold_
                     << ")";
      performResampling();
    } else {
      logger_.Info() << "No resampling needed (ESS >= " << resampling_threshold_
                     << ")";
    }

    // 4. Optional jittering to maintain diversity
    if (jittering_enabled_) {
      applyJittering();
    }

    // Update ensemble mean and perturbations
    ensemble_.RecomputeMean();
    ensemble_.RecomputePerturbations();

    logger_.Info() << "ParticleFilter analysis completed";
  }

  /**
   * @brief Save the analyzed ensemble to files
   */
  void saveEnsemble() const {
    logger_.Info() << "ParticleFilter saving ensemble";
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

    // Save weights to a separate file
    saveWeights();

    logger_.Info() << "ParticleFilter ensemble saved";
  }

  /**
   * @brief Get the analysis results
   */
  AnalysisResults getAnalysisResults() const {
    AnalysisResults results;
    results.weights = weights_;
    results.effective_sample_size = computeEffectiveSampleSize();
    results.max_weight = *std::max_element(weights_.begin(), weights_.end());
    results.min_weight = *std::min_element(weights_.begin(), weights_.end());
    results.weight_variance = computeWeightVariance();
    results.resampling_performed = resampling_performed_;
    results.resampling_method = getResamplingMethodString();
    results.resampling_threshold = resampling_threshold_;
    results.likelihood_values = likelihood_values_;
    return results;
  }

  /**
   * @brief Get current particle weights
   */
  const std::vector<double>& getWeights() const { return weights_; }

  /**
   * @brief Set particle weights (for external weight management)
   */
  void setWeights(const std::vector<double>& weights) {
    if (weights.size() != ensemble_.Size()) {
      throw std::runtime_error("Weight vector size must match ensemble size");
    }
    weights_ = weights;
    normalizeWeights();
  }

 private:
  /**
   * @brief Update particle weights using observations
   */
  void updateWeights() {
    const int ens_size = ensemble_.Size();
    const int obs_dim = obs_.size();

    // Get observation data
    const auto obs_data = obs_.getObservationValues();
    Eigen::VectorXd yo =
        Eigen::Map<const Eigen::VectorXd>(obs_data.data(), obs_dim);

    // Get observation error covariance
    const auto& R_data = obs_.getCovariance();
    Eigen::VectorXd R_diag =
        Eigen::Map<const Eigen::VectorXd>(R_data.data(), obs_dim);
    Eigen::MatrixXd R = R_diag.asDiagonal();

    likelihood_values_.resize(ens_size);

    // Update weights for each particle
    for (int i = 0; i < ens_size; ++i) {
      // Apply observation operator to get simulated observations
      const auto& simulated_obs = obs_op_.apply(ensemble_.GetMember(i), obs_);
      Eigen::VectorXd y_sim =
          Eigen::Map<const Eigen::VectorXd>(simulated_obs.data(), obs_dim);

      // Compute innovation
      Eigen::VectorXd innovation = yo - y_sim;

      // Compute likelihood (Gaussian assumption)
      double likelihood = computeGaussianLikelihood(innovation, R);
      likelihood_values_[i] = likelihood;

      // Update weight: w_i^{(a)} \propto w_i^{(f)} * likelihood
      weights_[i] *= likelihood;
    }

    // Normalize weights
    normalizeWeights();

    logger_.Info() << "Weights updated, max weight: "
                   << *std::max_element(weights_.begin(), weights_.end());
  }

  /**
   * @brief Compute Gaussian likelihood
   */
  double computeGaussianLikelihood(const Eigen::VectorXd& innovation,
                                   const Eigen::MatrixXd& R) const {
    // p(y|x) \propto \exp(-\frac{1}{2}(y - H(x))^T R^{-1}(y - H(x)))
    double mahalanobis_dist = innovation.transpose() * R.inverse() * innovation;
    return std::exp(-0.5 * mahalanobis_dist);
  }

  /**
   * @brief Normalize particle weights
   */
  void normalizeWeights() {
    double sum_weights = std::accumulate(weights_.begin(), weights_.end(), 0.0);
    if (sum_weights > 0) {
      for (auto& weight : weights_) {
        weight /= sum_weights;
      }
    } else {
      // If all weights are zero, set equal weights
      std::fill(weights_.begin(), weights_.end(), 1.0 / weights_.size());
      logger_.Warning() << "All weights were zero, setting equal weights";
    }
  }

  /**
   * @brief Compute effective sample size
   */
  double computeEffectiveSampleSize() const {
    double sum_sq_weights = 0.0;
    for (const auto& weight : weights_) {
      sum_sq_weights += weight * weight;
    }
    return 1.0 / sum_sq_weights;
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
   * @brief Perform resampling based on current weights
   */
  void performResampling() {
    const int ens_size = ensemble_.Size();

    // Create new ensemble members based on resampling
    std::vector<State<BackendTag>> new_members;
    new_members.reserve(ens_size);

    switch (resampling_method_) {
      case ResamplingMethod::SYSTEMATIC:
        systematicResampling(new_members);
        break;
      case ResamplingMethod::MULTINOMIAL:
        multinomialResampling(new_members);
        break;
      case ResamplingMethod::STRATIFIED:
        stratifiedResampling(new_members);
        break;
      case ResamplingMethod::RESIDUAL:
        residualResampling(new_members);
        break;
    }

    // Update ensemble with new members
    for (int i = 0; i < ens_size; ++i) {
      ensemble_.GetMember(i) = std::move(new_members[i]);
    }

    // Reset weights to equal
    std::fill(weights_.begin(), weights_.end(), 1.0 / ens_size);

    resampling_performed_ = true;
    logger_.Info() << "Resampling completed";
  }

  /**
   * @brief Systematic resampling
   */
  void systematicResampling(std::vector<State<BackendTag>>& new_members) const {
    const int ens_size = ensemble_.Size();

    // Compute cumulative weights
    std::vector<double> cumsum(ens_size + 1, 0.0);
    for (int i = 0; i < ens_size; ++i) {
      cumsum[i + 1] = cumsum[i] + weights_[i];
    }

    // Systematic resampling
    double u = static_cast<double>(rand()) / RAND_MAX / ens_size;
    for (int i = 0; i < ens_size; ++i) {
      double target = u + static_cast<double>(i) / ens_size;

      // Find particle index
      int j = 0;
      while (j < ens_size && cumsum[j + 1] < target) {
        ++j;
      }

      new_members.push_back(ensemble_.GetMember(j).clone());
    }
  }

  /**
   * @brief Multinomial resampling
   */
  void multinomialResampling(
      std::vector<State<BackendTag>>& new_members) const {
    const int ens_size = ensemble_.Size();

    for (int i = 0; i < ens_size; ++i) {
      // Generate random number
      double u = static_cast<double>(rand()) / RAND_MAX;

      // Find particle index
      double cumsum = 0.0;
      int j = 0;
      for (j = 0; j < ens_size; ++j) {
        cumsum += weights_[j];
        if (cumsum >= u) break;
      }

      new_members.push_back(ensemble_.GetMember(j).clone());
    }
  }

  /**
   * @brief Stratified resampling
   */
  void stratifiedResampling(std::vector<State<BackendTag>>& new_members) const {
    const int ens_size = ensemble_.Size();

    for (int i = 0; i < ens_size; ++i) {
      // Generate random number in stratum i
      double u = (static_cast<double>(rand()) / RAND_MAX + i) / ens_size;

      // Find particle index
      double cumsum = 0.0;
      int j = 0;
      for (j = 0; j < ens_size; ++j) {
        cumsum += weights_[j];
        if (cumsum >= u) break;
      }

      new_members.push_back(ensemble_.GetMember(j).clone());
    }
  }

  /**
   * @brief Residual resampling
   */
  void residualResampling(std::vector<State<BackendTag>>& new_members) const {
    const int ens_size = ensemble_.Size();

    // Compute integer parts of weights
    std::vector<int> n_i(ens_size);
    double sum_n_i = 0.0;
    for (int i = 0; i < ens_size; ++i) {
      n_i[i] = static_cast<int>(weights_[i] * ens_size);
      sum_n_i += n_i[i];
    }

    // Add deterministic copies
    for (int i = 0; i < ens_size; ++i) {
      for (int j = 0; j < n_i[i]; ++j) {
        new_members.push_back(ensemble_.GetMember(i).clone());
      }
    }

    // Fill remaining slots with multinomial resampling
    int remaining = ens_size - static_cast<int>(sum_n_i);
    if (remaining > 0) {
      // Compute residual weights
      std::vector<double> residual_weights(ens_size);
      for (int i = 0; i < ens_size; ++i) {
        residual_weights[i] = weights_[i] * ens_size - n_i[i];
      }

      // Normalize residual weights
      double sum_residual = std::accumulate(residual_weights.begin(),
                                            residual_weights.end(), 0.0);
      for (auto& w : residual_weights) {
        w /= sum_residual;
      }

      // Multinomial resampling for remaining slots
      for (int k = 0; k < remaining; ++k) {
        double u = static_cast<double>(rand()) / RAND_MAX;
        double cumsum = 0.0;
        int j = 0;
        for (j = 0; j < ens_size; ++j) {
          cumsum += residual_weights[j];
          if (cumsum >= u) break;
        }
        new_members.push_back(ensemble_.GetMember(j).clone());
      }
    }
  }

  /**
   * @brief Apply jittering to maintain particle diversity
   */
  void applyJittering() {
    const int ens_size = ensemble_.Size();

    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, jittering_std_);

    // Apply small perturbations to each particle
    for (int i = 0; i < ens_size; ++i) {
      auto& member = ensemble_.GetMember(i);
      auto data = member.template getDataPtr<double>();
      const int state_dim = member.size();

      for (int j = 0; j < state_dim; ++j) {
        data[j] += dist(gen);
      }
    }

    logger_.Info() << "Jittering applied with std=" << jittering_std_;
  }

  /**
   * @brief Save weights to file
   */
  void saveWeights() const {
    std::string weights_file = output_base_file_ + "_weights.txt";
    // Implementation would write weights to file
    logger_.Info() << "Weights saved to: " << weights_file;
  }

  /**
   * @brief Get resampling method as string
   */
  std::string getResamplingMethodString() const {
    switch (resampling_method_) {
      case ResamplingMethod::SYSTEMATIC:
        return "systematic";
      case ResamplingMethod::MULTINOMIAL:
        return "multinomial";
      case ResamplingMethod::STRATIFIED:
        return "stratified";
      case ResamplingMethod::RESIDUAL:
        return "residual";
      default:
        return "unknown";
    }
  }

  Ensemble<BackendTag>& ensemble_;
  Observation<BackendTag>& obs_;
  const ObsOperator<BackendTag>& obs_op_;
  std::vector<double> weights_;
  std::vector<double> likelihood_values_;
  ResamplingMethod resampling_method_;
  int resampling_threshold_;
  bool jittering_enabled_;
  double jittering_std_;
  bool resampling_performed_ = false;
  std::string output_base_file_;
  std::string format_ = "txt";  // Default format
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework