#pragma once

#include <cmath>
#include <vector>

namespace metada::framework {

/**
 * @brief Class for calculating ensemble data assimilation metrics
 *
 * This class provides methods for calculating various metrics used to evaluate
 * the performance of ensemble data assimilation algorithms. The metrics
 * include:
 * - RMSE (Root Mean Square Error)
 * - Ensemble spread (standard deviation)
 * - Bias
 * - Correlation coefficient
 * - CRPS (Continuous Ranked Probability Score)
 *
 * @tparam T The floating-point type to use for calculations (float or double)
 */
template <typename T = double>
class Metrics {
 public:
  /**
   * @brief Calculate ensemble mean
   *
   * Formula:
   * \f[ \bar{x} = \frac{1}{N} \sum_{i=1}^N x_i \f]
   * where N is ensemble size and x_i is the i-th ensemble member
   *
   * @param ensemble_data Vector of ensemble member data
   * @param state_dim Dimension of the state vector
   * @param ens_size Number of ensemble members
   * @return Vector containing the mean at each point
   */
  static std::vector<T> calculateEnsembleMean(
      const std::vector<std::vector<T>>& ensemble_data, size_t state_dim,
      size_t ens_size) {
    std::vector<T> mean(state_dim, 0.0);
    for (const auto& member : ensemble_data) {
      for (size_t i = 0; i < state_dim; ++i) {
        mean[i] += member[i];
      }
    }
    for (auto& val : mean) {
      val /= ens_size;
    }
    return mean;
  }

  /**
   * @brief Calculate ensemble spread (standard deviation)
   *
   * Formula:
   * \f[ \sigma = \sqrt{\frac{1}{N-1} \sum_{i=1}^N (x_i - \bar{x})^2} \f]
   * where N is ensemble size, x_i is the i-th ensemble member, and
   * \f$\bar{x}\f$ is ensemble mean
   *
   * @param ensemble_data Vector of ensemble member data
   * @param mean Ensemble mean
   * @param state_dim Dimension of the state vector
   * @param ens_size Number of ensemble members
   * @return Vector containing the spread at each point
   */
  static std::vector<T> calculateEnsembleSpread(
      const std::vector<std::vector<T>>& ensemble_data,
      const std::vector<T>& mean, size_t state_dim, size_t ens_size) {
    std::vector<T> spread(state_dim, 0.0);
    for (const auto& member : ensemble_data) {
      for (size_t i = 0; i < state_dim; ++i) {
        T diff = member[i] - mean[i];
        spread[i] += diff * diff;
      }
    }
    for (auto& val : spread) {
      val = std::sqrt(val / (ens_size - 1));
    }
    return spread;
  }

  /**
   * @brief Calculate bias
   *
   * Formula:
   * \f[ \text{bias} = \frac{1}{n} \sum_{i=1}^n (\bar{x}_i - x^t_i) \f]
   * where n is state dimension, \f$\bar{x}\f$ is ensemble mean, and \f$x^t\f$
   * is true state
   *
   * @param mean Ensemble mean
   * @param truth True state
   * @param state_dim Dimension of the state vector
   * @return Average bias
   */
  static T calculateBias(const std::vector<T>& mean,
                         const std::vector<T>& truth, size_t state_dim) {
    T bias = 0.0;
    for (size_t i = 0; i < state_dim; ++i) {
      bias += mean[i] - truth[i];
    }
    return bias / state_dim;
  }

  /**
   * @brief Calculate correlation coefficient
   *
   * Formula:
   * \f[ r = \frac{\sum_{i=1}^n (\bar{x}_i - \bar{\bar{x}})(x^t_i - \bar{x^t})}
   *           {\sqrt{\sum_{i=1}^n (\bar{x}_i - \bar{\bar{x}})^2}
   *            \sqrt{\sum_{i=1}^n (x^t_i - \bar{x^t})^2}} \f]
   * where n is state dimension
   *
   * @param mean Ensemble mean
   * @param truth True state
   * @param state_dim Dimension of the state vector
   * @return Correlation coefficient
   */
  static T calculateCorrelation(const std::vector<T>& mean,
                                const std::vector<T>& truth, size_t state_dim) {
    T sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0;
    T sum_x2 = 0.0, sum_y2 = 0.0;

    for (size_t i = 0; i < state_dim; ++i) {
      sum_x += mean[i];
      sum_y += truth[i];
      sum_xy += mean[i] * truth[i];
      sum_x2 += mean[i] * mean[i];
      sum_y2 += truth[i] * truth[i];
    }

    T n = static_cast<T>(state_dim);
    T numerator = n * sum_xy - sum_x * sum_y;
    T denominator =
        std::sqrt((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y));

    return numerator / denominator;
  }

  /**
   * @brief Calculate Continuous Ranked Probability Score (CRPS)
   *
   * Formula:
   * \f[ \text{CRPS} = \frac{1}{n} \sum_{i=1}^n \int_{-\infty}^{\infty}
   *                    (F_i(y) - H(y-x^t_i))^2 dy \f]
   * where F_i is the empirical CDF of ensemble at point i, and H is Heaviside
   * function
   *
   * For discrete ensemble, this reduces to:
   * \f[ \text{CRPS} = \frac{1}{n} \sum_{i=1}^n \left[
   *     \frac{1}{N} \sum_{j=1}^N |y_{ij} - x^t_i| -
   *     \frac{1}{2N^2} \sum_{j=1}^N \sum_{k=1}^N |y_{ij} - y_{ik}|
   *     \right] \f]
   * where N is ensemble size, y_{ij} is j-th ensemble member at point i
   *
   * @param ensemble_data Vector of ensemble member data
   * @param truth True state
   * @param state_dim Dimension of the state vector
   * @param ens_size Number of ensemble members
   * @return CRPS value
   */
  static T calculateCRPS(const std::vector<std::vector<T>>& ensemble_data,
                         const std::vector<T>& truth, size_t state_dim,
                         size_t ens_size) {
    T crps = 0.0;
    for (size_t i = 0; i < state_dim; ++i) {
      T sum_diff = 0.0;
      for (size_t j = 0; j < ens_size; ++j) {
        for (size_t k = 0; k < ens_size; ++k) {
          sum_diff += std::abs(ensemble_data[j][i] - ensemble_data[k][i]);
        }
      }
      T sum_truth_diff = 0.0;
      for (size_t j = 0; j < ens_size; ++j) {
        sum_truth_diff += std::abs(ensemble_data[j][i] - truth[i]);
      }
      crps +=
          sum_truth_diff / ens_size - sum_diff / (2.0 * ens_size * ens_size);
    }
    return crps / state_dim;
  }

  /**
   * @brief Calculate Root Mean Square Error (RMSE)
   *
   * Formula:
   * \f[ \text{RMSE} = \sqrt{\frac{1}{n} \sum_{i=1}^n (\bar{x}_i - x^t_i)^2} \f]
   * where n is state dimension
   *
   * @param mean Ensemble mean
   * @param truth True state
   * @param state_dim Dimension of the state vector
   * @return RMSE value
   */
  static T calculateRMSE(const std::vector<T>& mean,
                         const std::vector<T>& truth, size_t state_dim) {
    T rmse = 0.0;
    for (size_t i = 0; i < state_dim; ++i) {
      T diff = mean[i] - truth[i];
      rmse += diff * diff;
    }
    return std::sqrt(rmse / state_dim);
  }

  /**
   * @brief Calculate average spread
   *
   * Formula:
   * \f[ \bar{\sigma} = \frac{1}{n} \sum_{i=1}^n \sigma_i \f]
   * where n is state dimension and \f$\sigma_i\f$ is spread at point i
   *
   * @param spread Vector of spread values
   * @param state_dim Dimension of the state vector
   * @return Average spread value
   */
  static T calculateAverageSpread(const std::vector<T>& spread,
                                  size_t state_dim) {
    T avg_spread = 0.0;
    for (T val : spread) {
      avg_spread += val;
    }
    return avg_spread / state_dim;
  }
};

}  // namespace metada::framework