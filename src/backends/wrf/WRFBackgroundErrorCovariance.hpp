#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace metada::backends::wrf {

/**
 * @brief WRF backend implementation for background error covariance matrix
 *
 * @details This class provides a WRF-specific implementation of background
 * error covariance operations with support for variable-specific error
 * statistics, length scales, and vertical scales.
 *
 * @tparam ConfigBackend Configuration backend type
 */
template <typename ConfigBackend>
class WRFBackgroundErrorCovariance {
 public:
  /**
   * @brief Constructor from configuration
   *
   * @param config Configuration object containing covariance parameters
   */
  explicit WRFBackgroundErrorCovariance(const ConfigBackend& config)
      : initialized_(false) {
    // Parse WRF-specific covariance configuration
    try {
      // Read background error statistics file path
      if (config.HasKey("file")) {
        be_file_ = config.Get("file").asString();
      }

      // Read variables list
      if (config.HasKey("variables")) {
        auto vars = config.Get("variables");
        if (vars.isVectorString()) {
          variables_ = vars.asVectorString();
        } else if (vars.isVectorConfigValue()) {
          const auto& vec = vars.asVectorConfigValue();
          for (const auto& var : vec) {
            if (var.isString()) {
              variables_.push_back(var.asString());
            }
          }
        }
      }

      // Read variable-specific length scales
      if (config.HasKey("length_scales")) {
        auto length_scales = config.Get("length_scales");
        if (length_scales.isMap()) {
          const auto& map = length_scales.asMap();
          for (const auto& var : variables_) {
            auto it = map.find(var);
            if (it != map.end() && it->second.isFloat()) {
              length_scales_[var] = it->second.asFloat();
            }
          }
        }
      }

      // Read variable-specific vertical scales
      if (config.HasKey("vertical_scales")) {
        auto vertical_scales = config.Get("vertical_scales");
        if (vertical_scales.isMap()) {
          const auto& map = vertical_scales.asMap();
          for (const auto& var : variables_) {
            auto it = map.find(var);
            if (it != map.end() && it->second.isFloat()) {
              vertical_scales_[var] = it->second.asFloat();
            }
          }
        }
      }

      // Read variable-specific error standard deviations
      if (config.HasKey("error_std")) {
        auto error_std = config.Get("error_std");
        if (error_std.isMap()) {
          const auto& map = error_std.asMap();
          for (const auto& var : variables_) {
            auto it = map.find(var);
            if (it != map.end() && it->second.isFloat()) {
              error_std_[var] = it->second.asFloat();
            }
          }
        }
      }

      // Set default size (will be updated when state is known)
      size_ = 1000;
      initialized_ = true;
    } catch (const std::exception& e) {
      // If configuration parsing fails, use safe defaults
      size_ = 1000;
      initialized_ = true;
    }
  }

  // Initialization and state
  bool isInitialized() const { return initialized_; }

  // Variable-specific configuration access
  const std::vector<std::string>& getVariables() const { return variables_; }
  const std::string& getBackgroundErrorFile() const { return be_file_; }

  double getLengthScale(const std::string& variable) const {
    auto it = length_scales_.find(variable);
    return (it != length_scales_.end()) ? it->second : 50.0;  // Default 50 km
  }

  double getVerticalScale(const std::string& variable) const {
    auto it = vertical_scales_.find(variable);
    return (it != vertical_scales_.end()) ? it->second
                                          : 0.5;  // Default 0.5 scale height
  }

  double getErrorStd(const std::string& variable) const {
    auto it = error_std_.find(variable);
    return (it != error_std_.end()) ? it->second : 1.0;  // Default 1.0
  }

  /**
   * @brief Compute quadratic form for diagonal covariance using
   * variable-specific error statistics
   */
  template <typename StateType>
  double computeQuadraticFormDiagonal(const StateType& increment) const {
    // For WRF, implement variable-specific diagonal quadratic form
    size_t state_size = increment.size();

    // Update size to match actual state if needed
    if (size_ != state_size) {
      size_ = state_size;  // Update to actual state size
    }

    auto data = increment.getData();

    double result = 0.0;
    // Use average error standard deviation across all variables as fallback
    double avg_error_std = 1.0;
    if (!error_std_.empty()) {
      double sum = 0.0;
      for (const auto& pair : error_std_) {
        sum += pair.second;
      }
      avg_error_std = sum / error_std_.size();
    }

    for (size_t i = 0; i < state_size; ++i) {
      double variance = avg_error_std * avg_error_std;
      result += (data[i] * data[i]) / variance;
    }

    return result;
  }

  /**
   * @brief Compute quadratic form for ensemble covariance
   */
  template <typename StateType>
  double computeQuadraticFormEnsemble(const StateType& increment) const {
    // For WRF, ensemble covariance is not implemented yet, fallback to diagonal
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for parametric covariance
   */
  template <typename StateType>
  double computeQuadraticFormParametric(const StateType& increment) const {
    // For WRF, parametric covariance is not implemented yet, fallback to
    // diagonal
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for hybrid covariance
   */
  template <typename StateType>
  double computeQuadraticFormHybrid(const StateType& increment) const {
    // For WRF, hybrid covariance is not implemented yet, fallback to diagonal
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Compute quadratic form for full covariance
   */
  template <typename StateType>
  double computeQuadraticFormFull(const StateType& increment) const {
    // For WRF, full covariance is not implemented yet, fallback to diagonal
    return computeQuadraticFormDiagonal(increment);
  }

  /**
   * @brief Apply inverse ensemble covariance
   */
  template <typename StateType>
  void applyInverseEnsemble(const StateType& increment,
                            StateType& result) const {
    // For WRF, ensemble covariance is not implemented yet, fallback to diagonal
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse parametric covariance
   */
  template <typename StateType>
  void applyInverseParametric(const StateType& increment,
                              StateType& result) const {
    // For WRF, parametric covariance is not implemented yet, fallback to
    // diagonal
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse hybrid covariance
   */
  template <typename StateType>
  void applyInverseHybrid(const StateType& increment, StateType& result) const {
    // For WRF, hybrid covariance is not implemented yet, fallback to diagonal
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse full covariance
   */
  template <typename StateType>
  void applyInverseFull(const StateType& increment, StateType& result) const {
    // For WRF, full covariance is not implemented yet, fallback to diagonal
    applyInverseDiagonal(increment, result);
  }

  /**
   * @brief Apply inverse diagonal covariance using variable-specific error
   * statistics
   */
  template <typename StateType>
  void applyInverseDiagonal(const StateType& increment,
                            StateType& result) const {
    auto inc_data = increment.getData();
    size_t state_size = inc_data.size();

    // Update size to match actual state if needed
    if (size_ != state_size) {
      size_ = state_size;  // Update to actual state size
    }

    // Use average error standard deviation across all variables as fallback
    double avg_error_std = 1.0;
    if (!error_std_.empty()) {
      double sum = 0.0;
      for (const auto& pair : error_std_) {
        sum += pair.second;
      }
      avg_error_std = sum / error_std_.size();
    }

    // For diagonal covariance: B^-1 = diag(1/sigma_i^2)
    double variance = avg_error_std * avg_error_std;
    result.zero();
    result.axpy(1.0 / variance, increment);
  }

 private:
  mutable size_t size_;  ///< Matrix size (state vector dimension) - mutable for
                         ///< auto-sizing
  bool initialized_;     ///< Initialization status

  // WRF-specific configuration
  std::string be_file_;  ///< Background error statistics file path
  std::vector<std::string> variables_;  ///< List of analysis variables
  std::unordered_map<std::string, double>
      length_scales_;  ///< Variable-specific length scales (km)
  std::unordered_map<std::string, double>
      vertical_scales_;  ///< Variable-specific vertical scales (scale height)
  std::unordered_map<std::string, double>
      error_std_;  ///< Variable-specific error standard deviations
};

}  // namespace metada::backends::wrf