#pragma once
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "Location.hpp"
#include "Logger.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "PointObservation.hpp"
#include "ProgressBar.hpp"
#include "State.hpp"

namespace metada::framework {

/**
 * @brief Local Ensemble Transform Kalman Filter (LETKF) implementation.
 *
 * This class implements the LETKF algorithm for ensemble-based data
 * assimilation with localization. It updates each grid point using only
 * observations within a local region, reducing computational cost and
 * improving localization.
 *
 * The LETKF algorithm follows these steps for each grid point:
 * 1. Find observations within local radius
 * 2. Build local forecast quantities for the grid point
 * 3. Compute local innovation using nearby observations
 * 4. Perform local ETKF analysis
 * 5. Update the grid point with local analysis
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
   * @param config Configuration object containing localization parameters.
   */
  LETKF(Ensemble<BackendTag>& ensemble, Observation<BackendTag>& obs,
        const ObsOperator<BackendTag>& obs_op, const Config<BackendTag>& config)
      : ensemble_(ensemble),
        obs_(obs),
        obs_op_(obs_op),
        inflation_(config.Get("inflation").asFloat()),
        localization_radius_(config.Get("localization_radius").asFloat()),
        output_base_file_(config.Get("output_base_file").asString()),
        format_(config.Get("format").asString()) {
    logger_.Info() << "LETKF constructed with radius " << localization_radius_;
  }

  /**
   * @brief Perform the LETKF analysis step, updating the ensemble point by
   * point.
   */
  void Analyse() {
    logger_.Info() << "LETKF analysis started";
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int ens_size = ensemble_.Size();

    // Get observation data and locations
    const auto obs_data = obs_.template getData<std::vector<double>>();
    VectorXd yo = Eigen::Map<const VectorXd>(obs_data.data(), obs_.size());

    std::vector<Location> obs_locations;
    for (const auto& obs_point : obs_) {
      obs_locations.push_back(obs_point.location);
    }

    // Retrieve geometry from the first ensemble member's state
    const auto* geometry = ensemble_.GetMember(0).geometry();
    if (!geometry) {
      throw std::runtime_error("Geometry pointer is null in LETKF");
    }

    // Update each grid point locally using geometry iterator
    auto geometry_iter = geometry->begin();
    auto geometry_end = geometry->end();

    // Calculate total number of grid points for progress tracking
    size_t total_points = std::distance(geometry_iter, geometry_end);
    size_t current_point = 0;
    size_t progress_interval = ProgressBar::CalculateInterval(total_points);

    // Configure progress bar for LETKF
    ProgressBar::Config config;
    config.prefix = "LETKF progress: ";

    logger_.Info() << "Starting LETKF grid point updates for " << total_points
                   << " points";

    for (const auto& grid_point : *geometry) {
      updateGridPoint(grid_point, yo, obs_locations, ens_size);

      ++current_point;

      // Log progress at intervals using ProgressBar
      if (ProgressBar::ShouldLog(current_point, total_points,
                                 progress_interval)) {
        std::string progress_bar =
            ProgressBar::Create(current_point, total_points, config);
        logger_.Info() << progress_bar;
      }
    }

    // Compute analysis mean using Ensemble's Mean() method
    ensemble_.RecomputeMean();

    logger_.Info() << "LETKF analysis completed";
  }

  /**
   * @brief Save the analyzed ensemble to files
   * @param base_filename Base filename for saving (without extension)
   */
  void saveEnsemble() const {
    logger_.Info() << "LETKF saving ensemble";
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
    logger_.Info() << "LETKF ensemble saved";
  }

 private:
  /**
   * @brief Update a single grid point using local observations
   * @param grid_point The grid point object to update
   * @param yo Observation vector
   * @param obs_locations Vector of observation locations
   * @param ens_size Ensemble size
   */
  void updateGridPoint(const Location& grid_point, const Eigen::VectorXd& yo,
                       const std::vector<Location>& obs_locations,
                       int ens_size) {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // Find local observations using Location::distance_to
    std::vector<int> local_obs_indices;
    for (size_t i = 0; i < obs_locations.size(); ++i) {
      double distance = grid_point.distance_to(obs_locations[i]);
      if (distance <= localization_radius_) {
        local_obs_indices.push_back(i);
      }
    }

    if (local_obs_indices.empty()) {
      // No local observations found - apply inflation to maintain ensemble
      // spread

      // Get current ensemble values at this grid point
      VectorXd xb_local(ens_size);
      for (int i = 0; i < ens_size; ++i) {
        xb_local(i) = ensemble_.GetMember(i).at(grid_point);
      }

      // Compute ensemble mean and perturbations
      double xb_mean = xb_local.mean();
      VectorXd xb_pert = xb_local.array() - xb_mean;

      // Apply inflation to perturbations
      xb_pert *= std::sqrt(inflation_);

      // Update ensemble members with inflated values
      for (int i = 0; i < ens_size; ++i) {
        ensemble_.GetMember(i).at(grid_point) = xb_mean + xb_pert(i);
      }

      return;
    }

    // Build local quantities
    int local_obs_dim = local_obs_indices.size();
    VectorXd yo_local(local_obs_dim);
    MatrixXd H_local(local_obs_dim, ens_size);

    for (int i = 0; i < local_obs_dim; ++i) {
      yo_local(i) = yo(local_obs_indices[i]);
      // Use observation operator to apply H to each ensemble member
      for (int j = 0; j < ens_size; ++j) {
        // Apply observation operator to get observation space values
        const auto& obs_data = obs_op_.apply(ensemble_.GetMember(j), obs_);
        // Extract the value for the specific observation location
        H_local(i, j) = obs_data[local_obs_indices[i]];
      }
    }

    // Compute local forecast mean and anomalies
    VectorXd yb_mean = H_local.rowwise().mean();
    MatrixXd Yb_pert = H_local.colwise() - yb_mean;
    VectorXd d_local = yo_local - yb_mean;

    // Local ETKF analysis
    MatrixXd R_local = MatrixXd::Identity(local_obs_dim, local_obs_dim);
    MatrixXd Pa_local =
        (Yb_pert.transpose() * R_local.inverse() * Yb_pert +
         (ens_size - 1) * MatrixXd::Identity(ens_size, ens_size))
            .inverse();
    Pa_local *= inflation_;

    VectorXd wa_local =
        Pa_local * Yb_pert.transpose() * R_local.inverse() * d_local;
    MatrixXd Pa_sqrt = Pa_local.llt().matrixL();
    MatrixXd Wa_local = std::sqrt(ens_size - 1) * Pa_sqrt;

    // Update grid point
    VectorXd xb_local(ens_size);
    for (int i = 0; i < ens_size; ++i) {
      xb_local(i) = ensemble_.GetMember(i).at(grid_point);
    }

    VectorXd xb_mean_local = xb_local.mean() * VectorXd::Ones(ens_size);
    VectorXd xb_pert_local = xb_local - xb_mean_local;

    VectorXd xa_mean_local =
        xb_mean_local + xb_pert_local.asDiagonal() * wa_local;
    MatrixXd Xa_pert_local = xb_pert_local.asDiagonal() * Wa_local;
    VectorXd xa_local = xa_mean_local + Xa_pert_local.rowwise().sum();

    for (int i = 0; i < ens_size; ++i) {
      ensemble_.GetMember(i).at(grid_point) = xa_local(i);
    }
  }

  Ensemble<BackendTag>& ensemble_;
  Observation<BackendTag>& obs_;
  const ObsOperator<BackendTag>& obs_op_;
  double inflation_;
  double localization_radius_;
  std::string output_base_file_ = "analysis";
  std::string format_ = "nc";  // Default format
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework