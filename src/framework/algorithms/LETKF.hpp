#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

#include "Config.hpp"
#include "Ensemble.hpp"
#include "Logger.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "PointObservation.hpp"

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
        output_base_file_(config.Get("output_base_file").asString()) {
    logger_.Debug() << "LETKF constructed with radius " << localization_radius_;
  }

  /**
   * @brief Perform the LETKF analysis step, updating the ensemble point by
   * point.
   */
  void Analyse() {
    logger_.Debug() << "LETKF analysis started";
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int ens_size = ensemble_.Size();
    const int state_dim = ensemble_.GetMember(0).size();

    // Get ensemble data
    std::vector<VectorXd> member_data(ens_size);
    for (int i = 0; i < ens_size; ++i) {
      const auto data_ptr =
          ensemble_.GetMember(i).template getDataPtr<double>();
      member_data[i] = Eigen::Map<const VectorXd>(data_ptr, state_dim);
    }

    // Get observation data and locations
    const auto obs_data = obs_.template getData<std::vector<double>>();
    VectorXd yo = Eigen::Map<const VectorXd>(obs_data.data(), obs_.size());

    std::vector<ObservationLocation> obs_locations;
    for (const auto& obs_point : obs_) {
      obs_locations.push_back(obs_point.location);
    }

    // Retrieve geometry from the first ensemble member's state
    const auto* geometry = ensemble_.GetMember(0).geometry();
    if (!geometry) {
      throw std::runtime_error("Geometry pointer is null in LETKF");
    }

    // Update each grid point locally using geometry iterator
    for (const auto& grid_point : *geometry) {
      updateGridPoint(grid_point, member_data, yo, obs_locations, ens_size);
    }

    // Update ensemble members
    for (int i = 0; i < ens_size; ++i) {
      auto& member = ensemble_.GetMember(i);
      auto data_ptr = member.template getDataPtr<double>();
      Eigen::Map<VectorXd>(data_ptr, state_dim) = member_data[i];
    }

    logger_.Debug() << "LETKF analysis completed";
  }

  /**
   * @brief Save the analyzed ensemble to files
   * @param base_filename Base filename for saving (without extension)
   */
  void saveEnsemble() const {
    logger_.Debug() << "LETKF saving ensemble";
    const int ens_size = ensemble_.Size();

    // Save individual ensemble members
    for (int i = 0; i < ens_size; ++i) {
      const auto& member = ensemble_.GetMember(i);
      member.saveToFile(output_base_file_ + "_member_" + std::to_string(i) +
                        ".txt");
    }
    logger_.Debug() << "LETKF ensemble saved";
  }

 private:
  /**
   * @brief Update a single grid point using local observations
   * @param grid_point The grid point object to update
   * @param member_data Vector of ensemble member data
   * @param yo Observation vector
   * @param obs_locations Vector of observation locations
   * @param ens_size Ensemble size
   */
  template <typename GridPointType>
  void updateGridPoint(const GridPointType& grid_point,
                       std::vector<Eigen::VectorXd>& member_data,
                       const Eigen::VectorXd& yo,
                       const std::vector<ObservationLocation>& obs_locations,
                       int ens_size) {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    int x = grid_point.first;
    int y = grid_point.second;

    // Calculate linear index from 2D coordinates
    // This assumes row-major ordering: index = y * nx + x
    // We need to get the grid dimensions from the geometry
    const auto* geometry = ensemble_.GetMember(0).geometry();
    if (!geometry) {
      throw std::runtime_error("Geometry pointer is null in updateGridPoint");
    }

    // Calculate grid dimensions from geometry size
    // For a 2D grid, size = nx * ny, so we can estimate nx from size
    // This is a reasonable approximation for most regular grids
    size_t total_size = geometry->size();
    int nx = static_cast<int>(std::sqrt(total_size));  // Approximate nx
    int grid_index = y * nx + x;

    // Ensure grid_index is within bounds
    if (grid_index >= static_cast<int>(total_size)) {
      logger_.Warning() << "Grid index " << grid_index
                        << " out of bounds for size " << total_size;
      return;
    }

    // Find local observations using 2D Euclidean distance
    std::vector<int> local_obs_indices;
    for (size_t i = 0; i < obs_locations.size(); ++i) {
      double obs_x = obs_locations[i].longitude;
      double obs_y = obs_locations[i].latitude;
      double distance =
          std::sqrt((x - obs_x) * (x - obs_x) + (y - obs_y) * (y - obs_y));
      if (distance <= localization_radius_) {
        local_obs_indices.push_back(i);
      }
    }

    if (local_obs_indices.empty()) return;

    // Build local quantities
    int local_obs_dim = local_obs_indices.size();
    VectorXd yo_local(local_obs_dim);
    MatrixXd H_local(local_obs_dim, ens_size);

    for (int i = 0; i < local_obs_dim; ++i) {
      yo_local(i) = yo(local_obs_indices[i]);
      // Simplified local observation operator
      for (int j = 0; j < ens_size; ++j) {
        H_local(i, j) = member_data[j](local_obs_indices[i]);
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
      xb_local(i) = member_data[i](grid_index);  // Use grid_index
    }

    VectorXd xb_mean_local = xb_local.mean() * VectorXd::Ones(ens_size);
    VectorXd xb_pert_local = xb_local - xb_mean_local;

    VectorXd xa_mean_local =
        xb_mean_local + xb_pert_local.asDiagonal() * wa_local;
    MatrixXd Xa_pert_local = xb_pert_local.asDiagonal() * Wa_local;
    VectorXd xa_local = xa_mean_local + Xa_pert_local.rowwise().sum();

    for (int i = 0; i < ens_size; ++i) {
      member_data[i](grid_index) = xa_local(i);  // Use grid_index
    }
  }

  Ensemble<BackendTag>& ensemble_;
  Observation<BackendTag>& obs_;
  const ObsOperator<BackendTag>& obs_op_;
  double inflation_;
  double localization_radius_;
  std::string output_base_file_ = "analysis";
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework