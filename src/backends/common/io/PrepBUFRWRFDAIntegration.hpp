#pragma once

#include <string>
#include <vector>

#include "../../common/obsoperator/WRFDAObsOperator_c_api.h"
#include "PrepBUFRObsAdapter.hpp"

namespace metada::backends::common::io {

struct GridArrays {
  int nx{}, ny{}, nz{};
  const double* u{};
  const double* v{};
  const double* t{};
  const double* q{};
  const double* psfc{};
  const double* lats2d{};
  const double* lons2d{};
  const double* levels{};
};

class PrepBUFRWRFDAIntegration {
 public:
  static int applySurface(const std::string& token, const GridArrays& grid,
                          const SurfaceBatch& batch,
                          std::vector<double>& out_y) {
    out_y.resize(batch.lats.size());
    return wrfda_xtoy_apply_grid(
        token.c_str(), grid.nx, grid.ny, grid.nz, grid.u, grid.v, grid.t,
        grid.q, grid.psfc, grid.lats2d, grid.lons2d, grid.levels,
        static_cast<int>(batch.lats.size()), batch.lats.data(),
        batch.lons.data(), batch.levels.data(), out_y.data());
  }

  static int adjointSurface(const std::string& token, const GridArrays& grid,
                            const SurfaceBatch& batch,
                            const std::vector<double>& delta_y, double* inout_u,
                            double* inout_v, double* inout_t, double* inout_q,
                            double* inout_psfc) {
    return wrfda_xtoy_adjoint_grid(
        token.c_str(), grid.nx, grid.ny, grid.nz, delta_y.data(), grid.lats2d,
        grid.lons2d, grid.levels, static_cast<int>(batch.lats.size()),
        batch.lats.data(), batch.lons.data(), batch.levels.data(), inout_u,
        inout_v, inout_t, inout_q, inout_psfc);
  }

  static int applyProfiles(const std::string& token, const GridArrays& grid,
                           const ProfileBatch& batch,
                           std::vector<double>& out_y_flat) {
    out_y_flat.resize(batch.levels_flat.size());
    return wrfda_xtoy_apply_profiles(
        token.c_str(), grid.nx, grid.ny, grid.nz, grid.u, grid.v, grid.t,
        grid.q, grid.psfc, grid.lats2d, grid.lons2d, grid.levels,
        static_cast<int>(batch.counts.size()), batch.counts.data(),
        batch.lats.data(), batch.lons.data(), batch.levels_flat.data(),
        out_y_flat.data());
  }

  static int adjointProfiles(const std::string& token, const GridArrays& grid,
                             const ProfileBatch& batch,
                             const std::vector<double>& delta_y_flat,
                             double* inout_u, double* inout_v, double* inout_t,
                             double* inout_q, double* inout_psfc) {
    return wrfda_xtoy_adjoint_profiles(
        token.c_str(), grid.nx, grid.ny, grid.nz, delta_y_flat.data(),
        grid.lats2d, grid.lons2d, grid.levels,
        static_cast<int>(batch.counts.size()), batch.counts.data(),
        batch.lats.data(), batch.lons.data(), batch.levels_flat.data(), inout_u,
        inout_v, inout_t, inout_q, inout_psfc);
  }
};

}  // namespace metada::backends::common::io
