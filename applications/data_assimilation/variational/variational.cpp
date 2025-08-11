/**
 * @file variational.cpp
 * @brief Driver program for Variational Data Assimilation (4DVAR/3DVAR/FGAT)
 * @details This application implements variational data assimilation algorithms
 *          including 4DVAR, 3DVAR, and FGAT. It reads configuration from a
 * file, initializes required components like background state, observations,
 *          observation operators, model, and background error covariance, and
 *          performs the variational analysis.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *            Expected format: variational <config_file>
 * @return 0 on success, 1 on failure
 */

#include "Variational.hpp"

#include "ApplicationContext.hpp"
#include "BackgroundErrorCovariance.hpp"
#include "Config.hpp"
#include "Geometry.hpp"
#include "Model.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "State.hpp"
#include "WRFBackendTraits.hpp"
// PrepBUFR + WRFDA integration helpers
#include "../src/backends/common/io/BufrObsIO.hpp"
#include "../src/backends/common/io/PrepBUFRObsAdapter.hpp"
#include "../src/backends/common/io/PrepBUFRWRFDAIntegration.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::WRFBackendTag;

int main(int argc, char** argv) {
  try {
    // Validate command line arguments
    if (argc != 2) {
      std::cerr << "Usage: variational <config_file>" << std::endl;
      return 1;
    }

    // Initialize application context
    auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info() << "Variational data assimilation application starting...";

    // Get variational type from configuration
    std::string var_type = config.Get("variational_type").asString();
    logger.Info() << "Variational method: " << var_type;

    // Initialize geometry
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    // Optional: obtain obs families to assimilate from YAML (e.g., ["metar:t",
    // "synop:p"]) and log
    std::vector<std::string> obs_families;
    try {
      auto fam_val = config.Get("obs_families");
      if (fam_val.isVectorString()) {
        obs_families = fam_val.asVectorString();
      } else if (fam_val.isString()) {
        obs_families.push_back(fam_val.asString());
      }
    } catch (...) {
      // key not present; proceed with defaults
    }
    if (!obs_families.empty()) {
      logger.Info() << "Obs families from config:";
      for (const auto& f : obs_families) logger.Info() << "  - " << f;
    }

    // Initialize background state
    fwk::State<BackendTag> background(config.GetSubsection("background"),
                                      geometry);
    logger.Info() << "Background state: " << background;

    // Initialize model
    fwk::Model<BackendTag> model(config.GetSubsection("model"));

    // Initialize observations
    // For 3DVAR, we have a single observation time with multiple types
    std::vector<fwk::Observation<BackendTag>> observations;
    observations.emplace_back(config.GetSubsection("observations"));

    logger.Info() << "Loaded observations for " << var_type << " analysis";

    // Initialize observation operators
    // For 3DVAR, we need a single observation operator matching the observation
    // structure
    std::vector<fwk::ObsOperator<BackendTag>> obs_operators;
    obs_operators.emplace_back(config.GetSubsection("obs_operator"));

    logger.Info() << "Loaded observation operator for " << var_type
                  << " analysis";

    // Optional example: Read PrepBUFR, group, and call WRFDA H/HT per
    // configured families
    try {
      auto bufr_cfg = config.GetSubsection("prepbufr");
      using CfgBackend =
          metada::traits::BackendTraits<BackendTag>::ConfigBackend;
      metada::backends::io::BufrObsIO<CfgBackend> bufr(
          std::move(bufr_cfg.backend()));
      auto bufr_records = bufr.read();
      auto batches =
          metada::backends::common::io::PrepBUFRObsAdapter::groupByFamily(
              bufr_records);

      // Build simple grid arrays from geometry/background (placeholder example)
      const auto& simpleGeom = geometry.backend();
      const int nx = static_cast<int>(simpleGeom.x_dim());
      const int ny = static_cast<int>(simpleGeom.y_dim());
      int nz = 1;
      try {
        nz = config.Get("nz").asInt();
      } catch (...) {
        nz = 1;
      }
      const size_t nxy = static_cast<size_t>(nx) * ny;
      const size_t nxyz = nxy * nz;
      static std::vector<double> u(nxyz, 0.0), v(nxyz, 0.0), t(nxyz, 0.0),
          q(nxyz, 0.0);
      static std::vector<double> psfc(nxy, 0.0);
      static std::vector<double> lats2d(nxy, 0.0), lons2d(nxy, 0.0);
      static std::vector<double> levels(nz, 0.0);
      // Fill trivial lat/lon and levels for demo purposes
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          size_t idx = static_cast<size_t>(i) + static_cast<size_t>(j) * nx;
          lats2d[idx] = static_cast<double>(j);
          lons2d[idx] = static_cast<double>(i);
        }
      }
      for (int k = 0; k < nz; ++k)
        levels[static_cast<size_t>(k)] = static_cast<double>(k + 1);

      metada::backends::common::io::GridArrays grid{
          nx,           ny,       nz,          u.data(),      v.data(),
          t.data(),     q.data(), psfc.data(), lats2d.data(), lons2d.data(),
          levels.data()};

      if (!obs_families.empty()) {
        logger.Info() << "Calling WRFDA H/HT for configured families...";
        for (const auto& token : obs_families) {
          std::vector<double> y_out;
          if (token.rfind("airep:", 0) == 0) {
            metada::backends::common::io::PrepBUFRWRFDAIntegration::
                applyProfiles(token, grid, batches.airep, y_out);
            logger.Info() << "H(x) airep size=" << y_out.size();
          } else if (token.rfind("pilot:", 0) == 0) {
            metada::backends::common::io::PrepBUFRWRFDAIntegration::
                applyProfiles(token, grid, batches.pilot, y_out);
            logger.Info() << "H(x) pilot size=" << y_out.size();
          } else if (token.rfind("sound:", 0) == 0) {
            metada::backends::common::io::PrepBUFRWRFDAIntegration::
                applyProfiles(token, grid, batches.sound, y_out);
            logger.Info() << "H(x) sound size=" << y_out.size();
          } else {
            const metada::backends::common::io::SurfaceBatch* sb = nullptr;
            if (token.rfind("metar:", 0) == 0)
              sb = &batches.metar;
            else if (token.rfind("synop:", 0) == 0)
              sb = &batches.synop;
            else if (token.rfind("ships:", 0) == 0)
              sb = &batches.ships;
            else if (token.rfind("buoy:", 0) == 0)
              sb = &batches.buoy;
            else if (token.rfind("sonde_sfc:", 0) == 0)
              sb = &batches.sonde_sfc;
            if (sb && !sb->lats.empty()) {
              metada::backends::common::io::PrepBUFRWRFDAIntegration::
                  applySurface(token, grid, *sb, y_out);
              logger.Info()
                  << "H(x) surface " << token << " size=" << y_out.size();
            }
          }
          // Example adjoint: reuse y_out as delta_y
          if (!y_out.empty()) {
            if (token.rfind("airep:", 0) == 0) {
              metada::backends::common::io::PrepBUFRWRFDAIntegration::
                  adjointProfiles(token, grid, batches.airep, y_out, u.data(),
                                  v.data(), t.data(), q.data(), psfc.data());
            } else if (token.rfind("pilot:", 0) == 0) {
              metada::backends::common::io::PrepBUFRWRFDAIntegration::
                  adjointProfiles(token, grid, batches.pilot, y_out, u.data(),
                                  v.data(), t.data(), q.data(), psfc.data());
            } else if (token.rfind("sound:", 0) == 0) {
              metada::backends::common::io::PrepBUFRWRFDAIntegration::
                  adjointProfiles(token, grid, batches.sound, y_out, u.data(),
                                  v.data(), t.data(), q.data(), psfc.data());
            } else {
              const metada::backends::common::io::SurfaceBatch* sb = nullptr;
              if (token.rfind("metar:", 0) == 0)
                sb = &batches.metar;
              else if (token.rfind("synop:", 0) == 0)
                sb = &batches.synop;
              else if (token.rfind("ships:", 0) == 0)
                sb = &batches.ships;
              else if (token.rfind("buoy:", 0) == 0)
                sb = &batches.buoy;
              else if (token.rfind("sonde_sfc:", 0) == 0)
                sb = &batches.sonde_sfc;
              if (sb && !sb->lats.empty()) {
                metada::backends::common::io::PrepBUFRWRFDAIntegration::
                    adjointSurface(token, grid, *sb, y_out, u.data(), v.data(),
                                   t.data(), q.data(), psfc.data());
              }
            }
          }
        }
      }
    } catch (const std::exception& e) {
      logger.Warning() << "PrepBUFR/WRFDA example skipped: " << e.what();
    }

    // Initialize background error covariance
    fwk::BackgroundErrorCovariance<BackendTag> bg_error_cov(
        config.GetSubsection("background_covariance"));

    // Create variational algorithm instance
    fwk::Variational<BackendTag> variational(
        config, background, observations, obs_operators, model, bg_error_cov);

    logger.Info() << "Starting " << var_type << " analysis...";

    // Perform variational analysis
    auto results = variational.analyze();

    // Log results
    logger.Info() << "=== Analysis Results ===";
    logger.Info() << "Variational type: " << results.variational_type;
    logger.Info() << "Final cost: " << results.final_cost;
    logger.Info() << "Cost reduction: " << results.cost_reduction;
    logger.Info() << "Background cost: " << results.background_cost;
    logger.Info() << "Observation cost: " << results.observation_cost;
    logger.Info() << "Iterations: " << results.iterations;
    logger.Info() << "Converged: " << (results.converged ? "Yes" : "No");
    logger.Info() << "Convergence reason: " << results.convergence_reason;

    // Save analysis results
    variational.saveAnalysis(results);

    // Compute and display innovation statistics
    auto innovation_stats =
        variational.computeInnovationStatistics(results.analysis_state);

    logger.Info() << "=== Innovation Statistics ===";
    for (size_t i = 0; i < innovation_stats.size(); ++i) {
      logger.Info() << "Time " << i << ": mean=" << innovation_stats[i].first
                    << ", rms=" << innovation_stats[i].second;
    }

    // Optional: Perform gradient test if requested
    if (config.Get("perform_gradient_test").asBool()) {
      bool gradient_test_passed = variational.performGradientTest(background);
      logger.Info() << "Gradient test "
                    << (gradient_test_passed ? "PASSED" : "FAILED");
    }

    logger.Info()
        << "Variational data assimilation application completed successfully";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "Variational application failed: " << e.what() << std::endl;
    return 1;
  }
}