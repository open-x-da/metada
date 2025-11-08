#pragma once

#include <stdexcept>
#include <string>

#include "ControlVariableBackend.hpp"
#include "Config.hpp"
#include "Logger.hpp"

namespace metada::backends::wrf {

/**
 * @brief Control-variable backend that wraps WRFDA CV5 control space.
 *
 * @details Step 3.1 scaffolding: parses configuration and logs initialization,
 *          while delegating storage/transform logic to existing increment
 *          semantics. Later steps will replace the placeholders with calls to
 *          WRFDA's control-variable transforms.
 *
 * @tparam BackendTag Backend tag satisfying framework concepts (e.g.,
 *         traits::WRFBackendTag).
 */
template <typename BackendTag>
class WRFDACV5ControlVariableBackend
    : public framework::ControlVariableBackend<BackendTag> {
  using Base = framework::ControlVariableBackend<BackendTag>;

 public:
  using typename Base::GeometryBackendType;
  using typename Base::IncrementType;
  using typename Base::StateType;

  explicit WRFDACV5ControlVariableBackend(
      const framework::Config<BackendTag>& config)
      : wrfda_config_(extractControlConfig(config)),
        be_statistics_path_(parseStringOption("be_statistics")),
        cv_options_(parseIntOption("cv_options", 5)) {
    auto& logger = framework::Logger<BackendTag>::Instance();
    logger.Info() << "WRFDA CV5 control backend configured";
    logger.Info() << "  BE statistics: "
                  << (be_statistics_path_.empty() ? "<not provided>"
                                                  : be_statistics_path_);
    logger.Info() << "  cv_options: " << cv_options_;
    logger.Warning()
        << "WRFDA CV5 backend currently uses placeholder transforms (Step 3.1)";
  }

  ~WRFDACV5ControlVariableBackend() override = default;

  IncrementType createIncrement(
      const GeometryBackendType& geometry) const override {
    // Placeholder: reuse existing increment semantics. Later steps will allocate
    // CV5 control vectors and sync with WRFDA structures.
    return IncrementType::createFromGeometry(geometry);
  }

  void addIncrementToState(const IncrementType& increment,
                           StateType& state) const override {
    // Placeholder: direct addition. Future work will call WRFDA transforms
    // (da_transform_vtox / da_transfer_xatowrf).
    state += increment;

    if (!warned_add_increment_stub_) {
      auto& logger = framework::Logger<BackendTag>::Instance();
      logger.Warning()
          << "WRFDA CV5 backend addIncrementToState uses placeholder addition";
      warned_add_increment_stub_ = true;
    }
  }

  std::string name() const override { return "wrfda_cv5"; }

  const std::string& beStatisticsPath() const { return be_statistics_path_; }
  int cvOptions() const { return cv_options_; }

 private:
  static framework::Config<BackendTag> extractControlConfig(
      const framework::Config<BackendTag>& config) {
    if (!config.HasKey("wrfda_control")) {
      throw std::runtime_error(
          "wrfda_cv5 backend requires a 'wrfda_control' configuration "
          "section");
    }
    return config.GetSubsection("wrfda_control");
  }

  std::string parseStringOption(const std::string& key) const {
    if (wrfda_config_.HasKey(key)) {
      return wrfda_config_.Get(key).asString();
    }
    return {};
  }

  int parseIntOption(const std::string& key, int default_value) const {
    if (wrfda_config_.HasKey(key)) {
      return wrfda_config_.Get(key).asInt();
    }
    return default_value;
  }

  framework::Config<BackendTag> wrfda_config_;
  std::string be_statistics_path_;
  int cv_options_;
  mutable bool warned_add_increment_stub_ = false;
};

}  // namespace metada::backends::wrf


