#pragma once

#include "Config.hpp"
#include "ControlVariableBackend.hpp"
#include "ControlVariableBackendFactory.hpp"
#include "WRFDACV5ControlVariableBackend.hpp"

namespace metada::framework {

/**
 * @brief WRF-specific extension of ControlVariableBackendFactory
 *
 * @details This factory extends the generic ControlVariableBackendFactory
 * to support WRFDA CV5 control variable backend. It should be used instead
 * of the generic factory when WRF backend support is available.
 *
 * @tparam BackendTag Backend tag (typically WRFBackendTag)
 */
template <typename BackendTag>
class WRFDAControlVariableBackendFactory {
 public:
  using ConfigType = Config<BackendTag>;
  using BaseFactory = ControlVariableBackendFactory<BackendTag>;

  /**
   * @brief Determine backend kind from configuration
   *
   * @param config Configuration object
   * @return Backend kind
   */
  static ControlVariableBackendKind determineBackend(const ConfigType& config) {
    return BaseFactory::determineBackend(config);
  }

  /**
   * @brief Create a control variable backend (supports WRFDA CV5)
   *
   * @param kind Backend kind to create
   * @param config Configuration object
   * @return Unique pointer to the created backend
   */
  static std::unique_ptr<ControlVariableBackend<BackendTag>> createBackend(
      ControlVariableBackendKind kind, const ConfigType& config) {
    switch (kind) {
      case ControlVariableBackendKind::Identity:
        return BaseFactory::createBackend(kind, config);
      case ControlVariableBackendKind::WrfdaCv5:
        return std::make_unique<
            backends::wrf::WRFDACV5ControlVariableBackend<BackendTag>>(config);
      default:
        return BaseFactory::createBackend(kind, config);
    }
  }
};

}  // namespace metada::framework
