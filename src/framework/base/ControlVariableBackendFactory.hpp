#pragma once

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>

#include "Config.hpp"
#include "ControlVariableBackend.hpp"
#include "IdentityControlVariableBackend.hpp"
#include "backends/wrf/WRFDACV5ControlVariableBackend.hpp"

namespace metada::framework {

/**
 * @brief Supported control-variable backend implementations.
 */
enum class ControlVariableBackendKind {
  Identity,  ///< Direct grid%xa control vector (current default)
  WrfdaCv5   ///< WRFDA CV5 control vector (placeholder for future integration)
};

/**
 * @brief Convert backend kind to human-readable name.
 */
inline std::string controlBackendKindToString(ControlVariableBackendKind kind) {
  switch (kind) {
    case ControlVariableBackendKind::Identity:
      return "identity";
    case ControlVariableBackendKind::WrfdaCv5:
      return "wrfda_cv5";
    default:
      return "unknown";
  }
}

/**
 * @brief Factory utility for selecting and constructing control-variable
 *        backends.
 *
 * @details At this stage the factory routes all requests to the existing
 *          grid%xa-based backend while providing configuration plumbing for
 *          future backends (e.g., WRFDA CV5).
 */
template <typename BackendTag>
class ControlVariableBackendFactory {
 public:
  using ConfigType = Config<BackendTag>;

  static ControlVariableBackendKind determineBackend(const ConfigType& config) {
    if (!config.HasKey("control_variable_backend")) {
      return ControlVariableBackendKind::Identity;
    }

    auto backend_value = config.Get("control_variable_backend").asString();
    std::string normalized;
    normalized.resize(backend_value.size());
    std::transform(
        backend_value.begin(), backend_value.end(), normalized.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (normalized == "identity" || normalized == "grid_xa" ||
        normalized == "gridxa") {
      return ControlVariableBackendKind::Identity;
    }

    if (normalized == "wrfda_cv5" || normalized == "wrfda-cv5") {
      return ControlVariableBackendKind::WrfdaCv5;
    }

    throw std::runtime_error("Unsupported control_variable_backend: " +
                             backend_value);
  }

  static std::unique_ptr<ControlVariableBackend<BackendTag>> createBackend(
      ControlVariableBackendKind kind, const ConfigType& config) {
    switch (kind) {
      case ControlVariableBackendKind::Identity:
        return std::make_unique<IdentityControlVariableBackend<BackendTag>>();
      case ControlVariableBackendKind::WrfdaCv5:
        return std::make_unique<
            backends::wrf::WRFDACV5ControlVariableBackend<BackendTag>>(config);
      default:
        throw std::logic_error("Unknown control variable backend kind");
    }
  }
};

}  // namespace metada::framework
