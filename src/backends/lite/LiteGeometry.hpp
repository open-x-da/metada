#pragma once

#include <string>
#include <vector>

namespace metada::backends::lite {

/**
 * @brief Lite geometry backend for concrete testing
 */
class LiteGeometry {
 public:
  LiteGeometry() = default;

  // Constructor with config (for compatibility)
  template <typename ConfigBackend>
  LiteGeometry(const ConfigBackend& config) {
    // Initialize with default variables
    variable_names_ = {"temperature", "pressure", "humidity"};
  }

  // Variable information
  const std::vector<std::string>& getVariableNames() const {
    return variable_names_;
  }

  // Test helper methods
  void setVariables(const std::vector<std::string>& vars) {
    variable_names_ = vars;
  }

 private:
  std::vector<std::string> variable_names_;
};

}  // namespace metada::backends::lite