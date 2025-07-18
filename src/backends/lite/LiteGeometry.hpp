#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace metada::backends::lite {

/**
 * @brief Lite geometry backend for concrete testing
 *
 * This class is designed to comply with the GeometryBackendImpl concept.
 * - Deleted default constructor
 * - Deleted copy constructor and copy assignment operator
 * - Config constructor (templated for flexibility)
 * - STL-like container interface
 */
class LiteGeometry {
 public:
  // Type aliases for STL-like container
  using value_type = std::string;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using iterator = std::vector<std::string>::iterator;
  using const_iterator = std::vector<std::string>::const_iterator;

  // Deleted default constructor (concept compliance)
  LiteGeometry() = delete;

  // Deleted copy constructor and copy assignment operator (concept compliance)
  LiteGeometry(const LiteGeometry&) = delete;
  LiteGeometry& operator=(const LiteGeometry&) = delete;

  // Special copy-like constructor for concept compliance (not a general copy
  // constructor)
  LiteGeometry(const LiteGeometry& other, [[maybe_unused]] bool for_clone) {
    variable_names_ = other.variable_names_;
  }

  // Templated config constructor (for concept compliance and flexibility)
  template <typename ConfigBackend>
  explicit LiteGeometry([[maybe_unused]] const ConfigBackend& config) {
    variable_names_ = {"temperature", "pressure", "humidity"};
  }

  // Iterators
  iterator begin() { return variable_names_.begin(); }
  iterator end() { return variable_names_.end(); }
  const_iterator begin() const { return variable_names_.begin(); }
  const_iterator end() const { return variable_names_.end(); }
  const_iterator cbegin() const { return variable_names_.cbegin(); }
  const_iterator cend() const { return variable_names_.cend(); }

  // Size information
  size_type size() const { return variable_names_.size(); }
  bool empty() const { return variable_names_.empty(); }
  size_type max_size() const { return variable_names_.max_size(); }

  // Element access
  reference operator[](size_type idx) { return variable_names_[idx]; }
  const_reference operator[](size_type idx) const {
    return variable_names_[idx];
  }
  reference at(size_type idx) { return variable_names_.at(idx); }
  const_reference at(size_type idx) const { return variable_names_.at(idx); }
  reference front() { return variable_names_.front(); }
  const_reference front() const { return variable_names_.front(); }
  reference back() { return variable_names_.back(); }
  const_reference back() const { return variable_names_.back(); }

  // Concept-compliant clone method
  LiteGeometry clone() const { return LiteGeometry(*this, true); }

  // Variable information (legacy/test helper)
  const std::vector<std::string>& getVariableNames() const {
    return variable_names_;
  }
  void setVariables(const std::vector<std::string>& vars) {
    variable_names_ = vars;
  }

 private:
  std::vector<std::string> variable_names_;
};

}  // namespace metada::backends::lite