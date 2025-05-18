/**
 * @file MACOMGeometry.hpp
 * @brief MACOM geometry backend implementation (simplified)
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Forward declarations
namespace metada::backends::macom {
template <typename ConfigBackend>
class MACOMGeometryIterator;
template <typename ConfigBackend>
class MACOMGeometryConstIterator;
}  // namespace metada::backends::macom

namespace metada::backends::macom {

/**
 * @brief Simplified MACOM geometry backend implementation
 *
 * @details
 * This is a minimal implementation that provides the necessary interface
 * for the framework to operate, but with simplified functionality.
 */
template <typename ConfigBackend>
class MACOMGeometry {
 public:
  // Iterator type aliases
  using iterator = MACOMGeometryIterator<ConfigBackend>;
  using const_iterator = MACOMGeometryConstIterator<ConfigBackend>;

  // --- Deleted constructors and assignment operators ---
  MACOMGeometry() = delete;
  MACOMGeometry(const MACOMGeometry&) = delete;
  MACOMGeometry& operator=(const MACOMGeometry&) = delete;

  /**
   * @brief Constructor that takes a configuration backend
   *
   * @param config Configuration containing geometry settings
   */
  explicit MACOMGeometry(const ConfigBackend& config);

  /**
   * @brief Move constructor
   *
   * @param other MACOM geometry backend to move from
   */
  MACOMGeometry(MACOMGeometry&& other) noexcept;

  /**
   * @brief Move assignment operator
   *
   * @param other MACOM geometry backend to move from
   * @return Reference to this geometry after assignment
   */
  MACOMGeometry& operator=(MACOMGeometry&& other) noexcept;

  /**
   * @brief Destructor
   */
  ~MACOMGeometry() = default;

  /**
   * @brief Clone this geometry
   *
   * @return A new MACOMGeometry instance (by value)
   */
  MACOMGeometry clone() const;

  /**
   * @brief Get iterator to the beginning of the grid
   *
   * @return Iterator pointing to the first grid point
   */
  iterator begin();

  /**
   * @brief Get iterator to the end of the grid
   *
   * @return Iterator pointing past the last grid point
   */
  iterator end();

  /**
   * @brief Get const iterator to the beginning of the grid
   *
   * @return Const iterator pointing to the first grid point
   */
  const_iterator begin() const;

  /**
   * @brief Get const iterator to the end of the grid
   *
   * @return Const iterator pointing past the last grid point
   */
  const_iterator end() const;

  /**
   * @brief Get the total number of grid points
   *
   * @return Total number of grid points in the geometry
   */
  std::size_t totalGridSize() const;

  /**
   * @brief Check if geometry is properly initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const;

  /**
   * @brief Check if the geometry is periodic in X dimension
   *
   * @return Always false in this simplified implementation
   */
  bool isPeriodicX() const { return false; }

  /**
   * @brief Check if the geometry is periodic in Y dimension
   *
   * @return Always false in this simplified implementation
   */
  bool isPeriodicY() const { return false; }

  /**
   * @brief Check if the geometry is periodic in Z dimension
   *
   * @return Always false in this simplified implementation
   */
  bool isPeriodicZ() const { return false; }

  /**
   * @brief Dummy halo exchange implementation
   *
   * @param state The state (ignored in this simplified implementation)
   */
  template <typename StateBackend>
  void haloExchange([[maybe_unused]] StateBackend& state) {}

 private:
  /**
   * @brief Special constructor for cloning
   *
   * @param nx Grid dimension in X
   * @param ny Grid dimension in Y
   * @param nz Grid dimension in Z
   * @param init Initialization status
   */
  MACOMGeometry(std::size_t nx, std::size_t ny, std::size_t nz, bool init,
                const ConfigBackend* cfg);

  // Basic grid dimensions
  std::size_t nx_ = 10;  // Default X dimension
  std::size_t ny_ = 10;  // Default Y dimension
  std::size_t nz_ = 1;   // Default Z dimension
  bool initialized_ = false;
  const ConfigBackend* config_ptr_ =
      nullptr;  // Store a pointer to config for cloning

  // Friend declaration for iterator
  friend class MACOMGeometryIterator<ConfigBackend>;
  friend class MACOMGeometryConstIterator<ConfigBackend>;
};

// ConfigBackend constructor implementation
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend>::MACOMGeometry(const ConfigBackend& config)
    : initialized_(true), config_ptr_(&config) {
  // Simplified: Use defaults or try to read from config if keys exist
  try {
    if (config.HasKey("nx"))
      nx_ = static_cast<std::size_t>(config.Get("nx").asInt());
    if (config.HasKey("ny"))
      ny_ = static_cast<std::size_t>(config.Get("ny").asInt());
    if (config.HasKey("nz"))
      nz_ = static_cast<std::size_t>(config.Get("nz").asInt());
  } catch (const std::exception& e) {
    // Optional: Log warning or error if config keys are missing or invalid
    // For simplicity, we fall back to defaults if config reading fails.
  }
}

// Move constructor implementation
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend>::MACOMGeometry(MACOMGeometry&& other) noexcept
    : nx_(other.nx_),
      ny_(other.ny_),
      nz_(other.nz_),
      initialized_(other.initialized_),
      config_ptr_(other.config_ptr_) {
  other.initialized_ = false;  // Reset moved-from object
  other.nx_ = 0;
  other.ny_ = 0;
  other.nz_ = 0;
  other.config_ptr_ = nullptr;
}

// Move assignment operator implementation
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend>& MACOMGeometry<ConfigBackend>::operator=(
    MACOMGeometry&& other) noexcept {
  if (this != &other) {
    nx_ = other.nx_;
    ny_ = other.ny_;
    nz_ = other.nz_;
    initialized_ = other.initialized_;
    config_ptr_ = other.config_ptr_;

    other.initialized_ = false;  // Reset moved-from object
    other.nx_ = 0;
    other.ny_ = 0;
    other.nz_ = 0;
    other.config_ptr_ = nullptr;
  }
  return *this;
}

// Private constructor for cloning
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend>::MACOMGeometry(std::size_t nx, std::size_t ny,
                                            std::size_t nz, bool init,
                                            const ConfigBackend* cfg)
    : nx_(nx), ny_(ny), nz_(nz), initialized_(init), config_ptr_(cfg) {}

// Implementation of totalGridSize
template <typename ConfigBackend>
std::size_t MACOMGeometry<ConfigBackend>::totalGridSize() const {
  return nx_ * ny_ * nz_;
}

// isInitialized implementation
template <typename ConfigBackend>
bool MACOMGeometry<ConfigBackend>::isInitialized() const {
  return initialized_;
}

// 在类外添加这些实现：
template <typename ConfigBackend>
typename MACOMGeometry<ConfigBackend>::iterator
MACOMGeometry<ConfigBackend>::begin() {
  return iterator(this, 0);
}

template <typename ConfigBackend>
typename MACOMGeometry<ConfigBackend>::iterator
MACOMGeometry<ConfigBackend>::end() {
  return iterator(this, totalGridSize());
}

template <typename ConfigBackend>
typename MACOMGeometry<ConfigBackend>::const_iterator
MACOMGeometry<ConfigBackend>::begin() const {
  return const_iterator(this, 0);
}

template <typename ConfigBackend>
typename MACOMGeometry<ConfigBackend>::const_iterator
MACOMGeometry<ConfigBackend>::end() const {
  return const_iterator(this, totalGridSize());
}

// Implementation of clone
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend> MACOMGeometry<ConfigBackend>::clone() const {
  if (!config_ptr_) {
    throw std::runtime_error(
        "Cannot clone MACOMGeometry without a valid config pointer.");
  }
  // Use the private constructor to create a copy, passing the config pointer
  return MACOMGeometry<ConfigBackend>(nx_, ny_, nz_, initialized_, config_ptr_);
}

}  // namespace metada::backends::macom