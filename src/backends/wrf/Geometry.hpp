/**
 * @file Geometry.hpp
 * @brief WRF geometry backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <memory>
#include <netcdf>
#include <string>
#include <vector>
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>

namespace metada::backends::wrf {

/**
 * @brief WRF geometry backend implementation
 *
 * @details
 * This class implements a geometry backend for the WRF model. It reads
 * geographical data from a WRF NetCDF file and provides methods for
 * traversing the grid and managing boundary conditions.
 */
class Geometry {
 public:
  // Iterator classes (forward declarations)
  class iterator;
  class const_iterator;

  /**
   * @brief Default constructor is deleted
   */
  Geometry() = delete;

  /**
   * @brief Copy constructor is deleted
   */
  Geometry(const Geometry&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  Geometry& operator=(const Geometry&) = delete;

  /**
   * @brief Constructor that takes a configuration backend
   *
   * @param config Configuration containing WRF file path and timestamp
   */
  template <typename ConfigBackend>
  explicit Geometry(const ConfigBackend& config);

  /**
   * @brief Move constructor
   *
   * @param other WRF geometry backend to move from
   */
  Geometry(Geometry&& other) noexcept;

  /**
   * @brief Move assignment operator
   *
   * @param other WRF geometry backend to move from
   * @return Reference to this geometry after assignment
   */
  Geometry& operator=(Geometry&& other) noexcept;

  /**
   * @brief Destructor
   */
  ~Geometry() = default;

  /**
   * @brief Clone this geometry
   *
   * @return A new WRF geometry backend with the same state
   */
  Geometry clone() const;

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
   * @brief Check if the geometry is periodic in X dimension
   *
   * @return True if periodic in X, false otherwise
   */
  bool isPeriodicX() const;

  /**
   * @brief Check if the geometry is periodic in Y dimension
   *
   * @return True if periodic in Y, false otherwise
   */
  bool isPeriodicY() const;

  /**
   * @brief Check if the geometry is periodic in Z dimension
   *
   * @return True if periodic in Z, false otherwise
   */
  bool isPeriodicZ() const;

  /**
   * @brief Perform halo exchange on a state
   *
   * @param state The state on which to perform halo exchange
   */
  template <typename StateBackend>
  void haloExchange(StateBackend& state);

  /**
   * @brief Halo exchange implementation for the backend
   *
   * @param state_ptr Pointer to the state backend
   */
  void haloExchangeImpl(void* state_ptr);

  /**
   * @brief Check if geometry is properly initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const;

  // WRF specific methods

  /**
   * @brief Get the longitude array
   *
   * @return Reference to the longitude array
   */
  const xt::xarray<double>& getLongitude() const;

  /**
   * @brief Get the latitude array
   *
   * @return Reference to the latitude array
   */
  const xt::xarray<double>& getLatitude() const;

  /**
   * @brief Get the elevation array
   *
   * @return Reference to the elevation array
   */
  const xt::xarray<double>& getElevation() const;

 private:
  /**
   * @brief Load geometry data from the WRF NetCDF file
   *
   * @param filename Path to the WRF NetCDF file
   * @param timestamp Timestamp to read from the file
   */
  void loadGeometryData(const std::string& filename,
                        const std::string& timestamp);

  // WRF NetCDF file
  std::string wrfFilename_;
  std::string timestamp_;
  bool initialized_ = false;

  // Geographic data
  xt::xarray<double> longitude_;  // Longitude points
  xt::xarray<double> latitude_;   // Latitude points
  xt::xarray<double> elevation_;  // Elevation/terrain height

  // Grid dimensions
  std::size_t nx_ = 0;  // X dimension (west-east)
  std::size_t ny_ = 0;  // Y dimension (south-north)
  std::size_t nz_ = 0;  // Z dimension (vertical levels)

  // Periodicity flags
  bool periodicX_ = false;
  bool periodicY_ = false;
  bool periodicZ_ = false;

  // Friend declaration for iterator
  friend class Geometry::iterator;
};

// Template method implementation
template <typename StateBackend>
void Geometry::haloExchange(StateBackend& state) {
  haloExchangeImpl(static_cast<void*>(&state));
}

// Iterator method implementations
inline Geometry::iterator Geometry::begin() {
  return iterator(this, 0);
}

inline Geometry::iterator Geometry::end() {
  return iterator(this, totalGridSize());
}

inline Geometry::const_iterator Geometry::begin() const {
  return const_iterator(this, 0);
}

inline Geometry::const_iterator Geometry::end() const {
  return const_iterator(this, totalGridSize());
}

}  // namespace metada::backends::wrf

// Implementation of the WRF geometry iterator class
#include "WRFGeometryIterator.hpp"