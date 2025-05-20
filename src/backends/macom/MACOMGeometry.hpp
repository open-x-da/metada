/**
 * @file MACOMGeometry.hpp
 * @brief MACOM geometry backend implementation for hexagonal grid
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <netcdf>
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
 * @brief MACOM geometry backend implementation for hexagonal grid
 *
 * @details
 * This implementation handles MACOM's hexagonal grid structure and reads
 * grid data from NetCDF files. The grid data includes cell centers, edges,
 * areas and masks.
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
   * @return Always false in this implementation
   */
  bool isPeriodicX() const { return false; }

  /**
   * @brief Check if the geometry is periodic in Y dimension
   *
   * @return Always false in this implementation
   */
  bool isPeriodicY() const { return false; }

  /**
   * @brief Check if the geometry is periodic in Z dimension
   *
   * @return Always false in this implementation
   */
  bool isPeriodicZ() const { return false; }

  /**
   * @brief Get the number of grid points (nlpb)
   * @return nlpb_
   */

  template <typename StateBackend>
  void haloExchange([[maybe_unused]] StateBackend& state) {}

 private:
  /**
   * @brief Load grid dimensions from NetCDF file
   *
   * @param ncFile NetCDF file handle
   */
  void loadGridDimensions(netCDF::NcFile& ncFile);

  /**
   * @brief Load grid arrays from NetCDF file
   *
   * @param ncFile NetCDF file handle
   */
  void loadGridArrays(netCDF::NcFile& ncFile);

  /**
   * @brief Initialize grid from NetCDF file
   *
   * @param filename Path to the NetCDF grid file
   */
  void loadGeometryData(const std::string& filename);

  /**
   * @brief Special constructor for cloning
   *
   * @param nlpb Number of grid points
   * @param nk Number of vertical levels
   * @param init Initialization status
   */
  MACOMGeometry(std::size_t nlpb, std::size_t nk, bool init,
                const ConfigBackend* cfg)
      : nlpb_(nlpb), nk_(nk), initialized_(init), config_ptr_(cfg) {}

  // Grid dimensions
  std::size_t nlpb_ = 0;   // Number of grid points
  std::size_t nk_ = 0;     // Number of vertical levels
  std::size_t nkp1_ = 0;   // Number of vertical levels + 1
  std::size_t nl_ = 0;     // Total number of grid points
  std::size_t nlpbz_ = 0;  // Number of grid points for vorticity
  std::size_t nlz_ = 0;    // Total number of grid points for vorticity
  std::size_t nlbdy_ = 0;  // Number of boundary points
  std::size_t ni_ = 0;     // Number of iterations
  bool initialized_ = false;
  const ConfigBackend* config_ptr_ = nullptr;

  // Grid data
  std::vector<double> latC_;   // Latitude at cell centers
  std::vector<double> lonC_;   // Longitude at cell centers
  std::vector<double> dxC_;    // Grid spacing in x direction at cell centers
  std::vector<double> dyC_;    // Grid spacing in y direction at cell centers
  std::vector<double> dxW_;    // Grid spacing in x direction at west edges
  std::vector<double> dyW_;    // Grid spacing in y direction at west edges
  std::vector<double> dxS_;    // Grid spacing in x direction at south edges
  std::vector<double> dyS_;    // Grid spacing in y direction at south edges
  std::vector<double> rAc_;    // Cell areas
  std::vector<double> rAw_;    // Areas at west edges
  std::vector<double> rAs_;    // Areas at south edges
  std::vector<double> maskC_;  // Land/sea mask at cell centers
  std::vector<double> maskW_;  // Land/sea mask at west edges
  std::vector<double> maskS_;  // Land/sea mask at south edges
  std::vector<double> hFacC_;  // Cell thickness factors
  std::vector<double> hFacW_;  // West edge thickness factors
  std::vector<double> hFacS_;  // South edge thickness factors

  // Friend declaration for iterator
  friend class MACOMGeometryIterator<ConfigBackend>;
  friend class MACOMGeometryConstIterator<ConfigBackend>;
};

// ConfigBackend constructor implementation
template <typename ConfigBackend>
MACOMGeometry<ConfigBackend>::MACOMGeometry(const ConfigBackend& config)
    : initialized_(false) {
  std::string input_filename = config.Get("input_file").asString();
  if (input_filename.empty()) {
    throw std::runtime_error(
        "MACOM input grid file path not specified in configuration");
  }

  // Load geometry data from MACOM NetCDF file
  loadGeometryData(input_filename);
  initialized_ = true;
}

// Implementation of loadGridDimensions
template <typename ConfigBackend>
void MACOMGeometry<ConfigBackend>::loadGridDimensions(netCDF::NcFile& ncFile) {
  auto readDimension = [&ncFile](const std::string& name, std::size_t& value) {
    std::cout << "Attempting to read dimension: " << name << std::endl;
    auto dim = ncFile.getDim(name);
    if (dim.isNull()) {
      throw std::runtime_error("Dimension '" + name +
                               "' not found in grid file");
    }
    std::cout << "Found dimension: " << name << std::endl;
    value = dim.getSize();
    std::cout << "Read size for " << name << ": " << value << std::endl;
  };

  // Read all grid dimensions
  readDimension("nlpb", nlpb_);
  readDimension("nl", nl_);
  readDimension("nlpbz", nlpbz_);
  readDimension("nlz", nlz_);
  readDimension("nk", nk_);
  readDimension("nkp1", nkp1_);
  readDimension("ni", ni_);

  std::cout << "Loaded grid dimensions:" << std::endl;
  std::cout << "  nlpb=" << nlpb_ << ", nk=" << nk_ << ", nkp1=" << nkp1_
            << std::endl;
  std::cout << "  nl=" << nl_ << ", nlpbz=" << nlpbz_ << ", nlz=" << nlz_
            << std::endl;
  std::cout << "  ni=" << ni_ << std::endl;
}

// Implementation of loadGridArrays
template <typename ConfigBackend>
void MACOMGeometry<ConfigBackend>::loadGridArrays(netCDF::NcFile& ncFile) {
  // Resize vectors
  latC_.resize(nlpb_);
  lonC_.resize(nlpb_);
  dxC_.resize(nlpb_);
  dyC_.resize(nlpb_);
  dxW_.resize(nlpb_);
  dyW_.resize(nlpb_);
  dxS_.resize(nlpb_);
  dyS_.resize(nlpb_);
  rAc_.resize(nlpb_);
  rAw_.resize(nlpb_);
  rAs_.resize(nlpb_);
  maskC_.resize(nlpb_ * nk_);
  maskW_.resize(nlpb_ * nk_);
  maskS_.resize(nlpb_ * nk_);
  hFacC_.resize(nlpb_ * nk_);
  hFacW_.resize(nlpb_ * nk_);
  hFacS_.resize(nlpb_ * nk_);

  // Read variables
  auto readVar = [&ncFile](const std::string& name, std::vector<double>& data) {
    auto var = ncFile.getVar(name);
    if (var.isNull()) {
      throw std::runtime_error("Variable '" + name +
                               "' not found in grid file");
    }
    var.getVar(data.data());
  };

  // Read grid data
  readVar("lat_c", latC_);
  readVar("lon_c", lonC_);
  readVar("dxC", dxC_);
  readVar("dyC", dyC_);
  readVar("dxW", dxW_);
  readVar("dyW", dyW_);
  readVar("dxS", dxS_);
  readVar("dyS", dyS_);
  readVar("rAc", rAc_);
  readVar("rAw", rAw_);
  readVar("rAs", rAs_);
  readVar("maskC", maskC_);
  readVar("maskW", maskW_);
  readVar("maskS", maskS_);
  readVar("h0FacC", hFacC_);
  readVar("h0FacW", hFacW_);
  readVar("h0FacS", hFacS_);

  std::cout << "Loaded grid arrays for " << nlpb_ << " grid points and " << nk_
            << " vertical levels" << std::endl;
}

// Implementation of loadGeometryData
template <typename ConfigBackend>
void MACOMGeometry<ConfigBackend>::loadGeometryData(
    const std::string& filename) {
  try {
    std::cout << "Attempting to open NetCDF file: " << filename << std::endl;

    // Open the NetCDF file
    netCDF::NcFile ncFile(filename, netCDF::NcFile::read);

    if (ncFile.isNull()) {
      throw std::runtime_error("Failed to open NetCDF file: " + filename);
    }

    std::cout << "Successfully opened NetCDF file" << std::endl;

    // Step 1: Load grid dimensions
    loadGridDimensions(ncFile);

    // Step 2: Load grid arrays
    loadGridArrays(ncFile);

    // Close the file
    ncFile.close();

    initialized_ = true;

    std::cout << "Successfully initialized hexagonal grid from " << filename
              << std::endl;

  } catch (const netCDF::exceptions::NcException& e) {
    throw std::runtime_error("NetCDF error while reading grid file: " +
                             std::string(e.what()));
  } catch (const std::exception& e) {
    throw std::runtime_error("Error initializing grid: " +
                             std::string(e.what()));
  }
}

}  // namespace metada::backends::macom