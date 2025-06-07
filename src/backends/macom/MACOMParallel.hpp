/**
 * @file MACOMParallel.hpp
 * @brief MACOM parallel operations implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <memory>
#include <string>

#include "include/MACOMFortranInterface.hpp"
#include "include/MACOMlogging.hpp"

namespace metada::backends::macom {

/**
 * @brief MACOM parallel operations implementation
 *
 * @details
 * This class handles all parallel operations for the MACOM model,
 * including MPI initialization, finalization, and parallel execution.
 * Supports both C++ and Fortran-based MPI initialization.
 */
class MACOMParallel {
 public:
  /**
   * @brief Get the singleton instance
   *
   * @return Reference to the MACOMParallel instance
   */
  static MACOMParallel& getInstance() {
    static MACOMParallel instance;
    return instance;
  }

  /**
   * @brief Set MPI initialization mode
   *
   * @param use_fortran true to use Fortran MPI initialization, false for C++
   */
  void setFortranMode(bool use_fortran) {
    use_fortran_mpi_ = use_fortran;
    MACOM_LOG_INFO("MACOMParallel",
                   "MPI initialization mode set to: " +
                       std::string(use_fortran ? "Fortran" : "C++"));
  }

  /**
   * @brief Initialize parallel environment
   *
   * @param argc Command line argument count
   * @param argv Command line arguments
   * @return true if initialization successful, false otherwise
   */
  bool initialize(int argc, char** argv) {
#ifdef USE_MPI
    if (use_fortran_mpi_) {
      return initializeFortranMPI();
    } else {
      return initializeCppMPI(argc, argv);
    }
#else
    rank_ = 0;
    size_ = 1;
    initialized_ = true;
    return true;
#endif
  }

  /**
   * @brief Finalize parallel environment
   */
  void finalize() {
#ifdef USE_MPI
    if (initialized_) {
      int finalized;
      MPI_Finalized(&finalized);
      if (!finalized) {
        if (use_fortran_mpi_) {
          // Fortran finalization is handled by the Fortran code
          MACOM_LOG_INFO("MACOMParallel",
                         "Skipping MPI finalization in Fortran mode");
        } else {
          MPI_Finalize();
        }
      }
      initialized_ = false;
    }
#endif
  }

  /**
   * @brief Get current process rank
   *
   * @return Process rank
   */
  int getRank() const { return rank_; }

  /**
   * @brief Get total number of processes
   *
   * @return Number of processes
   */
  int getSize() const { return size_; }

  /**
   * @brief Check if running in parallel mode
   *
   * @return true if running with multiple processes
   */
  bool isParallel() const { return size_ > 1; }

  /**
   * @brief Check if this is the root process
   *
   * @return true if this is process 0
   */
  bool isRoot() const { return rank_ == 0; }

  /**
   * @brief Check if using Fortran MPI initialization
   *
   * @return true if using Fortran MPI initialization
   */
  bool isFortranMode() const { return use_fortran_mpi_; }

 private:
  // Fortran interface
  std::unique_ptr<MACOMFortranInterface> fortranInterface_;
  int mpi_rank_;

  // Private constructor for singleton
  MACOMParallel()
      : rank_(0), size_(1), initialized_(false), use_fortran_mpi_(false) {}

  // Delete copy constructor and assignment operator
  MACOMParallel(const MACOMParallel&) = delete;
  MACOMParallel& operator=(const MACOMParallel&) = delete;

  /**
   * @brief Initialize MPI using C++ interface
   */
  bool initializeCppMPI(int argc, char** argv) {
#ifdef USE_MPI
    // Check if we're running under mpiexec
    int initialized;
    MPI_Initialized(&initialized);

    if (!initialized) {
      // Not running under mpiexec, initialize MPI as a single process
      int provided;
      int result =
          MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
      if (result != MPI_SUCCESS) {
        MACOM_LOG_ERROR("MACOMParallel", "MPI_Init_thread failed");
        return false;
      }
    }

    // Get process info
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    if (rank_ == 0) {
      MACOM_LOG_INFO("MACOMParallel", "Running with " + std::to_string(size_) +
                                          " MPI processes (C++ mode)");
    }
    initialized_ = true;
    return true;
#else
    return false;
#endif
  }

  /**
   * @brief Initialize MPI using Fortran interface
   */
  bool initializeFortranMPI() {
#ifdef USE_MPI
    // Initialize MPI through Fortran interface
    if (!fortranInterface_) {
      fortranInterface_ = std::make_unique<MACOMFortranInterface>();
    }
    fortranInterface_->initializeMPI();

    mpi_rank_ = fortranInterface_->getRank();
    MACOM_LOG_INFO("MACOMModel", "MPI Initialized via Fortran. Model Rank: " +
                                     std::to_string(mpi_rank_));

    // Get process info
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    if (rank_ == 0) {
      MACOM_LOG_INFO("MACOMParallel", "Running with " + std::to_string(size_) +
                                          " MPI processes (Fortran mode)");
    }
    initialized_ = true;
    return true;
#else
    return false;
#endif
  }

  int rank_;              // Process rank
  int size_;              // Total number of processes
  bool initialized_;      // Initialization status
  bool use_fortran_mpi_;  // Whether to use Fortran MPI initialization
};

}  // namespace metada::backends::macom