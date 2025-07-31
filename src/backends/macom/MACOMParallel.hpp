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
    // MACOM_LOG_INFO("MACOMParallel",
    //                "MPI initialization mode set to: " +
    //                    std::string(use_fortran ? "Fortran" : "C++"));
  }

  /**
   * @brief Initialize the parallel environment
   * @param argc Command line argument count
   * @param argv Command line arguments
   * @return true if initialization successful, false otherwise
   */
  bool initialize([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {
#ifdef USE_MPI
    if (use_fortran_mpi_) {
      return initializeFortranMPI();
    } else {
      return initializeCppMPI(argc, argv);
    }
#else
    rank_ = 0;
    size_ = 1;
    io_procs_ = 0;
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
      if (use_fortran_mpi_) {
        finalizeFortranMPI();
      } else {
        finalizeCppMPI();
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
   * @brief Get number of I/O processes
   *
   * @return Number of I/O processes
   */
  int getIOProcs() const { return io_procs_; }

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

  /**
   * @brief Synchronize all processes using MPI_Barrier
   *
   * @details
   * This function ensures all processes reach this point before proceeding.
   * Useful for synchronizing processes before operations that should be
   * coordinated across all ranks.
   */
  void barrier() {
#ifdef USE_MPI
    if (initialized_) {
      // if (use_fortran_mpi_) {
      //   // Use Fortran interface for barrier
      //   if (fortranInterface_) {
      //     fortranInterface_->barrier();
      //   }
      // } else {
      // Use C++ MPI directly
      MPI_Barrier(MPI_COMM_WORLD);
      // }
    }
#endif
  }

  /**
   * @brief Configuration structure for MACOM model
   */
  struct ModelConfig {
    bool mitice_on{false};   // Whether mitice is enabled
    bool restart_in{false};  // Whether in restart mode
    bool assim_in{false};    // Whether in assimilation mode
    int nIter0{0};           // Initial iteration number
    int nIterMax{0};         // Maximum iteration number
  };

  /**
   * @brief Get number of compute processes
   *
   * @return Number of compute processes
   */
  int getCompProcs() const { return comp_procs_; }

  /**
   * @brief Initialize mitice components
   */
  void initializeMitice() {
    if (fortranInterface_) {
      fortranInterface_->initializeMitice();
    }
  }

  /**
   * @brief Get the model configuration
   *
   * @return const Config& Reference to the model configuration
   */
  const ModelConfig& getConfig() const { return modelconfig_; }

  /**
   * @brief Get the Fortran interface
   *
   * @return std::unique_ptr<MACOMFortranInterface>& Reference to the Fortran
   * interface
   */
  std::unique_ptr<MACOMFortranInterface>& getFortranInterface() {
    return fortranInterface_;
  }

 private:
  // Fortran interface
  std::unique_ptr<MACOMFortranInterface> fortranInterface_;

  // Private constructor for singleton
  MACOMParallel()
      : rank_(0),
        size_(1),
        io_procs_(0),
        comp_procs_(0),
        initialized_(false),
        use_fortran_mpi_(false) {}

  // Delete copy constructor and assignment operator
  MACOMParallel(const MACOMParallel&) = delete;
  MACOMParallel& operator=(const MACOMParallel&) = delete;

  /**
   * @brief Initialize MPI using C++ interface
   */
  bool initializeCppMPI([[maybe_unused]] int argc,
                        [[maybe_unused]] char** argv) {
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
    // Assuming I/O processes are the last quarter of total processes
    io_procs_ = size_ / 4;

    if (rank_ == 0) {
      MACOM_LOG_INFO("MACOMParallel", "Running with " + std::to_string(size_) +
                                          " MPI processes (C++ mode), " +
                                          std::to_string(io_procs_) +
                                          " I/O processes");
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

    io_procs_ = 2;

    fortranInterface_->initializeMPI(io_procs_);

    // mpi_rank_ = fortranInterface_->getRank(); // This line is removed
    MACOM_LOG_INFO("MACOMParallel",
                   "MPI Initialized via Fortran. Model Rank: " +
                       std::to_string(rank_));  // This line is updated

    // rank_ = mpi_rank_; // This line is removed
    // mpi_size_ = fortranInterface_->getSize(); // This line is removed
    size_ = fortranInterface_->getSize();  // This line is updated

    comp_procs_ = size_ - io_procs_;

    // // Get process info
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    // MPI_Comm_size(MPI_COMM_WORLD, &size_);
    // std::cout << "C++ rank: " << rank_ << std::endl;

    if (rank_ < comp_procs_) {
      // Read namelist and get config flags
      fortranInterface_->readNamelist();

      bool mitice, restart, assim;
      int init_iter, max_iter;
      fortranInterface_->getConfigFlags(mitice, restart, assim, init_iter,
                                        max_iter);

      // Store the config flags
      modelconfig_.mitice_on = mitice;
      modelconfig_.restart_in = restart;
      modelconfig_.assim_in = assim;
      modelconfig_.nIter0 = init_iter;
      modelconfig_.nIterMax = max_iter;
    }

    // if (rank_ == 0) {
    //   MACOM_LOG_INFO("MACOMParallel", "Running with " + std::to_string(size_)
    //   +
    //                                       " MPI processes (Fortran mode), " +
    //                                       std::to_string(io_procs_) +
    //                                       " I/O processes");
    // }
    initialized_ = true;
    return true;
#else
    return false;
#endif
  }

  /**
   * @brief Finalize MPI using C++ interface
   */
  void finalizeCppMPI() {
#ifdef USE_MPI
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized) {
      MACOM_LOG_INFO("MACOMParallel", "Finalizing MPI in C++ mode...");
      MPI_Finalize();
      MACOM_LOG_INFO("MACOMParallel", "MPI finalized in C++ mode");
    } else {
      MACOM_LOG_INFO("MACOMParallel", "MPI already finalized in C++ mode");
    }
#endif
  }

  /**
   * @brief Finalize MPI using Fortran interface
   */
  void finalizeFortranMPI() {
#ifdef USE_MPI
    if (fortranInterface_) {
      // mpi_rank_ = fortranInterface_->getRank(); // This line is removed

      if (rank_ == 0) {  // This line is updated
        MACOM_LOG_INFO("MACOMParallel", "Finalizing MPI in Fortran mode...");
      }
      fortranInterface_->finalizeMPI();

      // int finalized;
      // MPI_Finalized(&finalized);
      // if (!finalized) {
      //   MACOM_LOG_INFO("MACOMParallel", "MPI finalized in Fortran mode");
      // } else {
      //   MACOM_LOG_INFO("MACOMParallel",
      //                  "MPI already finalized in Fortran mode");
      // }
    }
#endif
  }

  int rank_;                 // Process rank
  int size_;                 // Total number of processes
  int io_procs_;             // Number of I/O processes
  int comp_procs_;           // Number of compute processes
  ModelConfig modelconfig_;  // Model configuration
  bool initialized_;         // Initialization status
  bool use_fortran_mpi_;     // Whether to use Fortran MPI initialization
};

}  // namespace metada::backends::macom