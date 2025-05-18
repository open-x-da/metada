/**
 * @file MACOMFortranInterface.cpp
 * @brief Implementation of the C++ interface to MACOM Fortran routines.
 * @ingroup backends_macom
 * @author Metada Framework Team
 */

#include "../include/MACOMFortranInterface.hpp"  // Adjust path as necessary

#include <iostream>  // For std::cout, std::cerr

namespace metada::backends::macom {

MACOMFortranInterface::MACOMFortranInterface(MPI_Comm comm)
    : mpi_comm_(comm),
      rank_(-1),
      mpi_initialized_by_this_instance_(false),
      model_components_initialized_(false) {
  // Convert the C MPI_Comm to a Fortran-compatible integer handle.
  // MPI_Comm_c2f is the standard way if mpi.h is correctly included and MPI
  // supports it.
  fortran_mpi_comm_ = MPI_Comm_c2f(mpi_comm_);
  std::cout << "[MACOMFortranInterface] Constructor: C MPI_Comm " << mpi_comm_
            << " converted to Fortran handle " << fortran_mpi_comm_
            << std::endl;
}

MACOMFortranInterface::~MACOMFortranInterface() {
  std::cout << "[MACOMFortranInterface] Destructor called." << std::endl;
  // Try to finalize model components if they were initialized
  if (model_components_initialized_) {
    try {
      finalizeModelComponents();
    } catch (const std::exception& e) {
      std::cerr << "[MACOMFortranInterface] Exception during "
                   "finalizeModelComponents in destructor: "
                << e.what() << std::endl;
    } catch (...) {
      std::cerr << "[MACOMFortranInterface] Unknown exception during "
                   "finalizeModelComponents in destructor."
                << std::endl;
    }
  }
  // Finalize MPI only if this instance initialized it.
  if (mpi_initialized_by_this_instance_) {
    try {
      finalizeMPI();
    } catch (const std::exception& e) {
      std::cerr << "[MACOMFortranInterface] Exception during finalizeMPI in "
                   "destructor: "
                << e.what() << std::endl;
    } catch (...) {
      std::cerr << "[MACOMFortranInterface] Unknown exception during "
                   "finalizeMPI in destructor."
                << std::endl;
    }
  }
}

void MACOMFortranInterface::initializeMPI() {
  if (mpi_initialized_by_this_instance_) {
    std::cout
        << "[MACOMFortranInterface] MPI already initialized by this instance."
        << std::endl;
    return;
  }
  std::cout << "[MACOMFortranInterface] Initializing MPI via Fortran..."
            << std::endl;
  c_macom_initialize_mpi(fortran_mpi_comm_);
  c_macom_get_mpi_rank(&rank_);  // Get rank after Fortran MPI init
  mpi_initialized_by_this_instance_ = true;
  std::cout << "[MACOMFortranInterface] MPI Initialized by Fortran. C++ knows "
               "rank as: "
            << rank_ << std::endl;
}

int MACOMFortranInterface::getRank() const {
  if (!mpi_initialized_by_this_instance_) {
    // Optionally, if MPI might be initialized externally and you still want
    // rank: int current_rank = -1; MPI_Comm_rank(mpi_comm_, &current_rank);
    // return current_rank;
    std::cerr << "[MACOMFortranInterface] Warning: getRank called but MPI not "
                 "initialized by this instance. Returning stored rank: "
              << rank_ << std::endl;
  }
  return rank_;
}

void MACOMFortranInterface::readNamelist() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "[MACOMFortranInterface] MPI must be initialized before reading "
        "namelist.");
  }
  std::cout << "[MACOMFortranInterface] Calling Fortran to read namelist..."
            << std::endl;
  c_macom_read_namelist();
  std::cout << "[MACOMFortranInterface] Fortran read namelist finished."
            << std::endl;
}

void MACOMFortranInterface::initializeModelComponents() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "[MACOMFortranInterface] MPI must be initialized before model "
        "components.");
  }
  std::cout << "[MACOMFortranInterface] Calling Fortran to initialize model "
               "components..."
            << std::endl;
  c_macom_initialize_model_components();
  model_components_initialized_ = true;
  std::cout
      << "[MACOMFortranInterface] Fortran initialize model components finished."
      << std::endl;
}

void MACOMFortranInterface::runModelStep(int current_iteration) {
  if (!model_components_initialized_) {
    throw std::runtime_error(
        "[MACOMFortranInterface] Model components must be initialized before "
        "running a step.");
  }
  std::cout << "[MACOMFortranInterface] Calling Fortran to run model step for "
               "iteration: "
            << current_iteration << std::endl;
  int status = 0;
  c_macom_run_model_step(current_iteration, &status);
  if (status != 0) {
    throw std::runtime_error(
        "[MACOMFortranInterface] Fortran model step reported an error. "
        "Status: " +
        std::to_string(status));
  }
  std::cout << "[MACOMFortranInterface] Fortran run model step finished for "
               "iteration: "
            << current_iteration << std::endl;
}

void MACOMFortranInterface::finalizeModelComponents() {
  if (!model_components_initialized_) {
    // std::cout << "[MACOMFortranInterface] Model components not initialized or
    // already finalized. Skipping Fortran call." << std::endl;
    return;
  }
  std::cout << "[MACOMFortranInterface] Calling Fortran to finalize model "
               "components..."
            << std::endl;
  c_macom_finalize_model_components();
  model_components_initialized_ = false;  // Mark as finalized
  std::cout
      << "[MACOMFortranInterface] Fortran finalize model components finished."
      << std::endl;
}

void MACOMFortranInterface::finalizeMPI() {
  if (!mpi_initialized_by_this_instance_) {
    // std::cout << "[MACOMFortranInterface] MPI not initialized by this
    // instance or already finalized. Skipping Fortran call." << std::endl;
    return;
  }
  std::cout << "[MACOMFortranInterface] Calling Fortran to finalize MPI..."
            << std::endl;
  c_macom_finalize_mpi();
  mpi_initialized_by_this_instance_ = false;  // Mark as finalized
  rank_ = -1;
  std::cout << "[MACOMFortranInterface] Fortran finalize MPI finished."
            << std::endl;
}

}  // namespace metada::backends::macom