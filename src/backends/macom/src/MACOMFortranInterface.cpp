/**
 * @file MACOMFortranInterface.cpp
 * @brief Implementation of the C++ interface to MACOM Fortran routines.
 * @ingroup backends_macom
 * @author Metada Framework Team
 */

#include "../include/MACOMFortranInterface.hpp"

#include <iostream>
#include <stdexcept>

namespace metada::backends::macom {

MACOMFortranInterface::MACOMFortranInterface()
    : rank_(0),
      mpi_initialized_by_this_instance_(false),
      model_components_initialized_(false) {
  std::cout << "[MACOMFortranInterface] Constructor called." << std::endl;
}

MACOMFortranInterface::~MACOMFortranInterface() {
  std::cout << "[MACOMFortranInterface] Destructor called." << std::endl;
  if (model_components_initialized_) {
    try {
      finalizeModelComponents();
    } catch (const std::exception& e) {
      std::cerr << "[MACOMFortranInterface] Exception during "
                   "finalizeModelComponents in destructor: "
                << e.what() << std::endl;
    }
  }
  if (mpi_initialized_by_this_instance_) {
    try {
      finalizeMPI();
    } catch (const std::exception& e) {
      std::cerr << "[MACOMFortranInterface] Exception during finalizeMPI in "
                   "destructor: "
                << e.what() << std::endl;
    }
  }
}

void MACOMFortranInterface::initializeMPI() {
  if (mpi_initialized_by_this_instance_) {
    std::cout << "[MACOMFortranInterface] Already initialized by this instance."
              << std::endl;
    return;
  }
  std::cout << "[MACOMFortranInterface] Initializing via Fortran..."
            << std::endl;
  c_macom_initialize_mpi(0);  // Pass 0 as dummy communicator
  c_macom_get_mpi_rank(&rank_);
  mpi_initialized_by_this_instance_ = true;
  std::cout << "[MACOMFortranInterface] Initialized by Fortran. Rank: " << rank_
            << std::endl;
}

int MACOMFortranInterface::getRank() const {
  return rank_;
}

void MACOMFortranInterface::readNamelist() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "[MACOMFortranInterface] Must be initialized before reading "
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
        "[MACOMFortranInterface] Must be initialized before model "
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
    return;
  }
  std::cout << "[MACOMFortranInterface] Calling Fortran to finalize model "
               "components..."
            << std::endl;
  c_macom_finalize_model_components();
  model_components_initialized_ = false;
  std::cout
      << "[MACOMFortranInterface] Fortran finalize model components finished."
      << std::endl;
}

void MACOMFortranInterface::finalizeMPI() {
  if (!mpi_initialized_by_this_instance_) {
    return;
  }
  std::cout << "[MACOMFortranInterface] Calling Fortran to finalize..."
            << std::endl;
  c_macom_finalize_mpi();
  mpi_initialized_by_this_instance_ = false;
  rank_ = 0;
  std::cout << "[MACOMFortranInterface] Fortran finalize finished."
            << std::endl;
}

}  // namespace metada::backends::macom