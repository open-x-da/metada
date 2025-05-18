/**
 * @file MACOMFortranInterface.hpp
 * @brief C++ interface to MACOM Fortran routines.
 * @ingroup backends_macom
 * @author Metada Framework Team
 */

#pragma once

#include <mpi.h>

#include <stdexcept>  // For std::runtime_error
#include <string>

// Forward declare the C functions from the Fortran wrapper
extern "C" {
void c_macom_initialize_mpi(int comm_c);
void c_macom_get_mpi_rank(int* rank);
void c_macom_read_namelist();
void c_macom_initialize_model_components();
void c_macom_run_model_step(int current_iter_c, int* model_status_c);
void c_macom_finalize_model_components();
void c_macom_finalize_mpi();
}

namespace metada::backends::macom {

class MACOMFortranInterface {
 public:
  /**
   * @brief Constructor.
   * @param comm The MPI communicator to be used (typically MPI_COMM_WORLD).
   */
  explicit MACOMFortranInterface(MPI_Comm comm = MPI_COMM_WORLD);

  /**
   * @brief Destructor.
   * Ensures MPI is finalized if initialized by this interface.
   */
  ~MACOMFortranInterface();

  // Delete copy and move constructors/assignments for simplicity
  MACOMFortranInterface(const MACOMFortranInterface&) = delete;
  MACOMFortranInterface& operator=(const MACOMFortranInterface&) = delete;
  MACOMFortranInterface(MACOMFortranInterface&&) = delete;
  MACOMFortranInterface& operator=(MACOMFortranInterface&&) = delete;

  /**
   * @brief Initializes the MPI environment through Fortran.
   */
  void initializeMPI();

  /**
   * @brief Gets the MPI rank of the current process from Fortran.
   * @return The MPI rank.
   */
  int getRank() const;

  /**
   * @brief Calls the Fortran routine to read the namelist.
   */
  void readNamelist();

  /**
   * @brief Calls the Fortran routine to initialize core model components.
   */
  void initializeModelComponents();

  /**
   * @brief Calls the Fortran routine to run a model step (or full loop).
   * @param current_iteration The current iteration number to pass to Fortran.
   * @throws std::runtime_error if the Fortran model step indicates an error.
   */
  void runModelStep(int current_iteration);

  /**
   * @brief Calls the Fortran routine to finalize model components.
   */
  void finalizeModelComponents();

  /**
   * @brief Finalizes the MPI environment through Fortran.
   */
  void finalizeMPI();

  /**
   * @brief Checks if the MPI environment was initialized by this interface.
   */
  bool isMPIInitialized() const { return mpi_initialized_by_this_instance_; }

 private:
  MPI_Comm mpi_comm_;
  int fortran_mpi_comm_;  // MPI_Comm converted for Fortran (usually an int)
  int rank_;
  bool mpi_initialized_by_this_instance_;
  bool model_components_initialized_;
};  // class MACOMFortranInterface

}  // namespace metada::backends::macom