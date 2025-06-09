/**
 * @file MACOMFortranInterface.hpp
 * @brief C++ interface to MACOM Fortran routines.
 * @ingroup backends_macom
 * @author Metada Framework Team
 */

#pragma once

#include <mpi.h>

// Forward declare the C functions from the Fortran wrapper
extern "C" {
void c_macom_initialize_mpi(int comm_c);
void c_macom_get_mpi_rank(int* rank);
void c_macom_get_mpi_size(int* size);
void c_macom_read_namelist();
void c_macom_initialize_mitice();
void c_get_macom_config_flags(bool* mitice, bool* restart, bool* assim,
                              int* init_iter, int* max_iter);
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
   */
  explicit MACOMFortranInterface();

  /**
   * @brief Destructor.
   */
  ~MACOMFortranInterface();

  // Delete copy and move constructors/assignments for simplicity
  MACOMFortranInterface(const MACOMFortranInterface&) = delete;
  MACOMFortranInterface& operator=(const MACOMFortranInterface&) = delete;
  MACOMFortranInterface(MACOMFortranInterface&&) = delete;
  MACOMFortranInterface& operator=(MACOMFortranInterface&&) = delete;

  /**
   * @brief Initializes the environment through Fortran.
   * @param io_procs Number of I/O processes
   */
  void initializeMPI(int io_procs);

  /**
   * @brief Gets the rank of the current process from Fortran.
   * @return The rank.
   */
  int getRank() const;

  /**
   * @brief Gets the size of the current process from Fortran.
   * @return The size.
   */
  int getSize() const;

  /**
   * @brief Gets the number of compute processes.
   * @return The number of compute processes.
   */
  int getCompProcs() const;

  /**
   * @brief Gets the Fortran communicator.
   * @return The Fortran communicator.
   */
  int getFortranComm() const;

  /**
   * @brief Calls the Fortran routine to read the namelist.
   */
  void readNamelist();

  /**
   * @brief Gets the configuration flags from Fortran.
   * @param mitice Output parameter for mitice flag
   * @param restart Output parameter for restart flag
   * @param assim Output parameter for assim flag
   * @param init_iter Output parameter for initial iteration
   * @param max_iter Output parameter for maximum iteration
   */
  void getConfigFlags(bool& mitice, bool& restart, bool& assim, int& init_iter,
                      int& max_iter);

  /**
   * @brief Initialize mitice components
   */
  void initializeMitice();

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
   * @brief Finalizes the environment through Fortran.
   */
  void finalizeMPI();

  /**
   * @brief Synchronize all processes using Fortran MPI_Barrier
   */
  void barrier();

  /**
   * @brief Checks if the environment was initialized by this interface.
   */
  bool isMPIInitialized() const { return mpi_initialized_by_this_instance_; }

 private:
  MPI_Comm mpi_comm_;  // C++ MPI communicator
  int f_comm_;         // Fortran format communicator
  int rank_;           // Current process rank
  int size_;           // Total number of processes
  int io_procs_;       // Number of I/O processes
  int comp_procs_;     // Number of compute processes
  bool mpi_initialized_by_this_instance_;
  bool model_components_initialized_;
};  // class MACOMFortranInterface

}  // namespace metada::backends::macom