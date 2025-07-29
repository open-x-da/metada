/**
 * @file MACOMFortranInterface.hpp
 * @brief C++ interface to MACOM Fortran routines.
 * @ingroup backends_macom
 * @author Metada Framework Team
 */

#pragma once

// MACOM_MODE_ENABLED is defined by CMake
#if MACOM_MODE_ENABLED
#include <mpi.h>
#else
// Define a dummy MPI_Comm type when MPI is not available
typedef void* MPI_Comm;
#define MPI_COMM_NULL nullptr
#define MPI_COMM_WORLD nullptr
#endif

// Forward declare the C functions from the Fortran wrapper
extern "C" {
#if MACOM_MODE_ENABLED
void c_macom_initialize_mpi(int comm_c);
void c_macom_get_mpi_rank(int* rank);
void c_macom_get_mpi_size(int* size);
void c_macom_read_namelist();
void c_macom_initialize_mitice();
void c_get_macom_config_flags(bool* mitice, bool* restart, bool* assim,
                              int* init_iter, int* max_iter);
void c_macom_finalize_mpi();
void c_macom_mpi_send_info_comp_to_io();
void c_macom_misc_run_info_open();
void c_macom_init_csp();
void c_macom_mitice_init_all();
void c_macom_restart_and_assim();
void c_macom_run_csp_step();
void c_macom_csp_io_main();
void c_macom_finalize_mitice();
#endif
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
   * @brief Checks if the environment was initialized by this interface.
   * @return true if MPI was initialized by this instance, false otherwise
   */
  bool isMPIInitialized() const { return mpi_initialized_by_this_instance_; }

  /**
   * @brief Synchronize all processes using MPI_Barrier
   */
  void barrier();

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
   * @brief Initialize mitice components (no-op for region mode)
   */
  void initializeMitice();

  /**
   * @brief Finalizes the environment through Fortran.
   */
  void finalizeMPI();

  /**
   * @brief Send information from compute processes to I/O processes
   */
  void sendInfoToIO();

  /**
   * @brief Open misc run info
   */
  void openMiscRunInfo();

  /**
   * @brief Initialize CSP
   */
  void initCSP();

  /**
   * @brief Initialize mitice components (no-op for region mode)
   */
  void miticeInitAll();

  /**
   * @brief Handle restart and assimilation
   */
  void restartAndAssim();

  /**
   * @brief Run CSP step
   */
  void runCspStep();

  /**
   * @brief Run CSP I/O main
   */
  void CspIoMain();

  /**
   * @brief Finalize mitice components (for global mode)
   */
  void finalizeMitice();

 private:
  MPI_Comm mpi_comm_;  // C++ MPI communicator
  int f_comm_;         // Fortran format communicator
  int rank_;           // Current process rank
  int size_;           // Total number of processes
  int io_procs_;       // Number of I/O processes
  int comp_procs_;     // Number of compute processes
  bool mpi_initialized_by_this_instance_;
  bool model_components_initialized_;
};

}  // namespace metada::backends::macom