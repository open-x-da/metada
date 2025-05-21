/**
 * @file MACOMFortranInterface.hpp
 * @brief C++ interface to MACOM Fortran routines.
 * @ingroup backends_macom
 * @author Metada Framework Team
 */

#pragma once

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
   */
  void initializeMPI();

  /**
   * @brief Gets the rank of the current process from Fortran.
   * @return The rank.
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
   * @brief Finalizes the environment through Fortran.
   */
  void finalizeMPI();

  /**
   * @brief Checks if the environment was initialized by this interface.
   */
  bool isMPIInitialized() const { return mpi_initialized_by_this_instance_; }

 private:
  int rank_;
  bool mpi_initialized_by_this_instance_;
  bool model_components_initialized_;
};  // class MACOMFortranInterface

}  // namespace metada::backends::macom