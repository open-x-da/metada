#include "FortranInterface.hpp"
#include "MaComState.hpp"
#include <iostream>
#include <cstdlib>
#include <string>

namespace applications {
namespace macom {

// ==================================================================
// I. Construction and Setup
// ==================================================================

FortranInterface::FortranInterface(MPI_Comm comm, int io_procs)
    : mpi_comm_(comm), io_procs_(io_procs), rank_(-1),
      mitice_on_(false), restart_in_(false), assim_in_(false),
      init_iter_(0), max_iter_(0),
      data_copy_enabled_(true) {
    
    // Set environment variable to tell Fortran the number of IO processors
    std::string io_procs_str = std::to_string(io_procs_);
#ifdef _WIN32
    // Windows 环境使用 _putenv_s
    _putenv_s("MACOM_IO_PROCS", io_procs_str.c_str());
#else
    // Unix/Linux 环境使用 setenv
    setenv("MACOM_IO_PROCS", io_procs_str.c_str(), 1);
#endif

    // Convert communicator to Fortran format
    f_comm_ = MPI_Comm_c2f(mpi_comm_);
}

// ==================================================================
// II. MPI Process Related Functions
// ==================================================================

// Initialize MPI processes
void FortranInterface::initializeMPI() {
    // Use saved Fortran format communicator to call Fortran interface
    c_mpi_process_init(f_comm_);
    
    // Get current process rank
    c_get_mpi_rank(&rank_);
}

// Get current process rank
int FortranInterface::getRank() const {
    return rank_;
}

// Check if this is the main process
bool FortranInterface::isMainProcess() const {
    return rank_ == 0;
}

// Check if this is a compute process
bool FortranInterface::isComputeProcess() const {
    return c_is_compute_process();
}

// ==================================================================
// III. Configuration and Initialization Functions
// ==================================================================

// Read configuration
void FortranInterface::readNamelist() {
    c_read_namelist();
    
    // Immediately sync flags after reading config
    syncConfigFlags();
}

// ==================================================================
// IV. Sea Ice Module Functions
// ==================================================================

void FortranInterface::initMiticeRunInfo() {
    c_init_mitice_run_info();
}

void FortranInterface::readMiticeParams() {
    c_read_mitice_params();
}

void FortranInterface::allocateMitice() {
    c_allocate_mitice();
}

void FortranInterface::initMiticeFixed() {
    c_init_mitice_fixed();
}

void FortranInterface::updateGpuSeaiceParameters() {
    c_update_gpu_seaice_parameters();
}

void FortranInterface::initMiticeVars() {
    c_init_mitice_vars();
}

void FortranInterface::finalizeMitice() {
    c_finalize_mitice();
}

// ==================================================================
// V. Data Exchange Functions
// ==================================================================

void FortranInterface::sendInfoToIO() {
    c_send_info_to_io();
}

void FortranInterface::sendSecondInfoToIO() {
    c_send_second_info_to_io();
}

void FortranInterface::sendOutput() {
    c_send_output();
}

// ==================================================================
// VI. Computation Functions
// ==================================================================

// CSP related
void FortranInterface::initCSP() {
    c_init_csp();
}

void FortranInterface::initCspAsm() {
    c_init_csp_asm();
}

// Execute single CSP computation
void FortranInterface::callCSP(int iter) {
    c_call_csp(iter);
}

// Execute CSP loop computation
void FortranInterface::callCSP2() {
    c_call_csp2();
}

// ==================================================================
// VII. IO Process Functions
// ==================================================================

void FortranInterface::initIO() {
    c_init_io();
}

void FortranInterface::runIOLoop() {
    // Call mpi_csp_io_main function in Fortran
    c_csp_io_main();
}

// ==================================================================
// VIII. Restart Related Functions
// ==================================================================

void FortranInterface::readRestart() {
    c_read_restart();
}

// ==================================================================
// IX. Run Information Functions
// ==================================================================

void FortranInterface::openMiscRunInfo() {
    c_open_misc_run_info();
}

void FortranInterface::closeMiscRunInfo() {
    c_close_misc_run_info();
}

// ==================================================================
// X. GPU Related Functions
// ==================================================================

void FortranInterface::initGPU() {
    c_init_gpu();
}

// ==================================================================
// XI. Status Query Functions
// ==================================================================

bool FortranInterface::isMiticeEnabled() const {
    return c_is_mitice_enabled();
}

bool FortranInterface::isRestartEnabled() const {
    return c_is_restart_enabled();
}

bool FortranInterface::isAssimEnabled() const {
    return c_is_assim_enabled();
}

// ==================================================================
// XII. Cleanup Functions
// ==================================================================

void FortranInterface::finalizeMPI() {
    c_finalize_mpi();
}

// ==================================================================
// XIII. Helper Functions
// ==================================================================

// Set iteration value
void FortranInterface::setIterationValue(int iter) {
    // In some cases we might need to pass the current iteration value to Fortran
    // This function is currently empty, but might need implementation in the actual project
}

// Sync configuration flags
void FortranInterface::syncConfigFlags() {
    c_get_config_flags(&mitice_on_, &restart_in_, &assim_in_, 
                     &init_iter_, &max_iter_);
    
    // if (isMainProcess()) {
    //     std::cout << "Configuration synchronized: "
    //             << "mitice=" << (mitice_on_ ? "ON" : "OFF")
    //             << ", restart=" << (restart_in_ ? "ON" : "OFF")
    //             << ", assim=" << (assim_in_ ? "ON" : "OFF")
    //             << ", iter=" << init_iter_ << "-" << max_iter_
    //             << std::endl;
    // }
}

}} // namespace applications::macom