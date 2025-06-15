/**
 * @file MACOMFortranInterface.cpp
 * @brief Implementation of the C++ interface to MACOM Fortran routines.
 * @ingroup backends_macom
 * @author Metada Framework Team
 */

#include "../include/MACOMFortranInterface.hpp"

#include <mpi.h>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

namespace metada::backends::macom {

// Internal logging functions to avoid direct dependency on MACOMlogging.hpp
namespace {
enum class LogLevel { DEBUG, INFO, WARNING, ERROR };

void logMessage(LogLevel level, const std::string& component,
                const std::string& message) {
  std::string prefix;
  switch (level) {
    case LogLevel::DEBUG:
      prefix = "[DEBUG]";
      break;
    case LogLevel::INFO:
      prefix = "[INFO]";
      break;
    case LogLevel::WARNING:
      prefix = "[WARNING]";
      break;
    case LogLevel::ERROR:
      prefix = "[ERROR]";
      break;
  }

  std::stringstream ss;
  ss << prefix << " [" << component << "] " << message;

  if (level == LogLevel::ERROR) {
    std::cerr << ss.str() << std::endl;
  } else {
    std::cout << ss.str() << std::endl;
  }
}

// void logDebug(const std::string& component, const std::string& message) {
//   logMessage(LogLevel::DEBUG, component, message);
// }

void logInfo(const std::string& component, const std::string& message) {
  logMessage(LogLevel::INFO, component, message);
}

// void logWarning(const std::string& component, const std::string& message) {
//   logMessage(LogLevel::WARNING, component, message);
// }

void logError(const std::string& component, const std::string& message) {
  logMessage(LogLevel::ERROR, component, message);
}
}  // namespace

MACOMFortranInterface::MACOMFortranInterface()
    : rank_(0),
      size_(1),
      io_procs_(0),
      comp_procs_(1),
      mpi_initialized_by_this_instance_(false),
      model_components_initialized_(false) {
  mpi_comm_ = MPI_COMM_NULL;
  f_comm_ = 0;
}

MACOMFortranInterface::~MACOMFortranInterface() {
  // logInfo("MACOMFortranInterface", "Destructor called");
  try {
    finalizeMPI();
  } catch (const std::exception& e) {
    logError(
        "MACOMFortranInterface",
        "Exception during finalizeMPI in destructor: " + std::string(e.what()));
  }
}

void MACOMFortranInterface::initializeMPI(int io_procs) {
  if (mpi_initialized_by_this_instance_) {
    logInfo("MACOMFortranInterface",
            "MPI already initialized by this instance");
    return;
  }

  // Set environment variable to tell Fortran the number of IO processors
  std::string io_procs_str = std::to_string(io_procs);
#ifdef _WIN32
  // Windows environment uses _putenv_s
  _putenv_s("MACOM_IO_PROCS", io_procs_str.c_str());
#else
  // Unix/Linux environment uses setenv
  setenv("MACOM_IO_PROCS", io_procs_str.c_str(), 1);
#endif

  // Get the C++ MPI communicator
  mpi_comm_ = MPI_COMM_WORLD;

  // // Convert to C integer for Fortran
  // f_comm_ = MPI_Comm_c2f(mpi_comm_);

  // logInfo("MACOMFortranInterface", "Initializing MPI via Fortran...");
  // Call Fortran MPI initialization with the converted communicator
  c_macom_initialize_mpi(f_comm_);

  // Get rank from Fortran
  c_macom_get_mpi_rank(&rank_);
  // logInfo("MACOMFortranInterface",
  //         "MPI initialized by Fortran. Rank: " + std::to_string(rank_) +
  //             ", I/O Procs: " + io_procs_str);

  // Get size from Fortran
  c_macom_get_mpi_size(&size_);

  // Set I/O processes and compute processes
  io_procs_ = io_procs;
  comp_procs_ = size_ - io_procs_;

  mpi_initialized_by_this_instance_ = true;
}

int MACOMFortranInterface::getRank() const {
  return rank_;
}

int MACOMFortranInterface::getSize() const {
  return size_;
}

int MACOMFortranInterface::getFortranComm() const {
  return f_comm_;
}

void MACOMFortranInterface::readNamelist() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error("MPI must be initialized before reading namelist");
  }
  // logInfo("MACOMFortranInterface", "Calling Fortran to read namelist...");
  c_macom_read_namelist();
  // logInfo("MACOMFortranInterface", "Fortran read namelist finished");
}

void MACOMFortranInterface::finalizeMPI() {
  if (!mpi_initialized_by_this_instance_) {
    return;
  }

  // logInfo("MACOMFortranInterface", "Calling Fortran to finalize MPI...");
  c_macom_finalize_mpi();

  mpi_comm_ = MPI_COMM_NULL;
  f_comm_ = 0;

  mpi_initialized_by_this_instance_ = false;
  rank_ = 0;
  // logInfo("MACOMFortranInterface", "Fortran MPI finalize finished");
}

void MACOMFortranInterface::barrier() {
  if (!mpi_initialized_by_this_instance_) {
    return;
  }
  MPI_Barrier(mpi_comm_);
}

void MACOMFortranInterface::getConfigFlags(bool& mitice, bool& restart,
                                           bool& assim, int& init_iter,
                                           int& max_iter) {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "MPI must be initialized before getting config flags");
  }
  // logInfo("MACOMFortranInterface", "Calling Fortran to get config flags...");
  c_get_macom_config_flags(&mitice, &restart, &assim, &init_iter, &max_iter);
  // logInfo("MACOMFortranInterface", "Fortran get config flags finished");
}

void MACOMFortranInterface::initializeMitice() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "MPI must be initialized before initializing mitice components");
  }
  // logInfo("MACOMFortranInterface", "Calling Fortran to initialize mitice
  // components...");
  c_macom_initialize_mitice();
  // logInfo("MACOMFortranInterface", "Fortran initialize mitice components
  // finished");
}

int MACOMFortranInterface::getCompProcs() const {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "MPI must be initialized before getting compute processes");
  }
  return comp_procs_;
}

void MACOMFortranInterface::sendInfoToIO() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "MPI must be initialized before sending info to IO");
  }
  // logInfo("MACOMFortranInterface", "Calling Fortran to send info to IO...");
  c_macom_mpi_send_info_comp_to_io();
  // logInfo("MACOMFortranInterface", "Fortran send info to IO finished");
}

void MACOMFortranInterface::openMiscRunInfo() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "MPI must be initialized before opening misc run info");
  }
  // logInfo("MACOMFortranInterface", "Calling Fortran to open misc run
  // info...");
  c_macom_misc_run_info_open();
  // logInfo("MACOMFortranInterface", "Fortran open misc run info finished");
}

void MACOMFortranInterface::initCSP() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error("MPI must be initialized before initializing CSP");
  }
  logInfo("MACOMFortranInterface", "Calling Fortran to initialize CSP...");
  c_macom_init_csp();
  logInfo("MACOMFortranInterface", "Fortran CSP initialization finished");
}

void MACOMFortranInterface::miticeInitAll() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "MPI must be initialized before initializing mitice (all)");
  }
  logInfo("MACOMFortranInterface",
          "Calling Fortran to initialize all mitice components...");
  c_macom_mitice_init_all();
  logInfo("MACOMFortranInterface",
          "Fortran mitice full initialization finished");
}

void MACOMFortranInterface::restartAndAssim() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "MPI must be initialized before restart/assimilation");
  }
  logInfo("MACOMFortranInterface",
          "Calling Fortran to handle restart and assimilation...");
  c_macom_restart_and_assim();
  logInfo("MACOMFortranInterface", "Fortran restart/assimilation finished");
}

void MACOMFortranInterface::runCspStep() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error("MPI must be initialized before running CSP step");
  }
  logInfo("MACOMFortranInterface", "Calling Fortran to run CSP step...");
  c_macom_run_csp_step();
  logInfo("MACOMFortranInterface", "Fortran CSP step finished");
}

void MACOMFortranInterface::CspIoMain() {
  logInfo("MACOMFortranInterface", "Calling Fortran IO process main...");
  c_macom_csp_io_main();
  logInfo("MACOMFortranInterface", "Fortran IO process main finished");
}

void MACOMFortranInterface::finalizeMitice() {
  logInfo("MACOMFortranInterface", "Calling Fortran to finalize mitice...");
  c_macom_finalize_mitice();
  logInfo("MACOMFortranInterface", "Fortran mitice finalize finished");
}

}  // namespace metada::backends::macom