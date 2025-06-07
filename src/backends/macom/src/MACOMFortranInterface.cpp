/**
 * @file MACOMFortranInterface.cpp
 * @brief Implementation of the C++ interface to MACOM Fortran routines.
 * @ingroup backends_macom
 * @author Metada Framework Team
 */

#include "../include/MACOMFortranInterface.hpp"

#include <mpi.h>

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
      io_procs_(0),
      mpi_initialized_by_this_instance_(false),
      model_components_initialized_(false) {
  mpi_comm_ = MPI_COMM_NULL;
  f_comm_ = 0;
}

MACOMFortranInterface::~MACOMFortranInterface() {
  logInfo("MACOMFortranInterface", "Destructor called");
  if (model_components_initialized_) {
    try {
      finalizeModelComponents();
    } catch (const std::exception& e) {
      logError("MACOMFortranInterface",
               "Exception during finalizeModelComponents in destructor: " +
                   std::string(e.what()));
    }
  }
  if (mpi_initialized_by_this_instance_) {
    try {
      finalizeMPI();
    } catch (const std::exception& e) {
      logError("MACOMFortranInterface",
               "Exception during finalizeMPI in destructor: " +
                   std::string(e.what()));
    }
  }
}

void MACOMFortranInterface::initializeMPI() {
  if (mpi_initialized_by_this_instance_) {
    logInfo("MACOMFortranInterface",
            "MPI already initialized by this instance");
    return;
  }

  // Get the C++ MPI communicator
  mpi_comm_ = MPI_COMM_WORLD;

  // Convert to C integer for Fortran
  f_comm_ = MPI_Comm_c2f(mpi_comm_);

  logInfo("MACOMFortranInterface", "Initializing MPI via Fortran...");
  // Call Fortran MPI initialization with the converted communicator
  c_macom_initialize_mpi(f_comm_);

  // Get rank from Fortran
  c_macom_get_mpi_rank(&rank_);

  // Get total number of processes
  int total_procs;
  MPI_Comm_size(mpi_comm_, &total_procs);
  // Assuming I/O processes are the last quarter of total processes
  io_procs_ = total_procs / 4;

  mpi_initialized_by_this_instance_ = true;
  logInfo("MACOMFortranInterface",
          "MPI initialized by Fortran. Rank: " + std::to_string(rank_) +
              ", I/O Procs: " + std::to_string(io_procs_));
}

int MACOMFortranInterface::getRank() const {
  return rank_;
}

int MACOMFortranInterface::getIOProcs() const {
  return io_procs_;
}

int MACOMFortranInterface::getFortranComm() const {
  return f_comm_;
}

void MACOMFortranInterface::readNamelist() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error("MPI must be initialized before reading namelist");
  }
  logInfo("MACOMFortranInterface", "Calling Fortran to read namelist...");
  c_macom_read_namelist();
  logInfo("MACOMFortranInterface", "Fortran read namelist finished");
}

void MACOMFortranInterface::initializeModelComponents() {
  if (!mpi_initialized_by_this_instance_) {
    throw std::runtime_error(
        "MPI must be initialized before initializing model components");
  }
  logInfo("MACOMFortranInterface",
          "Calling Fortran to initialize model components...");
  c_macom_initialize_model_components();
  model_components_initialized_ = true;
  logInfo("MACOMFortranInterface",
          "Fortran initialize model components finished");
}

void MACOMFortranInterface::runModelStep(int current_iteration) {
  if (!model_components_initialized_) {
    throw std::runtime_error(
        "Model components must be initialized before running a step");
  }
  logInfo("MACOMFortranInterface",
          "Calling Fortran to run model step for iteration: " +
              std::to_string(current_iteration));
  int status = 0;
  c_macom_run_model_step(current_iteration, &status);
  if (status != 0) {
    throw std::runtime_error("Fortran model step reported an error. Status: " +
                             std::to_string(status));
  }
  logInfo("MACOMFortranInterface",
          "Fortran run model step finished for iteration: " +
              std::to_string(current_iteration));
}

void MACOMFortranInterface::finalizeModelComponents() {
  if (!model_components_initialized_) {
    return;
  }
  logInfo("MACOMFortranInterface",
          "Calling Fortran to finalize model components...");
  c_macom_finalize_model_components();
  model_components_initialized_ = false;
  logInfo("MACOMFortranInterface",
          "Fortran finalize model components finished");
}

void MACOMFortranInterface::finalizeMPI() {
  if (!mpi_initialized_by_this_instance_) {
    return;
  }

  logInfo("MACOMFortranInterface", "Calling Fortran to finalize MPI...");
  c_macom_finalize_mpi();

  mpi_comm_ = MPI_COMM_NULL;
  f_comm_ = 0;

  mpi_initialized_by_this_instance_ = false;
  rank_ = 0;
  io_procs_ = 0;
  logInfo("MACOMFortranInterface", "Fortran MPI finalize finished");
}

}  // namespace metada::backends::macom