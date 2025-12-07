#pragma once

#ifdef METADA_USE_MPI
#include <mpi.h>
#endif

#include <cstddef>

namespace metada::framework::parallel {

/**
 * @brief MPI context for managing MPI communicator and process information
 *
 * @details This class provides a RAII wrapper for MPI initialization and
 * communicator management. It supports both MPI and serial execution:
 * - With MPI: Initializes MPI if needed, manages communicator
 * - Without MPI: Returns num_procs=1, my_rank=0 (serial execution)
 *
 * The class is non-copyable but movable, following RAII principles.
 *
 * Example usage:
 * @code
 * MPIContext mpi_ctx;
 * if (mpi_ctx.isRoot()) {
 *   // Root process logic
 * }
 * int global_sum = mpi_ctx.communication().allReduceSum(local_value);
 * @endcode
 */
class MPIContext {
 public:
  /**
   * @brief Initialize MPI context
   * @details If MPI is not initialized, initializes it.
   *          If already initialized, uses existing communicator.
   *          For serial execution (no MPI), sets num_procs=1, my_rank=0.
   */
  MPIContext();

#ifdef METADA_USE_MPI
  /**
   * @brief Initialize with specific communicator
   * @param comm MPI communicator (default: MPI_COMM_WORLD)
   */
  explicit MPIContext(MPI_Comm comm);
#endif

  /**
   * @brief Destructor
   * @details Finalizes MPI if this instance initialized it.
   */
  ~MPIContext();

  // Non-copyable, movable
  MPIContext(const MPIContext&) = delete;
  MPIContext& operator=(const MPIContext&) = delete;
  MPIContext(MPIContext&&) noexcept;
  MPIContext& operator=(MPIContext&&) noexcept;

  /**
   * @brief Get number of processes
   * @return Total number of MPI processes (1 for serial execution)
   */
  [[nodiscard]] int numProcs() const noexcept;

  /**
   * @brief Get current process rank
   * @return Process rank (0-based, 0 for serial execution)
   */
  [[nodiscard]] int myRank() const noexcept;

  /**
   * @brief Check if this is the root process
   * @return true if rank == 0
   */
  [[nodiscard]] bool isRoot() const noexcept;

#ifdef METADA_USE_MPI
  /**
   * @brief Get MPI communicator
   * @return MPI communicator
   */
  [[nodiscard]] MPI_Comm communicator() const noexcept;
#endif

  /**
   * @brief Synchronize all processes
   * @details Calls MPI_Barrier if MPI is available, no-op otherwise
   */
  void barrier() const;

 private:
#ifdef METADA_USE_MPI
  MPI_Comm comm_;
  bool owns_comm_;  // Whether we created a duplicate communicator
#else
  // Dummy type for serial execution
  struct DummyComm {};
  DummyComm comm_;
  bool owns_comm_;
#endif
  int num_procs_;
  int my_rank_;
  bool initialized_mpi_;  // Whether we initialized MPI
};

}  // namespace metada::framework::parallel
