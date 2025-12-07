#include "MPIContext.hpp"

#include <stdexcept>

namespace metada::framework::parallel {

MPIContext::MPIContext()
#ifdef METADA_USE_MPI
    : comm_(MPI_COMM_WORLD),
      owns_comm_(false),
      num_procs_(1),
      my_rank_(0),
      initialized_mpi_(false)
#else
    : comm_(MPIContext::DummyComm{}),
      owns_comm_(false),
      num_procs_(1),
      my_rank_(0),
      initialized_mpi_(false)
#endif
{
#ifdef METADA_USE_MPI
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);

  if (!mpi_initialized) {
    int provided;
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SINGLE, &provided);
    initialized_mpi_ = true;
  }

  MPI_Comm_size(comm_, &num_procs_);
  MPI_Comm_rank(comm_, &my_rank_);
#else
  // Serial execution: num_procs=1, my_rank=0 (already set in initializer list)
#endif
}

#ifdef METADA_USE_MPI
MPIContext::MPIContext(MPI_Comm comm)
    : comm_(comm),
      owns_comm_(false),
      num_procs_(1),
      my_rank_(0),
      initialized_mpi_(false) {
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);

  if (!mpi_initialized) {
    int provided;
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SINGLE, &provided);
    initialized_mpi_ = true;
  }

  MPI_Comm_size(comm_, &num_procs_);
  MPI_Comm_rank(comm_, &my_rank_);
}
#endif

MPIContext::~MPIContext() {
#ifdef METADA_USE_MPI
  if (initialized_mpi_) {
    int mpi_finalized = 0;
    MPI_Finalized(&mpi_finalized);
    if (!mpi_finalized) {
      MPI_Finalize();
    }
  }

  if (owns_comm_ && comm_ != MPI_COMM_NULL) {
    MPI_Comm_free(&comm_);
  }
#endif
}

MPIContext::MPIContext(MPIContext&& other) noexcept
#ifdef METADA_USE_MPI
    : comm_(other.comm_),
      owns_comm_(other.owns_comm_),
      num_procs_(other.num_procs_),
      my_rank_(other.my_rank_),
      initialized_mpi_(other.initialized_mpi_)
#else
    : comm_(other.comm_),
      owns_comm_(other.owns_comm_),
      num_procs_(other.num_procs_),
      my_rank_(other.my_rank_),
      initialized_mpi_(other.initialized_mpi_)
#endif
{
#ifdef METADA_USE_MPI
  other.comm_ = MPI_COMM_NULL;
  other.owns_comm_ = false;
#endif
  other.num_procs_ = 1;
  other.my_rank_ = 0;
  other.initialized_mpi_ = false;
}

MPIContext& MPIContext::operator=(MPIContext&& other) noexcept {
  if (this != &other) {
#ifdef METADA_USE_MPI
    if (owns_comm_ && comm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&comm_);
    }
    comm_ = other.comm_;
    owns_comm_ = other.owns_comm_;
    other.comm_ = MPI_COMM_NULL;
    other.owns_comm_ = false;
#endif
    num_procs_ = other.num_procs_;
    my_rank_ = other.my_rank_;
    initialized_mpi_ = other.initialized_mpi_;

    other.num_procs_ = 1;
    other.my_rank_ = 0;
    other.initialized_mpi_ = false;
  }
  return *this;
}

int MPIContext::numProcs() const noexcept {
  return num_procs_;
}

int MPIContext::myRank() const noexcept {
  return my_rank_;
}

bool MPIContext::isRoot() const noexcept {
  return my_rank_ == 0;
}

#ifdef METADA_USE_MPI
MPI_Comm MPIContext::communicator() const noexcept {
  return comm_;
}
#endif

void MPIContext::barrier() const {
#ifdef METADA_USE_MPI
  MPI_Barrier(comm_);
#else
  // No-op for serial execution
#endif
}

}  // namespace metada::framework::parallel
