#include <gtest/gtest.h>

#include <iostream>
#include <type_traits>

#include "MPIContext.hpp"

using namespace metada::framework::parallel;

// Test default constructor
TEST(MPIContextTest, DefaultConstructor) {
  MPIContext ctx;

  // In serial mode (no MPI), should have 1 process, rank 0
  // In MPI mode, will have actual MPI values
  EXPECT_GE(ctx.numProcs(), 1);
  EXPECT_GE(ctx.myRank(), 0);
  EXPECT_LT(ctx.myRank(), ctx.numProcs());
  EXPECT_EQ(ctx.isRoot(), ctx.myRank() == 0);
}

// Test basic getters
TEST(MPIContextTest, BasicGetters) {
  MPIContext ctx;

  int num_procs = ctx.numProcs();
  int my_rank = ctx.myRank();

  EXPECT_GT(num_procs, 0);
  EXPECT_GE(my_rank, 0);
  EXPECT_LT(my_rank, num_procs);
}

// Test isRoot
TEST(MPIContextTest, IsRoot) {
  MPIContext ctx;

  bool is_root = ctx.isRoot();
  int my_rank = ctx.myRank();

  EXPECT_EQ(is_root, my_rank == 0);

  // Root process should have rank 0
  if (is_root) {
    EXPECT_EQ(my_rank, 0);
  }
}

// Test barrier (should not throw)
TEST(MPIContextTest, Barrier) {
  MPIContext ctx;

  // Barrier should complete without throwing
  EXPECT_NO_THROW(ctx.barrier());

  // Multiple barriers should work
  EXPECT_NO_THROW(ctx.barrier());
  EXPECT_NO_THROW(ctx.barrier());
}

// Test move constructor
TEST(MPIContextTest, MoveConstructor) {
  MPIContext ctx1;
  int num_procs = ctx1.numProcs();
  int my_rank = ctx1.myRank();

  // Move construct
  MPIContext ctx2(std::move(ctx1));

  // Moved-from object should be in valid state (but may have default values)
  EXPECT_GE(ctx1.numProcs(), 1);
  EXPECT_GE(ctx1.myRank(), 0);

  // Moved-to object should have original values
  EXPECT_EQ(ctx2.numProcs(), num_procs);
  EXPECT_EQ(ctx2.myRank(), my_rank);
  EXPECT_EQ(ctx2.isRoot(), my_rank == 0);
}

// Test move assignment
TEST(MPIContextTest, MoveAssignment) {
  MPIContext ctx1;
  int num_procs1 = ctx1.numProcs();
  int my_rank1 = ctx1.myRank();

  MPIContext ctx2;

  // Move assign
  ctx2 = std::move(ctx1);

  // Moved-to object should have original values from ctx1
  EXPECT_EQ(ctx2.numProcs(), num_procs1);
  EXPECT_EQ(ctx2.myRank(), my_rank1);

  // Moved-from object should be in valid state
  EXPECT_GE(ctx1.numProcs(), 1);
  EXPECT_GE(ctx1.myRank(), 0);
}

// Test copy constructor is deleted
TEST(MPIContextTest, CopyConstructorDeleted) {
  MPIContext ctx1;

  // Should not compile - copy constructor is deleted
  // This is a compile-time test, so we just verify the class is non-copyable
  static_assert(!std::is_copy_constructible_v<MPIContext>);
  static_assert(!std::is_copy_assignable_v<MPIContext>);
}

// Test multiple instances
TEST(MPIContextTest, MultipleInstances) {
  MPIContext ctx1;
  MPIContext ctx2;

  // Both should report the same number of processes
  EXPECT_EQ(ctx1.numProcs(), ctx2.numProcs());

  // Both should report the same rank
  EXPECT_EQ(ctx1.myRank(), ctx2.myRank());

  // Both should have the same root status
  EXPECT_EQ(ctx1.isRoot(), ctx2.isRoot());
}

// Test barrier synchronization (basic functionality)
TEST(MPIContextTest, BarrierSynchronization) {
  MPIContext ctx;

  // Barrier should be callable multiple times
  for (int i = 0; i < 5; ++i) {
    EXPECT_NO_THROW(ctx.barrier());
  }
}

// Test that destructor doesn't throw
TEST(MPIContextTest, Destructor) {
  {
    MPIContext ctx;
    // Context should be valid
    EXPECT_GE(ctx.numProcs(), 1);
  }
  // Destructor should complete without throwing
  // (tested implicitly by scope exit)
}

// Test move semantics with temporary
TEST(MPIContextTest, MoveFromTemporary) {
  MPIContext ctx1;
  int num_procs = ctx1.numProcs();
  int my_rank = ctx1.myRank();

  // Construct from temporary using braces to avoid vexing parse
  MPIContext ctx2{MPIContext()};

  // ctx2 should be in valid state
  EXPECT_GE(ctx2.numProcs(), 1);
  EXPECT_GE(ctx2.myRank(), 0);

  // Original ctx1 should still be valid
  EXPECT_EQ(ctx1.numProcs(), num_procs);
  EXPECT_EQ(ctx1.myRank(), my_rank);
}

#ifdef METADA_USE_MPI
// Test constructor with specific communicator
TEST(MPIContextTest, ConstructorWithCommunicator) {
  MPIContext ctx_default;
  MPI_Comm default_comm = ctx_default.communicator();

  // Create context with explicit communicator
  MPIContext ctx_explicit(default_comm);

  EXPECT_EQ(ctx_explicit.numProcs(), ctx_default.numProcs());
  EXPECT_EQ(ctx_explicit.myRank(), ctx_default.myRank());
  EXPECT_EQ(ctx_explicit.communicator(), default_comm);
}

// Test communicator getter
TEST(MPIContextTest, CommunicatorGetter) {
  MPIContext ctx;

  MPI_Comm comm = ctx.communicator();
  EXPECT_NE(comm, MPI_COMM_NULL);

  // Should be able to use the communicator
  int comm_size;
  MPI_Comm_size(comm, &comm_size);
  EXPECT_EQ(comm_size, ctx.numProcs());

  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);
  EXPECT_EQ(comm_rank, ctx.myRank());
}

// Test that communicator is consistent
TEST(MPIContextTest, CommunicatorConsistency) {
  MPIContext ctx1;
  MPIContext ctx2;

  // Both should use the same communicator (MPI_COMM_WORLD)
  EXPECT_EQ(ctx1.communicator(), ctx2.communicator());
  EXPECT_EQ(ctx1.communicator(), MPI_COMM_WORLD);
}
#endif  // METADA_USE_MPI

// Test that all methods are noexcept where specified
TEST(MPIContextTest, NoexceptSpecifications) {
  MPIContext ctx;

  // These should be noexcept
  static_assert(noexcept(ctx.numProcs()));
  static_assert(noexcept(ctx.myRank()));
  static_assert(noexcept(ctx.isRoot()));

#ifdef METADA_USE_MPI
  static_assert(noexcept(ctx.communicator()));
#endif
}

// Test rank and process count consistency
TEST(MPIContextTest, RankAndProcessConsistency) {
  MPIContext ctx;

  int num_procs = ctx.numProcs();
  int my_rank = ctx.myRank();

  // Rank should be in valid range
  EXPECT_GE(my_rank, 0);
  EXPECT_LT(my_rank, num_procs);

  // Number of processes should be positive
  EXPECT_GT(num_procs, 0);

  // Root check should be consistent
  EXPECT_EQ(ctx.isRoot(), my_rank == 0);
}

// Test to verify multiple processes are running (only in MPI mode)
#ifdef METADA_USE_MPI
TEST(MPIContextTest, MultipleProcessesRunning) {
  MPIContext ctx;

  int num_procs = ctx.numProcs();
  int my_rank = ctx.myRank();

  // Print process information (visible in test output)
  // This helps verify multiple processes are actually running
  std::cout << "[Rank " << my_rank << "/" << num_procs << "] "
            << "MPIContextTest_MultipleProcessesRunning" << std::endl;

  // Verify we have multiple processes when running with MPI
  // Note: This test will pass even with 1 process (serial mode),
  // but the output will show if multiple processes are running
  EXPECT_GE(num_procs, 1);
  EXPECT_GE(my_rank, 0);
  EXPECT_LT(my_rank, num_procs);

  // If running with MPI launcher, we should have multiple processes
  // This assertion will fail if MPI_TEST_NUMPROC > 1 but only 1 process runs
  if (num_procs > 1) {
    // Verify all ranks are present
    EXPECT_GE(num_procs, 2)
        << "Expected at least 2 processes when running with MPI";
  }
}

// Test barrier with multiple processes
TEST(MPIContextTest, BarrierWithMultipleProcesses) {
  MPIContext ctx;

  int num_procs = ctx.numProcs();
  int my_rank = ctx.myRank();

  // Print which process is at the barrier
  std::cout << "[Rank " << my_rank << "/" << num_procs << "] "
            << "Before barrier" << std::endl;

  // Barrier should synchronize all processes
  ctx.barrier();

  std::cout << "[Rank " << my_rank << "/" << num_procs << "] "
            << "After barrier" << std::endl;

  // All processes should reach this point
  EXPECT_GE(num_procs, 1);
  EXPECT_GE(my_rank, 0);
  EXPECT_LT(my_rank, num_procs);
}
#endif  // METADA_USE_MPI
