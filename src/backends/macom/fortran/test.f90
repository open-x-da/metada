module mpi_initest
  use mod_mpi_variables
  use mpi
  implicit none
contains

  subroutine init_mpi

    implicit none
    INTEGER :: rank, size
    integer :: ierr

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    
    write(*,*) "rank = ", rank

  end subroutine init_mpi

end module mpi_initest