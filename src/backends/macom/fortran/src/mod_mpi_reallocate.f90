module mod_mpi_reallocate
  ! Contains routines to re-allocate allocatable arrays
   use mod_misc_basic
   use mod_mpi_variables
   implicit none

   private

   public :: reallocate

   interface reallocate
      module procedure reallocate_d1D
      module procedure reallocate_d2D
      module procedure reallocate_d3D
      module procedure reallocate_i1D
      module procedure reallocate_i2D
      module procedure reallocate_i3D
      module procedure reallocate_id1D
      module procedure reallocate_id2D
      module procedure reallocate_id3D
      module procedure reallocate_i21D
      module procedure reallocate_i22D
   end interface

contains

   ! reallocate array of 1 dimension for wp
   subroutine reallocate_d1D(this, n)
      real(wp), allocatable, intent(inout) :: this(:) !! 1D array
      integer, intent(in) :: n !! New allocation size
      real(wp), allocatable :: tmp(:)
      integer :: n0, nTmp

      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = size(this)
      if (n == n0) return ! Don't reallocate the same size

      allocate (tmp(n))
      nTmp = min(n, n0)
      tmp(1:nTmp) = this(1:nTmp)
      deallocate (this)
      call move_alloc(from = tmp, to = this)
   end subroutine
   
   ! reallocate array of 1 dimension for wp
   subroutine reallocate_d2D(this, n)
      real(wp), allocatable, intent(inout) :: this(:, :) !! 2D array
      integer, intent(in) :: n(2) !! New allocation shape
      real(wp), allocatable :: tmp(:, :)
      integer :: n0(2), nTmp(2)
      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = shape(this)
      if (all(n == n0)) return ! Don't reallocate the same size
      allocate (tmp(n(1), n(2)))
      nTmp = [min(n(1), n0(1)), min(n(2), n0(2))]
      tmp(1:nTmp(1), 1:nTmp(2)) = this(1:nTmp(1), 1:nTmp(2))
      deallocate (this)
      call move_alloc(from = tmp, to = this)
   end subroutine

   ! reallocate array of 3 dimension for wp
   subroutine reallocate_d3D(this, n)
      real(wp), allocatable, intent(inout) :: this(:, :, :) !! 3D array
      integer, intent(in) :: n(3) !! New allocation shape
      real(wp), allocatable :: tmp(:, :, :)
      integer :: n0(3), nTmp(3)
      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = shape(this)
      if (all(n == n0)) return ! Don't reallocate the same size
      allocate (tmp(n(1), n(2), n(3)))
      nTmp(1) = min(n(1), n0(1))
      nTmp(2) = min(n(2), n0(2))
      nTmp(3) = min(n(3), n0(3))
      tmp(1:nTmp(1), 1:nTmp(2), 1:nTmp(3)) = this(1:nTmp(1), 1:nTmp(2), 1:nTmp(3))
      deallocate (this)
      call move_alloc(from=tmp, to=this)
   end subroutine
   
   ! reallocate array of 1 dimension for integer
   subroutine reallocate_i1D(this, n)
      integer, allocatable, intent(inout) :: this(:) !! 1D array
      integer, intent(in) :: n !! New allocation size
      integer, allocatable :: tmp(:)
      integer :: n0, nTmp

      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = size(this)
      if (n == n0) return ! Don't reallocate the same size

      allocate (tmp(n))
      tmp = 0
      nTmp = min(n, n0)
      tmp(1:nTmp) = this(1:nTmp)
      deallocate (this)
      call move_alloc(from=tmp, to=this)
   end subroutine
   
   ! reallocate array of 2 dimension for integer
   subroutine reallocate_i2D(this, n)
      integer, allocatable, intent(inout) :: this(:, :) !! 2D array
      integer, intent(in) :: n(2) !! New allocation shape
      integer, allocatable :: tmp(:, :)
      integer :: n0(2), nTmp(2)
      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = shape(this)
      if (all(n == n0)) return ! Don't reallocate the same size
      allocate (tmp(n(1), n(2)))
      nTmp = [min(n(1), n0(1)), min(n(2), n0(2))]
      tmp(1:nTmp(1), 1:nTmp(2)) = this(1:nTmp(1), 1:nTmp(2))
      deallocate (this)
      call move_alloc(from=tmp, to=this)
   end subroutine
   
   ! reallocate array of 3 dimension for integer
   subroutine reallocate_i3D(this, n)
      integer, allocatable, intent(inout) :: this(:, :, :) !! 3D array
      integer, intent(in) :: n(3) !! New allocation shape
      integer, allocatable :: tmp(:, :, :)
      integer :: n0(3), nTmp(3)
      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = shape(this)
      if (all(n == n0)) return ! Don't reallocate the same size
      allocate (tmp(n(1), n(2), n(3)))
      nTmp = [min(n(1), n0(1)), min(n(2), n0(2)), min(n(3), n0(3))]
      tmp(1:nTmp(1), 1:nTmp(2), 1:nTmp(3)) = this(1:nTmp(1), 1:nTmp(2), 1:nTmp(3))
      deallocate (this)
      call move_alloc(from=tmp, to=this)
   end subroutine
   
   ! reallocate array of 1 dimension for DI
   subroutine reallocate_id1D(this, n)
      integer(DI), allocatable, intent(inout) :: this(:) !! 1D array
      integer, intent(in) :: n !! New allocation size
      integer(DI), allocatable :: tmp(:)
      integer :: n0, nTmp

      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = size(this)
      if (n == n0) return ! Don't reallocate the same size

      allocate (tmp(n))
      nTmp = min(n, n0)
      tmp(1:nTmp) = this(1:nTmp)
      deallocate (this)
      call move_alloc(from=tmp, to=this)
   end subroutine
   
   ! reallocate array of 1 dimension for DI
   subroutine reallocate_id2D(this, n)
      integer(DI), allocatable, intent(inout) :: this(:, :) !! 2D array
      integer, intent(in) :: n(2) !! New allocation shape
      integer(DI), allocatable :: tmp(:, :)
      integer :: n0(2), nTmp(2)
      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = shape(this)
      if (all(n == n0)) return ! Don't reallocate the same size
      allocate (tmp(n(1), n(2)))
      tmp = 0
      nTmp = [min(n(1), n0(1)), min(n(2), n0(2))]
      tmp(1:nTmp(1), 1:nTmp(2)) = this(1:nTmp(1), 1:nTmp(2))
      deallocate (this)
      call move_alloc(from=tmp, to=this)
   end subroutine
   
   ! reallocate array of 3 dimension for DI
   subroutine reallocate_id3D(this, n)
      integer(DI), allocatable, intent(inout) :: this(:, :, :) !! 3D array
      integer, intent(in) :: n(3) !! New allocation shape
      integer(DI), allocatable :: tmp(:, :, :)
      integer :: n0(3), nTmp(3)
      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = shape(this)
      if (all(n == n0)) return ! Don't reallocate the same size
      allocate (tmp(n(1), n(2), n(3)))
      tmp = 0
      nTmp = [min(n(1), n0(1)), min(n(2), n0(2)), min(n(3), n0(3))]
      tmp(1:nTmp(1), 1:nTmp(2), 1:nTmp(3)) = this(1:nTmp(1), 1:nTmp(2), 1:nTmp(3))
      deallocate (this)
      call move_alloc(from=tmp, to=this)
   end subroutine

   ! reallocate array of 1 dimension for i2
   subroutine reallocate_i21D(this, n)
      integer(i2), allocatable, intent(inout) :: this(:) !! 1D array
      integer, intent(in) :: n !! New allocation size
      integer(i2), allocatable :: tmp(:)
      integer :: n0, nTmp

      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = size(this)
      if (n == n0) return ! Don't reallocate the same size

      allocate (tmp(n))
      nTmp = min(n, n0)
      tmp(1:nTmp) = this(1:nTmp)
      deallocate (this)
      call move_alloc(from=tmp, to=this)
   end subroutine

   ! reallocate array of 2 dimension for i2
   SUBROUTINE reallocate_i22D(this, n)

      integer(i2), allocatable, intent(inout) :: this(:, :) !! 2D array
      integer, intent(in) :: n(2) !! New allocation shape
      integer(i2), allocatable :: tmp(:, :)
      integer :: n0(2), nTmp(2)
      if (.not. allocated(this)) THEN
         WRITE (*, *) "Cannot reallocate an unallocated array"
         STOP
      END IF
      n0 = shape(this)
      if (all(n == n0)) return ! Don't reallocate the same size
      allocate (tmp(n(1), n(2)))
      tmp = 0
      nTmp = [min(n(1), n0(1)), min(n(2), n0(2))]
      tmp(1:nTmp(1), 1:nTmp(2)) = this(1:nTmp(1), 1:nTmp(2))
      deallocate (this)
      call move_alloc(from=tmp, to=this)

   END SUBROUTINE reallocate_i22D

end module mod_mpi_reallocate
