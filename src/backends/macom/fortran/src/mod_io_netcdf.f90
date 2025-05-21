module mod_io_netcdf
   use netcdf
   use mod_misc_basic
   use mod_csp_basic
   use mod_mpi_variables
!   use mod_misc
   implicit none

   interface netcdf_read
      module procedure netcdf_read_dimension
      module procedure netcdf_read_array_1d_i2
      module procedure netcdf_read_array_2d_i2
      module procedure netcdf_read_array_3d_i2
      module procedure netcdf_read_array_1d_i4
      module procedure netcdf_read_array_2d_i4
      module procedure netcdf_read_array_3d_i4
      module procedure netcdf_read_scalar_i8
      module procedure netcdf_read_scalar_wp
      module procedure netcdf_read_array_1d_wp
      module procedure netcdf_read_array_2d_wp
      module procedure netcdf_read_array_3d_wp
   end interface

   interface netcdf_write
      module procedure netcdf_write_array_0d_dp
      module procedure netcdf_write_array_1d_wp
      module procedure netcdf_write_array_2d_wp
   end interface

contains

!==============================================================================
   subroutine netcdf_check(status, filename, varname, errorFlag)
!==============================================================================
      implicit none
      integer, intent(in) :: status
      integer, intent(inout), optional :: errorFlag
      character(*), intent(in) :: filename
      character(*), intent(in), optional :: varname

      if (status /= nf90_noerr) then
         write (run_out_unit, *) 'ERROR! : '//trim(nf90_strerror(status))
         if (present(errorFlag)) then
            errorFlag = errorFlag + status
         else
            write (run_out_unit, *) filename, " var/dim/file name is ", varname
            write (run_out_unit, *) "*** nmefc mcom io failed ***"
            stop
         end if
      end if

   end subroutine netcdf_check

!==============================================================================
   subroutine netcdf_read_dimension(filename, dimname, dimlen)
!==============================================================================
      implicit none
      integer :: ncid, dimid
      integer, intent(inout) :: dimlen
      character(*), intent(in) :: filename, dimname

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_dimid(ncid, dimname, dimid), trim(filename), varname=trim(dimname))
      call netcdf_check(nf90_inquire_dimension(ncid, dimid, len=dimlen), trim(filename), varname=trim(dimname))
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
      write (run_out_unit, *) 'inquire '//trim(dimname)//' from '//trim(filename)
   end subroutine netcdf_read_dimension

!==============================================================================
   subroutine netcdf_read_array_1d_i2(filename, varname, var, d1, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1
      integer(i2), intent(inout) :: var(d1)
      integer(i1), allocatable, dimension(:) :: var_i1
      integer(i2), allocatable, dimension(:) :: var_i2
      integer(i4), allocatable, dimension(:) :: var_i4
      real(sp), allocatable, dimension(:) :: var_sp
      real(dp), allocatable, dimension(:) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1))
         var = int(var_i1, i2)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i2, &
                                           start=(/1, ts/), count=(/d1, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i2), trim(filename), varname=trim(varname))
         end if
         var = var_i2
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i4, &
                                           start=(/1, ts/), count=(/d1, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i4), trim(filename), varname=trim(varname))
         end if
         var = int(var_i4, i2)
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, ts/), count=(/d1, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = int(var_sp, i2)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1))
         var = int(var_dp, i2)
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_1d_i2

!==============================================================================
   subroutine netcdf_read_array_1d_i4(filename, varname, var, d1, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1
      integer(i4), intent(inout) :: var(d1)
      integer(i1), allocatable, dimension(:) :: var_i1
      integer(i2), allocatable, dimension(:) :: var_i2
      integer(i4), allocatable, dimension(:) :: var_i4
      real(sp), allocatable, dimension(:) :: var_sp
      real(dp), allocatable, dimension(:) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1))
         var = int(var_i1, i4)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i2, &
                                           start=(/1, ts/), count=(/d1, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i2), trim(filename), varname=trim(varname))
         end if
         var = int(var_i2, i4)
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i4, &
                                           start=(/1, ts/), count=(/d1, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i4), trim(filename), varname=trim(varname))
         end if
         var = var_i4
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, ts/), count=(/d1, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = int(var_sp, i4)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1))
         var = int(var_dp, i4)
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_1d_i4

!==============================================================================
   subroutine netcdf_read_scalar_i4(filename, varname, var, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer(i4), intent(inout) :: var
      integer(i1) :: var_i1(1)
      integer(i2) :: var_i2(1)
      integer(i4) :: var_i4(1)
      real(sp) :: var_sp(1)
      real(dp) :: var_dp(1)

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         var = int(var_i1(1), i4)
      case (nf90_short)
         var = int(var_i2(1), i4)
      case (nf90_int)
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i4, &
                                           start=(/ts/), count=(/1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i4), trim(filename), varname=trim(varname))
         end if
         var = var_i4(1)
      case (nf90_float)
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/ts/), count=(/1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = int(var_sp(1), i4)
      case (nf90_double)
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_dp, &
                                           start=(/ts/), count=(/1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_dp), trim(filename), varname=trim(varname))
         end if
         var = int(var_dp(1), i4)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_scalar_i4

!==============================================================================
   subroutine netcdf_read_scalar_i8(filename, varname, var, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer(i8), intent(inout) :: var
      integer(i1) :: var_i1(1)
      integer(i2) :: var_i2(1)
      integer(i4) :: var_i4(1)
      real(sp) :: var_sp(1)
      real(dp) :: var_dp(1)

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         var = int(var_i1(1), i8)
      case (nf90_short)
         var = int(var_i2(1), i8)
      case (nf90_int)
         var = int(var_i4(1), i8)
      case (nf90_float)
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/ts/), count=(/1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = int(var_sp(1), i8)
      case (nf90_double)
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_dp, &
                                           start=(/ts/), count=(/1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_dp), trim(filename), varname=trim(varname))
         end if
         var = int(var_dp(1), i8)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_scalar_i8

!==============================================================================
   subroutine netcdf_read_scalar_wp(filename, varname, var, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      real(wp), intent(inout) :: var
      integer(i1) :: var_i1(1)
      integer(i2) :: var_i2(1)
      integer(i4) :: var_i4(1)
      real(sp) :: var_sp(1)
      real(dp) :: var_dp(1)

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         var = real(var_i1(1), wp)
      case (nf90_short)
         var = real(var_i2(1), wp)
      case (nf90_int)
         var = real(var_i4(1), wp)
      case (nf90_float)
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/ts/), count=(/1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = real(var_sp(1), wp)
      case (nf90_double)
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_dp, &
                                           start=(/ts/), count=(/1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_dp), trim(filename), varname=trim(varname))
         end if
         var = real(var_dp(1), wp)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_scalar_wp

!==============================================================================
   subroutine netcdf_read_array_1d_wp(filename, varname, var, d1, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1
      real(wp), intent(inout) :: var(d1)
      integer(i1), allocatable, dimension(:) :: var_i1
      integer(i2), allocatable, dimension(:) :: var_i2
      integer(i4), allocatable, dimension(:) :: var_i4
      real(sp), allocatable, dimension(:) :: var_sp
      real(dp), allocatable, dimension(:) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1))
         var = real(var_i1, wp)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1))
         var = real(var_i2, wp)
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1))
         var = real(var_i4, wp)
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, ts/), count=(/d1, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = real(var_sp, wp)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_dp, &
                                           start=(/1, ts/), count=(/d1, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_dp), trim(filename), varname=trim(varname))
         end if
         var = real(var_dp, wp)
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_1d_wp

!==============================================================================
   subroutine netcdf_read_array_2d_i2(filename, varname, var, d1, d2, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1, d2
      integer(i2), intent(inout) :: var(d1, d2)
      integer(i1), allocatable, dimension(:, :) :: var_i1
      integer(i2), allocatable, dimension(:, :) :: var_i2
      integer(i4), allocatable, dimension(:, :) :: var_i4
      real(sp), allocatable, dimension(:, :) :: var_sp
      real(dp), allocatable, dimension(:, :) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1, d2))
         var = int(var_i1, i2)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i2, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i2), trim(filename), varname=trim(varname))
         end if
         var = var_i2
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i4, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i4), trim(filename), varname=trim(varname))
         end if
         var = int(var_i4, i2)
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = int(var_sp, i2)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1, d2))
         var = int(var_dp, i2)
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_2d_i2

!==============================================================================
   subroutine netcdf_read_array_2d_i4(filename, varname, var, d1, d2, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1, d2
      integer(i4), intent(inout) :: var(d1, d2)
      integer(i1), allocatable, dimension(:, :) :: var_i1
      integer(i2), allocatable, dimension(:, :) :: var_i2
      integer(i4), allocatable, dimension(:, :) :: var_i4
      real(sp), allocatable, dimension(:, :) :: var_sp
      real(dp), allocatable, dimension(:, :) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1, d2))
         var = int(var_i1, i4)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i2, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i2), trim(filename), varname=trim(varname))
         end if
         var = int(var_i2, i4)
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i4, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i4), trim(filename), varname=trim(varname))
         end if
         var = var_i4
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = int(var_sp, i4)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1, d2))
         var = int(var_dp, i4)
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_2d_i4

!==============================================================================
   subroutine netcdf_read_array_2d_wp(filename, varname, var, d1, d2, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1, d2
      real(wp), intent(inout) :: var(d1, d2)
      integer(i1), allocatable, dimension(:, :) :: var_i1
      integer(i2), allocatable, dimension(:, :) :: var_i2
      integer(i4), allocatable, dimension(:, :) :: var_i4
      real(sp), allocatable, dimension(:, :) :: var_sp
      real(dp), allocatable, dimension(:, :) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1, d2))
         var = real(var_i1, wp)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i2, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i2), trim(filename), varname=trim(varname))
         end if
         var = real(var_i2, wp)
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1, d2))
         var = real(var_i4, wp)
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = real(var_sp, wp)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_dp, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_dp), trim(filename), varname=trim(varname))
         end if
         var = real(var_dp, wp)
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_2d_wp

!==============================================================================
   subroutine netcdf_read_array_3d_i2(filename, varname, var, d1, d2, d3, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1, d2, d3
      integer(i2), intent(inout) :: var(d1, d2, d3)
      integer(i1), allocatable, dimension(:, :, :) :: var_i1
      integer(i2), allocatable, dimension(:, :, :) :: var_i2
      integer(i4), allocatable, dimension(:, :, :) :: var_i4
      real(sp), allocatable, dimension(:, :, :) :: var_sp
      real(dp), allocatable, dimension(:, :, :) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1, d2, d3))
         var = int(var_i1, i2)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i2, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i2), trim(filename), varname=trim(varname))
         end if
         var = var_i2
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1, d2, d3))
         var = int(var_i4, i2)
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = int(var_sp, i2)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1, d2, d3))
         var = int(var_dp, i2)
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(varname))
   end subroutine netcdf_read_array_3d_i2

!==============================================================================
   subroutine netcdf_read_array_3d_i4(filename, varname, var, d1, d2, d3, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1, d2, d3
      integer(i4), intent(inout) :: var(d1, d2, d3)
      integer(i1), allocatable, dimension(:, :, :) :: var_i1
      integer(i2), allocatable, dimension(:, :, :) :: var_i2
      integer(i4), allocatable, dimension(:, :, :) :: var_i4
      real(sp), allocatable, dimension(:, :, :) :: var_sp
      real(dp), allocatable, dimension(:, :, :) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1, d2, d3))
         var = int(var_i1, i4)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i2, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i2), trim(filename), varname=trim(varname))
         end if
         var = int(var_i2, i4)
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i4, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i4), trim(filename), varname=trim(varname))
         end if
         var = var_i4
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = int(var_sp, i4)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1, d2, d3))
         var = int(var_dp, i4)
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_3d_i4

!==============================================================================
   subroutine netcdf_read_array_3d_wp(filename, varname, var, d1, d2, d3, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1, d2, d3
      real(wp), intent(inout) :: var(d1, d2, d3)
      integer(i1), allocatable, dimension(:, :, :) :: var_i1
      integer(i2), allocatable, dimension(:, :, :) :: var_i2
      integer(i4), allocatable, dimension(:, :, :) :: var_i4
      real(sp), allocatable, dimension(:, :, :) :: var_sp
      real(dp), allocatable, dimension(:, :, :) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      if (present(ts)) then
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)//' at time', ts
      else
         write (run_out_unit, *) 'read '//trim(varname)//' from '//trim(filename)
      end if
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1, d2, d3))
         var = real(var_i1, wp)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i2, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i2), trim(filename), varname=trim(varname))
         end if
         var = real(var_i2, wp)
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1, d2, d3))
         var = real(var_i4, wp)
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = real(var_sp, wp)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_dp, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_dp), trim(filename), varname=trim(varname))
         end if
         var = real(var_dp, wp)
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_3d_wp

!==============================================================================
   subroutine netcdf_write_array_0d_dp(filename, varname, var, fp, ts, nc_attr1)
!==============================================================================
      implicit none
      integer :: ncid, varid, dimid_t, status, status_var, status_def, status_dim
      integer, intent(in) :: ts
      character(*), intent(in) :: filename, varname, fp
      character(lc) :: time, date, zone, timestamp
      real(dp), intent(in) :: var
      integer(i1) :: var_i1(1)
      integer(i2) :: var_i2(1)
      integer(i4) :: var_i4(1)
      real(sp) :: var_sp(1)
      real(dp) :: var_dp(1)
      logical :: file_exist
      type(nc_attr) :: nc_attr1

      inquire (file=trim(filename), exist=file_exist)
      status = nf90_open(trim(filename), nf90_write, ncid)

      if (file_exist) then
         if (status == nf90_noerr) then
            status_var = nf90_inq_varid(ncid, trim(varname), varid)
            if (status_var == nf90_noerr) then
               continue
            else
               write (*, *) trim(varname)//' not exist in '//trim(filename)//', we create it now'
               status_def = nf90_redef(ncid)
               status_dim = nf90_inq_dimid(ncid, 'time', dimid_t)
               if (status_dim .ne. nf90_noerr) then
                  write (*, *) 'time not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
               end if
               select case (trim(fp))
               case ('i1')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                                 (/dimid_t/), varid), trim(filename))
               case ('i2')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                                 (/dimid_t/), varid), trim(filename))
               case ('i4')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                                 (/dimid_t/), varid), trim(filename))
               case ('sp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                                 (/dimid_t/), varid), trim(filename))
               case ('dp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                                 (/dimid_t/), varid), trim(filename))
               case default
                  write (*, *) 'ERROR! : define '//trim(varname)//' for ' &
                     //trim(filename)
                  write (*, *) 'unknown nf90 type = '//fp
                  write (*, *) 'support nf90 type is :'
                  write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
                  write (*, *) 'i4 for nf90_int, sp for nf90_float'
                  write (*, *) 'dp for nf90_double'
                  stop
               end select
               call netcdf_check(nf90_put_att(ncid, varid, "units", trim(nc_attr1%var_units)), trim(filename))
               call netcdf_check(nf90_enddef(ncid), trim(filename))
            end if
         else
            write (*, *) trim(filename)//' open '//trim(filename)//' failed!'
            stop
         end if
      else
         write (*, *) 'Create output file '//trim(filename)
         call netcdf_check(nf90_create(trim(filename), nf90_netcdf4, ncid), trim(filename), varname=trim(filename))
         call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
         call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                  "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), trim(filename))
         call date_and_time(date=date, time=time, zone=zone)
         timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                     time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//zone
         call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), trim(filename))
         select case (trim(fp))
         case ('i1')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                           (/dimid_t/), varid), trim(filename))
         case ('i2')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                           (/dimid_t/), varid), trim(filename))
         case ('i4')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                           (/dimid_t/), varid), trim(filename))
         case ('sp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                           (/dimid_t/), varid), trim(filename))
         case ('dp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                           (/dimid_t/), varid), trim(filename))
         case default
            write (*, *) 'ERROR! : define '//trim(varname)//' for ' &
               //trim(filename)
            write (*, *) 'unknown nf90 type = '//fp
            write (*, *) 'support nf90 type is :'
            write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
            write (*, *) 'i4 for nf90_int, sp for nf90_float'
            write (*, *) 'dp for nf90_double'
            stop
         end select
         call netcdf_check(nf90_put_att(ncid, varid, "units", trim(nc_attr1%var_units)), trim(filename))
         call netcdf_check(nf90_enddef(ncid), trim(filename))
      end if

      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename))
      ! call netcdf_check(nf90_inquire_variable(ncid, varid), trim(filename))
      write (run_out_unit, "(a,a8,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
         ' write ', trim(varname), ' for date: ', nowyearf, '-', nowmonf, '-', nowdayf, &
         ' ', nowhourf, ':', nowminf, ':', nowsecf, ' to file '//trim(filename)//' at ts = ', ts
      select case (trim(fp))
      case ('i1')
         var_i1(1) = int(var, i1)
      case ('i2')
         var_i2(1) = int(var, i2)
      case ('i4')
         var_i4(1) = int(var, i4)
      case ('sp')
         var_sp(1) = real(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_sp, &
                                        start=(/ts/), count=(/1/)), trim(filename))
      case ('dp')
         var_dp(1) = dble(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_dp, &
                                        start=(/ts/), count=(/1/)), trim(filename))
      case default
         write (run_out_unit, *) 'ERROR! : write '//trim(varname)//' for ' &
            //trim(filename)
         write (*, *) 'unknown nf90 type = '//fp
         write (*, *) 'support nf90 type is :'
         write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
         write (*, *) 'i4 for nf90_int, sp for nf90_float'
         write (*, *) 'dp for nf90_double'
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine netcdf_write_array_0d_dp

!==============================================================================
   subroutine netcdf_write_array_1d_wp(filename, varname, var, fp, d1, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, dimid_d1, dimid_t, status, status_var, status_def, status_dim
      integer, intent(in) :: ts
      character(*), intent(in) :: filename, varname, fp
      character(lc) :: time, date, zone, timestamp
      character(lc) :: dimname1
      integer, intent(in) :: d1
      real(wp), intent(in) :: var(d1)
      integer(i1), allocatable, dimension(:) :: var_i1
      integer(i2), allocatable, dimension(:) :: var_i2
      integer(i4), allocatable, dimension(:) :: var_i4
      real(sp), allocatable, dimension(:) :: var_sp
      real(dp), allocatable, dimension(:) :: var_dp
      logical :: file_exist

      if (d1 == mpi_total_nlpb) then
         dimname1 = 'nlpb'
      else if (d1 == mpi_total_nlpbz) then
         dimname1 = 'nlpbz'
      else if (d1 == nk) then
         dimname1 = 'nk'
      else if (d1 == nkp1) then
         dimname1 = 'nkp1'
      end if

      inquire (file = trim(filename), exist = file_exist)
      status = nf90_open(trim(filename), nf90_write, ncid)

      if (file_exist) then
         if (status == nf90_noerr) then
            status_var = nf90_inq_varid(ncid, trim(varname), varid)
            if (status_var == nf90_noerr) then
               continue
            else
               write (run_out_unit, *) trim(varname)//' not exist in '//trim(filename)//', we create it now'
               status_def = nf90_redef(ncid)
               status_dim = nf90_inq_dimid(ncid, trim(dimname1), dimid_d1)
               if (status_dim .ne. nf90_noerr) then
                  write (run_out_unit, *) trim(dimname1)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, 'time', dimid_t)
               if (status_dim .ne. nf90_noerr) then
                  write (run_out_unit, *) 'time not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
               end if
               select case (trim(fp))
               case ('i1')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                                 (/dimid_d1, dimid_t/), varid), trim(filename))
               case ('i2')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                                 (/dimid_d1, dimid_t/), varid), trim(filename))
               case ('i4')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                                 (/dimid_d1, dimid_t/), varid), trim(filename))
               case ('sp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                                 (/dimid_d1, dimid_t/), varid), trim(filename))
               case ('dp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                                 (/dimid_d1, dimid_t/), varid), trim(filename))
               case default
                  write (run_out_unit, *) 'ERROR! : define '//trim(varname)//' for ' &
                     //trim(filename)
                  write (run_out_unit, *) 'unknown nf90 type = '//fp
                  write (run_out_unit, *) 'support nf90 type is :'
                  write (run_out_unit, *) 'i1 for nf90_byte, i2 for nf90_short'
                  write (run_out_unit, *) 'i4 for nf90_int, sp for nf90_float'
                  write (run_out_unit, *) 'dp for nf90_double'
                  stop
               end select
               call netcdf_check(nf90_enddef(ncid), trim(filename))
            end if
         else
            write (run_out_unit, *) trim(filename)//' open '//trim(filename)//' failed!'
            stop
         end if
      else
         write (run_out_unit, *) 'Create output file '//trim(filename)
         call netcdf_check(nf90_create(trim(filename), nf90_netcdf4, ncid), trim(filename), varname=trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
         call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                  "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), trim(filename))
         call date_and_time(date=date, time=time, zone=zone)
         timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                     time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//zone
         call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), trim(filename))
         select case (trim(fp))
         case ('i1')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                           (/dimid_d1, dimid_t/), varid), trim(filename))
         case ('i2')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                           (/dimid_d1, dimid_t/), varid), trim(filename))
         case ('i4')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                           (/dimid_d1, dimid_t/), varid), trim(filename))
         case ('sp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                           (/dimid_d1, dimid_t/), varid), trim(filename))
         case ('dp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                           (/dimid_d1, dimid_t/), varid), trim(filename))
         case default
            write (run_out_unit, *) 'ERROR! : define '//trim(varname)//' for ' &
               //trim(filename)
            write (run_out_unit, *) 'unknown nf90 type = '//fp
            write (run_out_unit, *) 'support nf90 type is :'
            write (run_out_unit, *) 'i1 for nf90_byte, i2 for nf90_short'
            write (run_out_unit, *) 'i4 for nf90_int, sp for nf90_float'
            write (run_out_unit, *) 'dp for nf90_double'
            stop
         end select
         call netcdf_check(nf90_enddef(ncid), trim(filename))
      end if

      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename))
      ! call netcdf_check(nf90_inquire_variable(ncid, varid), trim(filename))
      write (run_out_unit, "(a,a15,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
         ' write ', trim(varname), ' for date: ', nowyearf, '-', nowmonf, '-', nowdayf, &
         ' ', nowhourf, ':', nowminf, ':', nowsecf, ' to file '//trim(filename)//' at ts = ', ts
      select case (trim(fp))
      case ('i1')
         allocate (var_i1(d1))
         var_i1 = int(var, i1)
         deallocate (var_i1)
      case ('i2')
         allocate (var_i2(d1))
         var_i2 = int(var, i2)
         deallocate (var_i2)
      case ('i4')
         allocate (var_i4(d1))
         var_i4 = int(var, i4)
         deallocate (var_i4)
      case ('sp')
         allocate (var_sp(d1))
         var_sp = real(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_sp, &
                                        start=(/1, ts/), count=(/d1, 1/)), trim(filename))
         deallocate (var_sp)
      case ('dp')
         allocate (var_dp(d1))
         var_dp = dble(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_dp, &
                                        start=(/1, ts/), count=(/d1, 1/)), trim(filename))
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : write '//trim(varname)//' for ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = '//fp
         write (run_out_unit, *) 'support nf90 type is :'
         write (run_out_unit, *) 'i1 for nf90_byte, i2 for nf90_short'
         write (run_out_unit, *) 'i4 for nf90_int, sp for nf90_float'
         write (run_out_unit, *) 'dp for nf90_double'
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine netcdf_write_array_1d_wp

!==============================================================================
   subroutine netcdf_write_array_2d_wp(filename, varname, var, fp, d1, d2, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, dimid_d1, dimid_d2, dimid_t, status, status_var, status_def, status_dim
      integer, intent(in) :: ts
      character(*), intent(in) :: filename, varname, fp
      character(lc) :: time, date, zone, timestamp
      character(lc) :: dimname1, dimname2
      integer, intent(in) :: d1, d2
      real(wp), intent(in) :: var(d1, d2)
      integer(i1), allocatable, dimension(:, :) :: var_i1
      integer(i2), allocatable, dimension(:, :) :: var_i2
      integer(i4), allocatable, dimension(:, :) :: var_i4
      real(sp), allocatable, dimension(:, :) :: var_sp
      real(dp), allocatable, dimension(:, :) :: var_dp
      logical :: file_exist

      if (d1 == mpi_total_nlpb) then
         dimname1 = 'nlpb'
      else if (d1 == mpi_total_nlpbz) then
         dimname1 = 'nlpbz'
      else if (d1 == nk) then
         dimname1 = 'nk'
      else if (d1 == nkp1) then
         dimname1 = 'nkp1'
      else
         write (run_out_unit, *) 'ERROR! : d1 not defined in model, now only support nlpb, nlpbz, nk or nkp1,', &
            ' but d1 = ', d1
         stop
      end if

      if (d2 == nk) then
         dimname2 = 'nk'
      else if (d2 == nkp1) then
         dimname2 = 'nkp1'
      else
         write (run_out_unit, *) 'ERROR! : d2 not defined in model, now only support nk or nkp1,', &
            ' but d2 = ', d2
         stop
      end if

      inquire (file=trim(filename), exist=file_exist)
      status = nf90_open(trim(filename), nf90_write, ncid)

      if (file_exist) then
         if (status == nf90_noerr) then
            status_var = nf90_inq_varid(ncid, trim(varname), varid)
            if (status_var == nf90_noerr) then
               continue
            else
               write (run_out_unit, *) trim(varname)//' not exist in '//trim(filename)//', we create it now'

               status_def = nf90_redef(ncid)
               status_dim = nf90_inq_dimid(ncid, trim(dimname1), dimid_d1)
               if (status_dim .ne. nf90_noerr) then
                  write (run_out_unit, *) trim(dimname1)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, trim(dimname2), dimid_d2)
               if (status_dim .ne. nf90_noerr) then
                  write (run_out_unit, *) trim(dimname2)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname2), d2, dimid_d2), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, 'time', dimid_t)
               if (status_dim .ne. nf90_noerr) then
                  write (run_out_unit, *) 'time not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
               end if

               select case (trim(fp))
               case ('i1')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('i2')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('i4')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('sp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('dp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case default
                  write (run_out_unit, *) 'ERROR! : define '//trim(varname)//' for ' &
                     //trim(filename)
                  write (run_out_unit, *) 'unknown nf90 type = '//fp
                  write (run_out_unit, *) 'support nf90 type is :'
                  write (run_out_unit, *) 'i1 for nf90_byte, i2 for nf90_short'
                  write (run_out_unit, *) 'i4 for nf90_int, sp for nf90_float'
                  write (run_out_unit, *) 'dp for nf90_double'
                  stop
               end select
               call netcdf_check(nf90_enddef(ncid), trim(filename))
            end if
         else
            write (run_out_unit, *) trim(filename)//' open '//trim(filename)//' failed!'
            stop
         end if
      else
         write (run_out_unit, *) 'Create output file '//trim(filename)
         call netcdf_check(nf90_create(trim(filename), nf90_netcdf4, ncid), trim(filename), varname=trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname2), d2, dimid_d2), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
         call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                  "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), trim(filename))
         call date_and_time(date=date, time=time, zone=zone)
         timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                     time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//zone
         call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), trim(filename))
         select case (trim(fp))
         case ('i1')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('i2')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('i4')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('sp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('dp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case default
            write (run_out_unit, *) 'ERROR! : define '//trim(varname)//' for ' &
               //trim(filename)
            write (run_out_unit, *) 'unknown nf90 type = '//fp
            write (run_out_unit, *) 'support nf90 type is :'
            write (run_out_unit, *) 'i1 for nf90_byte, i2 for nf90_short'
            write (run_out_unit, *) 'i4 for nf90_int, sp for nf90_float'
            write (run_out_unit, *) 'dp for nf90_double'
            stop
         end select
         call netcdf_check(nf90_enddef(ncid), trim(filename))
      end if

      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename))
      ! call netcdf_check(nf90_inquire_variable(ncid, varid), trim(filename))
      write (run_out_unit, "(a,a8,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
         ' write ', trim(varname), ' for date: ', nowyearf, '-', nowmonf, '-', nowdayf, &
         ' ', nowhourf, ':', nowminf, ':', nowsecf, ' to file '//trim(filename)//' at ts = ', ts
      select case (trim(fp))
      case ('i1')
         allocate (var_i1(d1, d2))
         var_i1 = int(var, i1)
         deallocate (var_i1)
      case ('i2')
         allocate (var_i2(d1, d2))
         var_i2 = int(var, i2)
         deallocate (var_i2)
      case ('i4')
         allocate (var_i4(d1, d2))
         var_i4 = int(var, i4)
         deallocate (var_i4)
      case ('sp')
         allocate (var_sp(d1, d2))
         var_sp = real(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_sp, &
                                        start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename))
         deallocate (var_sp)
      case ('dp')
         allocate (var_dp(d1, d2))
         var_dp = dble(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_dp, &
                                        start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename))
         deallocate (var_dp)
      case default
         write (run_out_unit, *) 'ERROR! : write '//trim(varname)//' for ' &
            //trim(filename)
         write (run_out_unit, *) 'unknown nf90 type = '//fp
         write (run_out_unit, *) 'support nf90 type is :'
         write (run_out_unit, *) 'i1 for nf90_byte, i2 for nf90_short'
         write (run_out_unit, *) 'i4 for nf90_int, sp for nf90_float'
         write (run_out_unit, *) 'dp for nf90_double'
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine netcdf_write_array_2d_wp

end module mod_io_netcdf
