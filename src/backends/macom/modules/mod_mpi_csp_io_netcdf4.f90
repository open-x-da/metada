module mod_mpi_csp_io_netcdf4
   use netcdf
   use mod_misc
   use mod_mpi_variables
   use mod_mpi_interfaces
   implicit none

   ! parallel writing netcdf files
   interface mpi_netcdf_write
      module procedure mpi_netcdf_write_array_1d_wp
      module procedure mpi_netcdf_write_array_2d_wp
   end interface

contains

   ! parallel writing wp data of 1 dimension
   subroutine mpi_netcdf_write_array_1d_wp(filename, varname, var, fp, d1, d1_p, &
                                           starts, counts, ts)
      implicit none
      integer :: ncid, varid, dimid_d1, dimid_t, status, status_var, status_def, &
                 status_dim
      integer, intent(in) :: ts
      character(*), intent(in) :: filename, varname, fp
      character(lc) :: time, date, zone, timestamp
      character(lc) :: dimname1
      integer, intent(in) :: d1, d1_p
      real(wp), intent(in) :: var(d1_p)
      integer(i1), allocatable, dimension(:) :: var_i1
      integer(i2), allocatable, dimension(:) :: var_i2
      integer(i4), allocatable, dimension(:) :: var_i4
      real(sp), allocatable, dimension(:) :: var_sp
      real(dp), allocatable, dimension(:) :: var_dp
      integer, intent(in) :: starts
      integer, intent(in) :: counts
      integer :: starts_ts(2)
      integer :: counts_ts(2)
      logical :: file_exist

      starts_ts(1) = starts
      starts_ts(2) = ts
      counts_ts(1) = counts
      counts_ts(2) = 1

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
      status = nf90_open(trim(filename), IOR(NF90_WRITE, NF90_MPIIO), ncid, &
                         comm = mpi_io_comm, info = MPI_INFO_NULL)

      if (file_exist) then
         if (status == nf90_noerr) then
            status_var = nf90_inq_varid(ncid, trim(varname), varid)
            if (status_var == nf90_noerr) then
               continue
            else
               if (mpi_rank == mpi_comp_procs) write (*, *) trim(varname)// &
                  ' not exist in '// trim(filename)//', we create it now'
               status_def = nf90_redef(ncid)
               status_dim = nf90_inq_dimid(ncid, trim(dimname1), dimid_d1)
               if (status_dim .ne. nf90_noerr) then
                  if (mpi_rank == mpi_comp_procs) write (*, *) trim(dimname1)// &
                     ' not exist in '// trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, 'time', dimid_t)
               if (status_dim .ne. nf90_noerr) then
                  if (mpi_rank == mpi_comp_procs) write (*, *) 'time not exist in '// &
                     trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
               end if
               select case (trim(fp))
               case ('i1')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                                 (/ dimid_d1, dimid_t /), varid), trim(filename))
               case ('i2')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                                 (/ dimid_d1, dimid_t /), varid), trim(filename))
               case ('i4')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                                 (/ dimid_d1, dimid_t /), varid), trim(filename))
               case ('sp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                                 (/ dimid_d1, dimid_t /), varid), trim(filename))
               case ('dp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                                 (/ dimid_d1, dimid_t /), varid), trim(filename))
               case default
                  if (mpi_rank == mpi_comp_procs) then
                     write (*, *) 'ERROR ! : define '// trim(varname)//' for '// trim(filename)
                     write (*, *) 'unknown nf90 type = '// fp
                     write (*, *) 'support nf90 type is :'
                     write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
                     write (*, *) 'i4 for nf90_int, sp for nf90_float'
                     write (*, *) 'dp for nf90_double'
                  end if
                  stop
               end select
               call netcdf_check(nf90_enddef(ncid), trim(filename))
            end if
         else
            if (mpi_rank == mpi_comp_procs) write (*, *) trim(filename)// &
               ' open '// trim(filename)//' failed !'
            stop
         end if
      else
         if (mpi_rank == mpi_comp_procs) write (*, *) 'Create output file '// &
            trim(filename)
         call netcdf_check(nf90_create(trim(filename), IOR(NF90_NETCDF4, NF90_MPIIO), &
                                       ncid, comm = mpi_io_comm, info = MPI_INFO_NULL), trim(filename), varname = trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
         call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                  "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), trim(filename))
         call date_and_time(date = date, time = time, zone = zone)
         timestamp = date(7:8)//"/"// date(5:6)//"/"// date(1:4)//" "// &
                     time(1:2)//":"// time(3:4)//":"// time(5:6)//" "// zone
         call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), trim(filename))
         select case (trim(fp))
         case ('i1')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                           (/ dimid_d1, dimid_t /), varid), trim(filename))
         case ('i2')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                           (/ dimid_d1, dimid_t /), varid), trim(filename))
         case ('i4')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                           (/ dimid_d1, dimid_t /), varid), trim(filename))
         case ('sp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                           (/ dimid_d1, dimid_t /), varid), trim(filename))
         case ('dp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                           (/ dimid_d1, dimid_t /), varid), trim(filename))
         case default
            if (mpi_rank == mpi_comp_procs) then
               write (*, *) 'ERROR ! : define '// trim(varname)//' for '// trim(filename)
               write (*, *) 'unknown nf90 type = '// fp
               write (*, *) 'support nf90 type is :'
               write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
               write (*, *) 'i4 for nf90_int, sp for nf90_float'
               write (*, *) 'dp for nf90_double'
            end if
            stop
         end select
         call netcdf_check(nf90_enddef(ncid), trim(filename))
      end if

      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename))
      call netcdf_check(nf90_var_par_access(ncid, varid, nf90_collective), trim(filename))
      if (mpi_rank == mpi_comp_procs) &
         write (*, "(a,a8,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
         ' write ', trim(varname), ' for date: ', nowyearf, '-', nowmonf, '-', nowdayf, &
         ' ', nowhourf, ':', nowminf, ':', nowsecf, ' to file '//trim(filename)//' at ts = ', ts
      select case (trim(fp))
      case ('i1')
         allocate (var_i1(d1_p))
         var_i1 = int(var, 1)
         deallocate (var_i1)
      case ('i2')
         allocate (var_i2(d1_p))
         var_i2 = int(var, 2)
         deallocate (var_i2)
      case ('i4')
         allocate (var_i4(d1_p))
         var_i4 = int(var, 4)
         deallocate (var_i4)
      case ('sp')
         allocate (var_sp(d1_p))
         var_sp = real(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_sp, starts_ts, counts_ts), &
                           trim(filename), varname = trim(varname))
         deallocate (var_sp)
      case ('dp')
         allocate (var_dp(d1_p))
         var_dp = dble(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_dp, starts_ts, counts_ts), &
                           trim(filename), varname = trim(varname))
         deallocate (var_dp)
      case default
         if (mpi_rank == mpi_comp_procs) then
            write (*, *) 'ERROR ! : write '//trim(varname)//' for '//trim(filename)
            write (*, *) 'unknown nf90 type = '//fp
            write (*, *) 'support nf90 type is :'
            write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
            write (*, *) 'i4 for nf90_int, sp for nf90_float'
            write (*, *) 'dp for nf90_double'
         end if
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine mpi_netcdf_write_array_1d_wp

   ! parallel writing wp data of 2 dimension
   subroutine mpi_netcdf_write_array_2d_wp(filename, varname, var, fp, d1, d2, &
                                           d1_p, d2_p, starts, counts, ts)
      implicit none
      integer :: ncid, varid, dimid_d1, dimid_d2, dimid_t, status, status_var, &
                 status_def, status_dim
      integer, intent(in) :: ts
      character(*), intent(in) :: filename, varname, fp
      character(lc) :: time, date, zone, timestamp
      character(lc) :: dimname1, dimname2
      integer, intent(in) :: d1, d2, d1_p, d2_p
      real(wp), intent(in) :: var(d1_p, d2_p)
      integer(i1), allocatable, dimension(:, :) :: var_i1
      integer(i2), allocatable, dimension(:, :) :: var_i2
      integer(i4), allocatable, dimension(:, :) :: var_i4
      real(sp), allocatable, dimension(:, :) :: var_sp
      real(dp), allocatable, dimension(:, :) :: var_dp
      integer, dimension(:), intent(in) :: starts
      integer, dimension(:), intent(in) :: counts
      integer :: starts_ts(3)
      integer :: counts_ts(3)
      logical :: file_exist

      starts_ts(1:2) = starts(1:2)
      starts_ts(3) = ts
      counts_ts(1:2) = counts(1:2)
      counts_ts(3) = 1

      if (d1 == mpi_total_nlpb) then
         dimname1 = 'nlpb'
      else if (d1 == mpi_total_nlpbz) then
         dimname1 = 'nlpbz'
      else if (d1 == nk) then
         dimname1 = 'nk'
      else if (d1 == nkp1) then
         dimname1 = 'nkp1'
      else
         write (*, *) &
            'ERROR ! : d1 not defined in model, now only support nlpb, nlpbz,', &
            ' nk or nkp1, but d1 = ', d1
         stop
      end if

      if (d2 == nk) then
         dimname2 = 'nk'
      else if (d2 == nkp1) then
         dimname2 = 'nkp1'
      else
         write (*, *) &
            'ERROR ! : d2 not defined in model, now only support nk or nkp1,', &
            ' but d2 = ', d2
         stop
      end if

      inquire (file = trim(filename), exist = file_exist)
      status = nf90_open(trim(filename), IOR(NF90_WRITE, NF90_MPIIO), ncid, &
                         comm = mpi_io_comm, info = MPI_INFO_NULL)

      if (file_exist) then
         if (status == nf90_noerr) then
            status_var = nf90_inq_varid(ncid, trim(varname), varid)
            if (status_var == nf90_noerr) then
               continue
            else
               if (mpi_rank == mpi_comp_procs) write (*, *) trim(varname)// &
                  ' not exist in '// trim(filename)//', we create it now'
               status_def = nf90_redef(ncid)
               status_dim = nf90_inq_dimid(ncid, trim(dimname1), dimid_d1)
               if (status_dim .ne. nf90_noerr) then
                  if (mpi_rank == mpi_comp_procs) write (*, *) trim(dimname1)// &
                     ' not exist in '// trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, &
                                                 dimid_d1), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, trim(dimname2), dimid_d2)
               if (status_dim .ne. nf90_noerr) then
                  if (mpi_rank == mpi_comp_procs) write (*, *) trim(dimname2)// &
                     ' not exist in '// trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname2), d2, &
                                                 dimid_d2), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, 'time', dimid_t)
               if (status_dim .ne. nf90_noerr) then
                  if (mpi_rank == mpi_comp_procs) write (*, *) 'time not exist in '// &
                     trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, &
                                                 dimid_t), trim(filename))
               end if

               select case (trim(fp))
               case ('i1')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                                 (/ dimid_d1, dimid_d2, dimid_t /), varid), trim(filename))
               case ('i2')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                                 (/ dimid_d1, dimid_d2, dimid_t /), varid), trim(filename))
               case ('i4')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                                 (/ dimid_d1, dimid_d2, dimid_t /), varid), trim(filename))
               case ('sp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                                 (/ dimid_d1, dimid_d2, dimid_t /), varid), trim(filename))
               case ('dp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                                 (/ dimid_d1, dimid_d2, dimid_t /), varid), trim(filename))
               case default
                  if (mpi_rank == mpi_comp_procs) then
                     write (*, *) 'ERROR ! : define '//trim(varname)//' for '// trim(filename)
                     write (*, *) 'unknown nf90 type = '//fp
                     write (*, *) 'support nf90 type is :'
                     write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
                     write (*, *) 'i4 for nf90_int, sp for nf90_float'
                     write (*, *) 'dp for nf90_double'
                  end if
                  stop
               end select
               call netcdf_check(nf90_enddef(ncid), trim(filename))
            end if
         else
            if (mpi_rank == mpi_comp_procs) write (*, *) trim(filename)// &
               ' open '//trim(filename)//' failed !'
            stop
         end if
      else
         if (mpi_rank == mpi_comp_procs) write (*, *) 'Create output file '// trim(filename)
         call netcdf_check(nf90_create(trim(filename), IOR(NF90_NETCDF4, NF90_MPIIO), &
                                       ncid, comm = mpi_io_comm, info = MPI_INFO_NULL), trim(filename), varname = trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname2), d2, dimid_d2), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
         call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                  "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), trim(filename))
         call date_and_time(date = date, time = time, zone = zone)
         timestamp = date(7:8)//"/"// date(5:6)//"/"// date(1:4)//" "// &
                     time(1:2)//":"// time(3:4)//":"// time(5:6)//" "// zone
         call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), trim(filename))
         select case (trim(fp))
         case ('i1')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                           (/ dimid_d1, dimid_d2, dimid_t /), varid), trim(filename))
         case ('i2')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                           (/ dimid_d1, dimid_d2, dimid_t /), varid), trim(filename))
         case ('i4')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                           (/ dimid_d1, dimid_d2, dimid_t /), varid), trim(filename))
         case ('sp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                           (/ dimid_d1, dimid_d2, dimid_t /), varid), trim(filename))
         case ('dp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                           (/ dimid_d1, dimid_d2, dimid_t /), varid), trim(filename))
         case default
            if (mpi_rank == mpi_comp_procs) then
               write (*, *) 'ERROR ! : define '// trim(varname)//' for '// trim(filename)
               write (*, *) 'unknown nf90 type = '// fp
               write (*, *) 'support nf90 type is :'
               write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
               write (*, *) 'i4 for nf90_int, sp for nf90_float'
               write (*, *) 'dp for nf90_double'
            end if
            stop
         end select
         call netcdf_check(nf90_enddef(ncid), trim(filename))
      end if

      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename))
      call netcdf_check(nf90_var_par_access(ncid, varid, nf90_collective), trim(filename))
      if (mpi_rank == mpi_comp_procs) &
         write (*, "(a,a8,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
         ' write ', trim(varname), ' for date: ', nowyearf, '-', nowmonf, '-', nowdayf, &
         ' ', nowhourf, ':', nowminf, ':', nowsecf, ' to file '// trim(filename)//' at ts = ', ts
      select case (trim(fp))
      case ('i1')
         allocate (var_i1(d1_p, d2_p))
         var_i1 = int(var, 1)
         deallocate (var_i1)
      case ('i2')
         allocate (var_i2(d1_p, d2_p))
         var_i2 = int(var, 2)
         deallocate (var_i2)
      case ('i4')
         allocate (var_i4(d1_p, d2_p))
         var_i4 = int(var, 4)
         deallocate (var_i4)
      case ('sp')
         allocate (var_sp(d1_p, d2_p))
         var_sp = real(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_sp, starts_ts, counts_ts), trim(filename))
         deallocate (var_sp)
      case ('dp')
         allocate (var_dp(d1_p, d2_p))
         var_dp = dble(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_dp, starts_ts, counts_ts), trim(filename))
         deallocate (var_dp)
      case default
         if (mpi_rank == mpi_comp_procs) then
            write (*, *) 'ERROR ! : write '// trim(varname)//' for '// trim(filename)
            write (*, *) 'unknown nf90 type = '// fp
            write (*, *) 'support nf90 type is :'
            write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
            write (*, *) 'i4 for nf90_int, sp for nf90_float'
            write (*, *) 'dp for nf90_double'
         end if
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine mpi_netcdf_write_array_2d_wp

end module mod_mpi_csp_io_netcdf4
