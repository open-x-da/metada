module mod_mpi_csp_io_bdy
   use mod_misc
   use mod_io_netcdf
   use mod_mpi_variables
   use mod_mpi_interfaces
   implicit none

contains

!==============================================================================
   subroutine mpi_csp_io_rdbdy
!==============================================================================
      implicit none
      integer :: bdy_ts, i, k
      integer :: bdyyear, bdymon, bdyday, bdyhour, bdymin, bdysec
      integer(i8) :: nowjulVEC, nowjulT, nowjulS, nowjulPbt
      character(lc) :: filebdy
      character(3) :: pbtname
      real(wp) :: dfjulvec, dfjult, dfjuls, dfjulpbt
      real(wp) :: fbvecr, fbtr, fbsr, fbpbtr
      logical :: readnxtVEC, readnxtT, readnxtS, readnxtPbt

      if (ln_bous) then
         pbtname = "ssh"
      else
         pbtname = "pbt"
      end if

      readnxtVEC = .false.
      readnxtT = .false.
      readnxtS = .false.
      readnxtPbt = .false.

      !if ((myIter .eq. 1) .and. (.not. (restart_in .and. date_res))) then
      if (myIter .eq. 1) then
         nxtjulVEC = nowjulian + int(fbvec, 8)
         nxtjulT = nowjulian + int(fbt, 8)
         nxtjulS = nowjulian + int(fbs, 8)
         nxtjulPbt = nowjulian + int(fbpbt, 8)

         call jul2greg(bdysec, bdymin, bdyhour, bdyday, bdymon, bdyyear, nowjulian)

         call mpi_csp_io_filebdy_name(nowjulian, bdy_ts, fbvec, nfbvec, 'u', bvecavg, filebdy, .false.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy u for date: ', bdyyear, '-', bdymon, '-', bdyday, ' ', bdyhour, ':', bdymin, ':', bdysec, &
               ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         ! 6 stands for exchanging boundary data,dimension size is nlbdy and nk
         CALL mpi_netcdf_read_exchange(filebdy, 'u', bdy_u(:,:,1), 6, optional_ts = bdy_ts)

         call mpi_csp_io_filebdy_name(nowjulian, bdy_ts, fbvec, nfbvec, 'v', bvecavg, filebdy, .false.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy v for date: ', bdyyear, '-', bdymon, '-', bdyday, ' ', bdyhour, ':', bdymin, ':', bdysec, &
               ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         CALL mpi_netcdf_read_exchange(filebdy, 'v', bdy_v(:,:,1), 6, optional_ts = bdy_ts)

         call mpi_csp_io_filebdy_name(nowjulian, bdy_ts, fbvec, nfbvec, 'w', bvecavg, filebdy, .false.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy w for date: ', bdyyear, '-', bdymon, '-', bdyday, ' ', bdyhour, ':', bdymin, ':', bdysec, &
               ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         CALL mpi_netcdf_read_exchange(filebdy, 'w', bdy_w(:,:,1), 6, optional_ts = bdy_ts)
         readnxtVEC = .true.

         call mpi_csp_io_filebdy_name(nowjulian, bdy_ts, fbt, nfbt, 't', btavg, filebdy, .false.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy t for date: ', bdyyear, '-', bdymon, '-', bdyday, ' ', bdyhour, ':', bdymin, ':', bdysec, &
               ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         CALL mpi_netcdf_read_exchange(filebdy, 't', bdy_t(:,:,1), 6, optional_ts = bdy_ts)
         readnxtT = .true.

         call mpi_csp_io_filebdy_name(nowjulian, bdy_ts, fbs, nfbs, 's', bsavg, filebdy, .false.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy s for date: ', bdyyear, '-', bdymon, '-', bdyday, ' ', bdyhour, ':', bdymin, ':', bdysec, &
               ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         CALL mpi_netcdf_read_exchange(filebdy, 's', bdy_s(:,:,1), 6, optional_ts = bdy_ts)
         readnxtS = .true.

         call mpi_csp_io_filebdy_name(nowjulian, bdy_ts, fbpbt, nfbpbt, pbtname, bpbtavg, filebdy, .false.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy pbt for date: ', bdyyear, '-', bdymon, '-', bdyday, ' ', bdyhour, ':', bdymin, ':', bdysec, &
               ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         ! 5 stands for exchanging boundary data,dimension size is nlbdy
         CALL mpi_netcdf_read_exchange(filebdy, pbtname, bdy_pbt(:,1), 5, optional_ts = bdy_ts)
         readnxtPbt = .true.

         !$ACC update device(bdy_u(:,:,1), bdy_v(:,:,1), bdy_w(:,:,1), &
         !$ACC               bdy_t(:,:,1), bdy_s(:,:,1), bdy_pbt(:,1))
      end if

      ! process 0 may not have boundary
      if(ln_bdy) then
         if (nowjulian .ge. nxtjulVEC) then
            if (mpi_check_bdy) then
               !$ACC kernels present(bdy_u, bdy_v, bdy_w)
               !$ACC loop collapse(2) independent
               do k = 1, nk
                  do i = 1, nlbdy
                     bdy_u(i,k,1) = bdy_u(i,k,2)
                     bdy_v(i,k,1) = bdy_v(i,k,2)
                     bdy_w(i,k,1) = bdy_w(i,k,2)
                  end do
               end do
               !$ACC end kernels
            end if
            readnxtVEC = .true.
            nxtjulVEC = nxtjulVEC + int(fbvec, 8)
         end if

         if (nowjulian .ge. nxtjulT) then
            if (mpi_check_bdy) then
               !$ACC kernels present(bdy_t)
               !$ACC loop collapse(2) independent
               do k = 1, nk
                  do i = 1, nlbdy
                     bdy_t(i,k,1) = bdy_t(i,k,2)
                  end do
               end do
               !$ACC end kernels
            end if
            readnxtT = .true.
            nxtjulT = nxtjulT + int(fbt, 8)
         end if

         if (nowjulian .ge. nxtjulS) then
            if (mpi_check_bdy) then
               !$ACC kernels present(bdy_s)
               !$ACC loop collapse(2) independent
               do k = 1, nk
                  do i = 1, nlbdy
                     bdy_s(i,k,1) = bdy_s(i,k,2)
                  end do
               end do
               !$ACC end kernels
            end if
            readnxtS = .true.
            nxtjulS = nxtjulS + int(fbs, 8)
         end if

         if (nowjulian .ge. nxtjulPbt) then
            if (mpi_check_bdy) then
               !$ACC kernels present(bdy_pbt)
               !$ACC loop independent
               do i = 1, nlbdy
                  bdy_pbt(i,1) = bdy_pbt(i,2)
               end do
               !$ACC end kernels
            end if
            readnxtPbt = .true.
            nxtjulPbt = nxtjulPbt + int(fbpbt, 8)
         end if
      end if

      if (readnxtVEC) then
         call jul2greg(bdysec, bdymin, bdyhour, bdyday, bdymon, bdyyear, nxtjulVEC)
         call mpi_csp_io_filebdy_name(nxtjulVEC, bdy_ts, fbvec, nfbvec, 'u', bvecavg, filebdy, .true.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy u for date: ', bdyyear, '-', bdymon, '-', bdyday, ' ', bdyhour, ':', bdymin, ':', bdysec, &
               ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         CALL mpi_netcdf_read_exchange(filebdy, 'u', bdy_u(:,:,2), 6, optional_ts = bdy_ts)

         call mpi_csp_io_filebdy_name(nxtjulVEC, bdy_ts, fbvec, nfbvec, 'v', bvecavg, filebdy, .true.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy v for date: ', bdyyear, '-', bdymon, '-', bdyday, ' ', bdyhour, ':', bdymin, ':', bdysec, &
               ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         CALL mpi_netcdf_read_exchange(filebdy, 'v', bdy_v(:,:,2), 6, optional_ts = bdy_ts)

         call mpi_csp_io_filebdy_name(nxtjulVEC, bdy_ts, fbvec, nfbvec, 'w', bvecavg, filebdy, .true.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy w for date: ', bdyyear, '-', bdymon, '-', bdyday, ' ', bdyhour, ':', bdymin, ':', bdysec, &
               ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         CALL mpi_netcdf_read_exchange(filebdy, 'w', bdy_w(:,:,2), 6, optional_ts = bdy_ts)
         if(ln_bdy .and. mpi_check_bdy) then
            fbvecr = real(fbvec, wp)
            !$ACC update device(bdy_u(:,:,2), bdy_v(:,:,2), bdy_w(:,:,2))
            !$ACC kernels present(bdy_u, bdy_v, bdy_w), copyin(fbvecr)
            !$ACC loop collapse(2) independent
            do k = 1, nk
               do i = 1, nlbdy
                  bdy_u(i,k,3) = (bdy_u(i,k,2) - bdy_u(i,k,1))/ fbvecr
                  bdy_v(i,k,3) = (bdy_v(i,k,2) - bdy_v(i,k,1))/ fbvecr
                  bdy_w(i,k,3) = (bdy_w(i,k,2) - bdy_w(i,k,1))/ fbvecr
               end do
            end do
            !$ACC end kernels
         end if
      end if

      if (readnxtT) then
         call jul2greg(bdysec, bdymin, bdyhour, bdyday, bdymon, bdyyear, nxtjulT)
         call mpi_csp_io_filebdy_name(nxtjulT, bdy_ts, fbt, nfbt, 't', btavg, filebdy, .true.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy t for date: ', bdyyear, '-', bdymon, '-', bdyday, &
               ' ', bdyhour, ':', bdymin, ':', bdysec, ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         CALL mpi_netcdf_read_exchange(filebdy, 't', bdy_t(:,:,2), 6, optional_ts = bdy_ts)
         if(ln_bdy .and. mpi_check_bdy) then
            fbtr = real(fbt, wp)
            !$ACC update device(bdy_t(:,:,2))
            !$ACC kernels present(bdy_t), copyin(fbtr)
            !$ACC loop collapse(2) independent
            do k = 1, nk
               do i = 1, nlbdy
                  bdy_t(i,k,3) = (bdy_t(i,k,2) - bdy_t(i,k,1))/ fbtr
               end do
            end do
            !$ACC end kernels
         end if
      end if

      if (readnxtS) then
         call jul2greg(bdysec, bdymin, bdyhour, bdyday, bdymon, bdyyear, nxtjulS)
         call mpi_csp_io_filebdy_name(nxtjulS, bdy_ts, fbs, nfbs, 's', bsavg, filebdy, .true.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy s for date: ', bdyyear, '-', bdymon, '-', bdyday, &
               ' ', bdyhour, ':', bdymin, ':', bdysec, ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         CALL mpi_netcdf_read_exchange(filebdy, 's', bdy_s(:,:,2), 6, optional_ts = bdy_ts)
         if(ln_bdy .and. mpi_check_bdy) then
            fbsr = real(fbs, wp)
            !$ACC update device(bdy_s(:,:,2))
            !$ACC kernels present(bdy_s), copyin(fbsr)
            !$ACC loop collapse(2) independent
            do k = 1, nk
               do i = 1, nlbdy
                  bdy_s(i,k,3) = (bdy_s(i,k,2) - bdy_s(i,k,1))/ fbsr
               end do
            end do
            !$ACC end kernels
         end if
      end if

      if (readnxtPbt) then
         call jul2greg(bdysec, bdymin, bdyhour, bdyday, bdymon, bdyyear, nxtjulPbt)
         call mpi_csp_io_filebdy_name(nxtjulPbt, bdy_ts, fbpbt, nfbpbt, pbtname, bpbtavg, filebdy, .true.)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read bdy pbt for date: ', bdyyear, '-', bdymon, '-', bdyday, &
               ' ', bdyhour, ':', bdymin, ':', bdysec, ' from file '// trim(filebdy)//' at ts = ', bdy_ts
         end if
         CALL mpi_netcdf_read_exchange(filebdy, pbtname, bdy_pbt(:,2), 5, optional_ts = bdy_ts)
         if(ln_bdy .and. mpi_check_bdy) then
            fbpbtr = real(fbpbt, wp)
            !$ACC update device(bdy_pbt(:,2))
            !$ACC kernels present(bdy_pbt), copyin(fbpbtr)
            !$ACC loop independent
            do i = 1, nlbdy
               bdy_pbt(i,3) = (bdy_pbt(i,2) - bdy_pbt(i,1))/ fbpbtr
            end do
            !$ACC end kernels
         end if
      end if

      if(ln_bdy .and. mpi_check_bdy) then
         dfjulvec = real(fbvec - int(nxtjulVec - nowjulian, 4) + int(dTtracer, 4), wp)
         dfjult = real(fbt - int(nxtjulT - nowjulian, 4) + int(dTtracer, 4), wp)
         dfjuls = real(fbs - int(nxtjulS - nowjulian, 4) + int(dTtracer, 4), wp)
         dfjulpbt = real(fbpbt - int(nxtjulPbt - nowjulian, 4) + int(dTtracer, 4), wp)

         !$ACC kernels present(bdy_u, bdy_v, bdy_w, bdy_t, bdy_s, bdy_pbt), &
         !$ACC          copyin(dfjulvec, dfjult, dfjuls, dfjulpbt)
         !$ACC loop collapse(2) independent
         do k = 1, nk
            do i = 1, nlbdy
               bdy_u(i,k,4) = bdy_u(i,k,1) + bdy_u(i,k,3)* dfjulvec
               bdy_v(i,k,4) = bdy_v(i,k,1) + bdy_v(i,k,3)* dfjulvec
               bdy_w(i,k,4) = bdy_w(i,k,1) + bdy_w(i,k,3)* dfjulvec
               bdy_t(i,k,4) = bdy_t(i,k,1) + bdy_t(i,k,3)* dfjult
               bdy_s(i,k,4) = bdy_s(i,k,1) + bdy_s(i,k,3)* dfjuls
            end do
         end do

         !$ACC loop
         do i = 1, nlbdy
            bdy_pbt(i,4) = bdy_pbt(i,1) + bdy_pbt(i,3)* dfjulpbt
         end do
         !$ACC end kernels
      end if

   end subroutine mpi_csp_io_rdbdy

!==============================================================================
   subroutine mpi_csp_io_filebdy_name(bdyjulian, bdyts, bdyff, bdynff, bdykind, bdyavg, bdyfilename, returnjulian)
!==============================================================================
      implicit none
      integer(i8), intent(inout) :: bdyjulian
      integer, intent(inout) :: bdyts
      integer, intent(in) :: bdyff
      character(*), intent(in) :: bdynff, bdykind
      character(*), intent(inout) :: bdyfilename
      logical, intent(in) :: bdyavg
      character(lc) :: fileprefix, domid_str, bdyyear_str, bdymon_str, bdyday_str
      integer :: bdyyear, bdymon, bdyday, bdyhour, bdymin, bdysec
      integer(i8) :: julian_st, julian_file, julian_en
      integer :: djulian, ff_cycle_last
      logical, intent(in) :: returnjulian

      write (domid_str, "(i2.2)") idom
      fileprefix = trim(fbdy_dir)//'nmefc_macom_reg_bdy_'

      ! align bdyjulian to file record time
      call jul2greg(bdysec, bdymin, bdyhour, bdyday, bdymon, bdyyear, bdyjulian)
      select case (trim(bdynff))
      case ('yearly')
         call greg2jul(0, 0, 0, 1, 1, bdyyear, julian_st)
         call greg2jul(0, 0, 0, 1, 1, bdyyear-1, julian_en)
         ff_cycle_last = (julian_st - julian_en)/bdyff
      case ('monthly')
         call greg2jul(0, 0, 0, 1, bdymon, bdyyear, julian_st)
         bdymon = bdymon - 1
         if (bdymon .eq. 0) then
            bdyyear = bdyyear - 1
            bdymon = 12
         end if
         call greg2jul(0, 0, 0, 1, bdymon, bdyyear, julian_en)
         ff_cycle_last = (julian_st - julian_en)/bdyff
      case ('daily')
         call greg2jul(0, 0, 0, bdyday, bdymon, bdyyear, julian_st)
         ff_cycle_last = 86400/bdyff
      case default
         write (run_out_unit, *) 'ERROR! : bdy file frequency kind not exist, bdynff for ', &
            trim(bdykind), ' is ', trim(bdynff)
         write (run_out_unit, *) 'model only support daily, monthly, yearly now, check 1 '
         stop
      end select

      if (bdyavg) then
         julian_st = julian_st + int(bdyff/2, 8)
      end if
      djulian = bdyjulian - julian_st
      if (djulian .ge. 0) then
         bdyts = int(djulian/bdyff) + 1
         julian_file = julian_st + (bdyts - 1)*int(bdyff, 8)
      else
         julian_file = julian_st - int(bdyff, 8)
         bdyts = ff_cycle_last
      end if
      if (returnjulian) bdyjulian = julian_file

      call jul2greg(bdysec, bdymin, bdyhour, bdyday, bdymon, bdyyear, julian_file)
      write (bdyyear_str, "(i4.4)") bdyyear
      write (bdymon_str, "(i2.2)") bdymon
      write (bdyday_str, "(i2.2)") bdyday

      select case (trim(bdynff))
      case ('yearly')
         bdyfilename = trim(fileprefix)//trim(bdykind)//'_d'//trim(domid_str)// &
                       '_Y'//trim(bdyyear_str)//'.nc'
      case ('monthly')
         bdyfilename = trim(fileprefix)//trim(bdykind)//'_d'//trim(domid_str)// &
                       '_Y'//trim(bdyyear_str)//'M'//trim(bdymon_str)//'.nc'
      case ('daily')
         bdyfilename = trim(fileprefix)//trim(bdykind)//'_d'//trim(domid_str)// &
                       '_Y'//trim(bdyyear_str)//'M'//trim(bdymon_str)//'D'//trim(bdyday_str)//'.nc'
      case default
         write (run_out_unit, *) 'ERROR! : bdy file frequency kind not exist, bdynff for ', &
            trim(bdykind), ' is ', trim(bdynff)
         write (run_out_unit, *) 'model only support daily, monthly, yearly now, check 2 '
         stop
      end select

   end subroutine mpi_csp_io_filebdy_name

end module mod_mpi_csp_io_bdy
