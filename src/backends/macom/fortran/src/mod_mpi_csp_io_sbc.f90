module mod_mpi_csp_io_sbc
   use mod_misc
   use mod_io_netcdf
   use mod_mpi_variables
   use mod_mpi_interfaces
   implicit none

contains

!==============================================================================
   subroutine mpi_csp_io_rdforce
!==============================================================================
      implicit none
      integer :: sbc_ts, i
      integer :: sbcyear, sbcmon, sbcday, sbchour, sbcmin, sbcsec
      integer(i8) :: nowjuluv, nowjultq, nowjulrad, nowjulprc, nowjulslp, &
                     nowjulrunoff, nowjulsfrs, nowjulsfrt
      character(lc) :: filesbc
      real(wp), dimension(8) :: dfjul
      logical :: readnxtuv, readnxttq, readnxtrad, readnxtprc, readnxtslp, &
                 readnxtrunoff, readnxtsfrs, readnxtsfrt

      readnxtuv = .false.
      readnxttq = .false.
      readnxtrad = .false.
      readnxtprc = .false.
      readnxtslp = .false.
      readnxtrunoff = .false.
      readnxtsfrs = .false.
      readnxtsfrt = .false.

      if ((myIter .eq. 1) .and. (.not. (restart_in .and. date_res))) then
         if (sbc_cli) then
            nowjuluv = nowjulian - int(dTtracer, 8)
            if (fuvavg) then
               nxtjuluv = nowjulian - int(dTtracer, 8) + int(ffuv / 2, 8)
            else
               nxtjuluv = nowjulian - int(dTtracer, 8) + int(ffuv, 8)
            end if

            nowjultq = nowjulian - int(dTtracer, 8)
            if (ftqavg) then
               nxtjultq = nowjulian - int(dTtracer, 8) + int(fftq / 2, 8)
            else
               nxtjultq = nowjulian - int(dTtracer, 8) + int(fftq, 8)
            end if

            nowjulrad = nowjulian - int(dTtracer, 8)
            if (fradavg) then
               nxtjulrad = nowjulian - int(dTtracer, 8) + int(ffrad / 2, 8)
            else
               nxtjulrad = nowjulian - int(dTtracer, 8) + int(ffrad, 8)
            end if

            nowjulprc = nowjulian - int(dTtracer, 8)
            if (fprcavg) then
               nxtjulprc = nowjulian - int(dTtracer, 8) + int(ffprc / 2, 8)
            else
               nxtjulprc = nowjulian - int(dTtracer, 8) + int(ffprc, 8)
            end if

            nowjulslp = nowjulian - int(dTtracer, 8)
            if (fslpavg) then
               nxtjulslp = nowjulian - int(dTtracer, 8) + int(ffslp / 2, 8)
            else
               nxtjulslp = nowjulian - int(dTtracer, 8) + int(ffslp, 8)
            end if

            nowjulrunoff = nowjulian - int(dTtracer, 8)
            if (frunoffavg) then
               nxtjulrunoff = nowjulian - int(dTtracer, 8) + int(ffrunoff / 2, 8)
            else
               nxtjulrunoff = nowjulian - int(dTtracer, 8) + int(ffrunoff, 8)
            end if

            nowjulsfrs = nowjulian - int(dTtracer, 8)
            if (fsfrsavg) then
               nxtjulsfrs = nowjulian - int(dTtracer, 8) + int(ffsfrs / 2, 8)
            else
               nxtjulsfrs = nowjulian - int(dTtracer, 8) + int(ffsfrs, 8)
            end if

            nowjulsfrt = nowjulian - int(dTtracer, 8)
            if (fsfrtavg) then
               nxtjulsfrt = nowjulian - int(dTtracer, 8) + int(ffsfrt / 2, 8)
            else
               nxtjulsfrt = nowjulian - int(dTtracer, 8) + int(ffsfrt, 8)
            end if
         else
            if (fuvavg) then
               nowjuluv = nowjulian - int(dTtracer, 8) - int(ffuv / 2, 8)
               nxtjuluv = nowjulian - int(dTtracer, 8) + int(ffuv / 2, 8)
            else
               nowjuluv = nowjulian - int(dTtracer, 8)
               nxtjuluv = nowjulian - int(dTtracer, 8) + int(ffuv, 8)
            end if

            if (ftqavg) then
               nowjultq = nowjulian - int(dTtracer, 8) - int(fftq / 2, 8)
               nxtjultq = nowjulian - int(dTtracer, 8) + int(fftq / 2, 8)
            else
               nowjultq = nowjulian - int(dTtracer, 8)
               nxtjultq = nowjulian - int(dTtracer, 8) + int(fftq, 8)
            end if

            if (fradavg) then
               nowjulrad = nowjulian - int(dTtracer, 8) - int(ffrad / 2, 8)
               nxtjulrad = nowjulian - int(dTtracer, 8) + int(ffrad / 2, 8)
            else
               nowjulrad = nowjulian - int(dTtracer, 8)
               nxtjulrad = nowjulian - int(dTtracer, 8) + int(ffrad, 8)
            end if

            if (fprcavg) then
               nowjulprc = nowjulian - int(dTtracer, 8) - int(ffprc / 2, 8)
               nxtjulprc = nowjulian - int(dTtracer, 8) + int(ffprc / 2, 8)
            else
               nowjulprc = nowjulian - int(dTtracer, 8)
               nxtjulprc = nowjulian - int(dTtracer, 8) + int(ffprc, 8)
            end if

            if (fslpavg) then
               nowjulslp = nowjulian - int(dTtracer, 8) - int(ffslp / 2, 8)
               nxtjulslp = nowjulian - int(dTtracer, 8) + int(ffslp / 2, 8)
            else
               nowjulslp = nowjulian - int(dTtracer, 8)
               nxtjulslp = nowjulian - int(dTtracer, 8) + int(ffslp, 8)
            end if

            if (frunoffavg) then
               nowjulrunoff = nowjulian - int(dTtracer, 8) - int(ffrunoff / 2, 8)
               nxtjulrunoff = nowjulian - int(dTtracer, 8) + int(ffrunoff / 2, 8)
            else
               nowjulrunoff = nowjulian - int(dTtracer, 8)
               nxtjulrunoff = nowjulian - int(dTtracer, 8) + int(ffrunoff, 8)
            end if

            if (fsfrsavg) then
               nowjulsfrs = nowjulian - int(dTtracer, 8) - int(ffsfrs / 2, 8)
               nxtjulsfrs = nowjulian - int(dTtracer, 8) + int(ffsfrs / 2, 8)
            else
               nowjulsfrs = nowjulian - int(dTtracer, 8)
               nxtjulsfrs = nowjulian - int(dTtracer, 8) + int(ffsfrs, 8)
            end if

            if (fsfrtavg) then
               nowjulsfrt = nowjulian - int(dTtracer, 8) - int(ffsfrt / 2, 8)
               nxtjulsfrt = nowjulian - int(dTtracer, 8) + int(ffsfrt / 2, 8)
            else
               nowjulsfrt = nowjulian - int(dTtracer, 8)
               nxtjulsfrt = nowjulian - int(dTtracer, 8) + int(ffsfrt, 8)
            end if
         end if

         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nowjuluv)
         call mpi_csp_io_filesbc_name(nowjuluv, sbc_ts, ffuv, nffuv, 'sbcu', fuvavg, filesbc)

# if !defined (NOWIND)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc u for date: ', sbcyear, '-', sbcmon, '-', sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
               ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'sbcu', u10(1:nlpb, 1), 1, optional_ts = sbc_ts)

         call mpi_csp_io_filesbc_name(nowjuluv, sbc_ts, ffuv, nffuv, 'sbcv', fuvavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc v for date: ', sbcyear, '-', sbcmon, '-', sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
               ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'sbcv', v10(1:nlpb, 1), 1, optional_ts = sbc_ts)
         readnxtuv = .true.
# endif

# if !defined (BARMOD)
         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nowjultq)
         call mpi_csp_io_filesbc_name(nowjultq, sbc_ts, fftq, nfftq, 'sbct', ftqavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc t for date: ', sbcyear, '-', sbcmon, '-', sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
               ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'sbct', t10(1:nlpb, 1), 1,optional_ts = sbc_ts)

         call mpi_csp_io_filesbc_name(nowjultq, sbc_ts, fftq, nfftq, 'sbcq', ftqavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc q for date: ', sbcyear, '-', sbcmon, '-', sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
               ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'sbcq', q10(1:nlpb, 1), 1,optional_ts = sbc_ts)
         readnxttq = .true.

         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nowjulrad)
         call mpi_csp_io_filesbc_name(nowjulrad, sbc_ts, ffrad, nffrad, 'swdn', fradavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc swdn for date: ', sbcyear, '-', sbcmon, '-', sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
               ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'swdn', swdn(1:nlpb, 1), 1,optional_ts = sbc_ts)

         call mpi_csp_io_filesbc_name(nowjulrad, sbc_ts, ffrad, nffrad, 'lwdn', fradavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc lwdn for date: ', sbcyear, '-', sbcmon, '-', sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
               ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'lwdn', lwdn(1:nlpb, 1), 1,optional_ts = sbc_ts)
         readnxtrad = .true.

         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nowjulprc)
         call mpi_csp_io_filesbc_name(nowjulprc, sbc_ts, ffprc, nffprc, 'prec', fprcavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc prec for date: ', sbcyear, '-', sbcmon, '-', sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
               ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'prec', prec(1:nlpb, 1), 1,optional_ts = sbc_ts)

         call mpi_csp_io_filesbc_name(nowjulprc, sbc_ts, ffprc, nffprc, 'snow', fprcavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc snow for date: ', sbcyear, '-', sbcmon, '-', sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
               ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'snow', snow(1:nlpb, 1), 1,optional_ts = sbc_ts)
         readnxtprc = .true.

         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nowjulrunoff)
         call mpi_csp_io_filesbc_name(nowjulrunoff, sbc_ts, ffrunoff, nffrunoff, 'runoff', frunoffavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc runoff for date: ', sbcyear, '-', sbcmon, '-', sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
               ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'runoff', runoff(1:nlpb, 1), 1,optional_ts = sbc_ts)
         readnxtrunoff = .true.

         if (gammas .gt. 0.0_wp) then
            call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nowjulsfrs)
            call mpi_csp_io_filesbc_name(nowjulsfrs, sbc_ts, ffsfrs, nffsfrs, 'sfrs', &
                                         fsfrsavg, filesbc)
            if (mpi_rank == 0) then
               write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
                  ' read surface restore s for date: ', sbcyear, '-', sbcmon, '-', &
                  sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
                  ' from file '// trim(filesbc)//' at ts = ', sbc_ts
            end if
            CALL mpi_netcdf_read_exchange(filesbc, 'sfrs', sfrs(1:nlpb, 1), 1,optional_ts = sbc_ts)
            readnxtsfrs = .true.
         end if

         if (gammat .gt. 0.0_wp) then
            call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nowjulsfrt)
            call mpi_csp_io_filesbc_name(nowjulsfrt, sbc_ts, ffsfrt, nffsfrt, 'sfrt', &
                                         fsfrtavg, filesbc)
            if (mpi_rank == 0) then
               write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
                  ' read surface restore t for date: ', sbcyear, '-', sbcmon, '-', &
                  sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
                  ' from file '// trim(filesbc)//' at ts = ', sbc_ts
            end if
            CALL mpi_netcdf_read_exchange(filesbc, 'sfrt', sfrt(1:nlpb, 1), 1,optional_ts = sbc_ts)
            readnxtsfrt = .true.
         end if
# endif

# if !defined (NOSLP)
         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nowjulslp)
         call mpi_csp_io_filesbc_name(nowjulslp, sbc_ts, ffslp, nffslp, 'slp', fslpavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc slp for date: ', sbcyear, '-', sbcmon, '-', sbcday, ' ', sbchour, ':', sbcmin, ':', sbcsec, &
               ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'slp', slp(1:nlpb, 1), 1,optional_ts = sbc_ts)
         readnxtslp = .true.
# endif

         ! ---- OPENACC: update force field at 1st time step
         !$ACC update device(u10(:,1),v10(:,1),t10(:,1),q10(:,1),swdn(:,1),  &
         !$ACC        lwdn(:,1),prec(:,1),snow(:,1),slp(:,1),runoff(:,1),sfrs(:,1),sfrt(:,1))
         ! ----
      end if

      !YY,ZY: SNOW and PREC can not be A NEGATIVE VALUE, RIGHT?? TO RULE OUT ANY UNEXPECTED VALUES:
      !$acc kernels loop present(prec,snow)
      do i = 1, nlpb
         if(prec(i,1) .lt. 0.0_wp) prec(i,1) = 0.0_wp
         if(snow(i,1) .lt. 0.0_wp) snow(i,1) = 0.0_wp
      end do
      !$acc end kernels
      
      if ((nowjulian - int(dTtracer, 8)) .ge. nxtjuluv) then
         nxtjuluv = nxtjuluv + int(ffuv, 8)
         !$acc kernels loop present(u10,v10)
         do i = 1, nlpb
            u10(i, 1) = u10(i, 2)
            v10(i, 1) = v10(i, 2)
         end do
         !$acc end kernels
         readnxtuv = .true.
         if (nxtjuluv .ge. endjulian) then
            if (sbc_cli) then
               nxtjuluv = endjulian
               readnxtuv = .false.
               !$acc kernels loop present(u10,v10)
               do i = 1,nlpb
                  u10(i, 3) = 0.0_wp
                  v10(i, 3) = 0.0_wp
               end do
               !$acc end kernels
            end if
         end if
      end if

      if ((nowjulian - int(dTtracer, 8)) .ge. nxtjultq) then
         nxtjultq = nxtjultq + int(fftq, 8)
         !$acc kernels loop present(t10,q10)
         do i = 1,nlpb
            t10(i, 1) = t10(i, 2)
            q10(i, 1) = q10(i, 2)
         end do
         !$acc end kernels
         readnxttq = .true.
         if (nxtjultq .ge. endjulian) then
            if (sbc_cli) then
               nxtjultq = endjulian
               readnxttq = .false.
               !$acc kernels loop present(t10,q10)
               do i = 1,nlpb
                  t10(i, 3) = 0.0_wp
                  q10(i, 3) = 0.0_wp
               end do
               !$acc end kernels
            end if
         end if
      end if

      if ((nowjulian - int(dTtracer, 8)) .ge. nxtjulrad) then
         nxtjulrad = nxtjulrad + int(ffrad, 8)
         !$acc kernels loop present(swdn,lwdn)
         do i = 1,nlpb
            swdn(i, 1) = swdn(i, 2)
            lwdn(i, 1) = lwdn(i, 2)
         end do
         !$acc end kernels
         readnxtrad = .true.
         if (nxtjulrad .ge. endjulian) then
            if (sbc_cli) then
               nxtjulrad = endjulian
               readnxtrad = .false.
               !$acc kernels loop present(swdn,lwdn)
               do i = 1,nlpb
                  swdn(i, 3) = 0.0_wp
                  lwdn(i, 3) = 0.0_wp
               end do
               !$acc end kernels
            end if
         end if
      end if

      if ((nowjulian - int(dTtracer, 8)) .ge. nxtjulprc) then
         nxtjulprc = nxtjulprc + int(ffprc, 8)
         !$acc kernels loop present(prec,snow)
         do i = 1,nlpb
            prec(i, 1) = prec(i, 2)
            snow(i, 1) = snow(i, 2)
         end do
         !$acc end kernels
         readnxtprc = .true.
         if (nxtjulprc .ge. endjulian) then
            if (sbc_cli) then
               nxtjulprc = endjulian
               readnxtprc = .false.
               !$acc kernels loop present(prec,snow)
               do i = 1,nlpb
                  prec(i, 3) = 0.0_wp
                  snow(i, 3) = 0.0_wp
               end do
               !$acc end kernels
            end if
         end if
      end if

      if ((nowjulian - int(dTtracer, 8)) .ge. nxtjulslp) then
         nxtjulslp = nxtjulslp + int(ffslp, 8)
         !$acc kernels loop present(slp)
         do i = 1,nlpb
            slp(i, 1) = slp(i, 2)
         end do
         !$acc end kernels
         readnxtslp = .true.
         if (nxtjulslp .ge. endjulian) then
            if (sbc_cli) then
               nxtjulslp = endjulian
               readnxtslp = .false.
               !$acc kernels loop present(slp)
               do i = 1,nlpb
                  slp(i, 3) = 0.0_wp
               end do
               !$acc end kernels
            end if
         end if
      end if

      if ((nowjulian - int(dTtracer, 8)) .ge. nxtjulrunoff) then
         nxtjulrunoff = nxtjulrunoff + int(ffrunoff, 8)
         !$acc kernels loop present(runoff)
         do i = 1,nlpb
            runoff(i, 1) = runoff(i, 2)
         end do
         !$acc end kernels
         readnxtrunoff = .true.
         if (nxtjulrunoff .ge. endjulian) then
            if (sbc_cli) then
               nxtjulrunoff = endjulian
               readnxtrunoff = .false.
               !$acc kernels loop present(runoff)
               do i = 1,nlpb
                  runoff(i, 3) = 0.0_wp
               end do
               !$acc end kernels
            end if
         end if
      end if

      if (((nowjulian - int(dTtracer, 8)) .ge. nxtjulsfrs) .and. &
          (gammas .gt. 0.0_wp)) then
         nxtjulsfrs = nxtjulsfrs + int(ffsfrs, 8)
         !$acc kernels loop present(sfrs)
         do i = 1,nlpb
            sfrs(i, 1) = sfrs(i, 2)
         end do
         !$acc end kernels
         readnxtsfrs = .true.
         if (nxtjulsfrs .ge. endjulian) then
            if (sbc_cli) then
               nxtjulsfrs = endjulian
               readnxtsfrs = .false.
               !$acc kernels loop present(sfrs)
               do i = 1,nlpb
                  sfrs(i, 3) = 0.0_wp
               end do
               !$acc end kernels
            end if
         end if
      end if

      if (((nowjulian - int(dTtracer, 8)) .ge. nxtjulsfrt) .and. &
          (gammat .gt. 0.0_wp)) then
         nxtjulsfrt = nxtjulsfrt + int(ffsfrt, 8)
         !$acc kernels loop present(sfrt)
         do i = 1,nlpb
            sfrt(i, 1) = sfrt(i, 2)
         end do
         !$acc end kernels
         readnxtsfrt = .true.
         if (nxtjulsfrt .ge. endjulian) then
            if (sbc_cli) then
               nxtjulsfrt = endjulian
               readnxtsfrt = .false.
               !$acc kernels loop present(sfrt)
               do i = 1,nlpb
                  sfrt(i, 3) = 0.0_wp
               end do
               !$acc end kernels
            end if
         end if
      end if

# if !defined (NOWIND)
      if (readnxtuv) then
         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nxtjuluv)
         call mpi_csp_io_filesbc_name(nxtjuluv, sbc_ts, ffuv, nffuv, 'sbcu', fuvavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc u for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'sbcu', u10(1:nlpb, 2), 1,optional_ts = sbc_ts)

         call mpi_csp_io_filesbc_name(nxtjuluv, sbc_ts, ffuv, nffuv, 'sbcv', fuvavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc v for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'sbcv', v10(1:nlpb, 2), 1,optional_ts = sbc_ts)
         ! ----  OPENACC: update next force field to GPU
         !$acc update device(u10(:,2), v10(:,2))
         ! ----
         !$acc kernels loop present(u10,v10,ffuv)
         do i = 1, nlpb
            u10(i, 3) = (u10(i, 2) - u10(i, 1))/ real(ffuv, wp)
            v10(i, 3) = (v10(i, 2) - v10(i, 1))/ real(ffuv, wp)
         end do
         !$acc end kernels
      end if
# endif 

# if !defined (BARMOD)
      if (readnxttq) then
         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nxtjultq)
         call mpi_csp_io_filesbc_name(nxtjultq, sbc_ts, fftq, nfftq, 'sbct', ftqavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc t for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'sbct', t10(1:nlpb, 2), 1,optional_ts = sbc_ts)

         call mpi_csp_io_filesbc_name(nxtjultq, sbc_ts, fftq, nfftq, 'sbcq', ftqavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc q for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'sbcq', q10(1:nlpb, 2), 1,optional_ts = sbc_ts)
         ! ----  OPENACC: update next force field to GPU
         !$acc update device(t10(:,2), q10(:,2))
         !$acc kernels loop present(t10,q10,fftq)
         do i = 1, nlpb
            t10(i, 3) = (t10(i, 2) - t10(i, 1))/ real(fftq, wp)
            q10(i, 3) = (q10(i, 2) - q10(i, 1))/ real(fftq, wp)
         end do
         !$acc end kernels
      end if

      if (readnxtrad) then
         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nxtjulrad)
         call mpi_csp_io_filesbc_name(nxtjulrad, sbc_ts, ffrad, nffrad, 'swdn', fradavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc swdn for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'swdn', swdn(1:nlpb, 2), 1,optional_ts = sbc_ts)

         call mpi_csp_io_filesbc_name(nxtjulrad, sbc_ts, ffrad, nffrad, 'lwdn', fradavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc lwdn for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'lwdn', lwdn(1:nlpb, 2), 1,optional_ts = sbc_ts)
         ! ----  OPENACC: update next force field to GPU
         !$acc update device(swdn(:,2), lwdn(:,2))
         !$acc kernels loop present(swdn,lwdn,ffrad)
         do i = 1, nlpb
            swdn(i, 3) = (swdn(i, 2) - swdn(i, 1))/ real(ffrad, wp)
            lwdn(i, 3) = (lwdn(i, 2) - lwdn(i, 1))/ real(ffrad, wp)
         end do
         !$acc end kernels
      end if

      if (readnxtprc) then
         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nxtjulprc)
         call mpi_csp_io_filesbc_name(nxtjulprc, sbc_ts, ffprc, nffprc, 'prec', fprcavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc prec for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'prec', prec(1:nlpb, 2), 1,optional_ts = sbc_ts)

         call mpi_csp_io_filesbc_name(nxtjulprc, sbc_ts, ffprc, nffprc, 'snow', fprcavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc snow for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'snow', snow(1:nlpb, 2), 1,optional_ts = sbc_ts)
         ! ----  OPENACC: update next force field to GPU
         !$acc update device(prec(:,2),snow(:,2))

         !YY,ZY: SNOW and PREC can not be A NEGATIVE VALUE, RIGHT?? TO RULE OUT ANY UNEXPECTED VALUES:
         !$acc kernels loop present(prec,snow)
         do i = 1, nlpb
            if(prec(i,2) .lt. 0.0_wp) prec(i,2) = 0.0_wp
            if(snow(i,2) .lt. 0.0_wp) snow(i,2) = 0.0_wp
         end do
         !$acc end kernels

         !$acc kernels loop present(prec,snow,ffprc)
         do i = 1, nlpb
            prec(i, 3) = (prec(i, 2) - prec(i, 1))/ real(ffprc, wp)
            snow(i, 3) = (snow(i, 2) - snow(i, 1))/ real(ffprc, wp)
         end do
         !$acc end kernels
      end if

      if (readnxtrunoff) then
         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nxtjulrunoff)
         call mpi_csp_io_filesbc_name(nxtjulrunoff, sbc_ts, ffrunoff, nffrunoff, 'runoff', frunoffavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc runoff for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'runoff', runoff(1:nlpb, 2), 1,optional_ts = sbc_ts)
         ! ----  OPENACC: update next force field to GPU
         !$acc update device(runoff(:,2))
         !$acc kernels loop present(runoff,ffrunoff)
         do i = 1, nlpb
            runoff(i, 3) = (runoff(i, 2) - runoff(i, 1))/ real(ffrunoff, wp)
         end do
         !$acc end kernels
      end if

      if (readnxtsfrs) then
         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nxtjulsfrs)
         call mpi_csp_io_filesbc_name(nxtjulsfrs, sbc_ts, ffsfrs, nffsfrs, 'sfrs', fsfrsavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read surface restore s for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'sfrs', sfrs(1:nlpb, 2), 1,optional_ts = sbc_ts)
         ! ----  OPENACC: update next force field to GPU
         !$acc update device(sfrs(:,2))
         !$acc kernels loop present(sfrs,ffsfrs)
         do i = 1, nlpb
            sfrs(i, 3) = (sfrs(i, 2) - sfrs(i, 1))/ real(ffsfrs, wp)
         end do
         !$acc end kernels
      end if

      if (readnxtsfrt) then
         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nxtjulsfrt)
         call mpi_csp_io_filesbc_name(nxtjulsfrt, sbc_ts, ffsfrt, nffsfrt, 'sfrt', fsfrtavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read surface restore t for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'sfrt', sfrt(1:nlpb, 2), 1,optional_ts = sbc_ts)
         ! ----  OPENACC: update next force field to GPU
         !$acc update device(sfrt(:,2))
         !$acc kernels loop present(sfrt,ffsfrt)
         do i = 1, nlpb
            sfrt(i, 3) = (sfrt(i, 2) - sfrt(i, 1))/ real(ffsfrt, wp)
         end do
         !$acc end kernels
      end if
# endif 

# if !defined (NOSLP)
      if (readnxtslp) then
         call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, nxtjulslp)
         call mpi_csp_io_filesbc_name(nxtjulslp, sbc_ts, ffslp, nffslp, 'slp', fslpavg, filesbc)
         if (mpi_rank == 0) then
            write (run_out_unit, "(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i5)") &
               ' read sbc slp for date: ', sbcyear, '-', sbcmon, '-', sbcday, &
               ' ', sbchour, ':', sbcmin, ':', sbcsec, ' from file '// trim(filesbc)//' at ts = ', sbc_ts
         end if
         CALL mpi_netcdf_read_exchange(filesbc, 'slp', slp(1:nlpb, 2), 1,optional_ts = sbc_ts)
         ! ----  OPENACC: update next force field to GPU
         !$acc update device(slp(:,2))
         !$acc kernels loop present(slp,ffslp)
         do i = 1, nlpb
            slp(i, 3) = (slp(i, 2) - slp(i, 1))/ real(ffslp, wp)
         end do
         !$acc end kernels
      end if
# endif 

      dfjul(1) = real(ffuv - int(nxtjuluv - nowjulian, 4) + int(dTtracer, 4), wp)
      dfjul(2) = real(fftq - int(nxtjultq - nowjulian, 4) + int(dTtracer, 4), wp)
      dfjul(3) = real(ffrad - int(nxtjulrad - nowjulian, 4) + int(dTtracer, 4), wp)
      dfjul(4) = real(ffprc - int(nxtjulprc - nowjulian, 4) + int(dTtracer, 4), wp)
      dfjul(5) = real(ffslp - int(nxtjulslp - nowjulian, 4) + int(dTtracer, 4), wp)
      dfjul(6) = real(ffrunoff - int(nxtjulrunoff - nowjulian, 4) + int(dTtracer, 4), wp)
      dfjul(7) = real(ffsfrs - int(nxtjulsfrs - nowjulian, 4) + int(dTtracer, 4), wp)
      dfjul(8) = real(ffsfrt - int(nxtjulsfrt - nowjulian, 4) + int(dTtracer, 4), wp)
      !$acc kernels present(nk,tFld,u10,v10,t10,q10,swdn,lwdn,    &
      !$acc             prec,snow,slp,runoff,sfrs,sfrt,phi0surf), &
      !$acc         copyin(dfjul)
      !$acc loop
      do i = 1, nlpb
         u10(i, 4) = u10(i, 1) + u10(i, 3)* dfjul(1)
         v10(i, 4) = v10(i, 1) + v10(i, 3)* dfjul(1)
         t10(i, 4) = t10(i, 1) + t10(i, 3)* dfjul(2)
         q10(i, 4) = q10(i, 1) + q10(i, 3)* dfjul(2)
         swdn(i, 4) = swdn(i, 1) + swdn(i, 3)* dfjul(3)
         lwdn(i, 4) = lwdn(i, 1) + lwdn(i, 3)* dfjul(3)
         prec(i, 4) = prec(i, 1) + prec(i, 3)* dfjul(4)
         snow(i, 4) = snow(i, 1) + snow(i, 3)* dfjul(4)
         runoff(i, 4) = runoff(i, 1) + runoff(i, 3)* dfjul(6)

         if (gammas .gt. 0.0_wp) then
            sfrs(i, 4) = sfrs(i, 1) + sfrs(i, 3)* dfjul(7)
         end if

         if (gammat .gt. 0.0_wp) then
            sfrt(i, 4) = sfrt(i, 1) + sfrt(i, 3)* dfjul(8)
         end if

# if !defined (NOSLP)
         slp(i, 4) = slp(i, 1) + slp(i, 3)* dfjul(5)
         phi0surf(i) = slp(i, 4)
# endif
      end do
      !$acc end kernels

   end subroutine mpi_csp_io_rdforce

!==============================================================================
   subroutine mpi_csp_io_filesbc_name(sbcjulian, sbcts, sbcff, sbcnff, sbckind, sbcavg, sbcfilename)
!==============================================================================
      implicit none
      integer(i8), intent(in) :: sbcjulian
      integer, intent(inout) :: sbcts
      integer, intent(in) :: sbcff
      character(*), intent(in) :: sbcnff, sbckind
      character(*), intent(inout) :: sbcfilename
      logical, intent(in) :: sbcavg
      character(lc) :: fileprefix, domid_str, sbcyear_str, sbcmon_str, sbcday_str
      integer :: sbcyear, sbcmon, sbcday, sbchour, sbcmin, sbcsec
      integer(i8) :: julian_st
      integer :: djulian

      write (domid_str, "(i2.2)") idom
      fileprefix = trim(fsbc_dir)//'nmefc_macom_csp_sbc_'

      call jul2greg(sbcsec, sbcmin, sbchour, sbcday, sbcmon, sbcyear, sbcjulian)

      write (sbcyear_str, "(i4.4)") sbcyear
      write (sbcmon_str, "(i2.2)") sbcmon
      write (sbcday_str, "(i2.2)") sbcday

      select case (trim(sbcnff))
      case ('yearly')
         sbcfilename = trim(fileprefix)//trim(sbckind)//'_d'//trim(domid_str)// &
                       '_Y'//trim(sbcyear_str)//'.nc'

         call greg2jul(0, 0, 0, 1, 1, sbcyear, julian_st)
         djulian = int(sbcjulian - julian_st, 4)
         sbcts = int(djulian/sbcff, 4)
         if (sbcavg) then
            sbcts = sbcts + 1
         else
            if (mod(djulian, sbcff) .eq. 0) then
               sbcts = sbcts + 1
            else
               write (run_out_unit, *) 'ERROR! : when sbc is instant field, sbc'// &
                  ' julian must divide exactly by sbc frequency'
               write (run_out_unit, *) 'this error find for '//trim(sbckind)
               call MPI_Abort(MPI_COMM_WORLD, 911, mpi_err)
            end if
         end if

      case ('monthly')
         sbcfilename = trim(fileprefix)//trim(sbckind)//'_d'//trim(domid_str)// &
                       '_Y'//trim(sbcyear_str)//'M'//trim(sbcmon_str)//'.nc'

         call greg2jul(0, 0, 0, 1, sbcmon, sbcyear, julian_st)
         djulian = int(sbcjulian - julian_st, 4)
         sbcts = int(djulian/sbcff, 4)
         if (sbcavg) then
            sbcts = sbcts + 1
         else
            if (mod(djulian, sbcff) .eq. 0) then
               sbcts = sbcts + 1
            else
               write (run_out_unit, *) 'ERROR! : when sbc is instant field, sbc'// &
                  ' julian must divide exactly by sbc frequency'
               write (run_out_unit, *) 'this error find for '//trim(sbckind)
               call MPI_Abort(MPI_COMM_WORLD, 911, mpi_err)
            end if
         end if

      case ('daily')
         sbcfilename = trim(fileprefix)//trim(sbckind)//'_d'//trim(domid_str)// &
                       '_Y'//trim(sbcyear_str)//'M'//trim(sbcmon_str)//'D'//trim(sbcday_str)//'.nc'

         call greg2jul(0, 0, 0, sbcday, sbcmon, sbcyear, julian_st)
         djulian = int(sbcjulian - julian_st, 4)
         sbcts = int(djulian/sbcff, 4)
         if (sbcavg) then
            sbcts = sbcts + 1
         else
            if (mod(djulian, sbcff) .eq. 0) then
               sbcts = sbcts + 1
            else
               write (run_out_unit, *) 'ERROR! : when sbc is instant field, sbc'// &
                  ' julian must divide exactly by sbc frequency'
               write (run_out_unit, *) 'this error find for '//trim(sbckind)
               call MPI_Abort(MPI_COMM_WORLD, 911, mpi_err)
            end if
         end if

      case default
         write (run_out_unit, *) 'ERROR! : sbc file frequency kind not exist, sbcnff for ', &
            trim(sbckind), ' is ', trim(sbcnff)
         write (run_out_unit, *) 'model only support daily, monthly, yearly now '
         call MPI_Abort(MPI_COMM_WORLD, 911, mpi_err)
      end select

      if (sbc_cyc) then
         sbcfilename = trim(fileprefix)//trim(sbckind)//'_d'//trim(domid_str)//'.nc'
      end if

   end subroutine mpi_csp_io_filesbc_name

end module mod_mpi_csp_io_sbc
