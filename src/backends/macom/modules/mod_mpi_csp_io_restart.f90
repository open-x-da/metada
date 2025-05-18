module mod_mpi_csp_io_restart
   use mod_misc
   use mod_io_netcdf
   use mod_csp_force, only : EmPmRrevise, EmPmRcount
   use mod_mpi_variables
   use mod_mpi_interfaces
   use mod_mpi_csp_io_misc
   !YY use seaice modules
   use mitice_parameters
   use mitice_vars
   use mitice_io
#ifdef SEAICE_ITD
   use mitice_itd
#endif
   implicit none

   real(dp), save, public         :: nowjulian_dp
   real(dp), save, public         :: nxtjultq_dp
   real(dp), save, public         :: nxtjuluv_dp
   real(dp), save, public         :: nxtjulrad_dp
   real(dp), save, public         :: nxtjulprc_dp
   real(dp), save, public         :: nxtjulslp_dp
   real(dp), save, public         :: nxtjulrunoff_dp
   real(dp), save, public         :: nxtjulsfrs_dp
   real(dp), save, public         :: nxtjulsfrt_dp

contains

!==============================================================================
   subroutine mpi_csp_io_restart_read
!==============================================================================
      implicit none
      integer :: i, k
      character(lc) :: filename, varname, domid_str

      FirstMomU = .false.
      FirstMomV = .false.
      FirstTracer = .false.

      write (domid_str, "(i2.2)") idom
      filename = "nmefc_macom_restart_d"//trim(domid_str)//".nc"

      IF (mpi_rank == 0) THEN

         if (date_res) then
            call netcdf_read(trim(filename), 'nowjulian', nowjulian)
            write (run_out_unit, *) "Last time it stop at second", nowjulian, " from 1949-10-01:00:00:00"
            call jul2greg(nowsec, nowmin, nowhour, nowday, nowmon, nowyear, nowjulian)

            print *, 'restart read, nowjulian: ', nowjulian
            !run_out_unit
            write (run_out_unit, "(a,i12,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)") &
               " Model date will continue from last stop : ", nowjulian, ' ', nowyear, '-', &
               nowmon, '-', nowday, ' ', nowhour, ':', nowmin, ':', nowsec
            call netcdf_read(trim(filename), 'nxtjultq', nxtjultq)
            call netcdf_read(trim(filename), 'nxtjuluv', nxtjuluv)
            call netcdf_read(trim(filename), 'nxtjulrad', nxtjulrad)
            call netcdf_read(trim(filename), 'nxtjulprc', nxtjulprc)
            call netcdf_read(trim(filename), 'nxtjulslp', nxtjulslp)
            call netcdf_read(trim(filename), 'nxtjulrunoff', nxtjulrunoff)
            call netcdf_read(trim(filename), 'nxtjulsfrs', nxtjulsfrs)
            call netcdf_read(trim(filename), 'nxtjulsfrt', nxtjulsfrt)

            nowjulian_dp = dble(nowjulian)
            nxtjultq_dp = dble(nxtjultq)
            nxtjuluv_dp = dble(nxtjuluv)
            nxtjulrad_dp = dble(nxtjulrad)
            nxtjulprc_dp = dble(nxtjulprc)
            nxtjulslp_dp = dble(nxtjulslp)
            nxtjulrunoff_dp = dble(nxtjulrunoff)
            nxtjulsfrs_dp = dble(nxtjulsfrs)
            nxtjulsfrt_dp = dble(nxtjulsfrt)

            CALL MPI_BCAST(nowjulian_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjultq_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjuluv_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulrad_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulprc_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulslp_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulrunoff_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulsfrs_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulsfrt_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
         end if

         if (ln_EmPmRevise) then
            call netcdf_read_scalar_i4(trim(filename), 'EmPmRcount', EmPmRcount)
            call netcdf_read(trim(filename), 'EmPmRrevise', EmPmRrevise)
            CALL MPI_BCAST(EmPmRcount, 1, MPI_INTEGER, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(EmPmRrevise, 1, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         end if

         ! read restart for ocean thermodynamic module
         call netcdf_read(trim(filename), 'rhoLev', rhoLev, d1=nk, ts=1)
         CALL MPI_BCAST(rhoLev, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)

      ELSE

         if (date_res) then
            CALL MPI_BCAST(nowjulian_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjultq_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjuluv_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulrad_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulprc_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulslp_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulrunoff_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulsfrs_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(nxtjulsfrt_dp, 1, MPI_DOUBLE_PRECISION, mpi_root_comp, mpi_comp_comm, mpi_err)
            
            nowjulian = int(nowjulian_dp, i8)
            nxtjultq = int(nxtjultq_dp, i8)
            nxtjuluv = int(nxtjuluv_dp, i8)
            nxtjulrad = int(nxtjulrad_dp, i8)
            nxtjulprc = int(nxtjulprc_dp, i8)
            nxtjulslp = int(nxtjulslp_dp, i8)
            nxtjulrunoff = int(nxtjulrunoff_dp, i8)
            nxtjulsfrs = int(nxtjulsfrs_dp, i8)
            nxtjulsfrt = int(nxtjulsfrt_dp, i8)

            call jul2greg(nowsec, nowmin, nowhour, nowday, nowmon, nowyear, nowjulian)
         end if

         if (ln_EmPmRevise) then
            CALL MPI_BCAST(EmPmRcount, 1, MPI_INTEGER, mpi_root_comp, mpi_comp_comm, mpi_err)
            CALL MPI_BCAST(EmPmRrevise, 1, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         end if

         CALL MPI_BCAST(rhoLev, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         ! read restart for ocean dynamic module
       END IF
       
       !YY, nowtime = nowjulian - startTime, nowTime defines model starts at 0 second.
       startTime = nowjulian

       ! read restart for ocean thermodynamic module
       CALL mpi_netcdf_read_exchange(filename, 'tFld', tFld(1:nlpb,1:nk,1), 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'sFld', tFld(1:nlpb,1:nk,2), 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'gTracerPastT', gTracerPast(1:nlpb,1:nk,1), 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'gTracerPastS', gTracerPast(1:nlpb,1:nk,2), 3, optional_ts = 1)
       
       if (ln_EmPmRevise) CALL mpi_netcdf_read_exchange(filename, 'EmPmR_long', EmPmR_long, 1, optional_ts = 1)

      ! in 2D, 1 for nlpb and ni, 2 for nlpbz and ni ,3 for nlpb and nk, 4 for nlpbz and nk, 5 for nlpb and sbct
       CALL mpi_netcdf_read_exchange(filename, 'u10', u10, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'v10', v10, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 't10', t10, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'q10', q10, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'slp', slp, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'swdn', swdn, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'lwdn', lwdn, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'prec', prec, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'snow', snow, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'runoff', runoff, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'sfrs', sfrs, 5, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'sfrt', sfrt, 5, optional_ts = 1)

       ! read restart for ocean dynamic module
       CALL mpi_netcdf_read_exchange(filename, 'pbt_base', pbt_base, 1,optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'etaH', etaH, 1,optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'etaN', etaN, 1,optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'uFld', uFld, 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'vFld', vFld, 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'wFld', wFld(1:nlpb, 2:nkp1), 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'gUpast', gUpast, 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'gVpast', gVpast, 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'KappaRM', KappaRM(1:nlpb, 2:nkp1), 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'KappaRT', KappaRT(1:nlpb, 2:nkp1), 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'hmxl', hmxl(1:nlpb, 2:nkp1), 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'en', en(1:nlpb, 2:nkp1), 3, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'rStarFacC', rStarFacC, 1, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'rStarFacW', rStarFacW, 1, optional_ts = 1)
       CALL mpi_netcdf_read_exchange(filename, 'rStarFacS', rStarFacS, 1, optional_ts = 1)

         ! read restart if ln_bdy is ture
         if (ln_bdy) then
            CALL mpi_netcdf_read_exchange(filename, 'uFldBcl_past', uFldBcl_past, 3, optional_ts = 1)
            CALL mpi_netcdf_read_exchange(filename, 'vFldBcl_past', vFldBcl_past, 3, optional_ts = 1)
         end if

         !YY: read seaice hotstart if mitice_on is true
         ! note that SEAICE_multDim can be 1
         if (mitice_on) then
            call mpi_netcdf_read_exchange(filename, 'TICES', TICES, 9, optional_dim1=SEAICE_multDim, optional_ts = 1)
            !$acc update device(TICES)
#ifdef SEAICE_ITD
            if (SEAICE_multDim .GT. 1) then      !debug, YUAN 
              call mpi_netcdf_read_exchange(filename, 'AREAITD', AREAITD, 9, optional_dim1=SEAICE_multDim, optional_ts = 1)
              call mpi_netcdf_read_exchange(filename, 'HEFFITD', HEFFITD, 9, optional_dim1=SEAICE_multDim, optional_ts = 1)
              call mpi_netcdf_read_exchange(filename, 'HSNOWITD', HSNOWITD, 9, optional_dim1=SEAICE_multDim, optional_ts = 1)
              !$acc update device(AREAITD,HEFFITD,HSNOWITD)
!             update total ice area as well as mean ice and snow thickness
              call seaice_itd_sum


            else
              call mpi_netcdf_read_exchange(filename, 'AREA', AREA, 1,optional_ts = 1)
              call mpi_netcdf_read_exchange(filename, 'HEFF', HEFF, 1,optional_ts = 1)
              call mpi_netcdf_read_exchange(filename, 'HSNOW', HSNOW, 1,optional_ts = 1)
              !$acc update device(AREA,HEFF,HSNOW)
!             redistribute over categories, assuming a log-normal distribtuion
              call seaice_itd_pickup
            endif
#else
            call mpi_netcdf_read_exchange(filename, 'AREA', AREA, 1,optional_ts = 1)
            call mpi_netcdf_read_exchange(filename, 'HEFF', HEFF, 1,optional_ts = 1)
            call mpi_netcdf_read_exchange(filename, 'HSNOW', HSNOW, 1,optional_ts = 1)
            !$acc update device(AREA,HEFF,HSNOW)
#endif
#ifdef SEAICE_VARIABLE_SALINITY
            call mpi_netcdf_read_exchange(filename, 'HSALT', HSALT, 1,optional_ts = 1)
            !$acc update device(HSALT)
#endif
            call mpi_netcdf_read_exchange(filename, 'UICE', UICE, 1,optional_ts = 1)
            call mpi_netcdf_read_exchange(filename, 'VICE', VICE, 1,optional_ts = 1)
            !$acc update device(UICE,VICE)
            if (SEAICEuseEVPpickup) then     !default is true, set in seaice namelist
              call mpi_netcdf_read_exchange(filename, 'seaice_sigma1', seaice_sigma1, 1,optional_ts = 1)
              call mpi_netcdf_read_exchange(filename, 'seaice_sigma2', seaice_sigma2, 1,optional_ts = 1)
              call mpi_netcdf_read_exchange(filename, 'seaice_sigma12', seaice_sigma12, 1,optional_ts = 1)
              !$acc update device(seaice_sigma1,seaice_sigma2,seaice_sigma12)
            endif
!YY: if hot start, the first frame is what just read in
          IF (SEAICEwriteState) THEN     !default is false. set in init.
            write(prefix, '(A,i4.4,i2.2,i2.2,i2.2,i2.2,i2.2)')   &
              trim(adjustl(foldername))//'dump_',nowyear,nowmon,nowday,nowhour,nowmin,nowsec
            call seaice_io_dump(prefix,1)
            call seaice_io_init
            NxtDumpTime = int(SEAICE_dumpFreq, 8)   ! put in the end
          ENDIF
         endif    ! end for mitice_on 

      do k = 1, nk
         hFacC(1:nlpb, k) = h0FacC(1:nlpb, k)*rStarFacC(1:nlpb)
         hFacW(1:nlpb, k) = h0FacW(1:nlpb, k)*rStarFacW(1:nlpb)
         hFacS(1:nlpb, k) = h0FacS(1:nlpb, k)*rStarFacS(1:nlpb)
      end do

      do k = 1, nk
         do i = 1, nlpb
            if (hFacC(i, k) .ne. 0.0_wp) then
               recip_hFacC(i, k) = 1.0_wp/hFacC(i, k)
            else
               recip_hFacC(i, k) = 0.0_wp
            end if
            if (hFacW(i, k) .ne. 0.0_wp) then
               recip_hFacW(i, k) = 1.0_wp/hFacW(i, k)
            else
               recip_hFacW(i, k) = 0.0_wp
            end if
            if (hFacS(i, k) .ne. 0.0_wp) then
               recip_hFacS(i, k) = 1.0_wp/hFacS(i, k)
            else
               recip_hFacS(i, k) = 0.0_wp
            end if
         end do
      end do

      do k = 1, nk
         if (mpi_rank == 0) then
            if (ln_bous) then
               write (*, "(a,i4,a,f16.6,a,e15.6,a)") &
                  'depth of level ', k, ' is ', rC(k), &
                  ' m, mean density inverse is ', rhoLev(k), ' m3/kg'
            else
               write (*, "(a,i4,a,f16.6,a,e15.6,a)") &
                  'depth of level ', k, ' is ', rC(k), &
                  ' Pa, mean density inverse is ', rhoLev(k), ' m3/kg'
            end if
         end if
      end do

   end subroutine mpi_csp_io_restart_read

!==============================================================================
   subroutine mpi_csp_io_restart_main
!==============================================================================
      implicit none

      call mpi_csp_io_restart_recv_info

      call mpi_csp_io_restart_write

   end subroutine mpi_csp_io_restart_main

!==============================================================================
   subroutine mpi_csp_io_restart_send
!==============================================================================
      implicit none
      integer :: mpi_req_restart
      logical :: mpi_flag_restart
      integer :: i, sbct

      sbct = 4

      if (mpi_rank == 0) call mpi_csp_io_restart_send_info

      if (mpi_rank == 0) CALL MPI_SEND(rhoLev, nk, mpi_real_wp, mpi_root_comp_io, &
                                       mpi_comp_io_tag, mpi_comp_io_comm, mpi_err)

      mpi_flag_restart = .false.

      !$ACC update self(tFld, gTracerPast)

      CALL mpi_csp_io_send_output_2d(tFld(:, :, 1), mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(tFld(:, :, 2), mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(gTracerPast(:, :, 1), mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(gTracerPast(:, :, 2), mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)

      !$ACC update self(u10, v10, t10, q10, slp, swdn, lwdn, prec, snow, runoff, sfrs, sfrt, EmPmR_long)

      CALL mpi_csp_io_send_output_1d(EmPmR_long, mpi_comp_buf_1d_restart, &
                                     mpi_req_restart, mpi_flag_restart)

      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(u10(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(v10(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(t10(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(q10(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(slp(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(swdn(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(lwdn(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(prec(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(snow(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(runoff(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(sfrs(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do
      do i = 1, sbct
         CALL mpi_csp_io_send_output_1d(sfrt(:, i), mpi_comp_buf_1d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end do

      !$ACC update self(pbt_base, etaN, etaH, uFld, vFld, wFld, gUpast, gVpast, &
      !$ACC             KappaRM, KappaRT, hmxl, en, rStarFacC, rStarFacW, rStarFacS)

      CALL mpi_csp_io_send_output_1d(pbt_base, mpi_comp_buf_1d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_1d(etaN, mpi_comp_buf_1d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_1d(etaH, mpi_comp_buf_1d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(uFld, mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(vFld, mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(wFld(:, 2:nkp1), mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(gUpast, mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(gVpast, mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(KappaRM(:, 2:nkp1), mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(KappaRT(:, 2:nkp1), mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(hmxl(:, 2:nkp1), mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_2d(en(:, 2:nkp1), mpi_comp_buf_2d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_1d(rStarFacC, mpi_comp_buf_1d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_1d(rStarFacW, mpi_comp_buf_1d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      CALL mpi_csp_io_send_output_1d(rStarFacS, mpi_comp_buf_1d_restart, &
                                     mpi_req_restart, mpi_flag_restart)

      if (ln_bdy) then
         !$ACC update self(uFldBcl_past, vFldBcl_past)

         CALL mpi_csp_io_send_output_2d(uFldBcl_past, mpi_comp_buf_2d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
         CALL mpi_csp_io_send_output_2d(vFldBcl_past, mpi_comp_buf_2d_restart, &
                                        mpi_req_restart, mpi_flag_restart)
      end if


!YY,ZY: add seaice components hotstart
      if (mitice_on) then
        !$ACC update self(seaice_sigma1,seaice_sigma2,seaice_sigma12,  &
#ifdef SEAICE_VARIABLE_SALINITY
        !$ACC             HSALT,                                       &
#endif
#ifdef SEAICE_ITD
        !$ACC             AREAITD,HEFFITD,HSNOWITD,                    &
#endif 
        !$ACC            TICES,HEFF,AREA,HSNOW,UICE,VICE)

        do i = 1,SEAICE_multDim    !SEAICE_multDim can be 1 if SEAICE_ITD not defined
          CALL mpi_csp_io_send_output_1d(TICES(:,i), mpi_comp_buf_1d_restart, &
                                       mpi_req_restart, mpi_flag_restart)
        enddo
#ifdef SEAICE_ITD
        do i = 1,SEAICE_multDim    !SEAICE_multDim can be 1 if SEAICE_ITD not defined
          CALL mpi_csp_io_send_output_1d(AREAITD(:,i), mpi_comp_buf_1d_restart, &
                                       mpi_req_restart, mpi_flag_restart)
        enddo
        do i = 1,SEAICE_multDim
          CALL mpi_csp_io_send_output_1d(HEFFITD(:,i), mpi_comp_buf_1d_restart, &
                                       mpi_req_restart, mpi_flag_restart)
        enddo
        do i = 1,SEAICE_multDim
          CALL mpi_csp_io_send_output_1d(HSNOWITD(:,i), mpi_comp_buf_1d_restart, &
                                       mpi_req_restart, mpi_flag_restart)
        enddo
#else
        CALL mpi_csp_io_send_output_1d(AREA, mpi_comp_buf_1d_restart, &
                                    mpi_req_restart, mpi_flag_restart)
        CALL mpi_csp_io_send_output_1d(HEFF, mpi_comp_buf_1d_restart, &
                                    mpi_req_restart, mpi_flag_restart)
        CALL mpi_csp_io_send_output_1d(HSNOW, mpi_comp_buf_1d_restart, &
                                       mpi_req_restart, mpi_flag_restart)
#endif
#ifdef SEAICE_VARIABLE_SALINITY
        CALL mpi_csp_io_send_output_1d(HSALT, mpi_comp_buf_1d_restart, &
                                       mpi_req_restart, mpi_flag_restart)
#endif

        CALL mpi_csp_io_send_output_1d(UICE, mpi_comp_buf_1d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
        CALL mpi_csp_io_send_output_1d(VICE, mpi_comp_buf_1d_restart, &
                                     mpi_req_restart, mpi_flag_restart)
      !if (SEAICEuseEVP)
        CALL mpi_csp_io_send_output_1d(seaice_sigma1, mpi_comp_buf_1d_restart, &
                                      mpi_req_restart, mpi_flag_restart)
        CALL mpi_csp_io_send_output_1d(seaice_sigma2, mpi_comp_buf_1d_restart, &
                                      mpi_req_restart, mpi_flag_restart)
        CALL mpi_csp_io_send_output_1d(seaice_sigma12, mpi_comp_buf_1d_restart, &
                                      mpi_req_restart, mpi_flag_restart)
      end if     !endif for mitice_on

      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)

   end subroutine mpi_csp_io_restart_send

!==============================================================================
   SUBROUTINE mpi_csp_io_restart_send_info
!==============================================================================
      nxtjultq_dp = dble(nxtjultq)
      nxtjuluv_dp = dble(nxtjuluv)
      nxtjulrad_dp = dble(nxtjulrad)
      nxtjulprc_dp = dble(nxtjulprc)
      nxtjulslp_dp = dble(nxtjulslp)
      nxtjulrunoff_dp = dble(nxtjulrunoff)
      nxtjulsfrs_dp = dble(nxtjulsfrs)
      nxtjulsfrt_dp = dble(nxtjulsfrt)

      mpi_buf_size = max_domain*20*8
      ALLOCATE (mpi_buf(mpi_buf_size))
      mpi_position = 0
      CALL MPI_PACK(nxtjultq_dp, 1, MPI_DOUBLE, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nxtjuluv_dp, 1, MPI_DOUBLE, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nxtjulrad_dp, 1, MPI_DOUBLE, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nxtjulprc_dp, 1, MPI_DOUBLE, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nxtjulslp_dp, 1, MPI_DOUBLE, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nxtjulrunoff_dp, 1, MPI_DOUBLE, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nxtjulsfrs_dp, 1, MPI_DOUBLE, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nxtjulsfrt_dp, 1, MPI_DOUBLE, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(ln_bdy, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(mitice_on, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nlpb, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(EmPmRcount, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(EmPmRrevise, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      if (mitice_on) &
      CALL MPI_PACK(SEAICE_multDim, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_SEND(mpi_position, 1, MPI_INTEGER, mpi_root_comp_io, mpi_comp_io_tag, mpi_comp_io_comm, mpi_err)
      CALL MPI_SEND(mpi_buf, mpi_position, MPI_PACKED, mpi_root_comp_io, mpi_comp_io_tag, mpi_comp_io_comm, mpi_err)
      DEALLOCATE (mpi_buf)

   END SUBROUTINE mpi_csp_io_restart_send_info

!==============================================================================
   SUBROUTINE mpi_csp_io_restart_recv_info
!==============================================================================
      CALL MPI_RECV(mpi_buf_size, 1, MPI_INTEGER, mpi_sub_root_comp_io, &
                    mpi_comp_io_tag, mpi_comp_io_comm, mpi_comp_io_status, mpi_err)
      ALLOCATE (mpi_buf(mpi_buf_size))
      CALL MPI_RECV(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_sub_root_comp_io, &
                    mpi_comp_io_tag, mpi_comp_io_comm, mpi_comp_io_status, mpi_err)
      mpi_position = 0
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nxtjultq_dp, 1, MPI_DOUBLE, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nxtjuluv_dp, 1, MPI_DOUBLE, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nxtjulrad_dp, 1, MPI_DOUBLE, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nxtjulprc_dp, 1, MPI_DOUBLE, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nxtjulslp_dp, 1, MPI_DOUBLE, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nxtjulrunoff_dp, 1, MPI_DOUBLE, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nxtjulsfrs_dp, 1, MPI_DOUBLE, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nxtjulsfrt_dp, 1, MPI_DOUBLE, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ln_bdy, 1, MPI_LOGICAL, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, mitice_on, 1, MPI_LOGICAL, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nlpb, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, EmPmRcount, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, EmPmRrevise, 1, mpi_real_wp, mpi_comp_io_comm, mpi_err)
      if(mitice_on) &
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_multDim, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      DEALLOCATE (mpi_buf)
   END SUBROUTINE mpi_csp_io_restart_recv_info

!==============================================================================
   subroutine mpi_csp_io_restart_write
!==============================================================================
      implicit none
      integer :: i
      integer :: ncid, dimid_nlpb, dimid_nk, dimid_sbct, dimid_t
      integer :: rhoLev_id, tFld_id, sFld_id, gTracerPastT_id, gTracerPastS_id
      integer :: nowjulian_id, nxtjultq_id, nxtjuluv_id, nxtjulrad_id, nxtjulprc_id, &
                 nxtjulslp_id, nxtjulrunoff_id, nxtjulsfrs_id, nxtjulsfrt_id, &
                 EmPmRrevise_id, EmPmRcount_id, &
                 EmPmR_long_id, u10_id, v10_id, t10_id, q10_id, slp_id, swdn_id, &
                 lwdn_id, prec_id, snow_id, runoff_id, sfrs_id, sfrt_id
      integer :: pbt_base_id, etaN_id, etaH_id, uFld_id, vFld_id, wFld_id, &
                 gUpast_id, gVpast_id, KappaRM_id, KappaRT_id
      integer :: hmxl_id, en_id, rStarFacC_id, rStarFacW_id, rStarFacS_id
      integer :: uFldBcl_past_id, vFldBcl_past_id
      integer :: xtype, ts, sbct
      character(lc) :: filename, varname, time_str, domid_str
      character(lc) :: time, date, zone, timestamp, RestartDate
      !YY
      integer :: dimid_nITD      ! debug, YUAN: even SEAICE_nITD not defiined,keep it.
#ifdef SEAICE_ITD
      integer :: areaitd_id,heffitd_id,hsnowitd_id
#endif
#ifdef SEAICE_VARIABLE_SALINITY
      integer :: hsalt_id
#endif
      integer :: tices_id,area_id,heff_id,hsnow_id,uice_id,vice_id
      integer :: seaice_sigma1_id,seaice_sigma2_id,seaice_sigma12_id

      sbct = 4
      ts = 1

      write (RestartDate, "(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)") &
         nowyear, '-', nowmon, '-', nowday, ' ', nowhour, ':', nowmin, ':', nowsec
      write (time_str, "(i10.10)") myIterm*int(dTtracer)
      write (domid_str, "(i2.2)") idom
      filename = "nmefc_macom_restart_"//trim(time_str)//"_d"//trim(domid_str)//".nc"

      call netcdf_check(nf90_create(trim(filename), nf90_netcdf4, ncid), trim(filename))
      call netcdf_check(nf90_def_dim(ncid, 'nlpb', mpi_total_nlpb, dimid_nlpb), trim(filename), varname='nlpb')
      call netcdf_check(nf90_def_dim(ncid, 'nk', nk, dimid_nk), trim(filename), varname='nk')
      call netcdf_check(nf90_def_dim(ncid, 'sbct', sbct, dimid_sbct), trim(filename), varname='sbct')
      if (mitice_on) &
      !YY:define seaice thickness categories nITD:  SEAICE_multDim can be 1 if SEAICE_ITD not defined
        call netcdf_check(nf90_def_dim(ncid, 'nITD', SEAICE_multDim, dimid_nITD), trim(filename), varname='nITD')
      !
      call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename), varname='time')
      call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                     "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), &
                        trim(filename), varname='Creater')
      call date_and_time(date=date, time=time, zone=zone)
      timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                  time(1:2)//":"//time(3:4)//":"//time(5:6)//" "// &
                  zone
      call netcdf_check(nf90_put_att(ncid, nf90_global, "RestartDate", &
                                     trim(RestartDate)), trim(filename), varname='RestartDate')
      call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), &
                        trim(filename), varname='TimeStamp')
      if (wp == sp) then
         xtype = nf90_float
      else if (wp == dp) then
         xtype = nf90_double
      else
         write (*, *) 'ERROR! : restart only support sp or dp, but wp is ', wp, ' now'
      end if

      ! define restart for ocean thermodynamic module
      call netcdf_check(nf90_def_var(ncid, 'rhoLev', xtype, &
                                     (/dimid_nk, dimid_t/), rhoLev_id), trim(filename), varname='rhoLev')
      call netcdf_check(nf90_def_var(ncid, 'tFld', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), tFld_id), trim(filename), varname='tFld')
      call netcdf_check(nf90_def_var(ncid, 'sFld', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), sFld_id), trim(filename), varname='sFld')
      call netcdf_check(nf90_def_var(ncid, 'gTracerPastT', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), gTracerPastT_id), trim(filename), &
                        varname='gTracerPastT')
      call netcdf_check(nf90_def_var(ncid, 'gTracerPastS', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), gTracerPastS_id), trim(filename), &
                        varname='gTracerPastS')

      ! define restart for Atmosphere-Ocean flux module
      call netcdf_check(nf90_def_var(ncid, 'nowjulian', nf90_double, &
                                     nowjulian_id), trim(filename), varname='nowjulian')
      call netcdf_check(nf90_put_att(ncid, nowjulian_id, "units", &
                                     "seconds since 1949-10-01 00:00:00"), trim(filename))
      call netcdf_check(nf90_def_var(ncid, 'nxtjultq', nf90_double, &
                                     nxtjultq_id), trim(filename), varname='nxtjultq')
      call netcdf_check(nf90_put_att(ncid, nxtjultq_id, "units", &
                                     "seconds since 1949-10-01 00:00:00"), trim(filename))
      call netcdf_check(nf90_def_var(ncid, 'nxtjuluv', nf90_double, &
                                     nxtjuluv_id), trim(filename), varname='nxtjuluv')
      call netcdf_check(nf90_put_att(ncid, nxtjuluv_id, "units", &
                                     "seconds since 1949-10-01 00:00:00"), trim(filename))
      call netcdf_check(nf90_def_var(ncid, 'nxtjulrad', nf90_double, &
                                     nxtjulrad_id), trim(filename), varname='nxtjulrad')
      call netcdf_check(nf90_put_att(ncid, nxtjulrad_id, "units", &
                                     "seconds since 1949-10-01 00:00:00"), trim(filename))
      call netcdf_check(nf90_def_var(ncid, 'nxtjulprc', nf90_double, &
                                     nxtjulprc_id), trim(filename), varname='nxtjulprc')
      call netcdf_check(nf90_put_att(ncid, nxtjulprc_id, "units", &
                                     "seconds since 1949-10-01 00:00:00"), trim(filename))
      call netcdf_check(nf90_def_var(ncid, 'nxtjulslp', nf90_double, &
                                     nxtjulslp_id), trim(filename), varname='nxtjulslp')
      call netcdf_check(nf90_put_att(ncid, nxtjulslp_id, "units", &
                                     "seconds since 1949-10-01 00:00:00"), trim(filename))
      call netcdf_check(nf90_def_var(ncid, 'nxtjulrunoff', nf90_double, &
                                     nxtjulrunoff_id), trim(filename), varname='nxtjulrunoff')
      call netcdf_check(nf90_put_att(ncid, nxtjulrunoff_id, "units", &
                                     "seconds since 1949-10-01 00:00:00"), trim(filename))
      call netcdf_check(nf90_def_var(ncid, 'nxtjulsfrs', nf90_double, &
                                     nxtjulsfrs_id), trim(filename), varname='nxtjulsfrs')
      call netcdf_check(nf90_put_att(ncid, nxtjulsfrs_id, "units", &
                                     "seconds since 1949-10-01 00:00:00"), trim(filename))
      call netcdf_check(nf90_def_var(ncid, 'nxtjulsfrt', nf90_double, &
                                     nxtjulsfrt_id), trim(filename), varname='nxtjulsfrt')
      call netcdf_check(nf90_put_att(ncid, nxtjulsfrt_id, "units", &
                                     "seconds since 1949-10-01 00:00:00"), trim(filename))
      call netcdf_check(nf90_def_var(ncid, 'EmPmRcount', nf90_int, &
                                     EmPmRcount_id), trim(filename), varname='EmPmRcount')
      call netcdf_check(nf90_def_var(ncid, 'EmPmRrevise', xtype, &
                                     EmPmRrevise_id), trim(filename), varname='EmPmRrevise')
      call netcdf_check(nf90_def_var(ncid, 'EmPmR_long', xtype, &
                                     (/dimid_nlpb, dimid_t/), EmPmR_long_id), trim(filename), varname='EmPmR_long')
      call netcdf_check(nf90_def_var(ncid, 'u10', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), u10_id), trim(filename), varname='u10')
      call netcdf_check(nf90_def_var(ncid, 'v10', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), v10_id), trim(filename), varname='v10')
      call netcdf_check(nf90_def_var(ncid, 't10', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), t10_id), trim(filename), varname='t10')
      call netcdf_check(nf90_def_var(ncid, 'q10', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), q10_id), trim(filename), varname='q10')
      call netcdf_check(nf90_def_var(ncid, 'slp', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), slp_id), trim(filename), varname='slp')
      call netcdf_check(nf90_def_var(ncid, 'swdn', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), swdn_id), trim(filename), varname='swdn')
      call netcdf_check(nf90_def_var(ncid, 'lwdn', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), lwdn_id), trim(filename), varname='lwdn')
      call netcdf_check(nf90_def_var(ncid, 'prec', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), prec_id), trim(filename), varname='prec')
      call netcdf_check(nf90_def_var(ncid, 'snow', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), snow_id), trim(filename), varname='snow')
      call netcdf_check(nf90_def_var(ncid, 'runoff', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), runoff_id), trim(filename), varname='runoff')
      call netcdf_check(nf90_def_var(ncid, 'sfrs', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), sfrs_id), trim(filename), varname='sfrs')
      call netcdf_check(nf90_def_var(ncid, 'sfrt', xtype, &
                                     (/dimid_nlpb, dimid_sbct, dimid_t/), sfrt_id), trim(filename), varname='sfrt')
      !YY: define seaice vars
      if (mitice_on) then
        call netcdf_check(nf90_def_var(ncid, 'TICES', xtype, &
                         (/dimid_nlpb, dimid_nITD, dimid_t/), tices_id), trim(filename), varname='TICES')
#ifdef SEAICE_ITD
        call netcdf_check(nf90_def_var(ncid, 'AREAITD', xtype, &
                         (/dimid_nlpb, dimid_nITD, dimid_t/), areaitd_id), trim(filename), varname='AREAITD')
        call netcdf_check(nf90_def_var(ncid, 'HEFFITD', xtype, &
                         (/dimid_nlpb, dimid_nITD, dimid_t/), heffitd_id), trim(filename), varname='HEFFITD')
        call netcdf_check(nf90_def_var(ncid, 'HSNOWITD', xtype, &
                         (/dimid_nlpb, dimid_nITD, dimid_t/), hsnowitd_id), trim(filename), varname='HSNOWITD')
#else
        call netcdf_check(nf90_def_var(ncid, 'AREA', xtype, &
                         (/dimid_nlpb, dimid_t/), area_id), trim(filename), varname='AREA')
        call netcdf_check(nf90_def_var(ncid, 'HEFF', xtype, &
                         (/dimid_nlpb, dimid_t/), heff_id), trim(filename), varname='HEFF')
        call netcdf_check(nf90_def_var(ncid, 'HSNOW', xtype, &
                         (/dimid_nlpb, dimid_t/), hsnow_id), trim(filename), varname='HSNOW')
#endif
#ifdef SEAICE_VARIABLE_SALINITY
        call netcdf_check(nf90_def_var(ncid, 'HSALT', xtype, &
                         (/dimid_nlpb, dimid_t/), hsalt_id), trim(filename), varname='HSALT')
#endif
        call netcdf_check(nf90_def_var(ncid, 'UICE', xtype, &
                         (/dimid_nlpb, dimid_t/), uice_id), trim(filename), varname='UICE')
        call netcdf_check(nf90_def_var(ncid, 'VICE', xtype, &
                         (/dimid_nlpb, dimid_t/), vice_id), trim(filename), varname='VICE')
        call netcdf_check(nf90_def_var(ncid, 'seaice_sigma1', xtype, &
                         (/dimid_nlpb, dimid_t/), seaice_sigma1_id), trim(filename), varname='seaice_sigma1')
        call netcdf_check(nf90_def_var(ncid, 'seaice_sigma2', xtype, &
                         (/dimid_nlpb, dimid_t/), seaice_sigma2_id), trim(filename), varname='seaice_sigma2')
        call netcdf_check(nf90_def_var(ncid, 'seaice_sigma12', xtype, &
                         (/dimid_nlpb, dimid_t/), seaice_sigma12_id), trim(filename), varname='seaice_sigma12')
      end if     !endif for mitice_on

      ! define restart for ocean dynamic module
      call netcdf_check(nf90_def_var(ncid, 'pbt_base', xtype, &
                                     (/dimid_nlpb, dimid_t/), pbt_base_id), trim(filename), varname='pbt_base')
      call netcdf_check(nf90_def_var(ncid, 'etaN', xtype, &
                                     (/dimid_nlpb, dimid_t/), etaN_id), trim(filename), varname='etaN')
      call netcdf_check(nf90_def_var(ncid, 'etaH', xtype, &
                                     (/dimid_nlpb, dimid_t/), etaH_id), trim(filename), varname='etaH')
      call netcdf_check(nf90_def_var(ncid, 'uFld', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), uFld_id), trim(filename), varname='uFld')
      call netcdf_check(nf90_def_var(ncid, 'vFld', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), vFld_id), trim(filename), varname='vFld')
      call netcdf_check(nf90_def_var(ncid, 'wFld', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), wFld_id), trim(filename), varname='wFld')
      call netcdf_check(nf90_def_var(ncid, 'gUpast', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), gUpast_id), trim(filename), varname='gUpast')
      call netcdf_check(nf90_def_var(ncid, 'gVpast', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), gVpast_id), trim(filename), varname='gVpast')
      call netcdf_check(nf90_def_var(ncid, 'KappaRM', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), KappaRM_id), trim(filename), varname='KappaRM')
      call netcdf_check(nf90_def_var(ncid, 'KappaRT', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), KappaRT_id), trim(filename), varname='KappaRT')
      call netcdf_check(nf90_def_var(ncid, 'hmxl', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), hmxl_id), trim(filename), varname='hmxl')
      call netcdf_check(nf90_def_var(ncid, 'en', xtype, &
                                     (/dimid_nlpb, dimid_nk, dimid_t/), en_id), trim(filename), varname='en')
      call netcdf_check(nf90_def_var(ncid, 'rStarFacC', xtype, &
                                     (/dimid_nlpb, dimid_t/), rStarFacC_id), trim(filename), varname='rStarFacC')
      call netcdf_check(nf90_def_var(ncid, 'rStarFacW', xtype, &
                                     (/dimid_nlpb, dimid_t/), rStarFacW_id), trim(filename), varname='rStarFacW')
      call netcdf_check(nf90_def_var(ncid, 'rStarFacS', xtype, &
                                     (/dimid_nlpb, dimid_t/), rStarFacS_id), trim(filename), varname='rStarFacS')

      ! define restart for regional bdy use
      if (ln_bdy) then
         call netcdf_check(nf90_def_var(ncid, 'uFldBcl_past', xtype, &
              (/dimid_nlpb, dimid_nk, dimid_t/), uFldBcl_past_id), trim(filename), varname='uFldBcl_past')
         call netcdf_check(nf90_def_var(ncid, 'vFldBcl_past', xtype, &
              (/dimid_nlpb, dimid_nk, dimid_t/), vFldBcl_past_id), trim(filename), varname='vFldBcl_past')
      end if

      call netcdf_check(nf90_enddef(ncid), trim(filename))

      ! write restart for ocean thermodynamic module
      if (.not. allocated(mpi_io_rhoLev)) allocate (mpi_io_rhoLev(nk))
      CALL MPI_RECV(mpi_io_rhoLev, nk, mpi_real_wp, mpi_sub_root_comp_io, &
                    mpi_comp_io_tag, mpi_comp_io_comm, mpi_comp_io_status, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "rhoLev", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, rhoLev_id, mpi_io_rhoLev, &
           start=(/1, ts/), count=(/nk, 1/)), trim(filename), varname='rhoLev')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "tFld", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, tFld_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='tFld')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "sFld", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, sFld_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='sFld')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "gTracerPastT", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, gTracerPastT_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='gTracerPastT')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "gTracerPastS", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, gTracerPastS_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='gTracerPastS')

      ! write restart for Atmosphere-Ocean flux module
      call netcdf_check(nf90_put_var(ncid, nowjulian_id, dble(nowjulian)), trim(filename), varname='nowjulian')
      call netcdf_check(nf90_put_var(ncid, nxtjultq_id, nxtjultq_dp), trim(filename), varname='nxtjultq')
      call netcdf_check(nf90_put_var(ncid, nxtjuluv_id, nxtjuluv_dp), trim(filename), varname='nxtjuluv')
      call netcdf_check(nf90_put_var(ncid, nxtjulrad_id, nxtjulrad_dp), trim(filename), varname='nxtjulrad')
      call netcdf_check(nf90_put_var(ncid, nxtjulprc_id, nxtjulprc_dp), trim(filename), varname='nxtjulprc')
      call netcdf_check(nf90_put_var(ncid, nxtjulslp_id, nxtjulslp_dp), trim(filename), varname='nxtjulslp')
      call netcdf_check(nf90_put_var(ncid, nxtjulrunoff_id, nxtjulrunoff_dp), trim(filename), varname='nxtjulrunoff')
      call netcdf_check(nf90_put_var(ncid, nxtjulsfrs_id, nxtjulsfrs_dp), trim(filename), varname='nxtjulsfrs')
      call netcdf_check(nf90_put_var(ncid, nxtjulsfrt_id, nxtjulsfrt_dp), trim(filename), varname='nxtjulsfrt')
      call netcdf_check(nf90_put_var(ncid, EmPmRcount_id, EmPmRcount), trim(filename), varname='EmPmRcount')
      call netcdf_check(nf90_put_var(ncid, EmPmRrevise_id, EmPmRrevise), trim(filename), varname='EmPmRrevise')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
           mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
           mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      write (*, "(a,a13,a,a)") " recieve and write ", "EmPmR_long", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, EmPmR_long_id, mpi_io_buf_adjust_1d, &
           start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='EmPmR_long')

      if (.not. allocated(mpi_io_sbc)) allocate (mpi_io_sbc(mpi_total_nlpb, sbct))
      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "u10", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, u10_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='u10')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "v10", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, v10_id, mpi_io_sbc, &
                                     start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='v10')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "t10", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, t10_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='t10')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "q10", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, q10_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='q10')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "slp", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, slp_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='slp')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "swdn", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, swdn_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='swdn')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "lwdn", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, lwdn_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='lwdn')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "prec", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, prec_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='prec')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "snow", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, snow_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='snow')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "runoff", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, runoff_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='runoff')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "sfrs", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, sfrs_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='sfrs')

      do i = 1, sbct
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
              mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         mpi_io_sbc(:, i) = mpi_io_buf_adjust_1d(:)
      end do
      write (*, "(a,a13,a,a)") " recieve and write ", "sfrt", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, sfrt_id, mpi_io_sbc, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, sbct, 1/)), trim(filename), varname='sfrt')

!      ! write restart for ocean dynamic module
      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
           mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
           mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      write (*, "(a,a13,a,a)") " recieve and write ", "pbt_base", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, pbt_base_id, mpi_io_buf_adjust_1d, &
           start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='pbt_base')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
           mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
           mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      write (*, "(a,a13,a,a)") " recieve and write ", "etaN", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, etaN_id, mpi_io_buf_adjust_1d, &
           start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='etaN')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
           mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
           mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      write (*, "(a,a13,a,a)") " recieve and write ", "etaH", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, etaH_id, mpi_io_buf_adjust_1d, &
           start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='etaH')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "uFld", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, uFld_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='uFld')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "vFld", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, vFld_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='vFld')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "wFld", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, wFld_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='wFld')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "gUpast", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, gUpast_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='gUpast')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "gVpast", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, gVpast_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='gVpast')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "KappaRM", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, KappaRM_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='KappaRM')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "KappaRT", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, KappaRT_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='KappaRT')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "hmxl", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, hmxl_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='hmxl')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
           mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
           mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      write (*, "(a,a13,a,a)") " recieve and write ", "en", " in ", trim(filename)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
           mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
           mpi_nlpb_block_1d_displs_all)
      call netcdf_check(nf90_put_var(ncid, en_id, mpi_io_buf_adjust_2d, &
           start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='en')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
           mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
           mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      write (*, "(a,a13,a,a)") " recieve and write ", "rStarFacC", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, rStarFacC_id, mpi_io_buf_adjust_1d, &
           start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='rStarFacC')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
           mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
           mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      write (*, "(a,a13,a,a)") " recieve and write ", "rStarFacW", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, rStarFacW_id, mpi_io_buf_adjust_1d, &
           start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='rStarFacW')

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
           mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
      CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
           mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      write (*, "(a,a13,a,a)") " recieve and write ", "rStarFacS", " in ", trim(filename)
      call netcdf_check(nf90_put_var(ncid, rStarFacS_id, mpi_io_buf_adjust_1d, &
           start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='rStarFacS')

      if (ln_bdy) then
         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
              mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
              mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         write (*, "(a,a13,a,a)") " recieve and write ", "uFldBcl_past", " in ", trim(filename)
         CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
              mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
              mpi_nlpb_block_1d_displs_all)
         call netcdf_check(nf90_put_var(ncid, uFldBcl_past_id, mpi_io_buf_adjust_2d, &
              start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='uFldBcl_past')

         CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_restart, &
              mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
              mpi_comp_io_comm, mpi_req_restart)
         CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
         write (*, "(a,a13,a,a)") " recieve and write ", "vFldBcl_past", " in ", trim(filename)
         CALL mpi_2d_data_adjust(mpi_io_buf_2d_restart, mpi_io_buf_adjust_2d, &
              mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
              mpi_nlpb_block_1d_displs_all)
         call netcdf_check(nf90_put_var(ncid, vFldBcl_past_id, mpi_io_buf_adjust_2d, &
              start=(/1, 1, ts/), count=(/mpi_total_nlpb, nk, 1/)), trim(filename), varname='vFldBcl_past')
      end if

!YY: seaice vars write. Note SEAICE_multDim can be 1.
      if (mitice_on) then
        if (.not. allocated(mpi_io_nITD)) allocate (mpi_io_nITD(mpi_total_nlpb, SEAICE_multDim))
!TICES
        do i = 1, SEAICE_multDim
           CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                             mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, &
                                             mpi_comp_io_comm, mpi_req_restart)
           CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
           CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                                   mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all,         &
                                   mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
           mpi_io_nITD(:, i) = mpi_io_buf_adjust_1d(:)
        end do
        write (*, "(a,a13,a,a)") " recieve and write ", "TICES", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, tices_id, mpi_io_nITD, &
                          start=(/1, 1, ts/), count=(/mpi_total_nlpb, SEAICE_multDim, 1/)), trim(filename), varname='TICES')
#ifdef SEAICE_ITD
!AREAITD
        do i = 1, SEAICE_multDim
           CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                             mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
           CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
           CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                                   mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all,         &
                                   mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
           mpi_io_nITD(:, i) = mpi_io_buf_adjust_1d(:)
        end do
        write (*, "(a,a13,a,a)") " recieve and write ", "AREAITD", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, areaitd_id, mpi_io_nITD, &
                          start=(/1, 1, ts/), count=(/mpi_total_nlpb, SEAICE_multDim, 1/)), trim(filename), varname='AREAITD')
!HEFFITD
        do i = 1, SEAICE_multDim
           CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                             mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, &
                                             mpi_comp_io_comm, mpi_req_restart)
           CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
           CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                                   mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all,         &
                                   mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
           mpi_io_nITD(:, i) = mpi_io_buf_adjust_1d(:)
        end do
        write (*, "(a,a13,a,a)") " recieve and write ", "HEFFITD", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, heffitd_id, mpi_io_nITD, &
                          start=(/1, 1, ts/), count=(/mpi_total_nlpb, SEAICE_multDim, 1/)), trim(filename), varname='HEFFITD')
!HSNOWITD
        do i = 1, SEAICE_multDim
           CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                             mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, &
                                             mpi_comp_io_comm, mpi_req_restart)
           CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
           CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                                   mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all,         &
                                   mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
           mpi_io_nITD(:, i) = mpi_io_buf_adjust_1d(:)
        end do
        write (*, "(a,a13,a,a)") " recieve and write ", "HSNOWITD", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, hsnowitd_id, mpi_io_nITD, &
                          start=(/1, 1, ts/), count=(/mpi_total_nlpb, SEAICE_multDim, 1/)), trim(filename), varname='HSNOWITD')
#else
!AREA
        CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
        CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
        CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
        write (*, "(a,a13,a,a)") " recieve and write ", "AREA", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, area_id, mpi_io_buf_adjust_1d, &
                          start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='AREA')
!HEFF
        CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
        CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
        CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
        write (*, "(a,a13,a,a)") " recieve and write ", "HEFF", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, heff_id, mpi_io_buf_adjust_1d, &
                          start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='HEFF')
!HSNOW
        CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
        CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
        CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
        write (*, "(a,a13,a,a)") " recieve and write ", "HSNOW", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, hsnow_id, mpi_io_buf_adjust_1d, &
                          start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='HSNOW')
#endif

#ifdef SEAICE_VARIABLE_SALINITY
!HSALT
        CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_restart)
        CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
        CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
        write (*, "(a,a13,a,a)") " recieve and write ", "HSALT", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, hsalt_id, mpi_io_buf_adjust_1d, &
                          start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='HSALT')
#endif
!UICE
        CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_restart)
        CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
        CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
        write (*, "(a,a13,a,a)") " recieve and write ", "UICE", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, uice_id, mpi_io_buf_adjust_1d, &
                          start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='UICE')
!VICE
        CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
        CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
        CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all,           &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
        write (*, "(a,a13,a,a)") " recieve and write ", "VICE", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, vice_id, mpi_io_buf_adjust_1d, &
                          start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='VICE')
!seaice_sigma1/2/12
      !if (SEAICEuseEVP)
        CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
        CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
        CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all,           &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
        write (*, "(a,a13,a,a)") " recieve and write ", "seaice_sigma1", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, seaice_sigma1_id, mpi_io_buf_adjust_1d, &
                          start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='seaice_sigma1')

        CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
        CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
        CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all,           &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
        write (*, "(a,a13,a,a)") " recieve and write ", "seaice_sigma2", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, seaice_sigma2_id, mpi_io_buf_adjust_1d, &
                          start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='seaice_sigma2')

        CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_restart, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_restart)
        CALL MPI_WAIT(mpi_req_restart, mpi_stat, mpi_err)
        CALL mpi_1d_data_adjust(mpi_io_buf_1d_restart, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all,           &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
        write (*, "(a,a13,a,a)") " recieve and write ", "seaice_sigma12", " in ", trim(filename)
        call netcdf_check(nf90_put_var(ncid, seaice_sigma12_id, mpi_io_buf_adjust_1d, &
                          start=(/1, ts/), count=(/mpi_total_nlpb, 1/)), trim(filename), varname='seaice_sigma12')


      endif    !endif for mitice_on

      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine mpi_csp_io_restart_write

end module mod_mpi_csp_io_restart
