module mod_csp_init
   use mod_csp_misc
   use mod_io_netcdf
   use mod_csp_dyn_phi
   use mod_csp_thdyn_force
   use mod_csp_vdiff_coef_gls
   use mod_csp_dyn_forward
   use mod_csp_tide
   use mod_csp_force, only : csp_force_init, EmPmRrevise
# if defined (PTIDE)
   use mod_csp_dyn_ptide
# endif
   use mod_mpi_reallocate
   use mod_mpi_csp_io
   !YY
   use mitice_init
   use mitice_parameters, only : NxtDumpFile, dumpFileIntv
   implicit none

contains

!==============================================================================
   subroutine csp_init
!==============================================================================
      implicit none
      integer :: i, k, w, s, bj
      integer :: i1, i2, i3, i4
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4
      character(lc) :: filename, domid_str
      real(dp) :: rhotemp, globalArea_temp, temp
      real(wp) :: mpi_tmp
      real(wp), allocatable, dimension(:,:) :: RcolWS
      real(wp) :: RcolZOpen
      real(wp) :: lat_ini, lon_ini, t0, dis, disc

      if (wp == dp) then
         avoid0 = 1.0e-20_wp
      else
         avoid0 = 1.0e-12_wp
      end if

      if (ln_bous) then
         gravityRau0 = 1.0_wp
         gravityDynForce = r1_rau0
      else
         gravityRau0 = gravity*rau0
         gravityDynForce = gravity
      end if

      rn_emin = gravityRau0*gravityRau0*1.0e-7_wp    ! minimum value of tke [Pa^2/s^2]
      rn_emin0 = gravityRau0*gravityRau0*1.0e-4_wp    ! surface minimum value of tke [Pa^2/s^2]
      rn_psimin = gravityRau0*gravityRau0*1.0e-12_wp   ! minimum value of eps [Pa^2/s^3]

      call csp_init_allocate

      call csp_misc_rho_init

      call csp_thdyn_force_qsr_init

      write (domid_str, "(i2.2)") idom
      filename = "nmefc_macom_csp_grid_d"//trim(domid_str)//".nc"

      IF (mpi_rank == 0) THEN

         if ((rhoref) .and. (.not. ln_bous)) then
            call netcdf_read(filename, 'drC', drC, d1=nkp1)
            CALL MPI_BCAST(drC, nkp1, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
            call netcdf_read(filename, 'drF', drF, d1=nk)
            CALL MPI_BCAST(drF, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
            call netcdf_read(filename, 'rC', rC, d1=nk)
            CALL MPI_BCAST(rC, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
            call netcdf_read(filename, 'rC_z', rC_z, d1=nk)
            CALL MPI_BCAST(rC_z, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
            call netcdf_read(filename, 'rF', rF, d1=nkp1)
            CALL MPI_BCAST(rF, nkp1, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
            call netcdf_read(filename, 'rhoLev', rhoLev, d1=nk)
            CALL MPI_BCAST(rhoLev, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         else
            call netcdf_read(filename, 'drC_z', drC, d1=nkp1)
            CALL MPI_BCAST(drC, nkp1, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
            call netcdf_read(filename, 'drF_z', drF, d1=nk)
            CALL MPI_BCAST(drF, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
            call netcdf_read(filename, 'rC_z', rC, d1=nk)
            CALL MPI_BCAST(rC, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
            call netcdf_read(filename, 'rC_z', rC_z, d1=nk)
            CALL MPI_BCAST(rC_z, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
            call netcdf_read(filename, 'rF_z', rF, d1=nkp1)
            CALL MPI_BCAST(rF, nkp1, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         end if

         call netcdf_read(filename, 'rF_z', rF_z, d1=nkp1)
         CALL MPI_BCAST(rF_z, nkp1, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)

      ELSE

         CALL MPI_BCAST(drC, nkp1, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         CALL MPI_BCAST(drF, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         CALL MPI_BCAST(rC, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         CALL MPI_BCAST(rC_z, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         CALL MPI_BCAST(rF, nkp1, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         if ((rhoref) .and. (.not. ln_bous)) then
            CALL MPI_BCAST(rhoLev, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         end if
         CALL MPI_BCAST(rF_z, nkp1, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
      END IF

      if (mpi_rank == 0) then
         if ((rhoref) .and. (ln_bous)) then
            call netcdf_read(filename, 'rhoLev', rhoLev, d1=nk)
            call MPI_BCAST(rhoLev, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         end if
      else
         if ((rhoref) .and. (ln_bous)) then
            call MPI_BCAST(rhoLev, nk, mpi_real_wp, mpi_root_comp, mpi_comp_comm, mpi_err)
         end if
      end if

      ! find nksr, below which light can not penetration, about 400m in general
      do k = nk, 1, -1
         if (rC_z(k) <= 400.0_wp) nksr = k
      end do
      if (nksr == 1) nksr = 2
      IF (mpi_comp_rank == 0) write (run_out_unit, *) 'nksr                = ', nksr

      if (.not. rhoref)  then
         drF = drF * gravityRau0
         drC = drC * gravityRau0
         rF = rF * gravityRau0
         rC = rC * gravityRau0
      end if

      ! in 1D, 1 for nlpb,  2 for nlpbz
      CALL mpi_netcdf_read_exchange(filename, 'lat_c', latC, 1)
      CALL mpi_netcdf_read_exchange(filename, 'lon_c', lonC, 1)
      CALL mpi_netcdf_read_exchange(filename, 'dxC', dxC, 1)
      CALL mpi_netcdf_read_exchange(filename, 'dyC', dyC, 1)
      CALL mpi_netcdf_read_exchange(filename, 'dxW', dxW, 1)
      CALL mpi_netcdf_read_exchange(filename, 'dyW', dyW, 1)
      CALL mpi_netcdf_read_exchange(filename, 'dxS', dxS, 1)
      CALL mpi_netcdf_read_exchange(filename, 'dyS', dyS, 1)
      CALL mpi_netcdf_read_exchange(filename, 'rAc', rAc, 1)
      CALL mpi_netcdf_read_exchange(filename, 'rAs', rAs, 1)
      CALL mpi_netcdf_read_exchange(filename, 'rAw', rAw, 1)
      CALL mpi_netcdf_read_exchange(filename, 'tw', tw, 1)
      CALL mpi_netcdf_read_exchange(filename, 'ts', ts, 1)
      CALL mpi_netcdf_read_exchange(filename, 'te', te, 1)
      CALL mpi_netcdf_read_exchange(filename, 'tn', tn, 1)
      CALL mpi_netcdf_read_exchange(filename, 'tw2', tw2, 1)
      CALL mpi_netcdf_read_exchange(filename, 'ts2', ts2, 1)
      CALL mpi_netcdf_read_exchange(filename, 'te2', te2, 1)
      CALL mpi_netcdf_read_exchange(filename, 'tn2', tn2, 1)
      CALL mpi_netcdf_read_exchange(filename, 'auz1', auz1, 1)
      CALL mpi_netcdf_read_exchange(filename, 'auz2', auz2, 1)
      CALL mpi_netcdf_read_exchange(filename, 'auz3', auz3, 1)
      CALL mpi_netcdf_read_exchange(filename, 'auz4', auz4, 1)
      CALL mpi_netcdf_read_exchange(filename, 'auz5', auz5, 1)
      CALL mpi_netcdf_read_exchange(filename, 'auz6', auz6, 1)
      CALL mpi_netcdf_read_exchange(filename, 'avz1', avz1, 1)
      CALL mpi_netcdf_read_exchange(filename, 'avz2', avz2, 1)
      CALL mpi_netcdf_read_exchange(filename, 'avz3', avz3, 1)
      CALL mpi_netcdf_read_exchange(filename, 'avz4', avz4, 1)
      CALL mpi_netcdf_read_exchange(filename, 'avz5', avz5, 1)
      CALL mpi_netcdf_read_exchange(filename, 'avz6', avz6, 1)
      CALL mpi_netcdf_read_exchange(filename, 'auz1', auz1, 1)
      CALL mpi_netcdf_read_exchange(filename, 'BottomDragFac', BottomDragFac, 1)

      CALL mpi_netcdf_read_exchange(filename, 'zw', zw, 2)
      CALL mpi_netcdf_read_exchange(filename, 'zs', zs, 2)
      CALL mpi_netcdf_read_exchange(filename, 'ze', ze, 2)
      CALL mpi_netcdf_read_exchange(filename, 'zn', zn, 2)

      CALL mpi_netcdf_read_exchange(filename, 'dxZ', dxZ, 2)
      CALL mpi_netcdf_read_exchange(filename, 'dyZ', dyZ, 2)
      CALL mpi_netcdf_read_exchange(filename, 'rAz', rAz, 2)
      CALL mpi_netcdf_read_exchange(filename, 'corf', fCoriG, 2)

      CALL mpi_netcdf_read_exchange(filename, 'kSurfC', kSurfC, 1)
      CALL mpi_netcdf_read_exchange(filename, 'kSurfS', kSurfS, 1)
      CALL mpi_netcdf_read_exchange(filename, 'kSurfW', kSurfW, 1)
      CALL mpi_netcdf_read_exchange(filename, 'kTopC', kTopC, 1)
      CALL mpi_netcdf_read_exchange(filename, 'kTopS', kTopS, 1)
      CALL mpi_netcdf_read_exchange(filename, 'kTopW', kTopW, 1)

      ! in 2D, 1 for nlpb and ni, 2 for nlpbz and ni ,3 for nlpb and nk, 4 for nlpbz and nk
      CALL mpi_netcdf_read_exchange(filename, 'z1', z1, 2)
      CALL mpi_netcdf_read_exchange(filename, 'z2', z2, 2)
      CALL mpi_netcdf_read_exchange(filename, 'z3', z3, 2)
      CALL mpi_netcdf_read_exchange(filename, 'z4', z4, 2)

      CALL mpi_netcdf_read_exchange(filename, 'uw', uw, 1)
      CALL mpi_netcdf_read_exchange(filename, 'us', us, 1)
      CALL mpi_netcdf_read_exchange(filename, 'ue', ue, 1)
      CALL mpi_netcdf_read_exchange(filename, 'un', un, 1)
      CALL mpi_netcdf_read_exchange(filename, 'vw', vw, 1)
      CALL mpi_netcdf_read_exchange(filename, 'vs', vs, 1)
      CALL mpi_netcdf_read_exchange(filename, 've', ve, 1)
      CALL mpi_netcdf_read_exchange(filename, 'vn', vn, 1)
      CALL mpi_netcdf_read_exchange(filename, 'au1', au1, 1)
      CALL mpi_netcdf_read_exchange(filename, 'au2', au2, 1)
      CALL mpi_netcdf_read_exchange(filename, 'au3', au3, 1)
      CALL mpi_netcdf_read_exchange(filename, 'au4', au4, 1)
      CALL mpi_netcdf_read_exchange(filename, 'av1', av1, 1)
      CALL mpi_netcdf_read_exchange(filename, 'av2', av2, 1)
      CALL mpi_netcdf_read_exchange(filename, 'av3', av3, 1)
      CALL mpi_netcdf_read_exchange(filename, 'av4', av4, 1)

      CALL mpi_netcdf_read_exchange(filename, 'maskC', maskC, 3)
      CALL mpi_netcdf_read_exchange(filename, 'maskW', maskW, 3)
      CALL mpi_netcdf_read_exchange(filename, 'maskS', maskS, 3)
      CALL mpi_netcdf_read_exchange(filename, 'maskZ', maskZ, 4)

      ! set boundary information
      if (ln_bdy) then
         CALL mpi_netcdf_read_exchange(filename, 'maskBdy', maskBdy, 1)
         CALL mpi_netcdf_read_exchange(filename, 'normal_idx', normal_idx, 1)
         bdywidth = maxval(maskBdy)
      else
         bdywidth = 0
         maskBdy = 0
      end if
      !write(*,*) mpi_check_bdy, nlbdy, mpi_comp_rank, bdywidth

      ! here ln_bdy refer to boundary switch on global
      ! mpi_check_bdy will T only when this mpi_rank have boundary
      ! but for mpi_rank 0, it always T, since rank 0 will use for bdy data read
      if (ln_bdy .and. mpi_check_bdy) then
         allocate(bdy_pbt(nlbdy,4), bdy_uBar(nlbdy), bdy_vBar(nlbdy))
         allocate(bdy_uBcl(nlbdy,nk), bdy_vBcl(nlbdy,nk))
         allocate(bdy_t(nlbdy,nk,4), bdy_s(nlbdy,nk,4))
         allocate(bdy_u(nlbdy,nk,4), bdy_v(nlbdy,nk,4), bdy_w(nlbdy,nk,4))
      end if
      ! because there only have mpi all gather for full global array, so we have
      ! allocate uFldBcl_past and vFldBcl_past for restart output use
      if (ln_bdy) then
         allocate(uFldBcl_past(nlpb,nk), vFldBcl_past(nlpb,nk))
      end if

      CALL mpi_netcdf_read_exchange(filename, 'h0FacC', h0FacC, 3)
      CALL mpi_netcdf_read_exchange(filename, 'h0FacW', h0FacW, 3)
      CALL mpi_netcdf_read_exchange(filename, 'h0FacS', h0FacS, 3)
      CALL mpi_netcdf_read_exchange(filename, 'SideDragFac', SideDragFac, 3)

      SideDragFac = 1.0_wp
      do k = 1, nk
         do i = 1, nlpb
            if (maskC(i, k)*maskW(i,k)*maskS(i,k) .eq. 0) then
               SideDragFac(i, k) = 2.0_wp
            end if
         end do
      end do

      !! ---- OPENACC: UPDATE INDEX VARS TO GPU FOR MPI_ADJUST_INDEXES
      call GPU_INDEX_VARS_INIT
      !! ----
      ! change global index to local index in index arrays
      CALL mpi_adjust_indexes
      ! only for acc test
      ! ---- after value change on GPU card, it should update on CPU for safe
      !$acc update self(                                           &
      !$acc        tw, te, ts, tn,                                 &
      !$acc        tw2, te2, ts2, tn2, zw, ze, zs, zn,             &
      !$acc        uw, ue, us, un, vw, ve, vs, vn,                 &
      !$acc        z1, z2, z3, z4,                                 &
      !$acc        auz1, auz2, auz3, auz4, auz5, auz6,             &
      !$acc        avz1, avz2, avz3, avz4, avz5, avz6,             &
      !$acc        au1, au2, au3, au4, av1, av2, av3, av4)

      globalArea = 0.0d0
      do i = 1, loc_nlpb
         globalArea = globalArea + dble(rAc(i))* maskC(i, nk)
      end do

      CALL mpi_comp_dp_allsum(globalArea)
      globalArea = dble(real(globalArea, 4))

      hFacC = h0FacC
      hFacS = h0FacS
      hFacW = h0FacW
! -- set arrays which need calculate at start
      recip_drF = 1.0_wp / drF
      recip_drC = 1.0_wp / drC

      recip_dxC = 1.0_wp / dxC
      recip_dyC = 1.0_wp / dyC
      recip_dxZ = 1.0_wp / dxZ
      recip_dyZ = 1.0_wp / dyZ
      recip_dxW = 1.0_wp / dxW
      recip_dyW = 1.0_wp / dyW
      recip_dxS = 1.0_wp / dxS
      recip_dyS = 1.0_wp / dyS

      recip_rAc = 1.0_wp / rAc
      recip_rAs = 1.0_wp / rAs
      recip_rAw = 1.0_wp / rAw
      recip_rAz = 1.0_wp / rAz

      KappaRMCon = KappaRMCon * gravityRau0 * gravityRau0
      KappaRTCon = KappaRTCon * gravityRau0 * gravityRau0

      kappaRT = KappaRTCon
      kappaRS = KappaRTCon
      kappaRM = KappaRMCon
      KappaRTevd = KappaRTCon * 10000.0_wp

      lsmagt2 = ((dxW * dyS)/(dxW + dyS)/ rpi)** 2

      do k = 1, nk
         do i = 1, nlpb
            if (hFacC(i, k) .ne. 0) then
               recip_hFacC(i, k) = 1.0_wp/hFacC(i, k)
            else
               recip_hFacC(i, k) = 0.0_wp
            end if
            if (hFacS(i, k) .ne. 0) then
               recip_hFacS(i, k) = 1.0_wp/hFacS(i, k)
            else
               recip_hFacS(i, k) = 0.0_wp
            end if
            if (hFacW(i, k) .ne. 0) then
               recip_hFacW(i, k) = 1.0_wp/hFacW(i, k)
            else
               recip_hFacW(i, k) = 0.0_wp
            end if
         end do
      end do

      ! set constant parameter for laternal diffusion
      do i = 1, nlpbz
         i1 = z1(i, 1)
         i2 = z2(i, 1)
         lsmagf2(i) = ((dxS(i2)* dyW(i1))/(dxS(i2) + dyW(i1))/ rpi)** 2
      end do

      do i = 1, nlpb
         rSurfC(i) = rF(kSurfC(i))
         rLowC(i) = rF(kTopC(i))
         rSurfC_min(i) = rF(kSurfC(i)) - 0.3_wp*drC(kSurfC(i))

         RcolC(i) = 0.0_wp
         do k = 1, nk
            RcolC(i) = RcolC(i) + drF(k)*h0FacC(i,k)*maskC(i,k)
         end do

         if (maskC(i, nk) .eq. 0) then
            recip_RcolC(i) = 1.0_wp
            RcolC(i) = 1.0_wp
            rSurfC(i) = 2.0_wp
            rLowC(i) = 1.0_wp
         else
            recip_RcolC(i) = 1.0_wp/RcolC(i)
         end if
      end do

      refMaxDep = maxval(RcolC(1:loc_nlpb))
      CALL MPI_ALLREDUCE(refMaxDep, mpi_tmp, 1, mpi_real_wp, MPI_MAX, mpi_comp_comm, mpi_err)
      refMaxDep = mpi_tmp

      ! set RcolW and RcolS
      do i = 1, nlpb
         RcolW(i) = 0.0_wp
         do k = 1, nk
            RcolW(i) = RcolW(i) + drF(k)*h0FacW(i,k)*maskW(i,k)
         end do

         if (maskW(i, nk) .eq. 0) then
            recip_RcolW(i) = 1.0_wp
         else
            recip_RcolW(i) = 1.0_wp/RcolW(i)
         end if

         RcolS(i) = 0.0_wp
         do k = 1, nk
            RcolS(i) = RcolS(i) + drF(k)*h0FacS(i,k)*maskS(i,k)
         end do

         if (maskS(i, nk) .eq. 0) then
            recip_RcolS(i) = 1.0_wp
         else
            recip_RcolS(i) = 1.0_wp/RcolS(i)
         end if
      end do

      ! set RcolZ
      allocate(RcolWS(nlpb, 2))
      do i = 1, nlpb
         RcolWS(i, iu) = RcolW(i)
         RcolWS(i, iv) = RcolS(i)
      end do

      do i = 1, nlpbz
         i1 = z1(i, 1)
         i2 = z2(i, 1)
         i3 = z3(i, 1)
         i4 = z4(i, 1)
         j1 = z1(i, 2)
         j2 = z2(i, 2)
         j3 = z3(i, 2)
         j4 = z4(i, 2)
         k1 = z1(i, 3)
         k2 = z2(i, 3)
         k3 = z3(i, 3)
         k4 = z4(i, 3)

         RcolZOpen = min(RcolWS(i1, j1)*abs(k1), RcolWS(i2, j2)*abs(k2))
         RcolZOpen = min(RcolZOpen, RcolWS(i3, j3)*abs(k3))
         if (k4 .ne. 0) then
            RcolZOpen = min(RcolZOpen, RcolWS(i4, j4)*abs(k4))
         end if
         RcolZ(i) = RcolZOpen
      end do

      do i = 1, nlpbz
         if (RcolZ(i) .ne. 0) then
            recip_RcolZ(i) = 1.0_wp/RcolZ(i)
         else
            recip_RcolZ(i) = 0.0_wp
         end if
      end do
      deallocate(RcolWS)

      ! set label for barotropic process
      nkbar = nk+5

      nIterSurf = int(dTtracer/dTsurf)
      iterSurf_st = nIterSurf/2+1
      iterSurf_en = nIterSurf/2 + nIterSurf
      iterSurfs = iterSurf_en - iterSurf_st + 1

      rSurfW(1:nlpb) = rSurfC(1:nlpb)
      rSurfS(1:nlpb) = rSurfC(1:nlpb)
      rLowW(1:nlpb) = rLowC(1:nlpb)
      rLowS(1:nlpb) = rLowC(1:nlpb)

      if (ln_bous) then
         Bo_surf(1:nlpb) = gravity
      else
         Bo_surf(1:nlpb) = r1_rau0
      end if
      rBo_surf(1:nlpb) = 1.0_wp/Bo_surf(1:nlpb)

      rStarFacC(1:nlpb) = 1.0_wp
      rStarFacW(1:nlpb) = 1.0_wp
      rStarFacS(1:nlpb) = 1.0_wp

      do k = 1, nk
         tFld(1:nlpb, k, 1) = ((16.-12.*TANH((rC_z(k) - 400)/700)) &
                               *(-TANH((500.-rC_z(k))/150.) + 1.)/2. &
                               +(15.*(1.-TANH((rC_z(k) - 50.)/1500.)) &
                                 - 1.4*TANH((rC_z(k) - 100.)/100.) &
                                 + 7.*(1500.-rC_z(k))/1500.) &
                               *(-TANH((rC_z(k) - 500.)/150.) + 1.)/2.)

         tFld(1:nlpb, k, 2) = ((36.25 - 1.13*TANH((rC_z(k) - 305)/460)) &
                               *(-TANH((500.-rC_z(k))/150.) + 1.)/2 &
                               + (35.55 + 1.25*(5000.-rC_z(k))/5000. &
                                  -1.62*TANH((rC_z(k) - 60.)/650.) &
                                  + 0.2*TANH((rC_z(k) - 35.)/100.) &
                                  + 0.2*TANH((rC_z(k) - 1000.)/5000.)) &
                               *(-TANH((rC_z(k) - 500.)/150.) + 1.)/2)
      end do

      do k = 1, nk
         tFld(1:nlpb, k, 1) = (((7.5 - 1.0*abs(latC(1:nlpb))/30.)* &
                                (1.-tanh((rC_z(k) - 80.)/30.)) + 10.*(5000.-rC_z(k))/5000.))
      end do

      !tFld(:,:,1) = 10.0_wp
      !tFld(:,:,2) = 35.0_wp
      uFld(1:nlpb, :) = 0.0_wp
      vFld(1:nlpb, :) = 0.0_wp
      wFld(1:nlpb, :) = 0.0_wp
      etaH(1:nlpb) = 0.0_wp

      if (.not. restart_in) then
         filename = "nmefc_macom_csp_initinal_d"//trim(domid_str)//".nc"
      ! in 2D, 1 for nlpb and ni, 2 for nlpbz and ni ,3 for nlpb and nk, 4 for nlpbz and nk
         if (ln_bous) then
            CALL mpi_netcdf_read_exchange(filename, 'ssh', etaH, 1, optional_ts = 1)
            pbt_base = 0.0_wp
            ssh = etaH
         else
            CALL mpi_netcdf_read_exchange(filename, 'pbt', etaH, 1, optional_ts = 1)
            if (ln_pbt_base .eq. 1) then
               CALL mpi_netcdf_read_exchange('nmefc_macom_pbt_base.nc', 'pbt_base', pbt_base, 1)
            else
               pbt_base = 0.0_wp
            end if
         end if
         !$ACC update device(etaH, ssh, pbt_base)

         CALL mpi_netcdf_read_exchange(filename, 'u', uFld, 3, optional_ts = 1)
         CALL mpi_netcdf_read_exchange(filename, 'v', vFld, 3, optional_ts = 1)
         !! pay attention, the second dimension of wFld starts from 2
         CALL mpi_netcdf_read_exchange(filename, 'w', wFld(1:nlpb,2:nkp1), 3, optional_ts = 1)
         CALL mpi_netcdf_read_exchange(filename, 't', tFld(:,:,1), 3, optional_ts = 1)
         CALL mpi_netcdf_read_exchange(filename, 's', tFld(:,:,2), 3, optional_ts = 1)

         if (.not. rhoref) then
            !$ACC update device(nk,nlpb,rC,tFld,maskC)
            call csp_misc_rhoInSitu
            !$ACC update self(rhoInSitu)

            do k = 1, nk
               rhotemp = 0.0d0
               globalArea_temp = 0.0d0
               do i = 1,loc_nlpb
                  globalArea_temp = globalArea_temp + dble(rAc(i))*maskC(i, k)
                  rhotemp = rhotemp + dble((rhoInSitu(i, k) + rau0))*dble(rAc(i))*maskC(i, k)
               end do

               CALL mpi_comp_dp_allsum(rhotemp)
               rhotemp = dble(real(rhotemp, 4))

               CALL mpi_comp_dp_allsum(globalArea_temp)
               globalArea_temp = max(dble(real(globalArea_temp, 4)), 1.0d0)

               rhoLev(k) = rhotemp/globalArea_temp
            end do
         end if

         do k = 1, nk
            if (rhoLev(k) > 0.0_wp) then
               rhoLev(k) = 1.0_wp/rhoLev(k)
            else
               rhoLev(k) = r1_rau0
            end if
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

         ! calculate pbt_base from initial field
         if (ln_pbt_base .eq. 2) then
            !$ACC update device(nk,nlpb,rC,tFld,maskC,rhoLev,kSurfC,rSurfC,drF)
            call csp_misc_rhoInSitu
            !$ACC kernels present(phiHydC,phiHydF)
            !$ACC loop
            do i = 1, nlpb
               phiHydC(i) = 0.0_wp
               phiHydF(i) = 0.0_wp
            end do
            !$ACC end kernels

            do k = 1, nk
              call csp_dyn_phi_hyd(k)
            end do

            !$ACC kernels present(phiHydC, phiHydF, ssh, pbt_base)
            !$ACC loop independent
            do i = 1, nlpb
               ssh(i) = phiHydC(i)/gravity
               pbt_base(i) = - phiHydF(i) * rau0
               phiHydC(i) = 0.0_wp
               phiHydF(i) = 0.0_wp
            end do
            !$ACC end kernels
         end if

         !$ACC kernels present(etaH, pbt_base)
         !$ACC loop
         do i = 1, nlpb
            etaH(i) = etaH(i) - pbt_base(i)
         end do
         !$ACC end kernels

         !$ACC update self(ssh, pbt_base, etaH)

      end if

      if (ln_bdy .and. mpi_check_bdy) then
         do i = 1, nlbdy
            bj = mpi_idx_bdy(i)

            do k = 1, nk
               bdy_u(i, k, :) = uFld(bj, k)
               bdy_v(i, k, :) = vFld(bj, k)
               bdy_w(i, k, :) = wFld(bj, k)
               bdy_t(i, k, :) = tFld(bj, k, 1)
               bdy_s(i, k, :) = tFld(bj, k, 2)
            end do
            bdy_pbt(i, :) = etaH(bj)

            bdy_uBar(i) = 0.0_wp
            bdy_vBar(i) = 0.0_wp
            do k = 1, nk
               bdy_uBar(i) = bdy_uBar(i) + bdy_u(i, k, 4)*drF(k)*recip_RcolW(bj)
               bdy_vBar(i) = bdy_vBar(i) + bdy_v(i, k, 4)*drF(k)*recip_RcolS(bj)
            end do

            do k = 1, nk
               bdy_uBcl(i,k) = bdy_u(i, k, 4) - bdy_uBar(i)
               bdy_vBcl(i,k) = bdy_v(i, k, 4) - bdy_vBar(i)
            end do
         end do
      end if

      if (ln_bdy) then
         uBar(:) = 0.0_wp
         vBar(:) = 0.0_wp
         do k = 1, nk
            do i = 1, nlpb
               uBar(i) = uBar(i) + uFld(i, k)*drF(k)*maskW(i,k)*recip_RcolW(i)
               vBar(i) = vBar(i) + vFld(i, k)*drF(k)*maskS(i,k)*recip_RcolS(i)
            end do
         end do

         do k = 1, nk
            uFldBcl_past(:,k) = uFld(:,k) - uBar(:)
            vFldBcl_past(:,k) = vFld(:,k) - vBar(:)
         end do
      end if

      ! here nowjulian from call csp_init_allocate
      ! endjulian will use in sbc read, here may have bug when sbc_cli = .true., need revised in future !!!
      endjulian = nowjulian + int(nIterMax*dTtracer, 8)

      IF (mpi_rank == 0) THEN
         if (mod(86400, int(dTtracer)) .ne. 0) then
            write (run_out_unit, *) 'ERROR: seconds of one day (86400) must divide exactly by dTtracer'
            stop
         end if
         if (mod(int(nOutFreq), int(dTtracer)) .ne. 0) then
            write (run_out_unit, *) 'ERROR: nOutFreq must divide exactly by dTtracer'
            stop
         end if
         if (fuvavg) then
            if (mod(ffuv, 2) .ne. 0) then
               write (run_out_unit, *) 'ERROR: ffuv must divide exactly by 2 when fuvavg is true'
               stop
            end if
         end if
         if (ftqavg) then
            if (mod(fftq, 2) .ne. 0) then
               write (run_out_unit, *) 'ERROR: fftq must divide exactly by 2 when ftqavg is true'
               stop
            end if
         end if
         if (fradavg) then
            if (mod(ffrad, 2) .ne. 0) then
               write (run_out_unit, *) 'ERROR: ffrad must divide exactly by 2 when fradavg is true'
            end if
         end if
         if (fprcavg) then
            if (mod(ffprc, 2) .ne. 0) then
               write (run_out_unit, *) 'ERROR: ffprc must divide exactly by 2 when fprcavg is true'
            end if
         end if
         if (fslpavg) then
            if (mod(ffslp, 2) .ne. 0) then
               write (run_out_unit, *) 'ERROR: ffslp must divide exactly by 2 when fslpavg is true'
            end if
         end if
         if (frunoffavg) then
            if (mod(ffrunoff, 2) .ne. 0) then
               write (run_out_unit, *) 'ERROR: ffrunoff must divide exactly by 2 when frunoffavg is true'
            end if
         end if
      END IF

      do i = 1, nlpb
         u10(i, :) = tFld(i, nk, 1)*0.01_wp
         v10(i, :) = tFld(i, nk, 1)*0.01_wp
         t10(i, :) = tFld(i, nk, 1) + 273.0_wp
         q10(i, :) = max(tFld(i, nk, 1)*0.0001_wp, 0.0001_wp)
         swdn(i, :) = max(tFld(i, nk, 1)*20_wp, 1.0_wp)
         lwdn(i, :) = max(tFld(i, nk, 1)*10_wp, 100.0_wp)
         prec(i, :) = max(tFld(i, nk, 1)*0.00001_wp, 0.0001_wp)
         snow(i, :) = 0.0_wp
         slp(i,:) = 101325.0_wp

         u10(i, 3) = 0.0_wp
         v10(i, 3) = 0.0_wp
         t10(i, 3) = 0.0_wp
         q10(i, 3) = 0.0_wp
         swdn(i, 3) = 0.0_wp
         lwdn(i, 3) = 0.0_wp
         prec(i, 3) = 0.0_wp
         snow(i, 3) = 0.0_wp
         slp(i, 3) = 0.0_wp

# if defined (NOSLP)
         phi0surf(i) = slp(i, 4)
# endif

      end do

      CALL mpi_csp_io_comp_init

      if (abs(zuv - zqt) < 0.01_wp) l_zt_equal_zu = .true.    ! testing "zu == zt" is risky with double precision

      call csp_vdiff_coef_gls_init

      if (ln_tide) then
         !! add tide initial
         call csp_tide_init

         !$ACC update device(tide_ssh, tide_u, tide_v, tide_freq_use, &
         !$ACC               tide_amp_ssh, tide_pha_ssh, tide_u_ssh, tide_v_ssh, tide_f_ssh, &
         !$ACC               tide_amp_u, tide_pha_u, tide_u_curu, tide_v_curu, tide_f_curu, &
         !$ACC               tide_amp_v, tide_pha_v, tide_u_curv, tide_v_curv, tide_f_curv)
         call csp_tide_update
         !$ACC update self(tide_u, tide_v, tide_ssh)

         !do k = 1, nk
         !   uFld(:,k) = uFld(:,k) + tide_u(:)*maskW(:,k)
         !   vFld(:,k) = vFld(:,k) + tide_v(:)*maskS(:,k)
         !end do
         etaH(:) = etaH(:) + tide_ssh(:)*gravityRau0*maskC(:,nk)
         ssh(:) = ssh(:) + tide_ssh(:)
      end if

# if defined (PTIDE)
      call csp_dyn_ptide_init
# endif

      surf_out_ratio = nOutFreq/dTsurf
      mom_out_ratio = nOutFreq/dTmom
      tracer_out_ratio = nOutFreq/dTtracer

      !YY,ZY: fcori at C-point. seaice module using fCori at grid center, not at vorticity point
      !YY: upload fCori in subroutine GPU_SEAICE_VARS_UPDATE
      fCori = 2.0_wp*(rpi/43082.0_wp)*sin(latC*radPi)

      call csp_force_init
      
   end subroutine csp_init

!==============================================================================
   subroutine csp_init_allocate
!==============================================================================
      implicit none
      character(lc) :: filename, domid_str

      IF (mpi_comp_rank == 0) THEN
! -- set scalar parameters
         write (domid_str, "(i2.2)") idom
         filename = "nmefc_macom_csp_grid_d"// trim(domid_str)//".nc"
         ! nlpb is used for real number of grids in each process, mpi_total_nlpb is global nlpb

         call netcdf_read(trim(filename), 'nlpb', mpi_total_nlpb)
         call netcdf_read(trim(filename), 'nl', mpi_total_nl)
         ! nlpbz is used for real number of grids in each process, mpi_total_nlpbz is global nlpbz
         call netcdf_read(trim(filename), 'nlpbz', mpi_total_nlpbz)
         call netcdf_read(trim(filename), 'nlz', mpi_total_nlz)
         if (ln_bdy) then
            call netcdf_read(trim(filename), 'nlbdy', mpi_total_nlbdy)
         else
            mpi_total_nlbdy = 0
         end if
         call netcdf_read(trim(filename), 'nk', nk)
         call netcdf_read(trim(filename), 'nkp1', nkp1)
         call netcdf_read(trim(filename), 'ni', ni)

         CALL mpi_netcdf_root_bcast_first

      ELSE

         CALL mpi_netcdf_nonroot_bcast_first

      END IF

      CALL mpi_index_graph_init_and_send

      abFac = 0.5_wp + abEps

      nowyear = startyear
      nowmon = startmon
      nowday = startday
      nowhour = starthour
      nowmin = startmin
      nowsec = startsec

      call greg2jul(nowsec, nowmin, nowhour, nowday, nowmon, nowyear, nowjulian)
      !YY
      startTime = nowjulian


      IF (mpi_rank == 0) THEN
         write (run_out_unit, *) 'model run parameters for domain ', idom
         write (run_out_unit, *) 'Total domain number is ', ndomain
         if (ln_bous) then
            write (run_out_unit, *) 'model will run in volume conservation mode'
         else
            write (run_out_unit, *) 'model will run in mass conservation mode'
         end if
         write (run_out_unit, *) 'nl                  = ', mpi_total_nl
         write (run_out_unit, *) 'nlpb                = ', mpi_total_nlpb
         write (run_out_unit, *) 'nlz                 = ', mpi_total_nlz
         write (run_out_unit, *) 'nlpbz               = ', mpi_total_nlpbz
         write (run_out_unit, *) 'nk                  = ', nk
         write (run_out_unit, *) 'nkp1                = ', nkp1
         write (run_out_unit, *) 'ni                  = ', ni
         write (run_out_unit, *) 'nsteps              = ', nIterMax - nIter0 + 1
         write (run_out_unit, *) 'parent              = ', parent
         write (run_out_unit, *) 'nRatio              = ', nRatio
         write (run_out_unit, *) 'nIter0              = ', nIter0
         write (run_out_unit, *) 'nIterMax            = ', nIterMax
         write (run_out_unit, *) 'nOutFreq            = ', nOutFreq
         write (run_out_unit, *) 'fOutFreq            = ', trim(fOutFreq)
         write (run_out_unit, *) 'nResFreq            = ', nResFreq
         write (run_out_unit, *) 'startyear           = ', nowyear
         write (run_out_unit, *) 'startmon            = ', nowmon
         write (run_out_unit, *) 'startday            = ', nowday
         write (run_out_unit, *) 'starthour           = ', nowhour
         write (run_out_unit, *) 'startmin            = ', nowmin
         write (run_out_unit, *) 'startsec            = ', nowsec
         write (run_out_unit, *) 'seconds from 1949-10-01 00:00:00 is ', nowjulian
         write (run_out_unit, *) 'dTtracer            = ', dTtracer
         write (run_out_unit, *) 'dTmom               = ', dTmom
         write (run_out_unit, *) 'dTsurf              = ', dTsurf
         write (run_out_unit, *) 'sbc_cli             = ', sbc_cli
         write (run_out_unit, *) 'sbc_cyc             = ', sbc_cyc
         write (run_out_unit, *) 'ffuv                = ', ffuv
         write (run_out_unit, *) 'fftq                = ', fftq
         write (run_out_unit, *) 'ffrad               = ', ffrad
         write (run_out_unit, *) 'ffprc               = ', ffprc
         write (run_out_unit, *) 'ffslp               = ', ffslp
         write (run_out_unit, *) 'ffrunoff            = ', ffrunoff
         write (run_out_unit, *) 'nffuv               = ', trim(nffuv)
         write (run_out_unit, *) 'nfftq               = ', trim(nfftq)
         write (run_out_unit, *) 'nffrad              = ', trim(nffrad)
         write (run_out_unit, *) 'nffprc              = ', trim(nffprc)
         write (run_out_unit, *) 'nffslp              = ', trim(nffslp)
         write (run_out_unit, *) 'nffrunoff           = ', trim(nffrunoff)
         write (run_out_unit, *) 'fuvavg              = ', fuvavg
         write (run_out_unit, *) 'ftqavg              = ', ftqavg
         write (run_out_unit, *) 'fradavg             = ', fradavg
         write (run_out_unit, *) 'fprcavg             = ', fprcavg
         write (run_out_unit, *) 'fslpavg             = ', fslpavg
         write (run_out_unit, *) 'frunoffavg          = ', frunoffavg
         write (run_out_unit, *) 'zqt                 = ', zqt
         write (run_out_unit, *) 'zuv                 = ', zuv
         write (run_out_unit, *) 'fsbc_dir            = ', trim(fsbc_dir)
         write (run_out_unit, *) 'fout_dir            = ', trim(fout_dir)
         write (run_out_unit, *) 'restart_in          = ', restart_in
         write (run_out_unit, *) 'restart_out         = ', restart_out
         write (run_out_unit, *) 'date_res            = ', date_res

         write (run_out_unit, *) 'rhoConst            = ', rau0
         write (run_out_unit, *) 'rhoa                = ', rhoa
         write (run_out_unit, *) 'gravity             = ', gravity
         write (run_out_unit, *) 'ViscAhDCon          = ', ViscAhDCon
         write (run_out_unit, *) 'ViscAhZCon          = ', ViscAhZCon
         write (run_out_unit, *) 'ViscA4DCon          = ', ViscA4DCon
         write (run_out_unit, *) 'ViscA4ZCon          = ', ViscA4ZCon
         write (run_out_unit, *) 'harmonic            = ', harmonic
         write (run_out_unit, *) 'biharmonic          = ', biharmonic
         write (run_out_unit, *) 'bottomDrag          = ', BottomDrag
         write (run_out_unit, *) 'bottomDragMax       = ', BottomDragMax
         write (run_out_unit, *) 'KappaRMCon          = ', KappaRMCon
         write (run_out_unit, *) 'abEps               = ', abEps
         write (run_out_unit, *) 'rStarFacLow         = ', rStarFacLow
         write (run_out_unit, *) 'rStarFacUp          = ', rStarFacUp
         write (run_out_unit, *) 'ntracer             = ', ntracer
         write (run_out_unit, *) 'rn_abs              = ', rn_abs
         write (run_out_unit, *) 'rn_si0              = ', rn_si0
         write (run_out_unit, *) 'KappaRTCon          = ', KappaRTCon
         write (run_out_unit, *) 'epsmin              = ', epsmin
         write (run_out_unit, *) 'harmonicT           = ', harmonicT
         write (run_out_unit, *) 'biharmonicT         = ', biharmonicT
         write (run_out_unit, *) 'diffKhCon           = ', diffKhCon
         write (run_out_unit, *) 'diffK4Con           = ', diffK4Con

         if (ln_bous .and. (ln_pbt_base .ne. 0)) then
            write (run_out_unit, *) 'Error: pbt_base only activate when model run in mass conservation'
            stop
         end if
      END IF

! -- allocate arrays for geometric
      allocate (drF(nk), drC(nkp1))
      allocate (rF(nkp1), rC(nk), rC_z(nk))
      !YY, ZY: add rF_z ??
      allocate (rF_z(nkp1))
      allocate (fCoriG(nlpbz),fCori(nlpb))
      allocate (latC(nlpb), lonC(nlpb))
      allocate (dxC(nlpb), dxS(nlpb), dxW(nlpb), dxZ(nlpbz))
      allocate (dyC(nlpb), dyS(nlpb), dyW(nlpb), dyZ(nlpbz))
      allocate (rAc(nlpb), rAs(nlpb), rAw(nlpb), rAz(nlpbz))
      allocate (tw(nlpb), ts(nlpb), te(nlpb), tn(nlpb))
      allocate (zw(nlpbz), zs(nlpbz), ze(nlpbz), zn(nlpbz))
      allocate (tw2(nlpb), ts2(nlpb), te2(nlpb), tn2(nlpb))
      allocate (auz1(nlpb), auz2(nlpb), auz3(nlpb), auz4(nlpb), auz5(nlpb), auz6(nlpb))
      allocate (avz1(nlpb), avz2(nlpb), avz3(nlpb), avz4(nlpb), avz5(nlpb), avz6(nlpb))
      allocate (z1(nlpbz, ni), z2(nlpbz, ni), z3(nlpbz, ni), z4(nlpbz, ni))
      allocate (uw(nlpb, ni), us(nlpb, ni), ue(nlpb, ni), un(nlpb, ni))
      allocate (vw(nlpb, ni), vs(nlpb, ni), ve(nlpb, ni), vn(nlpb, ni))
      allocate (au1(nlpb, ni), au2(nlpb, ni), au3(nlpb, ni), au4(nlpb, ni))
      allocate (av1(nlpb, ni), av2(nlpb, ni), av3(nlpb, ni), av4(nlpb, ni))
      allocate (kSurfC(nlpb), kSurfW(nlpb), kSurfS(nlpb))
      allocate (kTopC(nlpb), kTopW(nlpb), kTopS(nlpb))
      allocate (hFacC(nlpb, nk), hFacS(nlpb, nk))
      allocate (hFacW(nlpb, nk))
      allocate (hFacZ(nlpbz), recip_hFacZ(nlpbz))
      allocate (h0FacC(nlpb, nk), h0FacS(nlpb, nk))
      allocate (h0FacW(nlpb, nk))
      allocate (maskC(nlpb, nk), maskS(nlpb, nk))
      allocate (maskW(nlpb, nk), maskZ(nlpbz, nk))
      allocate (maskBdy(nlpb), normal_idx(nlpb))
      allocate (recip_drF(nk), recip_drC(nkp1))
      allocate (recip_dxC(nlpb), recip_dxZ(nlpbz))
      allocate (recip_dyC(nlpb), recip_dyZ(nlpbz))
      allocate (recip_dxW(nlpb), recip_dxS(nlpb))
      allocate (recip_dyW(nlpb), recip_dyS(nlpb))
      allocate (recip_rAc(nlpb), recip_rAs(nlpb))
      allocate (recip_rAw(nlpb), recip_rAz(nlpbz))
      allocate (recip_RcolC(nlpb), recip_RcolZ(nlpbz))
      allocate (recip_RcolW(nlpb), recip_RcolS(nlpb))
      allocate (RcolC(nlpb), RcolW(nlpb), RcolS(nlpb), RcolZ(nlpbz))
      allocate (recip_hFacC(nlpb, nk), recip_hFacS(nlpb, nk))
      allocate (recip_hFacW(nlpb, nk))
! -- allocate array for sea surface force
      allocate (taum(nlpb), utau(nlpb), vtau(nlpb), utau_out(loc_nlpb), vtau_out(loc_nlpb))
      utau_out = 0.0_wp
      vtau_out = 0.0_wp
      taum = 0.0_wp
      utau = 0.0_wp
      vtau = 0.0_wp
      allocate (guExt(nlpb), gvExt(nlpb))
      allocate (phi0surf(nlpb))
      !YY,ZY see mitice_dynsolver
      phi0surf = 0.0_wp
      allocate (EmPmR(nlpb), EmPmR_out(loc_nlpb), EmPmR_long(loc_nlpb))
      EmPmR = 0.0_wp
      EmPmR_long = 0.0_wp
      EmPmR_out = 0.0_wp
      allocate (addMass(nlpb, nk))
      addMass = 0.0_wp
! --allocate array for time forward

! --allocate array for dynamics
      ! allocate array for baroclinic-pressure-velocity coupling
      allocate (omega3(nlpbz))
      omega3 = 0.0_wp
      allocate (vort3(nlpbz))
      vort3 = 0.0_wp
      allocate (hDiv(nlpb))
      hDiv = 0.0_wp
      allocate (KE(nlpb))
      KE = 0.0_wp
      allocate (dKEdx(nlpb), dKEdy(nlpb))
      dKEdx = 0.0_wp
      dKEdy = 0.0_wp
      allocate (uCorTerm(nlpb), vCorTerm(nlpb))
      uCorTerm = 0.0_wp
      vCorTerm = 0.0_wp
      allocate (phiHydC(nlpb), phiHydF(nlpb))
      phiHydC = 0.0_wp
      phiHydF = 0.0_wp
      allocate (dPhiHydX(nlpb), dPhiHydY(nlpb))
      dPhiHydX = 0.0_wp
      dPhiHydY = 0.0_wp
      allocate (rhoLev(nk))
      rhoLev = r1_rau0
      allocate (uShearTerm(nlpb), vShearTerm(nlpb))
      uShearTerm = 0.0_wp
      vShearTerm = 0.0_wp
      allocate (uFld(nlpb, nk), vFld(nlpb, nk), wFld(nlpb, nkp1))
      allocate (uFldBcl(nlpb, nk), vFldBcl(nlpb, nk))
      allocate (uFldBclDash(nlpb, nk), vFldBclDash(nlpb, nk))
      uFldBcl = 0.0_wp
      vFldBcl = 0.0_wp
      uFldBclDash = 0.0_wp
      vFldBclDash = 0.0_wp
      allocate (uFld_out(loc_nlpb, nk), vFld_out(loc_nlpb, nk), wFld_out(loc_nlpb, nkp1))
      uFld = 0.0_wp
      vFld = 0.0_wp
      wFld = 0.0_wp
      uFld_out = 0.0_wp
      vFld_out = 0.0_wp
      wFld_out = 0.0_wp
      allocate (gUnow(nlpb, nk), gVnow(nlpb, nk))
      gUnow = 0.0_wp
      gVnow = 0.0_wp
      allocate (gUpast(nlpb, nk), gVpast(nlpb, nk))
      gUpast = 0.0_wp
      gVpast = 0.0_wp
      allocate (uBar(nlpb), vBar(nlpb))
      allocate (uBarPast(nlpb), vBarPast(nlpb))
      allocate (uBarCouple(nlpb), vBarCouple(nlpb))
      allocate (gUBarnow(nlpb), gVBarnow(nlpb))
      gUBarnow = 0.0_wp
      gVBarnow = 0.0_wp
      allocate (gUBarpast(nlpb), gVBarpast(nlpb))
      gUBarpast = 0.0_wp
      gVBarpast = 0.0_wp

      ! allocate array for barotropic - pressure - velocity coupling
      allocate (etaN(nlpb), etaH(nlpb), etaH_out(loc_nlpb), pbt_base(nlpb))
      etaH_out = 0.0_wp
      pbt_base = 0.0_wp
      allocate (ssh(nlpb), ssh_out(loc_nlpb))
      ssh = 0.0_wp
      ssh_out = 0.0_wp
      allocate (dphiSurfX(nlpb), dphiSurfY(nlpb))
      dphiSurfX = 0.0_wp
      dphiSurfY = 0.0_wp
      allocate (Bo_surf(nlpb), rBo_surf(nlpb))
      allocate (rSurfC(nlpb), rSurfW(nlpb), rSurfS(nlpb))
      allocate (rLowC(nlpb), rLowW(nlpb), rLowS(nlpb))
      allocate (rSurfC_min(nlpb))
      allocate (rStarFacC(nlpb), rStarFacW(nlpb), rStarFacS(nlpb))
      rStarFacC = 1.0_wp
      rStarFacW = 1.0_wp
      rStarFacS = 1.0_wp
      allocate (rStarExpC(nlpb))
      rStarExpC = 1.0_wp
      ! allocate array for horizontal momentum dissipation
      allocate (del2u(nlpb), del2v(nlpb))
      del2u = 0.0_wp
      del2v = 0.0_wp
      allocate (dStar(nlpb), zStar(nlpbz))
      dStar = 0.0_wp
      zStar = 0.0_wp
      allocate (uDissip(nlpb), vDissip(nlpb))
      uDissip = 0.0_wp
      vDissip = 0.0_wp
      allocate (viscAh_D(nlpb), viscAh_Z(nlpbz))
      allocate (viscA4_D(nlpb), viscA4_Z(nlpbz))
      allocate (lsmagt2(nlpb), lsmagf2(nlpbz))
      allocate (KappaRM(nlpb, nkp1))
      allocate (KappaRM_out(loc_nlpb, nkp1))
      KappaRM_out = 0.0_wp
      ! allocate array for momentum drag
      allocate (SideDragFac(nlpb, nk), BottomDragFac(nlpb))
      BottomDragFac = 1.0_wp
      allocate (uDragTerms(nlpb), vDragTerms(nlpb))
      uDragTerms = 0.0_wp
      vDragTerms = 0.0_wp
      allocate (uDragBottom(nlpb), vDragBottom(nlpb))
      uDragBottom = 0.0_wp
      vDragBottom = 0.0_wp
! --allocate array for thermodynamic
      ! allocate array for tracer
      allocate (alphaRho(nlpb))
      allocate (rhoInSitu(nlpb, nk))
      rhoInSitu = 0.0_wp
      allocate (tFld(nlpb, nk, ntracer))
      allocate (tFld_out(loc_nlpb, nk, ntracer))
      tFld_out = 0.0_wp
      allocate (gTracer(nlpb, nk), tracer(nlpb, nk))
      allocate (gTracerPast(nlpb, nk, ntracer))
      allocate (FirstTracer(ntracer))
      FirstTracer = .true.
      allocate (diffKhu(nlpb, nk), diffKhv(nlpb, nk))
      allocate (diffK4u(nlpb, nk), diffK4v(nlpb, nk))
      ! allocate array for horizontal tracer diffusion

      ! allocate array for vertical tracer diffusion
      allocate (rn2(nlpb, nkp1))
      allocate (KappaRT(nlpb, nkp1), KappaRS(nlpb, nkp1))
      allocate (KappaRT_out(loc_nlpb, nkp1))
      KappaRT_out = 0.0_wp
      allocate (hmxl(nlpb, nkp1))
      hmxl = 0.05_wp*gravity*rau0
      allocate (en(nlpb, nkp1))
      en = rn_emin

      ! allocate array for tracer external force
      allocate (gtExt(nlpb, nk))
      allocate (qsr(nlpb), qns(nlpb), qsr_out(loc_nlpb), qns_out(loc_nlpb))
      gtExt = 0.0_wp
      qsr = 0.0_wp
      qns = 0.0_wp
      qsr_out = 0.0_wp
      qns_out = 0.0_wp
      !YY,ZY: add for sea ice. definitions refer to mod_csp_forcing
      !later should move this to the seaice module ??
      allocate (Qnet(nlpb),Qsw(nlpb),saltflux(nlpb))
      Qnet = 0.0_wp
      Qsw = 0.0_wp
      saltflux = 0.0_wp

      ! allocate array for force
      allocate (u10(nlpb, 4), v10(nlpb, 4))
      allocate (t10(nlpb, 4), q10(nlpb, 4))
      allocate (lwdn(nlpb, 4), swdn(nlpb, 4))
      allocate (prec(nlpb, 4), snow(nlpb, 4))
      allocate (slp(nlpb, 4), runoff(nlpb, 4))
      allocate (sfrt(nlpb, 4), sfrs(nlpb, 4))
      u10 = 0.0_wp
      v10 = 0.0_wp
      t10 = 0.0_wp
      q10 = 0.0_wp
      lwdn = 0.0_wp
      swdn = 0.0_wp
      prec = 0.0_wp
      snow = 0.0_wp
      slp = 0.0_wp
      runoff = 0.0_wp
      sfrt = 0.0_wp
      sfrs = 0.0_wp
      !YY,ZY: add zevap, which used to be a local var in mod_csp_force
      allocate(zevap(nlpb))
      zevap = 0.0_wp

      ! array for long time tracer diagnosis
      allocate(tFld_adv(loc_nlpb, nk), sFld_adv(loc_nlpb, nk), tFld_vdiff(loc_nlpb, nk), sFld_vdiff(loc_nlpb, nk))
      allocate(tFld_hdiff(loc_nlpb, nk), sFld_hdiff(loc_nlpb, nk), tFld_force(loc_nlpb, nk), sFld_force(loc_nlpb, nk))
      tFld_adv = 0.0_wp
      sFld_adv = 0.0_wp
      tFld_vdiff = 0.0_wp
      sFld_vdiff = 0.0_wp
      tFld_hdiff = 0.0_wp
      sFld_hdiff = 0.0_wp
      tFld_force = 0.0_wp
      sFld_force = 0.0_wp

   end subroutine csp_init_allocate

!==============================================================================
   subroutine GPU_INDEX_VARS_INIT
!==============================================================================
      ! -- array for geometric indexes
      !$ACC update device(                                         &
      !$ACC        tw, te, ts, tn,                                 &
      !$ACC        tw2, te2, ts2, tn2, zw, ze, zs, zn,             &
      !$ACC        uw, ue, us, un, vw, ve, vs, vn,                 &
      !$ACC        z1, z2, z3, z4,                                 &
      !$ACC        auz1, auz2, auz3, auz4, auz5, auz6,             &
      !$ACC        avz1, avz2, avz3, avz4, avz5, avz6,             &
      !$ACC        au1, au2, au3, au4, av1, av2, av3, av4)

      !SCALAR VARS REGARDING DIMENTSIONS
      !$ACC update device(loc_nlpb, nlpb, nlpbz, nl, nlz, nk, nkp1, nksr)
   end subroutine GPU_INDEX_VARS_INIT

!==============================================================================
   subroutine GPU_VARS_INIT
!==============================================================================
      !$ACC UPDATE DEVICE(                                                       &
      !$ACC drF, drC, rF, rC, rC_z, fCoriG, latC, lonC,                          &
      !$ACC dxC, dxZ, dxW, dxS, dyC, dyZ, dyW, dyS, rAc, rAs, rAw, rAz,          &
      !$ACC kSurfC, kSurfW, kSurfS, kTopC, kTopW, kTopS,                         &
      !$ACC hFacC, hFacS, hFacW, hFacZ, recip_hFacZ, h0FacC, h0FacS, h0FacW,     &
      !$ACC maskC, maskS, maskW, maskZ, maskBdy, normal_idx,                     &
      !$ACC recip_drF, recip_drC, recip_dxC, recip_dxZ, recip_dyC, recip_dyZ,    &
      !$ACC recip_dxW, recip_dxS, recip_dyW, recip_dyS, recip_rAc, recip_rAs,    &
      !$ACC recip_rAw, recip_rAz, recip_RcolC, recip_RcolZ,                      &
      !$ACC recip_RcolW, recip_RcolS, RcolC, RcolW, RcolS, RcolZ,                &
      !$ACC recip_hFacC, recip_hFacS, recip_hFacW,                               &
      !!  COPY dynamic array from CPU to GPU
      !$ACC omega3, vort3, hDiv, KE, dKEdx, dKEdy, uCorTerm, vCorTerm,           &
      !$ACC phiHydC, phiHydF, dPhiHydX, dPhiHydY, rhoLev, uShearTerm, vShearTerm,&
      !$ACC uFld, vFld, wFld, uFldBcl, vFldBcl, uFldBcl_past, vFldBcl_past,      &
      !$ACC uFldBclDash, vFldBclDash, uFld_out, vFld_out, wFld_out,              &
      !$ACC gUnow, gVnow, gUpast, gVpast, uBar, vBar, uBarPast, vBarPast,        &
      !$ACC uBarCouple, vBarCouple, gUBarnow, gVBarnow, gUBarpast, gVBarpast,    &
      !$ACC etaN, etaH, etaH_out, pbt_base, ssh, ssh_out, dphiSurfX, dphiSurfY,  &
      !$ACC Bo_surf, rBo_surf, rSurfC, rSurfW, rSurfS, rLowC, rLowW, rLowS,      &
      !$ACC rSurfC_min, rStarFacC, rStarFacW, rStarFacS, rStarExpC,              &
      !$ACC del2u, del2v, dStar, zStar, uDissip, vDissip, viscAh_D, viscAh_Z,    &
      !$ACC viscA4_D, viscA4_Z, lsmagt2, lsmagf2, KappaRM, KappaRM_out,          &
      !$ACC SideDragFac, BottomDragFac, uDragTerms, vDragTerms, uDragBottom,     &
      !$ACC vDragBottom, alphaRho, rhoInSitu, tFld, tFld_out, gTracer,           &
      !$ACC gTracerPast, tracer, FirstTracer, diffKhu, diffKhv, diffK4u,         &
      !$ACC diffK4v, rn2, KappaRT, KappaRS, KappaRT_out, KappaRX, hmxl, en,      &
      !!  COPY surface boundary force array from CPU to GPU
      !$ACC taum, utau, vtau, utau_out, vtau_out, guExt, gvExt, phi0surf,        &
      !$ACC EmPmR, EmPmR_long, EmPmR_out, addMass, u10, v10, t10, q10, lwdn, swdn, prec,     &
      !$ACC snow, slp, runoff, sfrt, sfrs, gtExt, qsr, qsr_out, qns, qns_out,    &
      !!  COPY lateral boundary force array from CPU to GPU
      !$ACC bdy_uBar, bdy_vBar, bdy_uBcl, bdy_vBcl, bdy_pbt, bdy_t, bdy_s,       &
      !$ACC bdy_u, bdy_v, bdy_w, tFld_adv, sFld_adv, tFld_vdiff, sFld_vdiff,     &
      !$ACC tFld_hdiff, sFld_hdiff, tFld_force, sFld_force)

      !YY,ZY: modification due to inclusion of sea ice module
      !$ACC update device(zevap)

      ! here may have a compiler BUG, so the insignificance code for avoid it
      !$ACC update self(etaH_out, ssh_out, u10, v10, t10)

      !SCALAR VARS
      !$ACC update device(                                                       &
      !$ACC ndomain, idom, parent, nRatio, nl, nlpb, loc_nlpb, nlz, nlpbz, nk,   &
      !$ACC nksr, nlbdy, bdywidth, loc_nlbdy, sbct, nkp1, ni, dTtracer, dTmom,   &
      !$ACC dTsurf, nIterSurf, iterSurf_st, iterSurf_en, iterSurfs, nIter0,      &
      !$ACC nIterMax, nOutFreq, fOutFreq, nResFreq, startyear, startmon,         &
      !$ACC startday, starthour, startmin, startsec, nowyearf, nowmonf,          &
      !$ACC nowdayf, nowhourf, nowminf, nowsecf, nowyear, nowmon, nowday,        &
      !$ACC nowhour, nowmin, nowsec, nowjulian, endjulian, dateNodalJulian,      &
      !$ACC restart_out, date_res, output_nmuber, ln_bdy, ln_tide, ln_bous,      &
      !$ACC ln_pbt_base, myIter, globalArea, abFac, FirstMomU, FirstMomV,        &
      !$ACC gammas, gammat, ffsfrs, ffsfrt, nffsfrs, nffsfrt, EmPmRrevise,       &
      !$ACC fsfrsavg, fsfrtavg, sbc_cli, sbc_cyc, ffuv, fftq, ffrad,             &
      !$ACC ffprc, ffslp, ffrunoff, nffuv, nfftq, nffrad, nffprc, nffslp,        &
      !$ACC nffrunoff, fuvavg, ftqavg, fradavg, fprcavg, fslpavg, frunoffavg, ln_qdew,  &
      !$ACC nxtjultq, nxtjuluv, nxtjulrad, nxtjulprc, nxtjulslp, nxtjulrunoff,   &
      !$ACC nxtjulsfrs, nxtjulsfrt, rhoref, fbvec, fbt, fbs, fbpbt, nxtjulVEC,   &
      !$ACC nxtjulT, nxtjulS, nxtjulPbt, nfbvec, nfbt, nfbs, nfbpbt, bvecavg,    &
      !$ACC btavg, bsavg, bpbtavg, zqt, zuv, ViscAhDCon, ViscAhZCon, ViscA4DCon, &
      !$ACC ViscA4ZCon, harmonic, biharmonic, BottomDrag, BottomDragMax,         &
      !$ACC KappaRMCon, abEps, rStarFacLow, rStarFacUp, refMaxDep, nkbar,        &
      !$ACC rn_emin, rn_emin0, rn_psimin, gravityRau0, gravityDynForce, ntracer, &
      !$ACC ntide, tide2d, tideDim, rn_abs, rn_si0, KappaRTCon, KappaRTevd,      &
      !$ACC harmonicT, biharmonicT, diffKhCon, diffK4Con, surf_out_ratio,        &
      !$ACC mom_out_ratio, tracer_out_ratio, avoid0, l_zt_equal_zu)

   end subroutine GPU_VARS_INIT

!==============================================================================
   subroutine csp_init_step
!==============================================================================
      implicit none
      integer :: i
      integer :: FileDate_prev

      !$acc kernels present(nlpb, bdy_pbt, phiHydC, phiHydF, uBarCouple, vBarCouple, ssh)
      !$acc loop
      do i = 1, nlpb
         phiHydC(i) = 0.0_wp
         phiHydF(i) = 0.0_wp

         uBarCouple(i) = 0.0_wp
         vBarCouple(i) = 0.0_wp

         ssh(i) = 0.0_wp
      end do

      if (ln_bdy .and. mpi_check_bdy) then
         !$ACC loop independent
         do i = 1, nlbdy
            bdy_pbt(i,4) = 0.0_wp
         end do
      end if

      !$acc end kernels

      nowjulian = nowjulian + int(dTtracer, 8)

!YY nowTime used in seaice module
      nowTime = nowjulian - startTime

      if (mitice_on) then
        select case (dumpFileIntv)
          case (1)     !generate new ice state dump file when new month comes
            FileDate_prev = nowmon
          case (2)     !new year
            FileDate_prev = nowyear
          case default
            WRITE (*, *) "wrong options for dumpFileIntv, default is 1 or 2"
          STOP
        end select
      endif

      call jul2greg(nowsec, nowmin, nowhour, nowday, nowmon, nowyear, nowjulian)

      if (mitice_on) then
      !YY: this logical is for seaice dump files which can be named sequentially with regard to month/year
        select case (dumpFileIntv)
          case (1)     !generate new ice state dump file when new month comes
            if (FileDate_prev.NE.nowmon) NxtDumpFile = .True.
          case (2)     !new year
            if (FileDate_prev.NE.nowyear) NxtDumpFile = .True.
          case default
            WRITE (*, *) "wrong options for dumpFileIntv, default is 1 or 2"
          STOP
        end select
      endif

   end subroutine csp_init_step

!==============================================================================
  subroutine csp_init_asm
!==============================================================================
     implicit none
     character(lc) :: filename, domid_str
     real(wp) :: work3(nlpb, nk), work2(nlpb)

     work2 = 0.0_wp
     work3 = 0.0_wp

     write (domid_str, "(i2.2)") idom
     filename = "nmefc_mcom_assim_gain_d"//trim(domid_str)//".nc"

     CALL mpi_netcdf_read_exchange(filename, 'pbt', work2, 1, optional_ts = 1)
     etaH(1:nlpb) = etaH(1:nlpb) + work2(1:nlpb)

     CALL mpi_netcdf_read_exchange(filename, 'u', work3, 3, optional_ts = 1)
     uFld(1:nlpb, 1:nk) = uFld(1:nlpb, 1:nk) + work3(1:nlpb, 1:nk)

     CALL mpi_netcdf_read_exchange(filename, 'v', work3, 3, optional_ts = 1)
     vFld(1:nlpb, 1:nk) = vFld(1:nlpb, 1:nk) + work3(1:nlpb, 1:nk)

     CALL mpi_netcdf_read_exchange(filename, 'w', work3, 3, optional_ts = 1)
     wFld(1:nlpb,2:nkp1) = wFld(1:nlpb,2:nkp1) + work3(1:nlpb,1:nk)*gravity*rau0

     CALL mpi_netcdf_read_exchange(filename, 't', work3, 3, optional_ts = 1)
     tFld(1:nlpb, 1:nk, 1) = tFld(1:nlpb, 1:nk, 1) + work3(1:nlpb, 1:nk)

     CALL mpi_netcdf_read_exchange(filename, 's', work3, 3, optional_ts = 1)
     tFld(1:nlpb, 1:nk, 2) = tFld(1:nlpb, 1:nk, 2) + work3(1:nlpb, 1:nk)

  end subroutine csp_init_asm

end module mod_csp_init
