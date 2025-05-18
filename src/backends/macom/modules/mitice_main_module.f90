module mitice
use mitice_parameters
use mitice_vars
use mitice_init
use mitice_dynamics
use mitice_advdiff
use mitice_thermodyn
use mod_csp_basic
use mod_mpi_interfaces
use mod_mpi_variables

   implicit none

contains
   subroutine mitice_main
!     *===========================================================*
!     | o Time stepping of a dynamic/thermodynamic sea ice model. |
!     *===========================================================*
      implicit none

      INTEGER i

      if (nowTime.EQ.NxtIceTime) then
!YY,ZY: confirm if ssh can be used for sea level tilt term in the momentum eq
!YY,ZY: and how it initiated in the mod_csp_init
        !call mpi_data_exchange(mpi_sendbuf_1d, ssh)    !use etaH instead, YY, ZY
        !call mpi_data_exchange(mpi_sendbuf_1d, etaH)   ! exchanged somewhere else
!YY: is netcdf_read_exchange do the halo exchange after reading from NC ?
       ! call mpi_data_exchange(mpi_sendbuf_1d, u10(:,4))
       ! call mpi_data_exchange(mpi_sendbuf_1d, v10(:,4))

! solve ice momentum equations and update ocean surface stress u/vtau
        call mitice_dynsolver

!YY: TO DO LATER
!--   Apply ice velocity open boundary conditions
!#ifdef ALLOW_OBCS
!      IF ( useOBCS ) CALL OBCS_ADJUST_UVICE( uice, vice, myThid )
!#endif /* ALLOW_OBCS */

!Next seaice calling time
        NxtIceTime = nowTime + int(SEAICE_deltaTdyn, 8)

      endif

      if (nowTime.EQ.NxtTHDTime) then
! =====
! temporal variations of seaice thickness distribute (ITD) (LHS of ice conservation eqn)
! is contributed by horizontal advection and diffusion u*dH/dx, ridging processes, and 
! thermodynamical freezing/melting rate.
! 3-step splitting algorithm used to account for these sources.  
! =====

! advection and diffusion
        IF ( SEAICEadvHeff .OR. SEAICEadvArea .OR. SEAICEadvSnow .OR. SEAICEadvSalt ) THEN
          call mitice_adv_diff( uice, vice )
        ENDIF

!     After advection, the sea ice variables may have unphysical values
!     e.g., < 0, that are regularized here. Concentration as a special case
!     may be > 1 in convergent motion and a ridging algorithm redistributes
!     the ice to limit the concentration to 1.

        call seaice_reg_ridge

!     thermodynamics growth
        IF ( usePW79thermodynamics ) THEN
          call seaice_growth
        ENDIF
        
        NxtTHDTime = nowTime + int(SEAICE_deltaTtherm, 8)

!--   Apply ice tracer open boundary conditions
!#ifdef ALLOW_OBCS
!# ifndef DISABLE_SEAICE_OBCS
!       IF ( useOBCS ) CALL OBCS_APPLY_SEAICE( myThid )
!# endif /* DISABLE_SEAICE_OBCS */
!#endif /* ALLOW_OBCS */

!--     Update halo regions
        call mpi_data_exchange(mpi_sendbuf_1d, HEFF)
        call mpi_data_exchange(mpi_sendbuf_1d, AREA)
        call mpi_data_exchange(mpi_sendbuf_1d, HSNOW)
#ifdef SEAICE_ITD
        call mpi_data_exchange(mpi_sendbuf_2d_itd_f, HEFFITD,SEAICE_multDim)   !2d data exchange
        call mpi_data_exchange(mpi_sendbuf_2d_itd_f, AREAITD,SEAICE_multDim)
        call mpi_data_exchange(mpi_sendbuf_2d_itd_f, HSNOWITD,SEAICE_multDim)
#endif

#ifdef SEAICE_VARIABLE_SALINITY
         call mpi_data_exchange(mpi_sendbuf_1d, HSALT)
#endif

      endif

!update forcing fields modified by sea ice module
!YY: there is a problem, all forcings are the previous ice time step TO DO LATER
!$acc kernels present(seaiceMaskU,seaiceMaskV,HEFFM,utau,vtau,utau_ice,vtau_ice,qns,qns_ice,qsr,qsr_ice,EmPmR,EmPmR_ice)
!$acc loop
       do i = 1,nlpb
         if(seaiceMaskU(i).gt.0.0_wp) utau(i) = utau_ice(i)
         if(seaiceMaskV(i).gt.0.0_wp) vtau(i) = vtau_ice(i)
         if(AREA(i).gt.SEAICE_area_floor .and. HEFFM(i).gt.0.0_wp) then
           qsr(i) = qsr_ice(i)
           qns(i) = qns_ice(i)
           EmPmR(i) = EmPmR_ice(i)
         endif
       enddo
!$acc end kernels

       !YY,ZY: qns and qsr is only used in single grid, necessary to exchange??
       !====
       !call mpi_data_exchange(mpi_sendbuf_1d, EmPmR)
       !call mpi_data_exchange(mpi_sendbuf_1d, saltFlux)
       !call mpi_data_exchange(mpi_sendbuf_1d, qns)
       !call mpi_data_exchange(mpi_sendbuf_1d, qsr)
       !IF ( useRealFreshWaterFlux ) call mpi_data_exchange(mpi_sendbuf_1d, sIceLoad)
       !====



   end subroutine mitice_main

end module mitice
