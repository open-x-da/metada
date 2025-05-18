module mod_csp_force
   use mod_misc_basic
   use mod_csp_basic
   use mod_csp_misc

   implicit none
   real(wp) :: rn_vfac = 1.0_wp                         ! multiplication factor for ice/ocean velocity in the calculation of wind stress (clem)
   real(wp) :: rn_pfac = 1.0_wp                         ! multiplication factor for precipitation
   real(wp) :: rn_efac = 1.0_wp                         ! multiplication factor for evaporation (clem)
   real(wp) :: zcoef_qsatw = 0.98_wp*640380.0_wp/rhoa
   integer, save, public          :: EmPmR_cycle = 1
   real(wp), save, public         :: EmPmRrevise = 0.0_wp
   integer, save, public          :: EmPmRcount = 0

   !$ACC declare copyin(rn_vfac, rn_pfac, rn_efac, zcoef_qsatw)
   !$ACC declare create(EmPmRrevise)

contains

!==============================================================================
   subroutine csp_force_init
!==============================================================================
      implicit none

      EmPmR_cycle = EmPmR_cycle*86400/int(dTtracer)

   end subroutine csp_force_init

!==============================================================================
   subroutine csp_force_core
!==============================================================================
      !!---------------------------------------------------------------------
      !! ** Come from NEMO
      !! ** Purpose :   provide the momentum, heat and freshwater fluxes at
      !!      the ocean surface at each time step
      !!
      !! ** Method  :   CORE bulk formulea for the ocean using atmospheric
      !!      fields read in sbc_read
      !!
      !! ** Outputs : - utau    : i-component of the stress at U-point  (N/m2)
      !!              - vtau    : j-component of the stress at V-point  (N/m2)
      !!              - taum    : Wind stress module at T-point         (N/m2)
      !!              - wndm    : Wind speed module at T-point          (m/s)
      !!              - qsr     : Solar heat flux over the ocean        (W/m2)
      !!              - qns     : Non Solar heat flux over the ocean    (W/m2)
      !!              - emp     : evaporation minus precipitation       (kg/m2/s)
      !!
      !!  ** Nota  :   sf has to be a dummy argument for AGRIF on NEC
      !!---------------------------------------------------------------------
      implicit none
      integer :: i
      real(wp) :: zwnd_i(nlpb)                     ! wind speed relative current at U-point (i-component) (m/s)
      real(wp) :: zwnd_j(nlpb)                     ! wind speed relative current at V-point (j-component) (m/s)
      real(wp) :: wndm(nlpb)                       ! Wind speed module at T-point          (m/s)
      real(wp) :: zqlw(nlpb), zqsb(nlpb)           ! long wave and sensible heat fluxes
      !real(wp) :: zqla(nlpb), zevap(nlpb)          ! latent heat fluxes and evaporation
      !YY,ZY: zevap is now defined as global var
      real(wp) :: zqla(nlpb)                       ! latent heat fluxes
      real(wp) :: zt_zu(nlpb)                      ! air temperature at wind speed height
      real(wp) :: zq_zu(nlpb)                      ! air spec. hum.  at wind speed height
      real(wp) :: zztmp                            ! local variable
      real(wp) :: zqsatw(nlpb)                     ! specific humidity at pst
      real(wp) :: Cd(nlpb)                         ! transfer coefficient for momentum      (tau)
      real(wp) :: Ch(nlpb)                         ! transfer coefficient for sensible heat (Q_sens)
      real(wp) :: Ce(nlpb)                         ! tansfert coefficient for evaporation   (Q_lat)
      real(wp) :: zst(nlpb)                        ! surface temperature in Kelvin
      real(wp) :: vapre                            ! vapor pressure in mb
      real(wp) :: slptmp                           ! surface pressure im mb
      real(dp) :: sshtemp, gloavg1, gloavg2


!YY, ZY: zevap removed from 'create', now it is a global var
      !$ACC data create(zwnd_i,zwnd_j,wndm,zqlw,zqsb,zqla,zt_zu,zq_zu,zqsatw,Cd,Ch,Ce,zst),   &
      !$ACC present(nlpb,nk,uFld,vFld,tFld,u10,v10,utau,vtau,taum,qns,EmPmR,zqt,zuv,t10,q10, &
      !$ACC         qsr,lwdn,swdn,slp,maskC,l_zt_equal_zu,prec,runoff,snow,ln_qdew,  &
      !$ACC         rn_vfac,rn_pfac,rn_efac,zcoef_qsatw,  &
      !$ACC         Qnet,Qsw,zevap)

      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         ! convert SST from Celcius to Kelvin (and set minimum value far above 0 K)
         zst(i) = tFld(i, nk, 1) + rt0

         zwnd_i(i) = u10(i, 4) - rn_vfac*uFld(i, nk)
         zwnd_j(i) = v10(i, 4) - rn_vfac*vFld(i, nk)

         ! ... scalar wind ( = | U10m - U_oce | ) at T-point (masked)
         ! This is a very rough version, we should consider improve in future
         !YY: wndm should be averaged over C-grid?
         wndm(i) = sqrt(zwnd_i(i)*zwnd_i(i) + zwnd_j(i)*zwnd_j(i))

         ! ----------------------------------------------------------------------------- !
         !      I   Radiative FLUXES                                                     !
         ! ----------------------------------------------------------------------------- !
         ! We should consider include the diurnal cycle in future
         ! Short Wave
         !YY,ZY: here I should use a if statement to distinguish between openwater and ice albedo
         qsr(i) = (1.0_wp - albo)*swdn(i, 4)*maskC(i, nk)

         ! Long  Wave
!YY: zqlw is downward?  in MIT, lwup is scaled by ocean/ice/snow emissivity
!YY: is it nesscary here to add ocean_emissivity
!YY: following is from MITGCM
!C     lwup = emissivity*stefanBoltzmann*Tsrf^4 + (1-emissivity)*lwdown
!C     the second terms is the reflected incoming long wave radiation
!C     so that the net upward long wave heat flux is:
! lwflux(I,J,bi,bj) = - lwdownLoc(I,J)*SEAICE_emissivity+ SEAICE_emissivity*SEAICE_boltzmann*tsurfLoc(I,J)**4
         zqlw(i) = (lwdn(i, 4) - Stef*zst(i)*zst(i)*zst(i)*zst(i))*maskC(i, nk)

         ! ----------------------------------------------------------------------------- !
         !     II    Turbulent FLUXES                                                    !
         ! ----------------------------------------------------------------------------- !
         ! ... specific humidity at SST and IST
         zqsatw(i) = zcoef_qsatw*exp(-5107.4_wp/zst(i))
      end do
      !$ACC end kernels

      if (ln_qdew) then
         ! convert dew point temperature (K) to specific humidity (kg/kg)
         !$ACC kernels
         !$ACC loop independent private(zztmp,slptmp,vapre)
         do i = 1, nlpb
            zztmp = q10(i,4) - rt0  ! dew in degree C
            slptmp = slp(i,4)*0.01_wp ! surface pressure in mb
            vapre = 6.112_wp*exp((17.67_wp*zztmp)/(zztmp + 243.5_wp))
            q10(i,4) = maskC(i,nk)*(0.622_wp*vapre)/(slptmp - (0.378_wp*vapre))
        end do
         !$ACC end kernels
      end if

      ! ... NCAR Bulk formulae, computation of Cd, Ch, Ce at T-point :
      call csp_force_turb_core_2z(zqt, zuv, zst, t10(1:nlpb, 4), zqsatw, q10(1:nlpb, 4), wndm, &
                                  Cd, Ch, Ce, zt_zu, zq_zu)

      !$ACC kernels
      !$ACC loop independent private(zztmp)
      do i = 1, nlpb
         ! ... tau module, i and j component
         ! Both wndm Cd in T point, we need change zwnd_i and zwnd_j from U V point into T point in future
         zztmp = rhoa*wndm(i)*Cd(i)
         taum(i) = zztmp*wndm(i)
         utau(i) = zztmp*zwnd_i(i)
         vtau(i) = zztmp*zwnd_j(i)

         ! ... add the HF tau contribution to the wind stress module?
         ! we can consider this in future
         ! if( lhftau ) then
         !   taum(:,:) = taum(:,:) + sf(jp_tdif)%fnow(:,:,1)
         ! end if

         ! ... utau, vtau at U- and V_points, resp.
         !     Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
         !     Note the use of max(tmask(i,j),tmask(i+1,j) is to mask tau over ice shelves
         ! Since we don't consider U V in different point, so zwnd_i and zwnd_j still in their point,
         ! in future we should consider this
         ! utau(ji) = 0.5_wp * ( 2.0_wp - umask(ji,nk) ) * ( zwnd_i(ji) + zwnd_i(jj  ) )   &
         !                      * max(tmask(ji,jj,1),tmask(ji+1,jj,nk))
         ! vtau(ji) = 0.5_wp * ( 2.0_wp - vmask(ji,nk) ) * ( zwnd_j(ji) + zwnd_j(jj+1) )   &
         !                      * max(tmask(ji,jj,1),tmask(ji,jj+1,nk))

         !  Turbulent fluxes over ocean
         ! -----------------------------
         !! q_air and t_air are (or "are almost") given at 10m (wind reference height)
         !! q_air and t_air are not given at 10m (wind reference height)
         ! Values of temp. and hum. adjusted to height of wind during bulk algorithm iteration must be used
         ! Even q_air and t_air are (or "are almost") given at 10m, we also use iteration, in this time zq_zu
         ! and zt_zu is same with q10 and t10
         zevap(i) = rn_efac*max(0.0_wp, rhoa*Ce(i)*(zqsatw(i) - zq_zu(i))*wndm(i))   ! Evaporation
         zqsb(i) = cpa*rhoa*Ch(i)*(zst(i) - zt_zu(i))*wndm(i)                        ! Sensible Heat
         zqla(i) = Lv*zevap(i)                                                       ! Latent Heat
      end do
      !$ACC end kernels


      ! ----------------------------------------------------------------------------- !
      !     III    Total FLUXES                                                       !
      ! ----------------------------------------------------------------------------- !
      !
      ! mass flux (evap. - precip. - runoff)

      ! force global average (prec + runoff) equal evap at every step
      !do i = 1, nlpb
      !   EmPmR(i) = (prec(i, 4)*rn_pfac + runoff(i, 4))*maskC(i, nk)
      !end do
      !
      ! sshtemp = 0.0d0
      ! do i = 1, loc_nlpb
      !   sshtemp = sshtemp + dble(EmPmR(i))* dble(rAc(i))* maskC(i, nk)
      ! end do
      ! CALL mpi_comp_dp_allsum(sshtemp)
      ! gloavg1 = sshtemp / globalArea
      ! if (mpi_rank == 0) write (*,"(a,i15,a,e14.7)") "avg of prec at date ", &
      !    nowjulian, " is ", gloavg1
      !
      !sshtemp = 0.0d0
      !do i = 1, nlpb
      !   sshtemp = sshtemp + dble(zevap(i))*dble(rAc(i))*maskC(i, nk)
      !end do
      !CALL mpi_comp_dp_allsum(sshtemp)
      !gloavg2 = sshtemp/globalArea
      !if (mpi_rank == 0) write (*,"(a,i15,a,e14.7,f10.5)") "avg of evap at date ", nowjulian, " is ", gloavg2, gloavg2/gloavg1
      !
      !do i = 1, nlpb
      !   EmPmR(i) = (zevap(i) - (prec(i, 4)*rn_pfac + runoff(i, 4))*gloavg2/gloavg1)*maskC(i, nk)
      !end do

      ! here we use wndm as rain for save memory purpose
      if (ln_rain) then
         !$ACC kernels
         !$ACC loop independent
         do i = 1, nlpb
            wndm(i) = prec(i, 4)
         end do
         !$ACC end kernels
      else
         !$ACC kernels
         !$ACC loop independent
         do i = 1, nlpb
            wndm(i) = prec(i, 4) - snow(i, 4)
         end do
         !$ACC end kernels
      end if

      if (mitice_on) then
        !$ACC kernels
        !$ACC loop independent
        do i = 1, nlpb
           EmPmR(i) = (zevap(i) - (snow(i, 4)*rn_pfac + wndm(i)*rn_pfac + runoff(i, 4)))*maskC(i, nk)
           qns(i) = zqlw(i) - zqsb(i) - zqla(i) &                            ! Downward Non Solar
                    - zevap(i)*tFld(i, nk, 1)*rcp &                          ! remove evap heat content at SST
                    + wndm(i)*rn_pfac*(t10(i, 4) - rt0)*rcp &                ! add liquid precip (rain) heat content at Tair
                    + snow(i, 4)*rn_pfac*(min(t10(i, 4), rt0_snow) - rt0)*cpic*maskC(i, nk)  ! add solid precip heat content at min(Tair,Tsnow)
              !YY: snow's unit is m/s??
!YY: Add vars Qnet and Qsw in MaCOM for seaice
!YY: Definition of Qnet and Qsw is the same as MITgcm, which is positive upward.
!----  YY: definition of MITgcm. Note EmPmP unit, converted by rhofreshwater
!C     EmPmR :: Net upward freshwater flux in kg/m2/s
!C              EmPmR = Evaporation - precipitation - runoff
!C              > 0 for increase in salt (ocean salinity)
!C              Typical range: -1e-4 < EmPmR < 1e-4
!C              Southwest C-grid tracer point
!C           NOTE: for backward compatibility EmPmRfile is specified in
!C                 m/s when using external_fields_load.F.  It is converted
!C                 to kg/m2/s by multiplying by rhoConstFresh.
!C
!C  saltFlux :: Net upward salt flux in g/kg.kg/m^2/s = g/m^2/s
!C              flux of Salt taken out of the ocean per time unit (second).
!C              Note: only used when salty sea-ice forms or melts.
!C              > 0 for decrease in SSS.
!C              Southwest C-grid tracer point
!C
!C     Qnet  :: Net upward surface heat flux (including shortwave) in W/m^2
!C              Qnet = latent + sensible + net longwave + net shortwave
!C              > 0 for decrease in theta (ocean cooling)
!C              Typical range: -250 < Qnet < 600
!C              Southwest C-grid tracer point
!C
!C     Qsw   :: Net upward shortwave radiation in W/m^2
!C              Qsw = - ( downward - ice and snow absorption - reflected )
!C              > 0 for decrease in theta (ocean cooling)
!C              Typical range: -350 < Qsw < 0
!C              Southwest C-grid tracer point
!----
           Qnet(i) = -qns(i) - qsr(i)        !including sw radiation (upward)
           Qsw (i) = -qsr(i)
        end do
        !$ACC end kernels
      else
        !$ACC kernels
        !$ACC loop independent
        do i = 1, nlpb
          EmPmR(i) = (zevap(i) - (snow(i, 4)*rn_pfac + wndm(i)*rn_pfac + runoff(i, 4)))*maskC(i, nk)
          qns(i) = zqlw(i) - zqsb(i) - zqla(i) &                            ! Downward Non Solar
                   - snow(i, 4)*rn_pfac*lfus &                              ! remove latent melting heat for solid precip
                   - zevap(i)*tFld(i, nk, 1)*rcp &                          ! remove evap heat content at SST
                   + wndm(i)*rn_pfac*(t10(i, 4) - rt0)*rcp &                ! add liquid precip (rain) heat content at Tair
                   + snow(i, 4)*rn_pfac*(min(t10(i, 4), rt0_snow) - rt0)*cpic*maskC(i, nk)  ! add solid precip heat content at min(Tair,Tsnow)
        enddo
        !$ACC end kernels
      endif
      !$ACC end data
   end subroutine csp_force_core

!==============================================================================
   subroutine csp_force_turb_core_2z(zt, zu, sst, T_zt, q_sat, q_zt, dU, Cd, Ch, Ce, T_zu, q_zu)
!==============================================================================
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_core  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!
      !! ** Method : Monin Obukhov Similarity Theory
      !!             + Large & Yeager (2004,2008) closure: CD_n10 = f(U_n10)
      !!
      !! ** References :   Large & Yeager, 2004 / Large & Yeager, 2008
      !!
      !! ** Last update: Laurent Brodeau, June 2014:
      !!    - handles both cases zt=zu and zt/=zu
      !!    - optimized: less 2D arrays allocated and less operations
      !!    - better first guess of stability by checking air-sea difference of virtual temperature
      !!       rather than temperature difference only...
      !!    - added function "cd_neutral_10m" that uses the improved parametrization of
      !!      Large & Yeager 2008. Drag-coefficient reduction for Cyclone conditions!
      !!    - using code-wide physical constants defined into "phycst.mod" rather than redifining them
      !!      => 'vkarmn' and 'grav'
      !!----------------------------------------------------------------------
      integer, parameter :: nb_itt = 5       ! number of itterations
      real(wp), intent(in) :: zt             ! height for T_zt and q_zt                    [m]
      real(wp), intent(in) :: zu             ! height for dU                               [m]
      real(wp), intent(in) :: sst(nlpb)      ! sea surface temperature                [Kelvin]
      real(wp), intent(in) :: T_zt(nlpb)     ! potential air temperature              [Kelvin]
      real(wp), intent(in) :: q_sat(nlpb)    ! sea surface specific humidity           [kg/kg]
      real(wp), intent(in) :: q_zt(nlpb)     ! specific air humidity                   [kg/kg]
      real(wp), intent(in) :: dU(nlpb)       ! relative wind module at zu                [m/s]
      real(wp), intent(out) :: Cd(nlpb)       ! transfer coefficient for momentum         (tau)
      real(wp), intent(out) :: Ch(nlpb)       ! transfer coefficient for sensible heat (Q_sens)
      real(wp), intent(out) :: Ce(nlpb)       ! transfert coefficient for evaporation   (Q_lat)
      real(wp), intent(out) :: T_zu(nlpb)     ! air temp. shifted at zu                     [K]
      real(wp), intent(out) :: q_zu(nlpb)     ! spec. hum.  shifted at zu               [kg/kg]
      !
      integer :: i, j_itt
      !
      real(wp) :: U_zu          ! relative wind at zu                            [m/s]
      real(wp) :: Ce_n10        ! 10m neutral latent coefficient
      real(wp) :: Ch_n10        ! 10m neutral sensible coefficient
      real(wp) :: sqrt_Cd_n10   ! root square of Cd_n10
      real(wp) :: sqrt_Cd       ! root square of Cd
      real(wp) :: zeta_u        ! stability parameter at height zu
      real(wp) :: zeta_t        ! stability parameter at height zt
      real(wp) :: zpsi_h_u, zpsi_m_u
      real(wp) :: ztmp0, ztmp1, ztmp2
      real(wp) :: stab          ! 1st stability test integer
      real(wp) :: psi_h
      real(wp) :: rgt33, zw6, X2, X, stabit
      !!----------------------------------------------------------------------

      !$ACC kernels present(nlpb,zt,zu,sst,T_zt,q_sat,q_zt,dU,Cd,Ch,Ce,T_zu,q_zu,l_zt_equal_zu)
      !$ACC loop independent private(U_zu,Ce_n10,Ch_n10,sqrt_Cd_n10,sqrt_Cd,   &
      !$ACC                          zeta_u,zeta_t,zpsi_h_u,zpsi_m_u,ztmp0,    &
      !$ACC                          ztmp1,ztmp2,stab,zw6,rgt33,X2,X,stabit,psi_h)
      do i = 1, nlpb
         ! relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s
         U_zu = max(0.5_wp, dU(i))

         !! First guess of stability:
         ! air-sea difference of virtual pot. temp. at zt
         ztmp0 = T_zt(i)*(1.0_wp + 0.608_wp*q_zt(i)) - sst(i)*(1.0_wp + 0.608_wp*q_sat(i))
         ! stab = 1 if dTv > 0  => STABLE, 0 if unstable
         stab = 0.5_wp + sign(0.5_wp, ztmp0)

         !! Neutral coefficients at 10m:
         ! here can add wave drag in future
         ! function Cd_n10
         zw6 = U_zu*U_zu*U_zu*U_zu*U_zu*U_zu
         rgt33 = 0.5_wp + sign(0.5_wp, (U_zu - 33.0_wp))   ! If zw10 < 100. => 0, else => 1
         ztmp0 = 1.0e-3_wp*( &
                    (1.0_wp - rgt33)*(2.7_wp/U_zu + 0.142_wp + U_zu/13.09_wp - 3.14807e-10_wp*zw6) & ! zw10< 33.
                    + rgt33*2.34_wp)                                                ! zw10 >= 33.
         ztmp0 = max(ztmp0, 1.0e-6_wp)

         sqrt_Cd_n10 = sqrt(ztmp0)
         Ce_n10 = 1.0e-3_wp*(34.6_wp*sqrt_Cd_n10)
         Ch_n10 = 1.0e-3_wp*sqrt_Cd_n10*(18.0_wp*stab + 32.7_wp*(1.0_wp - stab))

         !! Initializing transf. coeff. with their first guess neutral equivalents :
         Cd(i) = ztmp0
         Ce(i) = Ce_n10
         Ch(i) = Ch_n10
         sqrt_Cd = sqrt_Cd_n10

         !! Initializing values at z_u with z_t values:
         T_zu(i) = T_zt(i)
         q_zu(i) = q_zt(i)

         !!  * Now starting iteration loop
         !$ACC loop seq
         do j_itt = 1, nb_itt
            ! Updating air/sea differences
            ztmp1 = T_zu(i) - sst(i)
            ztmp2 = q_zu(i) - q_sat(i)

            ! Updating turbulent scales :   (L&Y 2004 eq. (7))
            ztmp1 = Ch(i)/sqrt_Cd*ztmp1   ! theta*
            ztmp2 = Ce(i)/sqrt_Cd*ztmp2   ! q*

            ! virtual potential temperature at zu
            ztmp0 = T_zu(i)*(1.0_wp + 0.608_wp*q_zu(i))

            ! Estimate the inverse of Monin-Obukov length (1/L) at height zu:
            ztmp0 = (vkarmn*gravity/ztmp0*(ztmp1*(1.0_wp + 0.608_wp*q_zu(i)) &
                                                 + 0.608_wp*T_zu(i)*ztmp2))/(Cd(i)*U_zu*U_zu)   ! ( Cd*U_zu*U_zu is U*^2 at zu)

            !! Stability parameters :
            zeta_u = zu*ztmp0
            zeta_u = sign(min(abs(zeta_u), 10.0_wp), zeta_u)

            ! function psi_h
            X2 = sqrt(abs(1.0_wp - 16.0_wp*zeta_u))
            X2 = max(X2, 1.0_wp)
            X = sqrt(X2)
            stabit = 0.5_wp + sign(0.5_wp, zeta_u)
            zpsi_h_u = -5.0_wp*zeta_u*stabit &                       ! Stable
                          + (1.0_wp - stabit)*(2.0_wp*log((1.0_wp + X2)*0.5_wp))   ! Unstable

            ! function psi_m
            X2 = sqrt(abs(1.0_wp - 16.0_wp*zeta_u))
            X2 = max(X2, 1.0_wp)
            X = sqrt(X2)
            stabit = 0.5_wp + sign(0.5_wp, zeta_u)
            zpsi_m_u = -5.0_wp*zeta_u*stabit &                         ! Stable
                          + (1.0_wp - stabit)*(2.0_wp*log((1.0_wp + X)*0.5_wp) &     ! Unstable
                                               + log((1.0_wp + X2)*0.5_wp) - 2.0_wp*atan(X) + rpi*0.5_wp)  ! Unstable

            !! Shifting temperature and humidity at zu (L&Y 2004 eq. (9b-9c))
            if (.not. l_zt_equal_zu) then
               zeta_t = zt*ztmp0
               zeta_t = sign(min(abs(zeta_t), 10.0_wp), zeta_t)

               ! function psi_h
               X2 = sqrt(abs(1.0_wp - 16.0_wp*zeta_t))
               X2 = max(X2, 1.0_wp)
               X = sqrt(X2)
               stabit = 0.5_wp + sign(0.5_wp, zeta_t)
               psi_h = -5.0_wp*zeta_t*stabit &                       ! Stable
                       + (1.0_wp - stabit)*(2.0_wp*log((1.0_wp + X2)*0.5_wp))   ! Unstable

               stab = log(zu/zt) - zpsi_h_u + psi_h  ! stab just used as temp array

               T_zu(i) = T_zt(i) + ztmp1/vkarmn*stab    ! ztmp1 is still theta*
               q_zu(i) = q_zt(i) + ztmp2/vkarmn*stab    ! ztmp2 is still q*
               q_zu(i) = max(0.0_wp, q_zu(i))
            end if

            ! surface wave case, need check in future
            ! sqrt_Cd = vkarmn / ( vkarmn / sqrt_Cd_n10 - zpsi_m_u )
            ! Cd      = sqrt_Cd * sqrt_Cd

            ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 eq. 9a)...
            !   In very rare low-wind conditions, the old way of estimating the
            !   neutral wind speed at 10m leads to a negative value that causes the code
            !   to crash. To prevent this a threshold of 0.25m/s is imposed.
            ztmp0 = max(0.25_wp, U_zu/(1.0_wp + sqrt_Cd_n10/vkarmn*(log(zu/10.0_wp) - zpsi_m_u))) !  U_n10

            ! function Cd_n10
            zw6 = ztmp0*ztmp0*ztmp0*ztmp0*ztmp0*ztmp0
            rgt33 = 0.5_wp + sign(0.5_wp, (ztmp0 - 33.0_wp))   ! If zw10 < 100. => 0, else => 1
            ztmp0 = 1.0e-3_wp*((1.0_wp - rgt33)*(2.7_wp/ztmp0 + 0.142_wp + &
                      ztmp0/13.09_wp - 3.14807e-10_wp*zw6) &    ! zw10< 33.
                                  + rgt33*2.34_wp)                 ! zw10 >= 33.
            ztmp0 = max(ztmp0, 1.0e-6_wp)

            sqrt_Cd_n10 = sqrt(ztmp0)

            Ce_n10 = 1.0e-3_wp*(34.6_wp*sqrt_Cd_n10)      ! L&Y 2004 eq. (6b)
            stab = 0.5_wp + sign(0.5_wp, zeta_u)       ! update stability
            Ch_n10 = 1.0e-3_wp*sqrt_Cd_n10*(18.0_wp*stab + &
                                                  32.7_wp*(1.0_wp - stab))     ! L&Y 2004 eq. (6c-6d)

            !! Update of transfer coefficients:
            ztmp1 = 1.0_wp + sqrt_Cd_n10/vkarmn*(log(zu/10.0_wp) - zpsi_m_u)   ! L&Y 2004 eq. (10a)
            Cd(i) = ztmp0/(ztmp1*ztmp1)
            sqrt_Cd = sqrt(Cd(i))

            !
            ztmp0 = (log(zu/10.0_wp) - zpsi_h_u)/vkarmn/sqrt_Cd_n10
            ztmp2 = sqrt_Cd/sqrt_Cd_n10
            ztmp1 = 1.0_wp + Ch_n10*ztmp0
            Ch(i) = Ch_n10*ztmp2/ztmp1  ! L&Y 2004 eq. (10b)
            !
            ztmp1 = 1.0_wp + Ce_n10*ztmp0
            Ce(i) = Ce_n10*ztmp2/ztmp1  ! L&Y 2004 eq. (10c)
            !
         end do ! end j_itt iteration loop

      end do
      !$ACC end kernels

   end subroutine csp_force_turb_core_2z

!==============================================================================
   subroutine csp_force_EmPmR_revise
!==============================================================================
      implicit none 
      integer :: i
      real(dp) :: reviseTemp
      integer(i8) :: JulianYearStart, JulianYearEnd, dJulian

      !$ACC kernels present(EmPmR, EmPmRrevise)
      !$ACC loop independent
      do i = 1, nlpb 
         EmPmR(i) = (EmPmR(i) + EmPmRrevise)*maskC(i, nk)
      end do
      !$ACC end kernels

      if (EmPmRcount .lt. EmPmR_cycle) then
         !$ACC kernels present(EmPmR_long, EmPmR)
         !$ACC loop independent
         do i = 1, loc_nlpb 
            EmPmR_long(i) = EmPmR_long(i) + EmPmR(i)
         end do
         !$ACC end kernels
         EmPmRcount = EmPmRcount + 1
      end if

      if (EmPmRcount .eq. EmPmR_cycle) then 
         !$ACC update self(EmPmR_long)
         reviseTemp = 0.0d0
         do i = 1, loc_nlpb 
            reviseTemp = reviseTemp + dble(EmPmR_long(i))*dble(rAc(i))*maskC(i,nk)
         end do
         CALL mpi_comp_dp_allsum(reviseTemp)
         reviseTemp = reviseTemp / globalArea / EmPmRcount
         EmPmRrevise = EmPmRrevise - reviseTemp
         write(*,"(i3,a,e23.16,a,i8,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)") &
             mpi_rank, " EmPmRrevise = ", EmPmRrevise, " at iter ", myIter, &
            ', date ', nowyear, '-', nowmon, '-', nowday, ' ', &
            nowhour, ':', nowmin, ':', nowsec

         !$ACC update device(EmPmRrevise)

         EmPmRcount = 0
         !$ACC kernels present(EmPmR_long)
         !$ACC loop independent
         do i = 1, nlpb 
            EmPmR_long(i) = 0.0_wp
         end do
         !$ACC end kernels
      end if

   end subroutine csp_force_EmPmR_revise
end module mod_csp_force
