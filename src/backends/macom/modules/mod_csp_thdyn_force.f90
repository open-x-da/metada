module mod_csp_thdyn_force
   use mod_csp_misc
   implicit none
   real(wp), save, public :: zrgb(4, 61)    ! Full tabulated attenuation coefficients for RGB absorption
   real(wp), save, public :: rkrgb(3, 61)   ! tabulated attenuation coefficients for RGB absorption

contains

!==============================================================================
   subroutine csp_thdyn_force_qsr_init
!==============================================================================
    !!----------------------------------------------------------------------
    !! ** Purpose :   Compute the temperature surface boundary condition trend of
    !!      flux through the interface, concentration/dilution effect
    !!      and add it to the general trend of temperature equations.
    !!----------------------------------------------------------------------
      implicit none

      !  Chlorophyll      !   Blue attenuation   !   Green attenuation  !  Red attenuation    !
      zrgb(1, 1) = 0.010; zrgb(2, 1) = 0.01618; zrgb(3, 1) = 0.07464; zrgb(4, 1) = 0.37807
      zrgb(1, 2) = 0.011; zrgb(2, 2) = 0.01654; zrgb(3, 2) = 0.07480; zrgb(4, 2) = 0.37823
      zrgb(1, 3) = 0.013; zrgb(2, 3) = 0.01693; zrgb(3, 3) = 0.07499; zrgb(4, 3) = 0.37840
      zrgb(1, 4) = 0.014; zrgb(2, 4) = 0.01736; zrgb(3, 4) = 0.07518; zrgb(4, 4) = 0.37859
      zrgb(1, 5) = 0.016; zrgb(2, 5) = 0.01782; zrgb(3, 5) = 0.07539; zrgb(4, 5) = 0.37879
      zrgb(1, 6) = 0.018; zrgb(2, 6) = 0.01831; zrgb(3, 6) = 0.07562; zrgb(4, 6) = 0.37900
      zrgb(1, 7) = 0.020; zrgb(2, 7) = 0.01885; zrgb(3, 7) = 0.07586; zrgb(4, 7) = 0.37923
      zrgb(1, 8) = 0.022; zrgb(2, 8) = 0.01943; zrgb(3, 8) = 0.07613; zrgb(4, 8) = 0.37948
      zrgb(1, 9) = 0.025; zrgb(2, 9) = 0.02005; zrgb(3, 9) = 0.07641; zrgb(4, 9) = 0.37976
      zrgb(1, 10) = 0.028; zrgb(2, 10) = 0.02073; zrgb(3, 10) = 0.07672; zrgb(4, 10) = 0.38005
      zrgb(1, 11) = 0.032; zrgb(2, 11) = 0.02146; zrgb(3, 11) = 0.07705; zrgb(4, 11) = 0.38036
      zrgb(1, 12) = 0.035; zrgb(2, 12) = 0.02224; zrgb(3, 12) = 0.07741; zrgb(4, 12) = 0.38070
      zrgb(1, 13) = 0.040; zrgb(2, 13) = 0.02310; zrgb(3, 13) = 0.07780; zrgb(4, 13) = 0.38107
      zrgb(1, 14) = 0.045; zrgb(2, 14) = 0.02402; zrgb(3, 14) = 0.07821; zrgb(4, 14) = 0.38146
      zrgb(1, 15) = 0.050; zrgb(2, 15) = 0.02501; zrgb(3, 15) = 0.07866; zrgb(4, 15) = 0.38189
      zrgb(1, 16) = 0.056; zrgb(2, 16) = 0.02608; zrgb(3, 16) = 0.07914; zrgb(4, 16) = 0.38235
      zrgb(1, 17) = 0.063; zrgb(2, 17) = 0.02724; zrgb(3, 17) = 0.07967; zrgb(4, 17) = 0.38285
      zrgb(1, 18) = 0.071; zrgb(2, 18) = 0.02849; zrgb(3, 18) = 0.08023; zrgb(4, 18) = 0.38338
      zrgb(1, 19) = 0.079; zrgb(2, 19) = 0.02984; zrgb(3, 19) = 0.08083; zrgb(4, 19) = 0.38396
      zrgb(1, 20) = 0.089; zrgb(2, 20) = 0.03131; zrgb(3, 20) = 0.08149; zrgb(4, 20) = 0.38458
      zrgb(1, 21) = 0.100; zrgb(2, 21) = 0.03288; zrgb(3, 21) = 0.08219; zrgb(4, 21) = 0.38526
      zrgb(1, 22) = 0.112; zrgb(2, 22) = 0.03459; zrgb(3, 22) = 0.08295; zrgb(4, 22) = 0.38598
      zrgb(1, 23) = 0.126; zrgb(2, 23) = 0.03643; zrgb(3, 23) = 0.08377; zrgb(4, 23) = 0.38676
      zrgb(1, 24) = 0.141; zrgb(2, 24) = 0.03842; zrgb(3, 24) = 0.08466; zrgb(4, 24) = 0.38761
      zrgb(1, 25) = 0.158; zrgb(2, 25) = 0.04057; zrgb(3, 25) = 0.08561; zrgb(4, 25) = 0.38852
      zrgb(1, 26) = 0.178; zrgb(2, 26) = 0.04289; zrgb(3, 26) = 0.08664; zrgb(4, 26) = 0.38950
      zrgb(1, 27) = 0.200; zrgb(2, 27) = 0.04540; zrgb(3, 27) = 0.08775; zrgb(4, 27) = 0.39056
      zrgb(1, 28) = 0.224; zrgb(2, 28) = 0.04811; zrgb(3, 28) = 0.08894; zrgb(4, 28) = 0.39171
      zrgb(1, 29) = 0.251; zrgb(2, 29) = 0.05103; zrgb(3, 29) = 0.09023; zrgb(4, 29) = 0.39294
      zrgb(1, 30) = 0.282; zrgb(2, 30) = 0.05420; zrgb(3, 30) = 0.09162; zrgb(4, 30) = 0.39428
      zrgb(1, 31) = 0.316; zrgb(2, 31) = 0.05761; zrgb(3, 31) = 0.09312; zrgb(4, 31) = 0.39572
      zrgb(1, 32) = 0.355; zrgb(2, 32) = 0.06130; zrgb(3, 32) = 0.09474; zrgb(4, 32) = 0.39727
      zrgb(1, 33) = 0.398; zrgb(2, 33) = 0.06529; zrgb(3, 33) = 0.09649; zrgb(4, 33) = 0.39894
      zrgb(1, 34) = 0.447; zrgb(2, 34) = 0.06959; zrgb(3, 34) = 0.09837; zrgb(4, 34) = 0.40075
      zrgb(1, 35) = 0.501; zrgb(2, 35) = 0.07424; zrgb(3, 35) = 0.10040; zrgb(4, 35) = 0.40270
      zrgb(1, 36) = 0.562; zrgb(2, 36) = 0.07927; zrgb(3, 36) = 0.10259; zrgb(4, 36) = 0.40480
      zrgb(1, 37) = 0.631; zrgb(2, 37) = 0.08470; zrgb(3, 37) = 0.10495; zrgb(4, 37) = 0.40707
      zrgb(1, 38) = 0.708; zrgb(2, 38) = 0.09056; zrgb(3, 38) = 0.10749; zrgb(4, 38) = 0.40952
      zrgb(1, 39) = 0.794; zrgb(2, 39) = 0.09690; zrgb(3, 39) = 0.11024; zrgb(4, 39) = 0.41216
      zrgb(1, 40) = 0.891; zrgb(2, 40) = 0.10374; zrgb(3, 40) = 0.11320; zrgb(4, 40) = 0.41502
      zrgb(1, 41) = 1.000; zrgb(2, 41) = 0.11114; zrgb(3, 41) = 0.11639; zrgb(4, 41) = 0.41809
      zrgb(1, 42) = 1.122; zrgb(2, 42) = 0.11912; zrgb(3, 42) = 0.11984; zrgb(4, 42) = 0.42142
      zrgb(1, 43) = 1.259; zrgb(2, 43) = 0.12775; zrgb(3, 43) = 0.12356; zrgb(4, 43) = 0.42500
      zrgb(1, 44) = 1.413; zrgb(2, 44) = 0.13707; zrgb(3, 44) = 0.12757; zrgb(4, 44) = 0.42887
      zrgb(1, 45) = 1.585; zrgb(2, 45) = 0.14715; zrgb(3, 45) = 0.13189; zrgb(4, 45) = 0.43304
      zrgb(1, 46) = 1.778; zrgb(2, 46) = 0.15803; zrgb(3, 46) = 0.13655; zrgb(4, 46) = 0.43754
      zrgb(1, 47) = 1.995; zrgb(2, 47) = 0.16978; zrgb(3, 47) = 0.14158; zrgb(4, 47) = 0.44240
      zrgb(1, 48) = 2.239; zrgb(2, 48) = 0.18248; zrgb(3, 48) = 0.14701; zrgb(4, 48) = 0.44765
      zrgb(1, 49) = 2.512; zrgb(2, 49) = 0.19620; zrgb(3, 49) = 0.15286; zrgb(4, 49) = 0.45331
      zrgb(1, 50) = 2.818; zrgb(2, 50) = 0.21102; zrgb(3, 50) = 0.15918; zrgb(4, 50) = 0.45942
      zrgb(1, 51) = 3.162; zrgb(2, 51) = 0.22703; zrgb(3, 51) = 0.16599; zrgb(4, 51) = 0.46601
      zrgb(1, 52) = 3.548; zrgb(2, 52) = 0.24433; zrgb(3, 52) = 0.17334; zrgb(4, 52) = 0.47313
      zrgb(1, 53) = 3.981; zrgb(2, 53) = 0.26301; zrgb(3, 53) = 0.18126; zrgb(4, 53) = 0.48080
      zrgb(1, 54) = 4.467; zrgb(2, 54) = 0.28320; zrgb(3, 54) = 0.18981; zrgb(4, 54) = 0.48909
      zrgb(1, 55) = 5.012; zrgb(2, 55) = 0.30502; zrgb(3, 55) = 0.19903; zrgb(4, 55) = 0.49803
      zrgb(1, 56) = 5.623; zrgb(2, 56) = 0.32858; zrgb(3, 56) = 0.20898; zrgb(4, 56) = 0.50768
      zrgb(1, 57) = 6.310; zrgb(2, 57) = 0.35404; zrgb(3, 57) = 0.21971; zrgb(4, 57) = 0.51810
      zrgb(1, 58) = 7.079; zrgb(2, 58) = 0.38154; zrgb(3, 58) = 0.23129; zrgb(4, 58) = 0.52934
      zrgb(1, 59) = 7.943; zrgb(2, 59) = 0.41125; zrgb(3, 59) = 0.24378; zrgb(4, 59) = 0.54147
      zrgb(1, 60) = 8.912; zrgb(2, 60) = 0.44336; zrgb(3, 60) = 0.25725; zrgb(4, 60) = 0.55457
      zrgb(1, 61) = 10.000; zrgb(2, 61) = 0.47804; zrgb(3, 61) = 0.27178; zrgb(4, 61) = 0.56870

      rkrgb(:, :) = zrgb(2:4, :)

   end subroutine csp_thdyn_force_qsr_init

!==============================================================================
   subroutine csp_thdyn_force_temp
!==============================================================================
    !!----------------------------------------------------------------------
    !! ** Purpose :   Compute the temperature surface boundary condition trend of
    !!      flux through the interface, concentration/dilution effect
    !!      and add it to the general trend of temperature equations.
    !!----------------------------------------------------------------------
      implicit none
      integer :: i, k
      real(wp) :: t_temp, r_temp1, r_temp2

      !$ACC kernels present(nlpb,nk,gTracer,qns,recip_drF,recip_hFacC,prec,snow,tracer,gravityDynForce)
      !$ACC loop
      do i = 1, nlpb
         gTracer(i, nk) = gTracer(i, nk) + qns(i)*gravityDynForce*recip_drF(nk)*recip_hFacC(i, nk)/rcp

         ! we assume rain have same temperature with sst, it is very rough, need improve in furture
         !gTracer(i, nk) = gTracer(i, nk) - prec(i, 4)*gravityDynForce*5.0_wp*recip_drF(nk)*recip_hFacC(i, nk)

         !gTracer(i, nk) = gTracer(i, nk) - snow(i, 4)*gravityDynForce*tracer(i, nk)*recip_drF(nk)*recip_hFacC(i, nk)
      end do
      !$ACC end kernels

      call csp_thdyn_force_qsr

      ! only for simulate ice effect
      if (gammat .gt. 0.0_wp) then
         r_temp1 = 1/(gammat*86400.0_wp)
         !$ACC kernels present(gammat,gTracer,sfrt,tFld,tracer),copyin(r_temp1)
         !$ACC loop
         do i = 1, nlpb
            if ((sfrt(i, 4) .lt. -1.0_wp) .or. (tFld(i, nk, 1) .lt. -2.0_wp)) then
               gTracer(i, nk) = gTracer(i, nk) + (sfrt(i, 4) - tracer(i, nk))*r_temp1
            end if
         end do
         !$ACC end kernels
      end if

      !r_temp2 = 1/dTtracer
      !!$ACC kernels present(dTtracer,gTracer,tracer),copyin(r_temp2)
      !!$ACC loop private(t_temp)
      !do i = 1, nlpb
      !   t_temp = tracer(i, nk) + gTracer(i, nk)*dTtracer
      !   if (t_temp .lt. -2.5_wp) then
      !      gTracer(i, nk) = (-2.5_wp - tracer(i, nk))*r_temp2
      !   end if
      !end do
      !!$ACC end kernels

# if defined (EXDIAG)
      !$ACC kernels present(dTtracer,gTracer,tFld_force)
      !$ACC loop collapse(2) independent
      do k = 1, nk
         do i = 1, loc_nlpb
            tFld_force(i, k) = tFld_force(i, k) + gTracer(i, k)*dTtracer
         end do
      end do
      !$ACC end kernels
# endif

   end subroutine csp_thdyn_force_temp

!==============================================================================
   subroutine csp_thdyn_force_salt
!==============================================================================
    !!----------------------------------------------------------------------
    !! ** Purpose :   Compute the tracer surface boundary condition trend of
    !!      flux through the interface, concentration/dilution effect
    !!      and add it to the general trend of tracer equations.
    !!----------------------------------------------------------------------
      implicit none
      integer :: i, k
      real(wp) :: t_temp, r_temp1, r_temp2

      ! note: EmPmR(:)*gravity*dTtracer is FreshWaterFlux in Pa unit,
      ! it better not larger than level thickness of top surface
      !$ACC kernels present(nlpb,nk,gTracer,EmPmR,tracer,recip_hFacC,drF,dTtracer,gravityDynForce)
      !$ACC loop
      do i = 1, nlpb
         gTracer(i, nk) = gTracer(i, nk) + EmPmR(i)*gravityDynForce*tracer(i, nk)*recip_hFacC(i, nk) &
                          /(drF(nk) - EmPmR(i)*gravityDynForce*dTtracer)
        !YY: SEAICE//
!  saltFlux :: Net upward salt flux in g/kg.kg/m^2/s = g/m^2/s
!  fromMIT     flux of Salt taken out of the ocean per time unit (second).
!              Note: only used when salty sea-ice forms or melts.
!              > 0 for decrease in SSS.
!YY: hsalt = saltflux*deltaTtherm
         gTracer(i, nk) = gTracer(i, nk) + saltflux(i)*gravityDynForce*recip_hFacC(i, nk) &
                          /(drF(nk) - EmPmR(i)*gravityDynForce*dTtracer)
      end do
      !$ACC end kernels

      ! only for simulate ice effect
      if (gammas .gt. 0.0_wp) then
         r_temp1 = 1/(gammas*86400.0_wp)
         !$ACC kernels present(gammas,sfrt,tFld,gTracer,sfrs,tracer),copyin(r_temp1)
         !$ACC loop
         do i = 1, nlpb
            if ((sfrt(i, 4) .lt. -1.0_wp) .or. (tFld(i, nk, 1) .lt. -2.0_wp)) then
               gTracer(i, nk) = gTracer(i, nk) + (sfrs(i, 4) - tracer(i, nk))*r_temp1
            end if
         end do
         !$ACC end kernels
      end if

      r_temp2 = 1/dTtracer
      !$ACC kernels present(dTtracer,gTracer,tracer),copyin(r_temp2)
      !$ACC loop private(t_temp)
      do i = 1, nlpb
         t_temp = tracer(i, nk) + gTracer(i, nk)*dTtracer
         if (t_temp .lt. 0.1_wp) then
            gTracer(i, nk) = (0.1_wp - tracer(i, nk))*r_temp2
         end if
      end do
      !$ACC end kernels
      
# if defined (EXDIAG)
      !$ACC kernels present(dTtracer,gTracer,sFld_force)
      !$ACC loop collapse(2) independent
      do k = 1, nk
         do i = 1, loc_nlpb
            sFld_force(i, k) = sFld_force(i, k) + gTracer(i, k)*dTtracer
         end do
      end do
      !$ACC end kernels
# endif

   end subroutine csp_thdyn_force_salt

!==============================================================================
   subroutine csp_thdyn_force_qsr
!==============================================================================
      !!----------------------------------------------------------------------
      !! ** Purpose :   Compute the temperature trend due to the solar radiation
      !!      penetration and add it to the general temperature trend.
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 3 wavebands (RGB) and a ratio rn_abs, the computation is only
      !!      done down to the level where I(k) < 1.e-15 W/m2. In addition,
      !!      the coefficients used for the computation are calculated one for once
      !!      as they depends on k only.
      !!
      !! Reference  : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!              Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!----------------------------------------------------------------------
      implicit none
      integer :: i, k, irgb
      real(wp) :: zchl, zcoef, gp1
      real(wp) :: zekb, zekg, zekr
      real(wp) :: ze0, ze1, ze2, ze3, zea1, zea2
      real(wp) :: zc0, zc1, zc2, zc3

      gp1 = 1.0_wp/gravityRau0

      ! constant chlorophyll
      ! in future we can try variable one use bio-model
      zchl = 0.05_wp
      irgb = nint(41.0_wp + 20.0_wp*log10(zchl) + 1.0e-15_wp)
      zekb = rkrgb(1, irgb)                       ! Separation in R-G-B depending of the chlorophyll
      zekg = rkrgb(2, irgb)
      zekr = rkrgb(3, irgb)

      zcoef = (1.-rn_abs)/3.0_wp                  ! equi-partition in R-G-B

      !$ACC kernels present(nlpb,nk,nksr,gTracer,qsr,rn_abs,rn_si0,drF,  &
      !$ACC         hFacC,maskC,recip_drF,recip_hFacC,gravityDynForce),  &
      !$ACC         copyin(zekb,zekg,zekr,zcoef,gp1)
      !$ACC loop private(ze0,ze1,ze2,ze3,zea1,zc0,zc1,zc2,zc3,zea2)
      do i = 1, nlpb
         ze0 = rn_abs*qsr(i)
         ze1 = zcoef*qsr(i)
         ze2 = zcoef*qsr(i)
         ze3 = zcoef*qsr(i)
         zea1 = qsr(i)
         !$ACC loop seq
         do k = nk - 1, nksr - 1, -1
            zc0 = ze0*exp(-drF(k + 1)*gp1*hFacC(i, k + 1)/rn_si0)
            zc1 = ze1*exp(-drF(k + 1)*gp1*hFacC(i, k + 1)*zekb)
            zc2 = ze2*exp(-drF(k + 1)*gp1*hFacC(i, k + 1)*zekg)
            zc3 = ze3*exp(-drF(k + 1)*gp1*hFacC(i, k + 1)*zekr)
            ze0 = zc0
            ze1 = zc1
            ze2 = zc2
            ze3 = zc3
            zea2 = (zc0 + zc1 + zc2 + zc3)*maskC(i, k)

            gTracer(i, k + 1) = gTracer(i, k + 1) + (zea1 - zea2)*gravityDynForce*recip_drF(k + 1)*recip_hFacC(i, k + 1)/rcp

            zea1 = zea2
         end do
      end do
      !$ACC end kernels

   end subroutine csp_thdyn_force_qsr

end module mod_csp_thdyn_force
