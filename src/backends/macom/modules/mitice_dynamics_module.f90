module mitice_dynamics
use mitice_itd
use mitice_parameters
use mitice_vars
use mod_csp_basic
use mod_csp_dyn_mom_vecinv
use mod_csp_dyn_vort
use mod_mpi_interfaces
#ifdef SeaiceDebug
use mod_mpi_variables
use mitice_debug
#endif

   implicit none


contains
   
   subroutine mitice_dynsolver

!  *===========================================================*
!  Ice dynamics using explicit EVP solver
!  *===========================================================*

!     === Local variables ===
!  TAUX   :: zonal      wind stress over seaice at U point
!  TAUY   :: meridional wind stress over seaice at V point
      INTEGER i
      real(wp) ::  phiSurf(nlpb)
      real(wp) :: TAUX(nlpb), TAUY(nlpb)
      real(wp) ::  mask_uice
      real(wp) :: heffTooHeavy, dzSurf  !surface layer thickness in meters
      integer :: w, s

#ifdef SeaiceDebug
      !real(wp) :: term13(nlpb), term14(nlpb), term15(nlpb)
      !real(wp) :: term16(nlpb), term17(nlpb)
      !real(wp) :: termAve
#endif

!YY: makesure phi0surf is right;
!$acc data present(tw,ts,seaiceMassC,seaiceMassU,seaiceMassV,AREA,HEFF,HSNOW,   &
!$acc              seaiceMaskU,seaiceMaskV,HEFFM,Bo_surf,etaH,phi0surf,sIceLoad,&
!$acc             FORCEX0,FORCEY0,recip_dxC,recip_dyC),                         &
!$acc      create(phiSurf,TAUX,TAUY)

!$acc kernels
!--   now set up mass per unit area and coriolis term
      IF ( SEAICEaddSnowMass ) THEN    !default is true
!$acc loop private(w,s)
        do i = 1,nlpb
          w = tw(i)
          s = ts(i)
          seaiceMassC(i)=SEAICE_rhoIce*HEFF(i) + SEAICE_rhoSnow*HSNOW(i)
          seaiceMassU(i)=SEAICE_rhoIce*0.5_wp*( HEFF(i) + HEFF(w) )  &
              + SEAICE_rhoSnow*0.5_wp*( HSNOW(i) + HSNOW(w) )
          seaiceMassV(i)=SEAICE_rhoIce*0.5_wp*( HEFF(i) + HEFF(s) )  &
              + SEAICE_rhoSnow*0.5_wp*( HSNOW(i) + HSNOW(s) )
        enddo
      ELSE
!$acc loop private(w,s)
        do i = 1,nlpb
          w = tw(i)
          s = ts(i)
          seaiceMassC(i)=SEAICE_rhoIce*HEFF(i)
          seaiceMassU(i)=SEAICE_rhoIce*0.5_wp*( HEFF(i) + HEFF(w) )
          seaiceMassV(i)=SEAICE_rhoIce*0.5_wp*( HEFF(i) + HEFF(s) )
        enddo
      ENDIF

!     dynamic masking of areas with no ice
!YY: for ice-water interface grid, uice can be blown up if computing denomU with 
!YY: extremely small values of areaW and areaS
!YY: so here for these points, put seaiceMaskU/V to zero
!$acc loop private(w,s,mask_uice)
      DO i= 1, nlpb
         w = tw(i)
         s = ts(i)
         seaiceMaskU(i)=0.5_wp*(AREA(i)+AREA(w))
         mask_uice=HEFFM(i)+HEFFM(w)
         IF ( (seaiceMaskU(i) .GT. SEAICE_area_floor) .AND. (mask_uice .GT. 1.5_wp) ) THEN
           seaiceMaskU(i) = 1.0_wp
         ELSE
           seaiceMaskU(i) = 0.0_wp
         ENDIF
         seaiceMaskV(i)=0.5_wp*(AREA(i)+AREA(s))
         mask_uice=HEFFM(i)+HEFFM(s)
         IF ( (seaiceMaskV(i) .GT. SEAICE_area_floor) .AND. (mask_uice .GT. 1.5_wp) ) THEN
          seaiceMaskV(i) = 1.0_wp
         ELSE
          seaiceMaskV(i) = 0.0_wp
         ENDIF
      ENDDO

!--   NOW SET UP FORCING FIELDS

! Initialise local variables
!$acc loop 
      do i=1,nlpb      
        TAUX(i) = 0.0_wp
        TAUY(i) = 0.0_wp
      enddo
!$acc end kernels

!--   atmopheric forcing fields (wind stress over ice: taux,tauy)
      call seaice_get_dynforcing(TAUX,TAUY)

!--   Compute surface pressure at z==0:
!-    use actual sea surface height for tilt computations
!$acc kernels loop
      do i = 1,nlpb
        phiSurf(i) = Bo_surf(i)*etaH(i)   !here Bo_surf is gravity for p-coord
      enddo
!$acc end kernels

!--   add atmospheric loading and Sea-Ice loading

!YY,ZY: in mod_csp_dyn_phi, phi0surf is removed from ssh. here phi0surf is added here again
!YY,ZY: phi0surf is not initiated or zero-ed. here I add phi0surf=0.0 at mod_csp_init
!YY,ZY: r1_rau0 is divided(convert to water-equivalent height) 
      IF ( useRealFreshWaterFlux ) THEN
!     Cutoff for iceload
      !dzSurf   = drF(kSurface)*hFacC(i,kSurface)/gravityRau0  !YY: donot consider hFacC now
      !dzSurf   = drF(kSurface)/gravityRau0
      !heffTooHeavy = dzSurf * 0.2_wp
!$acc kernels loop
        do i = 1,nlpb
          ! Sea Ice Load on the sea surface.
          !YY: here sIceLoad is recomputed, as hot start require initial value of sIceload
          !YY: in thermodyn module, sIceLoad is computed first.
          sIceLoad(i) = HEFF(I)*SEAICE_rhoIce + HSNOW(I)*SEAICE_rhoSnow
          !if ice is too heavy compared with the surface layer of water
          !sIceLoad(i)= MIN(HEFF(I)*SEAICE_rhoIce + HSNOW(I)*SEAICE_rhoSnow , heffTooHeavy*rau0)
          phiSurf(i) = phiSurf(i)                &
                      + (phi0surf(i)             & 
                      + sIceLoad(i)*gravity*sIceLoadFac) * r1_rau0   
#ifdef SeaiceDebug
          !term13(i) = phi0surf(i)*r1_rau0
          !term14(i) = sIceLoad(i)*gravity*sIceLoadFac * r1_rau0
          !term15(i) = Bo_surf(i)*etaH(i)
#endif
        enddo
!$acc end kernels
      ELSE
!$acc kernels loop
        do i = 1,nlpb
          phiSurf(i) = phiSurf(i) + phi0surf(i)*r1_rau0
        enddo
!$acc end kernels
      ENDIF

#ifdef SeaiceDebug
     ! termAve = sum(abs(term13))/nlpb
     ! if (mpi_rank==21) write(mitice_runmsg_unit,*) 'term13 ave: ', myIter, termAve
     ! termAve = sum(abs(term14))/nlpb
     ! if (mpi_rank==21) write(mitice_runmsg_unit,*) 'term14 ave: ', myIter, termAve
     ! termAve = sum(abs(term15))/nlpb
     ! if (mpi_rank==21) write(mitice_runmsg_unit,*) 'term15 ave: ', myIter, termAve
#endif

!--   basic forcing by wind stress
!$acc kernels
      IF ( SEAICEscaleSurfStress ) THEN
!$acc loop private(w,s)
        do i = 1,nlpb
          w = tw(i)
          s = ts(i)
          FORCEX0(i)=TAUX(i) * 0.5_wp*(AREA(i)+AREA(w))
          FORCEY0(i)=TAUY(i) * 0.5_wp*(AREA(i)+AREA(s))
        enddo
      ELSE
!$acc loop private(w,s)
        do i = 1,nlpb
          w = tw(i)
          s = ts(i)
          FORCEX0(i)=TAUX(i)
          FORCEY0(i)=TAUY(i)
        enddo
      ENDIF

      IF ( SEAICEuseTILT ) then
!--   add in tilt
!$acc loop private(w,s)
        do i = 1,nlpb
          w = tw(i)
          s = ts(i)
          FORCEX0(i)=FORCEX0(i)-seaiceMassU(i)*recip_dxC(i)*( phiSurf(i)-phiSurf(w) )
          FORCEY0(i)=FORCEY0(i)-seaiceMassV(i)*recip_dyC(i)*( phiSurf(i)-phiSurf(s) )
        enddo
      ENDIF
!$acc end kernels

      call seaice_calc_ice_strength

      IF ( SEAICEuseDYNAMICS ) THEN
#ifdef SEAICE_ALLOW_FREEDRIFT
        !IF ( SEAICEuseFREEDRIFT .OR. SEAICEuseEVP ) THEN   ! YY: use EVP for what?
        IF ( SEAICEuseFREEDRIFT) THEN
         call seaice_freedrift
        ENDIF
        IF ( SEAICEuseFREEDRIFT ) THEN
          do i = 1,nlpb
            uIce(i) = uIce_fd(i)
            vIce(i) = vIce_fd(i)
            stressDivergenceX(i) = 0.0_wp
            stressDivergenceY(i) = 0.0_wp
          enddo
        ENDIF
#endif /* SEAICE_ALLOW_FREEDRIFT */

        IF ( SEAICEuseEVP ) THEN
!       Elastic-Viscous-Plastic solver, following Hunke (2001)
          call seaice_evp
        ENDIF
      ENDIF    !ENDIF FOR SEAICEuseDYNAMICS

! Update ocean surface stress
!YY: taux/y is wind stress over ice, while utau is a MaCOM var defining windstress over open water
!YY: u/vtau are updated considering ice cover. 
      call seaice_ocean_stress (taux, tauy)
      
#ifdef SEAICE_ALLOW_DYNAMICS
      IF ( SEAICEuseDYNAMICS .AND. SEAICE_clipVelocities) THEN
! Put a cap on ice velocity
! limit velocity to 0.40 m s-1 to avoid potential CFL violations
! in open water areas (drift of zero thickness ice)
!$acc kernels loop
        do i = 1,nlpb
          uIce(i)=MAX(MIN(uIce(i),0.4_wp),-0.4_wp)
          vIce(i)=MAX(MIN(vIce(i),0.4_wp),-0.4_wp)
        enddo
!$acc end kernels
      ENDIF
#endif /* SEAICE_ALLOW_DYNAMICS */

!$acc end data

   end subroutine mitice_dynsolver



   subroutine seaice_get_dynforcing(TAUX,TAUY)

!     *===========================================================*
!     |   compute surface stress from atmopheric forcing fields
!     *===========================================================*
!YY: utau/vtau in MaCOM is wind stress over grid; taux/y here is wind stress over ice

      implicit none
!LOCAL VARS
      INTEGER i
      real(wp)  COSWIN, SINWIN
      real(wp) CDAIR  !local wind stress coefficient
      real(wp) AAA
      real(wp) uTmp, vTmp 
      real(wp) :: TAUX(nlpb), TAUY(nlpb)

!$acc data present(TAUX,TAUY,u10,v10,uIce,vIce,latC,SIMaskU,SIMaskV,    &
!$acc              SEAICE_EPS_SQ,SEAICE_EPS,SEAICE_rhoAir,         &
!$acc              SEAICE_drag,SEAICE_drag_south)

!--   introduce turning angle (default is zero)
      SINWIN=SIN(SEAICE_airTurnAngle*radPi)
      COSWIN=COS(SEAICE_airTurnAngle*radPi)

!--   NOW SET UP FORCING FIELDS

!!YY: utau/vtau in MaCOM is computed by relative wind (to ocean current uFld)
!!YY: taux/y's computation use relative wind (to sea ice velocity)

!$acc kernels copyin(SINWIN,COSWIN)
!$acc loop private(uTmp,vTmp,AAA,CDAIR)
      do i = 1,nlpb
        uTmp = u10(i,4) - uIce(i)    ! YY: suppose u10 at U/V gridpoint
        vTmp = v10(i,4) - vIce(i)
        
!--   Now compute ice surface stress
        AAA=uTmp**2+vTmp**2
        IF ( AAA .LE. SEAICE_EPS_SQ ) THEN
          AAA=SEAICE_EPS
        ELSE
          AAA=SQRT(AAA)
        ENDIF
        IF ( latC(i) .LT. ZERO ) THEN
         CDAIR = SEAICE_rhoAir*SEAICE_drag_south*AAA
        ELSE
         CDAIR = SEAICE_rhoAir*SEAICE_drag*AAA
        ENDIF
        taux(i) = ( CDAIR*(COSWIN*uTmp-SINWIN*vTmp) )*SIMaskU(i)
        tauy(i) = ( CDAIR*(SINWIN*uTmp+COSWIN*vTmp) )*SIMaskV(i)
      enddo
!$acc end kernels
!$acc end data

   end subroutine seaice_get_dynforcing


#ifdef SEAICE_ALLOW_FREEDRIFT
   subroutine seaice_freedrift
!     *==========================================================*
!     | o Solve ice approximate momentum equation analytically
!     *==========================================================*
      implicit none

!     === Local variables ===
      INTEGER i, kSrf
      integer w, s, uel, uei, uep, vnl, vni, vnp

      real(wp) tmpscal1,tmpscal2,tmpscal3,tmpscal4

      real(wp) taux_onIce_cntr, tauy_onIce_cntr, uvel_cntr, vvel_cntr
      real(wp) mIceCor, rhs_x, rhs_y, rhs_n, rhs_a, sol_n, sol_a

      real(wp)  :: uice_cntr(nlpb)
      real(wp)  :: vice_cntr(nlpb)

      real(wp) :: vecFldLoc(nlpb, 2)
      real(wp) :: ForceLoc (nlpb, 2)

      !YY: only for MaCOM
       kSrf = nk

! initialize fields:
      do i = 1,nlpb
        uice_fd(i) = 0.0_wp
        vice_fd(i) = 0.0_wp
        uice_cntr(i)=0.0_wp
        vice_cntr(i)=0.0_wp
      enddo

      do i = 1,nlpb
        ForceLoc(i,iu) = FORCEX0(i)
        ForceLoc(i,iv) = FORCEY0(i)
        vecFldLoc(i, iu) = uFld(i, kSrf)
        vecFldLoc(i, iv) = vFld(i, kSrf)
      enddo
      do i = 1,nlpb   !not include halo point
        uel = ue(i, 1)
        uei = ue(i, 2)
        uep = ue(i, 3)
        vnl = vn(i, 1)
        vni = vn(i, 2)
        vnp = vn(i, 3)

! preliminary computations:
! =========================
! air-ice stress at cell center
! FORCEX0 is defined at U-grid.
!YY,ZY: check uep necessary?
        taux_onIce_cntr=HALF * (ForceLoc(i,iu)+ForceLoc(uel,uei)*uep)
        tauy_onIce_cntr=HALF * (ForceLoc(i,iv)+ForceLoc(vnl,vni)*vnp)
! mass of ice per unit area (kg/m2) times coriolis f
        mIceCor=SEAICE_rhoIce*HEFF(i)*fCori(i)
! ocean velocity at cell center
        uvel_cntr=HALF*(vecFldLoc(i,iu)+vecFldLoc(uel,uei)*uep)
        vvel_cntr=HALF*(vecFldLoc(i,iv)+vecFldLoc(vnl,vni)*vnp)
! right hand side of free drift equation:
        rhs_x= -taux_onIce_cntr -mIceCor*vvel_cntr
        rhs_y= -tauy_onIce_cntr +mIceCor*uvel_cntr
! norm of angle of rhs
        tmpscal1=rhs_x*rhs_x + rhs_y*rhs_y
        IF ( tmpscal1.GT.ZERO ) THEN
         rhs_n=SQRT( rhs_x*rhs_x + rhs_y*rhs_y )
         rhs_a=ATAN2(rhs_y,rhs_x)
        ELSE
         rhs_n=0.0_wp
         rhs_a=0.0_wp
        ENDIF

! solve for norm:
        IF ( latC(i) .LT. ZERO ) THEN
         tmpscal1 = 1.0_wp /rau0/SEAICE_waterDrag_south
        ELSE
         tmpscal1 = 1.0_wp /rau0/SEAICE_waterDrag
        ENDIF
! polynomial coefficients
        tmpscal2= tmpscal1*tmpscal1*mIceCor*mIceCor
        tmpscal3= tmpscal1*tmpscal1*rhs_n*rhs_n
! discriminant
        tmpscal4=tmpscal2*tmpscal2+4.0_wp*tmpscal3
        IF ( tmpscal3.GT.ZERO ) THEN
         sol_n=SQRT(HALF*(SQRT(tmpscal4)-tmpscal2))
        ELSE
         sol_n=0.0_wp
        ENDIF

! solve for angle:
        IF ( latC(i) .LT. ZERO ) THEN
         tmpscal1 = SEAICE_waterDrag_south*rau0
        ELSE
         tmpscal1 = SEAICE_waterDrag*rau0
        ENDIF

        tmpscal2= tmpscal1*sol_n*sol_n
        tmpscal3= mIceCor*sol_n

        tmpscal4=tmpscal2*tmpscal2 + tmpscal3*tmpscal3
        IF ( tmpscal4.GT.ZERO ) THEN
         sol_a=rhs_a-ATAN2(tmpscal3,tmpscal2)
        ELSE
         sol_a=0.0_wp
        ENDIF

! compute uice, vice at cell center:
        uice_cntr(i)=uvel_cntr-sol_n*COS(sol_a)
        vice_cntr(i)=vvel_cntr-sol_n*SIN(sol_a)

      enddo

! interpolated to velocity points:

      call mpi_data_exchange(mpi_sendbuf_1d, uice_cntr)
      call mpi_data_exchange(mpi_sendbuf_1d, vice_cntr)

      do i = 1,nlpb    !not use halo point
        w = tw(i)
        s = ts(i)
        uice_fd(i)=HALF*(uice_cntr(w)+uice_cntr(i))
        vice_fd(i)=HALF*(vice_cntr(s)+vice_cntr(i))
      enddo

      call mpi_data_exchange(mpi_sendbuf_1d, uice_fd)
      call mpi_data_exchange(mpi_sendbuf_1d, vice_fd)

!     Apply masks (same/similar to seaice_evp.F)
      do i = 1,nlpb
        uIce_fd(i)=uIce_fd(i)*SIMaskU(i)   !SIMaskU depends only geometry
        vIce_fd(i)=vIce_fd(i)*SIMaskV(i)
      enddo

   end subroutine seaice_freedrift
#endif /* SEAICE_ALLOW_FREEDRIFT */


#ifdef SEAICE_ALLOW_DYNAMICS
   subroutine seaice_evp
!  *==========================================================*
!  | o Ice dynamics using an EVP solver following
!  |   E. C. Hunke and J. K. Dukowicz. An
!  |   Elastic-Viscous-Plastic Model for Sea Ice Dynamics,
!  |   J. Phys. Oceanogr., 27, 1849-1867 (1997).
!  *==========================================================*
!  | written by Martin Losch, March 2006
!  *==========================================================*

      implicit none

!     kSrf           :: vertical index of surface layer
!     nEVPstep       :: number of timesteps within the EVP solver
!     SIN/COSWAT     :: sine/cosine of turning angle
!     (recip_)ecc2   :: (one over) eccentricity squared
!     recip_evpAlpha :: 1/SEAICE_evpAlpha
!     recip_deltaT   :: 1/SEAICE_deltaTdyn
!     evpStarFac     :: 1 if SEAICEuseEVPstar = .true., 0 otherwise
!     betaFac        :: SEAICE_evpBeta/SEAICE_deltaTdyn=1/SEAICE_deltaTevp
!     betaFacP1      :: betaFac + evpStarFac/SEAICE_deltaTdyn
!     e11,e12,e22    :: components of strain rate tensor
!     seaice_div     :: divergence strain rates at C-points times P
!                       /divided by Delta minus 1
!     seaice_tension :: tension    strain rates at C-points times P
!                       /divided by Delta
!     seaice_shear   :: shear      strain rates, defined at Z-points times P
!                       /divided by Delta
!     sig11, sig22   :: sum and difference of diagonal terms of stress tensor
!     modification for adaptive alpha and beta
!               (see Kimmritz, Danilov, Losch 2015 for gamma << alpha beta)
!     EVPcFac        :: SEAICE_deltaTdyn*SEAICEaEVPcStar*(SEAICEaEVPcoeff*PI)**2
!                        with
!     SEAICEaEVPcStar:: multiple of stabilty factor: alpha*beta = cstar*gamma
!     SEAICEaEVPcoeff:: largest stabilized frequency according to
!                        gamma = zeta * (cfac/cellarea)*deltaT/m
!                                with   (cfac/cellarea) <= pi**2/cellarea
!     evpAlphaC/Z    :: alpha field on C points and on Z points
!                        := sqrt(cstar gamma)
!     evpBetaU/V     :: beta field on u and on v points
!                        := sqrt(cstar gamma)
!     evpAlphaMin    :: lower limit of alpha and beta, regularisation
!                     to prevent singularities of system matrix,
!                     e.g. when ice concentration is too low.
!     betaFacP1U/V   :: = betaFacP1 in standard case,
!                          with varying beta in the adaptive case
!                           on u and on v point
!     betaFacU/V     :: analog betaFacP1U/V

      integer i
      integer kSrf
      integer nEVPstep, iEVPstep
      integer nEVPstepMax

      real(wp) COSWAT, SINWAT
      real(wp) ecc2, recip_ecc2, recip_deltaT
      !real(wp) recip_evpAlpha
      real(wp) betaFacP1, betaFac, evpStarFac, evpRevFac, recip_evpRevFac

      real(wp)  ::  seaice_div(nlpb), seaice_tension(nlpb)    !at C grid
      !YUAN, ZY: seaice_shear should have dimension of nlpbz
      real(wp)  ::  seaice_shear(nlpbz)    !at Z grid
      real(wp)  ::  sig11(nlpb), sig22(nlpb)
!     fractional area at velocity points
      real(wp)  ::  areaW(nlpb), areaS(nlpb)
!     auxilliary variables
      real(wp)  ::  ep(nlpb), em(nlpb), e12Csq(nlpb)    !e12C at C-point
      real(wp)  ::  pressC(nlpb), zetaC(nlpb), deltaZ(nlpb)
#ifdef SEAICE_ALLOW_MOM_ADVECTION
!     tendency due to advection of momentum
      real(wp)  ::  gUmom(nlpb), gVmom(nlpb)
#endif /*  SEAICE_ALLOW_MOM_ADVECTION */
      real(wp) deltaCreg, deltaSq, deltaMinSq, tmp
      real(wp) etaDenC, zetaMaxC, etaDenZ, zetaMaxZ
      real(wp) zMaxZ, zMinZ, fac
      real(wp)  ::  denom1, denom2
      real(wp) sumNorm, denomU, denomV
      real(wp) locMaskU, locMaskV
      real(wp) EVPcFac
      real(wp)  ::  evpAlphaC(nlpb)
      real(wp)  ::  evpAlphaZ(nlpbz)
      real(wp)  ::  evpBetaU (nlpb)
      real(wp)  ::  evpBetaV (nlpb)
      real(wp)  betaFacP1U, betaFacP1V, betaFacU, betaFacV
      LOGICAL useAdaptiveEVP

      real(wp)  resTile
      real(wp)  resLoc
      !YY:note that acc declare clause may not function normally in a subroutine.
      real(wp), allocatable  ::  uIcePm1 (:), vIcePm1 (:)
      !$acc declare create(uIcePm1,vIcePm1)
      real(wp), allocatable  ::  sig11pm1(:), sig22pm1(:), sig12pm1(:)
      !$acc declare create(sig11pm1,sig22pm1,sig12pm1)
      LOGICAL printResidual

      integer :: w,s,e,n, x1,x2
      integer :: i1,i2,i3,i4,i5
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4
      integer :: p1, p2, p3, p4
      real(wp) :: dyS_dxW(nlpb,2),dxS_dyW(nlpb,2),dyW_dxS(nlpb,2),dxW_dyS(nlpb,2)
      real(wp) :: vecIceLoc(nlpb,2), vecFldLoc(nlpb,2)

#ifdef SeaiceDebug
     real(wp) :: velIceMax, velMax, mpi_tmp
     real(wp) :: velIceMag(nlpb), velMag(nlpb)
!     real(wp) :: term1(nlpb), term2(nlpb),term3(nlpb),term4(nlpb),term5(nlpb),term6(nlpb),term7(nlpb)
!     real(wp) :: term8(nlpb), term9(nlpb),term10(nlpb),term11(nlpb),term12(nlpb)
#endif

!$acc data present(AREA,rAc,tw,ts,uIce,vIce,seaice_sigma1,seaice_sigma2,seaice_sigma12, &
!$acc              zetaZ,deltaC,avz1,avz2,avz4,avz5,e11,e22,e12,recip_rAc,recip_rAz,rAz, &
!$acc              press0,tensileStrFac,seaiceMassC,z1,z2,z3,z4,dyS,dxW,dxS,dyW,    &
!$acc              uw,un,ve,vs,zn,ze, stressDivergenceX,stressDivergenceY,recip_rAw,&
!$acc              uFld,vFld,recip_rAs,au1,au2,au3,au4,seaiceMassU,seaiceMassV,     &
!$acc              FORCEX,FORCEX0,FORCEY,FORCEY0,DWATN,fCori,av1,av2,av3,av4,       &
!$acc              seaiceMaskU,seaiceMaskV,uIceNm1,vIceNm1),     &
!$acc      create(seaice_div,seaice_tension,seaice_shear,sig11,sig22,areaW,areaS,   &
!$acc             ep,em,e12Csq,pressC,zetaC,deltaZ,evpAlphaC,evpAlphaZ,evpBetaU,    &
!$acc             evpBetaV,dyS_dxW,dxS_dyW,dyW_dxS,dxW_dyS,vecIceLoc,vecFldLoc,     &
!$acc             evpRevFac,recip_evpRevFac)


!     set tuning parameters for adaptive EVP
      useAdaptiveEVP = .FALSE.
      IF ( SEAICEaEvpCoeff .NE. UNSET_RL ) useAdaptiveEVP = .TRUE.
      EVPcFac = 0.0_wp
      IF (useAdaptiveEVP) EVPcFac = SEAICE_deltaTdyn*SEAICEaEVPcStar    &
                     * (SEAICEaEvpCoeff * rpi)**2

!#ifdef ALLOW_SEAICE_EVP_RESIDUAL
!YY: if drop ifdef, then above arrays should be declared as allocatable, and allocate them here
      printResidual = debuglevel.GE.1
      !printResidual = .TRUE.
      if (printResidual) then
        allocate(uIcePm1 (nlpb), vIcePm1 (nlpb) )
        allocate( sig11pm1(nlpb), sig22pm1(nlpb), sig12pm1(nlpbz) )
      endif
!#endif /* ALLOW_SEAICE_EVP_RESIDUAL */

!     surface level for pressure coord
      kSrf = nk

!--   introduce turning angles
      SINWAT=SIN(SEAICE_waterTurnAngle*radPi)
      COSWAT=COS(SEAICE_waterTurnAngle*radPi)

!     abbreviation eccentricity squared
      ecc2=SEAICE_eccen**2
      recip_ecc2 = 0.0_wp
      IF ( ecc2 .NE. 0.0_wp ) recip_ecc2=ONE/ecc2
      deltaMinSq = SEAICE_deltaMin**2
!     copy number of internal time steps from previously defined parameter
      nEVPstep = SEAICEnEVPstarSteps
!     SEAICE_evpAlpha = 2. * SEAICE_evpTauRelax/SEAICE_deltaTevp
      denom1 = 1.0_wp / ( SEAICE_evpAlpha + 1.0_wp )
      denom2 = 1.0_wp / ( SEAICE_evpAlpha + ecc2 )
      recip_deltaT = 1.0_wp / SEAICE_deltaTdyn
      !recip_evpAlpha = 0.0_wp
      !IF ( SEAICE_evpAlpha .GT. 0.0_wp ) recip_evpAlpha = 1.0_wp / SEAICE_evpAlpha
      evpStarFac = 0.0_wp
      evpRevFac  = 0.0_wp
      recip_evpRevFac = 1.0_wp
      IF ( SEAICEuseEVPstar ) evpStarFac = 1.0_wp
      IF ( SEAICEuseEVPrev  ) THEN
!     the Bouillon et al. (2013) discretization in time has more explicit terms
        evpRevFac       = 1.0_wp
        evpStarFac      = 1.0_wp
        recip_evpRevFac = recip_ecc2
        denom1     = 1.0_wp / SEAICE_evpAlpha
        denom2     = denom1
      ENDIF
! local device vars updated
!$acc update device(evpRevFac,recip_evpRevFac)

!YY: the following betafac parameters commented. need test
      !betaFac    = SEAICE_evpBeta*recip_deltaT
      !betaFacU   = betaFac
      !betaFacV   = betaFac

      !betaFacP1  = betaFac + evpStarFac*recip_deltaT
      !betaFacP1U = betaFacP1
      !betaFacP1V = betaFacP1

      nEVPstepMax = nEVPstep

!$acc kernels
!$acc loop
      do i = 1,nlpb
!     use u/vIce as work arrays: copy previous time step to u/vIceNm1
        uIceNm1(i) = uIce(i)
        vIceNm1(i) = vIce(i)
!     initialise strain rates
        e11(i)   = 0.0_wp
        e22(i)   = 0.0_wp
!     initialise adative-EVP-specific fields
        evpAlphaC(i) = SEAICE_evpAlpha
        evpBetaU (i) = SEAICE_evpBeta
        evpBetaV (i) = SEAICE_evpBeta
      enddo
!     YUAN, ZY:for those defined at Z-point, dimension is nlpbz
!$acc loop
      do i = 1,nlpbz
        e12(i) = 0.0_wp
        evpAlphaZ(i) = SEAICE_evpAlpha
      enddo
#ifdef SeaiceDebug
     ! velIceMag = sqrt(uIce*uIce+vIce*vIce)
     ! velMag = sqrt(uFld(:,nk)*uFld(:,nk)+vFld(:,nk)*vFld(:,nk))
     ! velIceMax = maxval(velIceMag(1:loc_nlpb))
     ! velMax = maxval(velMag(1:loc_nlpb))
     ! CALL MPI_ALLREDUCE(velIceMax, mpi_tmp, 1, mpi_real_wp, MPI_MAX, mpi_comp_comm, mpi_err)
     ! velIceMax = mpi_tmp
     ! CALL MPI_ALLREDUCE(velMax, mpi_tmp, 1, mpi_real_wp, MPI_MAX, mpi_comp_comm, mpi_err)
     ! velMax = mpi_tmp
     ! !if(mpi_rank==21) then
     ! if(mpi_rank==0) then
     !    print *, 'myIter= ',myIter,mpi_rank, ' uIceMax: ', velIceMax, 'uMax', velMax
     !    write(mitice_runmsg_unit,*)'myIter= ',myIter,mpi_rank, ' uIceMax: ', velIceMax, 'uMax', velMax
     ! endif
#endif
!     initialise fractional areas at velocity points
      IF ( SEAICEscaleSurfStress ) THEN
!$acc loop private(w,s)
        do i = 1,nlpb
          w = tw(i)
          s = ts(i)
          areaW(i) = 0.5_wp*(AREA(i)+AREA(w))
          areaS(i) = 0.5_wp*(AREA(i)+AREA(s))
        enddo
      ELSE
!$acc loop
        do i = 1,nlpb
          areaW(i) = 1.0_wp
          areaS(i) = 1.0_wp
        enddo
      ENDIF
!$acc end kernels


!#ifdef SEAICE_ALLOW_CLIPZETA
!     damping constraint (Hunke, J.Comp.Phys.,2001)
      IF ( SEAICE_evpDampC .GT. 0.0_wp ) THEN
        fac = 0.25_wp * SEAICE_evpDampC * SEAICE_evpAlpha /SEAICE_deltaTevp
!$acc kernels loop copyin(fac)
        do i = 1,nlpb
          zMax (i)   = rAc(i) * fac
        enddo
!$acc end kernels
      ENDIF
!#endif  

!     start of the main time loop
      DO iEVPstep = 1, nEVPstepMax
       IF (iEVPstep.LE.nEVPstep) THEN

!     first calculate strain rates and bulk moduli/viscosities
        call seaice_calc_strainrates(uice, vice, e11, e22, e12, ievpstep)
!
!#ifdef ALLOW_SEAICE_EVP_RESIDUAL
!     save previous (p-1) iteration
        IF ( printResidual ) THEN
!$acc kernels present(sig11Pm1,sig22Pm1,sig12Pm1,uIcePm1,vIcePm1)
!$acc loop
          do i = 1,nlpb
            sig11Pm1(i) = seaice_sigma1(i)
            sig22Pm1(i) = seaice_sigma2(i)
            uIcePm1 (i) = uIce(i)
            vIcePm1 (i) = vIce(i)
          enddo
!$acc loop
          do i = 1,nlpbz
            sig12Pm1(i) = seaice_sigma12(i)
          enddo
!$acc end kernels
        ENDIF
!#endif /* ALLOW_SEAICE_EVP_RESIDUAL */

!$acc kernels
!$acc loop
        do i = 1,nlpb
          seaice_div    (i) = 0.0_wp
          seaice_tension(i) = 0.0_wp
          pressC        (i) = 0.0_wp
          e12Csq        (i) = 0.0_wp
          zetaC         (i) = 0.0_wp
          deltaC        (i) = 0.0_wp
          ep(i) = e11(i) + e22(i)    !principal strain
          em(i) = e11(i) - e22(i)
        enddo

        !YY, ZY: nlpbz
!$acc loop
        do i = 1,nlpbz
          seaice_shear(i) = 0.0_wp
          deltaZ      (i) = 0.0_wp
          zetaZ       (i) = 0.0_wp
        enddo
        
!     average strain rates to C points
        IF ( SEAICEetaZmethod .EQ. 0 ) THEN
!YY, ZY: Z-point averaged to C-point
!$acc loop private(i1,i2,i4,i5,tmp)
          do i = 1, nlpb
            i1 = avz1(i)
            i2 = avz2(i)
            i4 = avz4(i)
            i5 = avz5(i)
            tmp = 0.25_wp*( e12(i1) + e12(i2) + e12(i4) + e12(i5) )
            e12Csq(i) = tmp*tmp
          enddo

        ELSEIF ( SEAICEetaZmethod .EQ. 3 ) THEN
!$acc loop private(i1,i2,i4,i5)
          do i = 1, nlpb
            i1 = avz1(i)
            i2 = avz2(i)
            i4 = avz4(i)
            i5 = avz5(i)
!     area weighted average of the squares of e12 is more accurate
!     (and energy conserving) according to Bouillon et al. 2013, eq (11)
            e12Csq(i) = 0.25_wp * recip_rAc(i) *             &
                ( rAz(i1)*e12(i1)**2 + rAz(i2)*e12(i2)**2    &
                + rAz(i4)*e12(i4)**2 + rAz(i5)*e12(i5)**2 )
          enddo
        ENDIF
!$acc end kernels

!$acc kernels copyin(recip_ecc2,EVPcFac,deltaMinSq,ecc2)
!$acc loop private(deltaSq,deltaCreg)
        do i = 1,nlpb
          deltaSq = ep(i)**2 + recip_ecc2*em(i)**2 + recip_ecc2*4.0_wp*e12Csq(i)
          deltaC(i) = SQRT(deltaSq)
!     smooth regularization (without max-function) of delta for
!     better differentiability
          deltaCreg  = deltaC(i) + SEAICE_deltaMin
          zetaC(i) = HALF*( press0(i) * ( 1.0_wp + tensileStrFac(i) ) )/deltaCreg
        enddo

        IF ( useAdaptiveEVP ) THEN
!$acc loop  
          do i = 1,nlpb
!ML   I do not like these hidden regularisations, why do we need to
!ML   divide by mass?
            evpAlphaC(i) = SQRT(zetaC(i)                     &
                * EVPcFac / MAX(seaiceMassC(i), 1.0D-04)     &
                * recip_rAc(i) ) * HEFFM(i)
            evpAlphaC(i) = MAX(evpAlphaC(i),SEAICEaEVPalphaMin)
          enddo
        ENDIF
!     compute zetaZ and deltaZ by simple averaging (following
!     Bouillon et al., 2013)

        !averaging from C-grid to Z-grid
        !YY,ZY: at special vorticity point with only 3 neighboring grids, the
        !averaging below is not accurate.
        !Similar situation also exists below
!$acc loop private(i1,i2,i3,i4,sumNorm)
        do i = 1,nlpbz
          i1 = z1(i, 1)
          i2 = z2(i, 1)
          i3 = z3(i, 1)
          i4 = z4(i, 1)
          sumNorm = HEFFM(i1)+HEFFM(i3)+HEFFM(i4)+HEFFM( tw(i3) )
          IF ( sumNorm.GT.0.0_wp ) sumNorm = 1.0_wp / sumNorm
          zetaZ(i) = sumNorm*( zetaC(i1) + zetaC(i3) + zetaC(i4) + zetaC( tw(i3) ) )
          deltaZ(i) = sumNorm*( deltaC(i1) + deltaC(i3) + deltaC(i4) + deltaC(tw(i3)) )
        enddo

!#ifdef SEAICE_ALLOW_CLIPZETA
!!     regularize zeta if necessary
!        do i = 1,nlpb
!          w = tw(i)
!          s = ts(i)
!
!          zetaC(i)  = MAX(zMin(i),MIN(zMax(i) ,zetaC(i)))
!
!!     zMin, zMax are defined at C-points, make sure that they are not
!!     masked by boundaries/land points
!          zMaxZ = MAX( MAX(zMax(i),zMax(s)), MAX(zMax(w),zMax(ts(w))) )
!          zMinZ = MAX( MAX(zMin(i),zMin(s)), MAX(zMin(w),zMin(ts(w))) )
!          zetaZ(i) = MAX(zMinZ,MIN(zMaxZ,zetaZ(i)))
!        enddo
!#endif /* SEAICE_ALLOW_CLIPZETA */

!     recompute pressure
!YY: replacement pressure
!$acc loop
        do i = 1,nlpb
          pressC(i) = ( press0(i) * ( 1.0_wp - SEAICEpressReplFac )                    &
              + TWO*zetaC(i)*deltaC(i)*SEAICEpressReplFac /(1.0_wp + tensileStrFac(i)) &
              ) * (1.0_wp - tensileStrFac(i))
        enddo

        IF ( SEAICEuseTEM ) THEN
!$acc loop private(etaDenC,zetaMaxC)
          do i = 1,nlpb
            etaDenC   = em(i)**2 + 4.0_wp * e12Csq(i)
            etaDenC  = SQRT(MAX(deltaMinSq,etaDenC))
            zetaMaxC = ecc2*zetaC(i)*(deltaC(i)-ep(i))/etaDenC
            seaice_div(i) = ( 2.0_wp *zetaC(i)*ep(i) - pressC(i) ) * HEFFM(i)
            seaice_tension(i) = 2.0_wp*MIN(zetaC(i),zetaMaxC) * em(i) * HEFFM(i)
          enddo
          !YUAN, ZY: etaDenZ defined at Z-grid
!$acc loop private(i1,i2,i3,i4,sumNorm,etaDenZ,zetaMaxZ)
          do i = 1,nlpbz
            i1 = z1(i, 1)
            i2 = z2(i, 1)
            i3 = z3(i, 1)
            i4 = z4(i, 1)

            sumNorm = 0.25_wp
!     Averaging the squares gives more accurate viscous-plastic behavior
!     than squaring the averages
!YY,ZY: ep em at C-points, etaDenZ at Z-points
            etaDenZ = sumNorm * recip_rAz(i) *      &
                      ( rAc(  i1  ) * em(  i1  )**2   &
                      + rAc(  i3  ) * em(  i3  )**2   &
                      + rAc(  i4  ) * em(  i4  )**2   &
                      + rAc(tw(i3)) * em(tw(i3))**2 ) &
                      + 4.0_wp*e12(i)**2
            etaDenZ  = SQRT(MAX(deltaMinSq,etaDenZ))
            zetaMaxZ = ecc2*zetaZ(i) * ( deltaZ(i)  &
                - sumNorm * ( ep( i1 ) + ep(  i3  )   &
                            + ep( i4 ) + ep(tw(i3)) ) )/etaDenZ
            seaice_shear(i) = 2.0_wp*MIN(zetaZ(i),zetaMaxZ) * 2.0_wp*e12(i)
          enddo
        ELSE     ! YY: else for useTEM
!$acc loop
          do i = 1,nlpb
            seaice_div(i) = ( 2.0_wp *zetaC(i)*ep(i) - pressC(i) ) * HEFFM(i)
            seaice_tension(i) = 2.0_wp*zetaC(i) * em(i) * HEFFM(i)
          enddo
!$acc loop
          do i = 1,nlpbz
            seaice_shear(i) = 2.0_wp*zetaZ(i)*e12(i)
          enddo
        ENDIF    !YY: endif for useTEM
!$acc end kernels

!     first step stress equations
!
        IF ( useAdaptiveEVP ) THEN
!$acc kernels
!$acc loop private(denom1,denom2)
          do i = 1,nlpb
            denom1 = 1.0_wp / evpAlphaC(i)
            denom2 = denom1
!         sigma1 and sigma2 are computed on C points
            seaice_sigma1(i) =                                     &
              ( seaice_sigma1(i)*( evpAlphaC(i) - evpRevFac ) + seaice_div(i) ) &
              * denom1*HEFFM(i)
            seaice_sigma2(i) =                                     &
              ( seaice_sigma2(i)*( evpAlphaC(i) - evpRevFac ) + seaice_tension(i)*recip_evpRevFac ) &
              * denom2*HEFFM(i)
!     Code to avoid very small numbers that can degrade performance.
!     Many compilers can handle this more efficiently with the help of
!     a flag (copied from CICE after correspondence with Elizabeth Hunke)
            seaice_sigma1(i) = SIGN(MAX(                                  &
               ABS( seaice_sigma1(i) ), SEAICE_EPS ), seaice_sigma1(i) )
            seaice_sigma2(i) = SIGN(MAX(                                  &
               ABS( seaice_sigma2(i) ), SEAICE_EPS ), seaice_sigma2(i) )
!
!     recover sigma11 and sigma22
            sig11(i) = 0.5_wp*( seaice_sigma1(i)+seaice_sigma2(i) )
            sig22(i) = 0.5_wp*( seaice_sigma1(i)-seaice_sigma2(i) )
          enddo
!$acc end kernels
        ELSE
!$acc kernels copyin(denom1,denom2)
!$acc loop     ! YY: OPENACC: here may be a problem: denom1
          do i = 1,nlpb
!     sigma1 and sigma2 are computed on C points
            seaice_sigma1(i) =                                     &
              ( seaice_sigma1(i)*( evpAlphaC(i) - evpRevFac ) + seaice_div(i) ) &
              * denom1*HEFFM(i)
            seaice_sigma2(i) =                                     &
              ( seaice_sigma2(i)*( evpAlphaC(i) - evpRevFac ) + seaice_tension(i)*recip_evpRevFac ) &
              * denom2*HEFFM(i)
!     Code to avoid very small numbers that can degrade performance.
!     Many compilers can handle this more efficiently with the help of
!     a flag (copied from CICE after correspondence with Elizabeth Hunke)
            seaice_sigma1(i) = SIGN(MAX(                                  &
               ABS( seaice_sigma1(i) ), SEAICE_EPS ), seaice_sigma1(i) )
            seaice_sigma2(i) = SIGN(MAX(                                  &
               ABS( seaice_sigma2(i) ), SEAICE_EPS ), seaice_sigma2(i) )
!
!     recover sigma11 and sigma22
            sig11(i) = 0.5_wp*( seaice_sigma1(i)+seaice_sigma2(i) )
            sig22(i) = 0.5_wp*( seaice_sigma1(i)-seaice_sigma2(i) )
          enddo
!$acc end kernels
        ENDIF

!     sigma12 is computed on Z points
        IF ( useAdaptiveEVP ) THEN
!$acc kernels loop private(i1,i2,i3,i4,denom2)
          do i = 1,nlpbz
            i1 = z1(i, 1)
            i2 = z2(i, 1)
            i3 = z3(i, 1)
            i4 = z4(i, 1)
            evpAlphaZ(i) = 0.25_wp *                                        &
               ( evpAlphaC(i1)+evpAlphaC(i3)+evpAlphaC(i4)+evpAlphaC(tw(i3)) )
            denom2 = 1.0_wp /  evpAlphaZ(i)
            seaice_sigma12(i) =                                     &
              ( seaice_sigma12(i)*( evpAlphaZ(i) - evpRevFac ) + seaice_shear(i)*recip_evpRevFac) &
               * denom2
            seaice_sigma12(i) = SIGN(MAX( ABS( seaice_sigma12(i) ), SEAICE_EPS ), seaice_sigma12(i) )
          enddo
!$acc end kernels
        ELSE    ! YY: here denom2(i) is changed to a scalar
!$acc kernels copyin(denom2)
!$acc loop
          do i = 1,nlpbz
            seaice_sigma12(i) =                                     &
              ( seaice_sigma12(i)*( evpAlphaZ(i) - evpRevFac ) + seaice_shear(i)*recip_evpRevFac) &
               * denom2
            seaice_sigma12(i) = SIGN(MAX(ABS( seaice_sigma12(i) ), SEAICE_EPS ), seaice_sigma12(i) )
          enddo
!$acc end kernels
        ENDIF

!     compute divergence of stress tensor

!YY: sigma12 defined at Z-point, sig11 at C-point
!YY: divergence about u/v-point
!YY: formulated by ZY. need to test
!$acc kernels 
!$acc loop private(w,n,e,s)
        do i = 1,nlpb
           w = uw(i,1)
           dyS_dxW(i,1) = dyS(w)
           dyS_dxW(i,2) = dxW(w)
        
           n = un(i,1)
           dxS_dyW(i,1) = dxS(n)
           dxS_dyW(i,2) = dyW(n)
        
           e = ve(i,1)
           dyW_dxS(i,1) = dxS(e)
           dyW_dxS(i,2) = dyW(e)
        
           s = vs(i,1)
           dxW_dyS(i,1) = dyS(s)
           dxW_dyS(i,2) = dxW(s)
        end do
!$acc loop private(w,n,s,e,x1,x2)
        do i = 1,nlpb
          w = tw(i)
          n = zn(i)
          x1 = uw(i,2)
          x2 = un(i,2)
          stressDivergenceX(i) = ( sig11(i) * dyS(i) - sig11(w) * dyS_dxW(i,x1)   &
               + seaice_sigma12(n) * dxS_dyW(i,x2) - seaice_sigma12(i) * dxS(i)   &
               ) * recip_rAw(i)

          s = ts(i)
          e = ze(i)
          x1 = vs(i,2)
          x2 = ve(i,2)
          stressDivergenceY(i) = ( sig22(i) * dxW(i) - sig22(s) * dxW_dyS(i,x1)  &
               + seaice_sigma12(e) * dyW_dxS(i,x2) - seaice_sigma12(i) * dyW(i)  &
               ) * recip_rAs(i)
        end do
!$acc end kernels

!#ifdef ALLOW_SEAICE_EVP_RESIDUAL
        IF ( printResidual ) THEN
          resTile = 0.0_wp
!$acc kernels present(sig11pm1,sig22pm1,sig12pm1)
!$acc loop
          do i = 1,nlpb
            sig11Pm1(i) = seaice_sigma1(i)-sig11pm1(i)
            sig22Pm1(i) = seaice_sigma2(i)-sig22pm1(i)
            sig11Pm1(i) = evpAlphaC(i) * sig11Pm1(i)
            sig22Pm1(i) = evpAlphaC(i) * sig22Pm1(i)
          enddo
!$acc loop
          do i = 1,nlpbz
            sig12Pm1(i) = seaice_sigma12(i)-sig12Pm1(i)
            sig12Pm1(i) = evpAlphaZ(i) * sig12Pm1(i)
          enddo
          IF ( .NOT. SEAICEscaleSurfStress ) THEN
              !YY, ZY: remember sig12Pm1 dimension is nlpbz, TO DO LATER
!$acc loop reduction(+:resTile)
            do i = 1,nlpb
!     multiply with mask (concentration) to only count ice contributions
              resTile = resTile + AREA(i) *     &
                  ( sig11Pm1(i)*sig11Pm1(i)     &
                  + sig22Pm1(i)*sig22Pm1(i)     &
                  + sig12Pm1(i)*sig12Pm1(i) )
            enddo
          ELSE
!     in this case the scaling with AREA is already done
!$acc loop reduction(+:resTile)
            do i = 1,nlpb
              resTile = resTile                 &
                  + sig11Pm1(i)*sig11Pm1(i)     &
                  + sig22Pm1(i)*sig22Pm1(i)     &
                  + sig12Pm1(i)*sig12Pm1(i)
            enddo
          ENDIF
!$acc end kernels

          CALL mpi_comp_dp_allsum(resTile)
          resTile = SQRT(resTile)
          if(mpi_rank==0) write(*,*) 'SEAICE_EVP: iEVPstep, residual sigma = ', iEVPstep, resTile
        ENDIF
!#endif /* ALLOW_SEAICE_EVP_RESIDUAL */

!     set up rhs for stepping the velocity field

        call seaice_oceandrag_coeffs(uice, vice, HEFFM, DWATN, ievpstep)

#ifdef SEAICE_ALLOW_BOTTOMDRAG
        call seaice_bottomdrag_coeffs(uIce, vIce, HEFFM, &
            HEFF, AREA, CbotC) 
#endif /* SEAICE_ALLOW_BOTTOMDRAG */

!$acc kernels copyin(kSrf,COSWAT,SINWAT)
!$acc loop
        do i = 1,nlpb
          vecFldLoc(i, iu) = uFld(i,kSrf)
          vecFldLoc(i, iv) = vFld(i,kSrf)
          vecIceLoc(i, iu) = uIce(i)
          vecIceLoc(i, iv) = vIce(i)
        end do

!$acc loop private(w,j1, j2, j3, j4,k1, k2, k3, k4,p1, p2, p3, p4,locMaskU)
        do i = 1,nlpb
          w = tw(i)
          j1 = au1(i, 1)
          j2 = au2(i, 1)
          j3 = au3(i, 1)
          j4 = au4(i, 1)
          k1 = au1(i, 2)
          k2 = au2(i, 2)
          k3 = au3(i, 2)
          k4 = au4(i, 2)
          p1 = au1(i, 3)
          p2 = au2(i, 3)
          p3 = au3(i, 3)
          p4 = au4(i, 3)
!     over open water, all terms that contain sea ice mass drop out
!     the balance is determined by the atmosphere-ice and ice-ocean stress;
!     the staggering of uIce and vIce can cause stripes in the open ocean
!     solution when the turning angle is non-zero (SINWAT.NE.0);
!     we mask this term here in order to avoid the stripes and because
!     over open ocean, u/vIce do not advect anything, so that the associated
!     error is small and most likely only confined to the ice edge but may
!     propagate into the ice covered regions.
          locMaskU = seaiceMassU(i)
          IF ( locMaskU .NE. 0.0_wp ) locMaskU = 1.0_wp
!     set up anti symmetric drag force and add in ice ocean stress
!     ( remember to average to correct velocity points )
!YY,ZY: confirm again. forcex/y at U/V-point
          FORCEX(i)=FORCEX0(i)+                    &
           (                                       &
            0.5_wp * ( DWATN(i)+DWATN(w) ) *       &
            COSWAT * vecFldLoc(i,iu)               &
            - SIGN(SINWAT, fCori(i))* 0.5_wp *     &
            ( DWATN(i) * 0.5_wp *                  & 
             (vecFldLoc(j4,k4)*p4-vecIceLoc(j4,k4)*p4  &
             +vecFldLoc(j2,k2)*p2-vecIceLoc(j2,k2)*p2) & 
            + DWATN(w) * 0.5_wp *                      &
             (vecFldLoc(j3,k3)*p3-vecIceLoc(j3,k3)*p3  &
             +vecFldLoc(j1,k1)*p1-vecIceLoc(j1,k1)*p1) )*locMaskU  &
              ) * areaW(I)
!     coriols terms
          FORCEX(i) = FORCEX(i) + HALF*                &
              ( seaiceMassC(i) * fCori(i)              &
               *0.5_wp*( vecIceLoc(j4,k4)*p4+vecIceLoc(j2,k2)*p2 )   &
              + seaiceMassC(w) * fCori(w)              &
               *0.5_wp*( vecIceLoc(j3,k3)*p3+vecIceLoc(j1,k1)*p1 ) )
        enddo

!$acc loop private(s,j1, j2, j3, j4,k1, k2, k3, k4,p1, p2, p3, p4,locMaskV)
        do i = 1,nlpb
          s = ts(i)
          j1 = av1(i, 1)
          j2 = av2(i, 1)
          j3 = av3(i, 1)
          j4 = av4(i, 1)
          k1 = av1(i, 2)
          k2 = av2(i, 2)
          k3 = av3(i, 2)
          k4 = av4(i, 2)
          p1 = av1(i, 3)
          p2 = av2(i, 3)
          p3 = av3(i, 3)
          p4 = av4(i, 3)

          locMaskV = seaiceMassV(i)
          IF ( locMaskV .NE. 0.0_wp ) locMaskV = 1.0_wp

          FORCEY(i)=FORCEY0(i)+                                &
           (                                                   &
            0.5_wp * ( DWATN(i)+DWATN(s) ) *                   &
           COSWAT * vecFldLoc(i,iv)                            &
           + SIGN(SINWAT, fCori(i)) * 0.5_wp *                 &
           ( DWATN(i) * 0.5_wp *                               &
            (vecFldLoc(j1,k1)*p1-vecIceLoc(j1,k1)*p1           &
            +vecFldLoc(j3,k3)*p3-vecIceLoc(j3,k3)*p3)          &
           + DWATN(s) * 0.5_wp *                               &
            (vecFldLoc(j2,k2)*p2-vecIceLoc(j2,k2)*p2           &
            +vecFldLoc(j4,k4)*p4-vecIceLoc(j4,k4)*p4) )*locMaskV  & 
             ) * areaS(i)
          FORCEY(i) = FORCEY(i) - HALF*                        &
              ( seaiceMassC(i) * fCori(i)                      &
              *0.5_wp*( vecIceLoc(j1,k1)*p1+vecIceLoc(j3,k3)*p3 ) &
              + seaiceMassC(s) * fCori(s)                         &
              *0.5_wp*( vecIceLoc(j2,k2)*p2+vecIceLoc(j4,k4)*p4 ) )
        enddo
!$acc end kernels

#ifdef SEAICE_ALLOW_MOM_ADVECTION
!YY: ADVECTION TERMS ARE OMITTED IN EVP ALGORITHM
        IF ( SEAICEmomAdvection ) THEN
         do i = 1,nlpb
           gUmom(i) = 0.0_wp
           gVmom(i) = 0.0_wp
         enddo
         call seaice_mom_advection(uIce, vIce, gUmom, gVmom,iEVPstep, nEVPstep)
         do i = 1,nlpb
           FORCEX(i) = FORCEX(i) + gUmom(i)
           FORCEY(i) = FORCEY(i) + gVmom(i)
         enddo
        ENDIF
#endif /* SEAICE_ALLOW_MOM_ADVECTION */

!     step momentum equations with ice-ocean stress term treated implicitly

        IF ( useAdaptiveEVP ) THEN
!$acc kernels loop private(w,s)
          do i = 1,nlpb
            w = tw(i)
            s = ts(i)
!     compute and adjust parameters that are constant otherwise
            evpBetaU(i) = 0.5_wp*(evpAlphaC(w)+evpAlphaC(i))
            evpBetaV(i) = 0.5_wp*(evpAlphaC(s)+evpAlphaC(i))
          enddo
!$acc end kernels
        ENDIF

!YY: at U-point, average of C-point
!$acc kernels copyin(COSWAT,evpStarFac,recip_deltaT)
!$acc loop private(w,s,betaFacU,betaFacV,tmp,betaFacP1V,betaFacP1U,denomU,denomV)
        do i = 1,nlpb
          w = tw(i)
          s = ts(i)

          betaFacU   = evpBetaU(i)*recip_deltaT
          betaFacV   = evpBetaV(i)*recip_deltaT
          tmp        = evpStarFac*recip_deltaT
          betaFacP1V = betaFacV + tmp
          betaFacP1U = betaFacU + tmp
          denomU = seaiceMassU(i)*betaFacP1U     &
              + 0.5_wp*( DWATN(i) + DWATN(w) )   &
              * COSWAT * areaW(i)
          denomV = seaiceMassV(i)*betaFacP1V     &
              + 0.5_wp*( DWATN(i) + DWATN(s) )   &
              * COSWAT * areaS(i)
#ifdef SEAICE_ALLOW_BOTTOMDRAG
          denomU = denomU + areaW(I) * 0.5_wp*( CbotC(I) + CbotC(w) )
          denomV = denomV + areaS(I) * 0.5_wp*( CbotC(I) + CbotC(s) )
#endif /* SEAICE_ALLOW_BOTTOMDRAG */
          IF ( denomU .EQ. 0.0_wp ) denomU = 1.0_wp
          IF ( denomV .EQ. 0.0_wp ) denomV = 1.0_wp
          uIce(i) = seaiceMaskU(i) *                                  &
              ( seaiceMassU(i)*betaFacU * uIce(i)                     &
              + seaiceMassU(i)*recip_deltaT*evpStarFac * uIceNm1(i)   &
              + FORCEX(i)                                             &
              + stressDivergenceX(i) ) / denomU
          vIce(i) = seaiceMaskV(i) *                                  &
              ( seaiceMassV(i)*betaFacV * vIce(i)                     &
              + seaiceMassV(i)*recip_deltaT*evpStarFac * vIceNm1(i)   &
              + FORCEY(i)                                             &
              + stressDivergenceY(i) ) / denomV  
          !term1(i) = seaiceMassU(i)*betaFacU * uIce(i)
          !term2(i) = seaiceMassU(i)*recip_deltaT*evpStarFac * uIceNm1(i)
          !term3(i) = FORCEX(i)
          !term4(i) = stressDivergenceX(i)
         ! term5(i) = seaiceMaskU(i)
          !term6(i) = denomU
         ! term7(i) = term1(i)+term2(i)+term3(i)+term4(i)
          !term8(i) = seaiceMassU(i)*betaFacP1U

          !term9(i) = 0.5_wp*( DWATN(i) + DWATN(w) )* COSWAT * areaW(i)
          !term10(i) = areaW(I) * 0.5_wp*( CbotC(I) + CbotC(w) )

          !term10(i) = areaW(i)
          !term11(i) = seaiceMassU(i)
          !term12(i) = betaFacP1U
          !term1(i) = seaiceMassV(i)*betaFacV * vIce(i)
          !term2(i) = seaiceMassV(i)*recip_deltaT*evpStarFac * vIceNm1(i)
          !term3(i) = FORCEY(i)
          !term4(i) = stressDivergenceY(i)
          !term5(i) = seaiceMaskV(i)
          !term6(i) = denomV
          !term7(i) = term1(i)+term2(i)+term3(i)+term4(i)
        enddo
!$acc end kernels
#ifdef SeaiceDebug
      !if (myIter==4 .and. mpi_rank==14) write(mitice_runmsg_unit,*) 'term1: ', myIter, term1(2367)
      !if (myIter==4 .and. mpi_rank==14) write(mitice_runmsg_unit,*) 'term2: ', myIter, term2(2367)
      !if (myIter==4 .and. mpi_rank==14) write(mitice_runmsg_unit,*) 'term3: ', myIter, term3(2367)
      !if (myIter==4 .and. mpi_rank==14) write(mitice_runmsg_unit,*) 'term4: ', myIter, term4(2367)
     !! if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term5: ', myIter, maxval(abs(term5(1:loc_nlpb)))
      !if (myIter==4 .and. mpi_rank==14) write(mitice_runmsg_unit,*) 'denom: ', myIter, abs(term6(2367))
      !if (myIter==4 .and. mpi_rank==14) print *, 'term1: ', myIter, term1(2367)
      !if (myIter==4 .and. mpi_rank==14) print *, 'term2: ', myIter, term2(2367)
      !if (myIter==4 .and. mpi_rank==14) print *, 'term3: ', myIter, term3(2367)
      !if (myIter==4 .and. mpi_rank==14) print *, 'term4: ', myIter, term4(2367)
      !if (myIter==4 .and. mpi_rank==14) print *,  'denom: ', myIter, abs(term6(2367))
    !  if (iEVPstep == nEVPstep) then
    !   ! if ( mpi_rank==14) print *, 'term1: ', myIter, term1(2367)
    !   ! if ( mpi_rank==14) print *, 'term2: ', myIter, term2(2367)
    !   ! if ( mpi_rank==14) print *, 'term3: ', myIter, term3(2367)
    !   ! if ( mpi_rank==14) print *, 'term4: ', myIter, term4(2367)
    !   ! if ( mpi_rank==14) print *, 'denom: ', myIter, abs(term6(2367))
    !   ! if ( mpi_rank==14) print *, 'term8: ', myIter, term8(2367)

    !    print *, 'term9: ', myIter, maxval(term9)
    !    print *, 'term10: ', myIter, maxval(term10)

    !   ! if ( mpi_rank==14) print *, 'term11: ', myIter, term11(2367)
    !   ! if ( mpi_rank==14) print *, 'term12: ', myIter, term12(2367)
    !  endif
     ! if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term sum: ', myIter, abs(term7(2754))
     ! if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term8: ', myIter, abs(term8(2754))
     ! if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term9: ', myIter, abs(term9(2754))
     ! if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term10: ', myIter, abs(term10(2754))
     ! if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term11: ', myIter, abs(term11(2754))
     ! if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term12: ', myIter, abs(term12(2754))

      !if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term1: ', myIter, abs(term1(2754))
      !if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term2: ', myIter, abs(term2(2754))
      !if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term3: ', myIter, abs(term3(2754))
      !if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term4: ', myIter, abs(term4(2754))
      !if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term5: ', myIter, abs(term5(2754))
      !if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'denom: ', myIter, abs(term6(2754))
      !if (myIter==5 .and. mpi_rank==21) write(mitice_runmsg_unit,*) 'term sum: ', myIter, abs(term7(2754))
    !  velIceMag = sqrt(uIce*uIce+vIce*vIce)
    !  velIceMax = maxval(velIceMag(1:loc_nlpb))
     !if (myIter==4) write(mitice_runmsg_unit,*) 'uIceMax: ',myIter,mpi_rank, velIceMax,maxloc(abs(velIceMag))
     !if (myIter==4)  print *, 'uIceMax: ',myIter,mpi_rank, velIceMax,maxloc(abs(velIceMag))
      !if (iEVPstep == nEVPstep) then
        !if ( mpi_rank==14) print *, 'uIceMax: ',myIter,velIceMag(2367),maxloc(velIceMag)
      !endif
      !if (myIter==4 .and. mpi_rank==14) write(mitice_runmsg_unit,*) 'uIceMax: ',myIter,velIceMax,maxloc(velIceMag)
#endif


!YY: comment it for now
!        do i = 1,nlpb
!          w = tw(i)
!          s = ts(u)
!          
!          locMaskU = maskInC(i)*maskInC(w)
!          locMaskV = maskInC(i)*maskInC(s)
!          uIce(i) = uIce(i)*locMaskU + uIceNm1(i)*(ONE-locMaskU)
!          vIce(i) = vIce(i)*locMaskV + vIceNm1(i)*(ONE-locMaskV)
!        enddo

!$acc kernels loop private(w,s,locMaskU,locMaskV)
        do i = 1,nlpb
          w = tw(i)
          s = ts(i)
          
          locMaskU = HEFFM(i)*HEFFM(w)
          locMaskV = HEFFM(i)*HEFFM(s)
          uIce(i) = uIce(i)*locMaskU + uIceNm1(i)*(ONE-locMaskU)
          vIce(i) = vIce(i)*locMaskV + vIceNm1(i)*(ONE-locMaskV)
        enddo
!$acc end kernels

        call mpi_data_exchange(mpi_sendbuf_1d,uIce)
        call mpi_data_exchange(mpi_sendbuf_1d,vIce)


!#ifdef ALLOW_SEAICE_EVP_RESIDUAL
        IF ( printResidual ) THEN
          resTile = 0.0_wp
!$acc kernels present(uIcePm1,vIcePm1)
!$acc loop
          do i = 1,nlpb
            uIcePm1(i) = seaiceMaskU(i) * ( uIce(i)-uIcePm1(i) ) * evpBetaU(i)
            vIcePm1(i) = seaiceMaskV(i) * ( vIce(i)-vIcePm1(i) ) * evpBetaV(i)
          enddo
          IF ( .NOT. SEAICEscaleSurfStress ) THEN
!     multiply with mask (concentration) to only count ice contributions
!$acc loop reduction(+:resTile)
            do i = 1,nlpb
              resTile = resTile + AREA(i) * (uIcePm1(i)*uIcePm1(i)+vIcePm1(i)*vIcePm1(i))
            enddo
          ELSE
!$acc loop reduction(+:resTile)
            do i = 1,nlpb
              resTile = resTile + uIcePm1(i)*uIcePm1(i) + vIcePm1(i)*vIcePm1(i)
            enddo
          ENDIF
!$acc end kernels
          CALL mpi_comp_dp_allsum(resTile)
          resTile = SQRT(resTile)
          write(*,*) 'SEAICE_EVP: iEVPstep, residual U = ', iEVPstep, resTile        
        ENDIF

!#endif /* ALLOW_SEAICE_EVP_RESIDUAL */

       ENDIF   !ENDIF /(iEVPstep.LE.nEVPstep)/
      ENDDO    !ENDDO /iEVPstep = 1, nEVPstepMax/

!YY
      if (printResidual) then
        deallocate(uIcePm1, vIcePm1)
        deallocate( sig11pm1, sig22pm1, sig12pm1)
      endif
!$acc end data

   end subroutine seaice_evp
#endif /* SEAICE_ALLOW_DYNAMICS */

   subroutine seaice_oceandrag_coeffs(uIceLoc, vIceLoc, HEFFMLoc,CwatC,iStep)

!  *==========================================================*
!  | o Compute the drag coefficients for ice-ocean drag,
!  |   so that we can use the same code for different solvers
!  *==========================================================*
!  | written by Martin Losch, Oct 2012
!  *==========================================================*

      implicit none

!     iStep  :: current sub-time step iterate
      INTEGER iStep
!     u/vIceLoc :: local copies of the current ice velocity
      real(wp)  ::  uIceLoc(nlpb), vIceLoc(nlpb)
!     HEFFMLoc  :: local copy of land-sea masks
      real(wp)  ::  HEFFMLoc(nlpb) 
!     CwatC     :: drag coefficients
      real(wp)  ::  CwatC(nlpb)

      INTEGER i
      INTEGER kSrf
      real(wp) ::  tempVar, tempMin
      real(wp) ::  dragCoeff

      integer :: uel,uei,uep,vnl,vni,vnp
      real(wp)  ::  vecFldLoc(nlpb,2), vecIceLoc(nlpb,2) 
      integer(i2) :: maskWS(nlpb,2)

!$acc data present(uIceLoc, vIceLoc, HEFFMLoc,CwatC,uFld,vFld,maskW,maskS,latC),   &
!$acc      create(vecFldLoc,vecIceLoc,maskWS)

!!YY:  Pressure coords
      kSrf = nk

      tempMin = SEAICEdWatMin*SEAICEdWatMin

!!YY,ZY: check
!$acc kernels copyin(kSrf,tempMin)
!$acc loop
      do i = 1,nlpb
        vecFldLoc(i, iu) = uFld(i,kSrf)
        vecFldLoc(i, iv) = vFld(i,kSrf)
        vecIceLoc(i, iu) = uIceLoc(i)
        vecIceLoc(i, iv) = vIceLoc(i)
        maskWS(i,iu) = maskW(i,kSrf)
        maskWS(i,iv) = maskS(i,kSrf)
      end do
!$acc loop private(uel,uei,uep,vnl,vni,vnp,tempVar,dragCoeff)
      do i = 1,nlpb
!     non-linear water drag coefficients CwatC (DWATN)
!        tempVar = 0.25 _d 0*(
!     &         ( ( uIceLoc(i)-uVel(i,kSrf) )
!     &          *maskInW(i)
!     &          +( uIceLoc(i+1,j,bi,bj)-uVel(i+1,j,kSrf,bi,bj) )
!     &          *maskInW(i+1,j,bi,bj) )**2
!     &       + ( ( vIceLoc(i,j  ,bi,bj)-vVel(i,j  ,kSrf,bi,bj) )
!     &          *maskInS(i, j ,bi,bj)
!     &          +( vIceLoc(i,j+1,bi,bj)-vVel(i,j+1,kSrf,bi,bj) )
!     &          *maskInS(i,j+1,bi,bj) )**2 )
!!YY: why using maskinW/S. anyway, here I use maskW/S instead
!!YY: use maskW/S, not SIMaskU/V
!!YY: maskinW/S specify zero on and beyond OB, so for global domain,
        uel = ue(i, 1)
        uei = ue(i, 2)
        uep = ue(i, 3)
        vnl = vn(i, 1)
        vni = vn(i, 2)
        vnp = vn(i, 3)
        tempVar = 0.25_wp * (                                        &
              ( ( vecIceLoc(i,iu)-vecFldLoc(i,iu) )*maskWS(i,iu)     &
               +( vecIceLoc(uel,uei)*uep - vecFldLoc(uel,uei)*uep )*maskWS(uel,uei) )**2  &
            + ( ( vecIceLoc(i,iv)-vecFldLoc(i,iv) ) *maskWS(i,iv)    &
               +( vecIceLoc(vnl,vni)*vnp-vecFldLoc(vnl,vni)*vnp )*maskWS(vnl,vni) )**2 )
        IF ( latC(i) .LT. ZERO ) THEN
         dragCoeff = SEAICE_waterDrag_south*rau0
        ELSE
         dragCoeff = SEAICE_waterDrag      *rau0
        ENDIF
        CwatC(i) = SEAICEdWatMin
        IF ( dragCoeff*dragCoeff * tempVar .GT. tempMin )    &
          CwatC(i) = dragCoeff*SQRT(tempVar)
        CwatC(i) = CwatC(i) * HEFFMLoc(i)
      enddo
!$acc end kernels
!$acc end data

   end subroutine seaice_oceandrag_coeffs


#ifdef SEAICE_ALLOW_MOM_ADVECTION
   subroutine seaice_mom_advection(uIceLoc, vIceLoc, gU, gV, iEVPstep, nEVPstep)
!     *==========================================================*
!     | o Form the advection of sea ice momentum to be added to  |
!     |   the right hand-side of the momentum equation.          |
!     *==========================================================*
      implicit none

!     gU      :: advection tendency (all explicit terms), u component
!     gV      :: advection tendency (all explicit terms), v component
      real(wp) ::  uIceLoc(nlpb)
      real(wp) ::  vIceLoc(nlpb)
      real(wp) ::  gU(nlpb)
      real(wp) ::  gV(nlpb)

!     == Local variables ==
      real(wp) ::  uCf(nlpb)
      real(wp) ::  vCf(nlpb)
      real(wp) ::  hFacZtmp   (nlpbz)
      real(wp) ::  recip_hFacZtmp (nlpbz)   !YY: declaration needed?? global vars in MaCOM
      !real(wp) ::  KE      (nlpb)
      real(wp) ::  vort3   (nlpbz)

      INTEGER i,kSrf
      LOGICAL vorticityFlag
#ifdef SeaiceDebug
      integer iEVPstep, nEVPstep
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--   Initialise intermediate terms
      do i = 1,nlpb
        uCf  (i) = 0.0_wp
        vCf  (i) = 0.0_wp
        gU   (i) = 0.0_wp
        gV   (i) = 0.0_wp
      enddo
      do i = 1,nlpbz
        vort3(i) = 0.0_wp
      enddo

      kSrf = nk       !YY: surface layer
!--   Calculate open water fraction at vorticity points
      !CALL MOM_CALC_HFACZ(bi,bj,k,hFacZ,r_hFacZ,myThid)
      !YY,ZY: hFacZ/recip_hFacZ are updated with surface layer,  OK?
      !YY: otherwise make a local copy
      call mitice_dyn_mom_update_hFacZ(kSrf,hFacZtmp, recip_hFacZtmp)    

!--   Bernoulli term
      !CALL MOM_CALC_KE(bi,bj,k,SEAICEselectKEscheme,uFld,vFld,KE,myThid)
      call mitice_calc_KE(kSrf,uIceLoc,vIceLoc,uCf,vCf)

      do i = 1,nlpb
        gU(i) = gU(i)+uCf(i)
        gV(i) = gV(i)+vCf(i)
      enddo

      !CALL MOM_CALC_RELVORT3(bi,bj,k,uFld,vFld,hFacZ,vort3,myThid)
      call mitice_dyn_hrelvort(uIceLoc,vIceLoc,vort3)
       
     ! do i = 1,nlpbz
     !    vort3(i) = vort3(i)*maskZ(i,kSrf)
     ! enddo

!--   Horizontal advection of relative (or absolute) vorticity
!YY: here use relative vorticity, so use vort3
      call mitice_dyn_vort_een(kSrf,uIceLoc,vIceLoc,recip_hFacZtmp, vort3,uCf,vCf)

      do i = 1,nlpb
        gU(i) = gU(i)+uCf(i)
        gV(i) = gV(i)+vCf(i)
      enddo

!!YY: here we comment below, with no consideration on OB
!!--   Set du/dt & dv/dt on boundaries to zero
!!     apply masks for interior (important when we have open boundaries)
!      DO j=jMin,jMax
!       DO i=iMin,iMax
!        gU(i,j) = gU(i,j)*maskInW(i,j,bi,bj)
!        gV(i,j) = gV(i,j)*maskInS(i,j,bi,bj)
!       ENDDO
!      ENDDO

   end subroutine seaice_mom_advection

!==============================================================================
   subroutine mitice_dyn_hrelvort(uIcetmp,vIcetmp, Vort3Loc)
!==============================================================================
      implicit none
      integer :: i
      integer :: i1, i2, i3, i4
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4
      real(wp), intent(in) :: uIcetmp(nlpb), vIcetmp(nlpb)
      real(wp), intent(out) :: Vort3Loc(nlpbz)
      real(wp) :: VecFld(nlpb, 2), dlc(nlpb, 2)

      !$ACC kernels present(uFldLoc,vFldLoc,Vort3Loc,dxC,dyC,z1,z2,z3,z4,recip_rAz),   &
      !$ACC         create(VecFld,dlc)
      !$ACC loop independent
      do i = 1, nlpb
         VecFld(i, iu) = uIcetmp(i)
         VecFld(i, iv) = vIcetmp(i)
         dlc(i, iu) = dxC(i)
         dlc(i, iv) = dyC(i)
      end do

      !$ACC loop independent private(i1, i2, i3, i4,j1, j2, j3, j4,k1, k2, k3, k4)
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
         Vort3Loc(i) = (-VecFld(i1, j1)*dlc(i1, j1)*k1 + VecFld(i2, j2)*dlc(i2, j2)*k2 &
                        + VecFld(i3, j3)*dlc(i3, j3)*k3 - VecFld(i4, j4)*dlc(i4, j4)*k4)*recip_rAz(i)

      end do

      !$ACC end kernels

   end subroutine mitice_dyn_hrelvort

!==============================================================================
   subroutine mitice_calc_KE(k,uIcetmp,vIcetmp,dKEdx,dKEdy)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i, w, s, uel, uei, vnl, vni
      real(wp) :: uIcetmp(nlpb), vIcetmp(nlpb)
      real(wp) :: dKEdx(nlpb), dKEdy(nlpb)
      real(wp) :: vecFldLoc(nlpb, 2), rAws(nlpb, 2), hFacWS(nlpb, 2)
      real(wp) :: KE  (nlpb)

      !$ACC data present(uIcetmp,vIcetmp,rAw,rAs,hFacW,hFacS,ue,vn,tw,ts,  &
      !$ACC              KE,dKEdx,dKEdy,   &
      !$ACC              recip_rAc,recip_hFacC,recip_dxC,recip_dyC,maskW,maskS),  &
      !$ACC      create(vecFldLoc,rAws,hFacWS)
      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         KE   (i) = 0.0_wp    !initialize
         vecFldLoc(i, iu) = uIcetmp(i)
         vecFldLoc(i, iv) = vIcetmp(i)
         rAws(i, iu) = rAw(i)
         rAws(i, iv) = rAs(i)
         hFacWS(i, iu) = hFacW(i, k)
         hFacWS(i, iv) = hFacS(i, k)
      end do

      ! KE is define on C grid same with T S
      !$ACC loop independent private(uel,uei,vnl,vni)
      do i = 1, nlpb
         uel = ue(i, 1)
         uei = ue(i, 2)
         vnl = vn(i, 1)
         vni = vn(i, 2)
         KE(i) = 0.25_wp*( &
                 vecFldLoc(i, iu)*vecFldLoc(i, iu)*rAws(i, iu)*hFacWS(i, iu) &
                 + vecFldLoc(uel, uei)*vecFldLoc(uel, uei)*rAws(uel, uei)*hFacWS(uel, uei) &
                 + vecFldLoc(i, iv)*vecFldLoc(i, iv)*rAws(i, iv)*hFacWS(i, iv) &
                 + vecFldLoc(vnl, vni)*vecFldLoc(vnl, vni)*rAws(vnl, vni)*hFacWS(vnl, vni) &
                 )*recip_rAc(i)*recip_hFacC(i, k)
      end do
      !$ACC end kernels

      !$ACC kernels
      !$ACC loop private(w,s)
      do i = 1, nlpb
         w = tw(i)
         !YY: I replaced maskW with SImasku
         !dKEdx(i) = -recip_dxC(i)*(KE(i) - KE(w))*maskW(i, k)
         dKEdx(i) = -recip_dxC(i)*(KE(i) - KE(w))*SIMaskU(i)

         s = ts(i)
         !dKEdy(i) = -recip_dyC(i)*(KE(i) - KE(s))*maskS(i,k)
         dKEdy(i) = -recip_dyC(i)*(KE(i) - KE(s))*SIMaskV(i)
      end do
      !$ACC end kernels
      !$ACC end data

   end subroutine mitice_calc_KE


   subroutine mitice_dyn_vort_een(k,uIcetmp,vIcetmp,recip_hFacZtmp,omega3,uCorT,vCorT)

      implicit none
      integer, intent(in) :: k
      real(wp),intent(in) :: uIcetmp(nlpb), vIcetmp(nlpb)
      real(wp),intent(in) :: omega3(nlpb)
      real(wp), intent(in) :: recip_hFacZtmp(nlpbz)
      real(wp),intent(out)  :: uCorT(nlpb), vCorT(nlpb)
      integer :: i
      integer :: i1, i2, i3, i4, i5, i6
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4
      integer :: p1, p2, p3, p4
      real(wp) :: oneThird, oneFour
      real(wp) :: vort3uw, vort3ue, vort3dw, vort3de
      real(wp) :: VecFld(nlpb,2), dlz(nlpbz,2), hFacWS(nlpb,2)
      real(wp) :: maskW_s, maskS_s
      integer(4) :: maskWS(nlpb,2)

      !$ACC kernels present(uIcetmp,vIcetmp,hFacW,hFacS,dyZ,dxZ,  &
      !$ACC         auz1,auz2,auz3,auz4,auz5,auz6,au1,au2,au3,au4,  &
      !$ACC         recip_hFacZ,omega3,recip_dxC,maskW,  &
      !$ACC         avz1,avz2,avz3,avz4,avz5,avz6,av1,av2,av3,av4,  &
      !$ACC         recip_dyC,maskS,maskZ),  &
      !$ACC         copyin(k),  &
      !$ACC         create(VecFld,dlz,hFacWS,maskWS)

      !$ACC loop independent
      do i = 1, nlpb
         VecFld(i, iu) = uIcetmp(i)
         VecFld(i, iv) = vIcetmp(i)
         hFacWS(i, iu) = hFacW(i, k)
         hFacWS(i, iv) = hFacS(i, k)
         !maskWS(i, iu) = maskW(i, k)
         !maskWS(i, iv) = maskS(i, k)
         maskWS(i, iu) = SIMaskU(i)
         maskWS(i, iv) = SIMaskV(i)
      end do

      !$ACC loop independent
      do i = 1, nlpbz
         dlz(i, iu) = dyZ(i)
         dlz(i, iv) = dxZ(i)
      end do

      !$ACC loop independent private(i1, i2, i3, i4, i5, i6, j1, j2, j3, j4,  &
      !$ACC                  k1, k2, k3, k4, p1, p2, p3, p4, oneThird, oneFour, &
      !$ACC                  vort3uw, vort3dw, vort3ue, vort3de)
      do i = 1, nlpb
         i1 = auz1(i)
         i2 = auz2(i)
         i3 = auz3(i)
         i4 = auz4(i)
         i5 = auz5(i)
         i6 = auz6(i)
         j1 = au1(i, 1)
         j2 = au2(i, 1)
         j3 = au3(i, 1)
         j4 = au4(i, 1)
         k1 = au1(i, 2)
         k2 = au2(i, 2)
         k3 = au3(i, 2)
         k4 = au4(i, 2)
         p1 = au1(i, 3)
         p2 = au2(i, 3)
         p3 = au3(i, 3)
         p4 = au4(i, 3)

         oneThird = maskZ(i1,k)+maskZ(i2,k)+maskZ(i5,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3uw = (recip_hFacZtmp(i1)*omega3(i1) &
                    + recip_hFacZtmp(i2)*omega3(i2) &
                    + recip_hFacZtmp(i5)*omega3(i5)) &
                   *oneThird*VecFld(j1, k1)*p1*dlz(j1, k1)*hFacWS(j1, k1)

         oneThird = maskZ(i3,k)+maskZ(i2,k)+maskZ(i5,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3ue = (recip_hFacZtmp(i3)*omega3(i3) &
                    + recip_hFacZtmp(i2)*omega3(i2) &
                    + recip_hFacZtmp(i5)*omega3(i5)) &
                   *oneThird*VecFld(j2, k2)*p2*dlz(j2, k2)*hFacWS(j2, k2)

         oneThird = maskZ(i4,k)+maskZ(i5,k)+maskZ(i2,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3dw = (recip_hFacZtmp(i4)*omega3(i4) &
                    + recip_hFacZtmp(i5)*omega3(i5) &
                    + recip_hFacZtmp(i2)*omega3(i2)) &
                   *oneThird*VecFld(j3, k3)*p3*dlz(j3, k3)*hFacWS(j3, k3)

         oneThird = maskZ(i6,k)+maskZ(i5,k)+maskZ(i2,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3de = (recip_hFacZtmp(i6)*omega3(i6) &
                    + recip_hFacZtmp(i5)*omega3(i5) &
                    + recip_hFacZtmp(i2)*omega3(i2)) &
                   *oneThird*VecFld(j4, k4)*p4*dlz(j4, k4)*hFacWS(j4, k4)

         oneFour = maskWS(j1,k1) + maskWS(j2,k2) + maskWS(j3,k3) + maskWS(j4,k4)
         if (oneFour .eq. 0.0_wp) then
            oneFour = 0.0_wp
         else
            oneFour = 1.0_wp/oneFour
         end if
         !uCorT(i) = +(vort3uw + vort3ue + vort3dw + vort3de)*oneFour*recip_dxC(i)*maskW(i, k)
         uCorT(i) = +(vort3uw + vort3ue + vort3dw + vort3de)*oneFour*recip_dxC(i)*SIMaskU(i)
      end do

      !$ACC loop independent private(i1, i2, i3, i4, i5, i6, j1, j2, j3, j4,  &
      !$ACC                  k1, k2, k3, k4, p1, p2, p3, p4, oneThird, oneFour, &
      !$ACC                  vort3uw, vort3dw, vort3ue, vort3de)
      do i = 1, nlpb
         i1 = avz1(i)
         i2 = avz2(i)
         i3 = avz3(i)
         i4 = avz4(i)
         i5 = avz5(i)
         i6 = avz6(i)
         j1 = av1(i, 1)
         j2 = av2(i, 1)
         j3 = av3(i, 1)
         j4 = av4(i, 1)
         k1 = av1(i, 2)
         k2 = av2(i, 2)
         k3 = av3(i, 2)
         k4 = av4(i, 2)
         p1 = av1(i, 3)
         p2 = av2(i, 3)
         p3 = av3(i, 3)
         p4 = av4(i, 3)

         oneThird = maskZ(i1,k)+maskZ(i2,k)+maskZ(i5,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3uw = (recip_hFacZtmp(i1)*omega3(i1) &
                    + recip_hFacZtmp(i2)*omega3(i2) &
                    + recip_hFacZtmp(i5)*omega3(i5)) &
                   *oneThird*VecFld(j1, k1)*p1*dlz(j1, k1)*hFacWS(j1, k1)

         oneThird = maskZ(i3,k)+maskZ(i2,k)+maskZ(i5,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3dw = (recip_hFacZtmp(i3)*omega3(i3) &
                    + recip_hFacZtmp(i2)*omega3(i2) &
                    + recip_hFacZtmp(i5)*omega3(i5)) &
                   *oneThird*VecFld(j2, k2)*p2*dlz(j2, k2)*hFacWS(j2, k2)

         oneThird = maskZ(i4,k)+maskZ(i5,k)+maskZ(i2,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3ue = (recip_hFacZtmp(i4)*omega3(i4) &
                    + recip_hFacZtmp(i5)*omega3(i5) &
                    + recip_hFacZtmp(i2)*omega3(i2)) &
                   *oneThird*VecFld(j3, k3)*p3*dlz(j3, k3)*hFacWS(j3, k3)

         oneThird = maskZ(i6,k)+maskZ(i5,k)+maskZ(i2,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3de = (recip_hFacZtmp(i6)*omega3(i6) &
                    + recip_hFacZtmp(i5)*omega3(i5) &
                    + recip_hFacZtmp(i2)*omega3(i2)) &
                   *oneThird*VecFld(j4, k4)*p4*dlz(j4, k4)*hFacWS(j4, k4)

         oneFour = maskWS(j1,k1) + maskWS(j2,k2) + maskWS(j3,k3) + maskWS(j4,k4)
         if (oneFour .eq. 0.0_wp) then
            oneFour = 0.0_wp
         else
            oneFour = 1.0_wp/oneFour
         end if
         !vCorT(i) = -(vort3uw + vort3dw + vort3ue + vort3de)*oneFour*recip_dyC(i)*maskS(i, k)
         vCorT(i) = -(vort3uw + vort3dw + vort3ue + vort3de)*oneFour*recip_dyC(i)*SIMaskV(i)
      end do

      !$ACC end kernels

   end subroutine mitice_dyn_vort_een

!==============================================================================
   subroutine mitice_dyn_mom_update_hFacZ(k,hFacZtmp,recip_hFacZtmp)
!==============================================================================
      implicit none
      integer :: i, k
      integer :: i1, i2, i3, i4
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4
      real(wp) :: hFacZOpen, hFacWS(nlpb, 2)
      real(wp) :: hFacZtmp(nlpbz), recip_hFacZtmp(nlpbz)

      !$ACC kernels present(nlpbz,hFacZ,recip_hFacZ,hFacW,hFacS,z1,z2,z3,z4),  &
      !$ACC         copyin(k),  &
      !$ACC         create(hFacWS)
      !$ACC loop
      do i = 1, nlpbz
         hFacZtmp(i) = 0.0_wp
         recip_hFacZtmp(i) = 0.0_wp
      end do

      !$ACC loop
      do i = 1, nlpb
         hFacWS(i, iu) = hFacW(i, k)
         hFacWS(i, iv) = hFacS(i, k)
      end do

      !$ACC loop private(i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,k4,hFacZOpen)
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
         hFacZOpen = min(hFacWS(i1, j1)*abs(k1), hFacWS(i2, j2)*abs(k2))
         hFacZOpen = min(hFacZOpen, hFacWS(i3, j3)*abs(k3))
         if (k4 .ne. 0) then
            hFacZOpen = min(hFacZOpen, hFacWS(i4, j4)*abs(k4))
         end if
         hFacZtmp(i) = hFacZOpen
      end do

      !$ACC loop
      do i = 1, nlpbz
         if (hFacZtmp(i) .ne. 0) then
            recip_hFacZtmp(i) = 1.0_wp/hFacZtmp(i)
         else
            recip_hFacZtmp(i) = 0.0_wp
         end if
      end do
      !$ACC end kernels

   end subroutine mitice_dyn_mom_update_hFacZ
#endif /* SEAICE_ALLOW_MOM_ADVECTION */

   subroutine seaice_ocean_stress(windTauX, windTauY)
!     *==========================================================*
!     | o Calculate ocean surface stresses                       |
!     *==========================================================*
!YY: update ocean surface wind stress due to update of sea ice area

      implicit none

!     windTauX :: X-direction wind stress over seaice at U point
!     windTauY :: Y-direction wind stress over seaice at V point
      real(wp) :: windTauX(nlpb), windTauY(nlpb)

!     kSrf      :: vertical index of surface layer
      INTEGER i
      INTEGER kSrf
      real(wp)  COSWAT
      real(wp)  SINWAT
      real(wp)  fuIceLoc, fvIceLoc
      real(wp)  areaW, areaS
      real(wp) :: vecIceLoc(nlpb,2), vecFldLoc(nlpb,2)
      integer :: w, s
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4
      integer :: p1, p2, p3, p4

!$acc data present(windTauX, windTauY,AREA,utau,vtau,uFld,vFld,uIce,vIce,  &
!$acc              stressDivergenceX,stressDivergenceY,au1,au2,au3,au4,    &
!$acc              fCori,DWATN,av1,av2,av3,av4,utau_ice,vtau_ice),   &
!$acc      create(vecIceLoc,vecFldLoc)

!     surrface level
      kSrf = nk

!     introduce turning angle (default is zero)
      SINWAT=SIN(SEAICE_waterTurnAngle*radPi)
      COSWAT=COS(SEAICE_waterTurnAngle*radPi)

!$acc kernels copyin(KSrf,SINWAT,COSWAT)
      !YY: default value is false
      IF ( useHB87StressCoupling ) THEN

!     use an intergral over ice and ocean surface layer to define
!     surface stresses on ocean following Hibler and Bryan (1987, JPO)

!$acc loop private(w,s,areaW,areaS)
        do i = 1,nlpb
          w = tw(i)
          s = ts(i)
!     average wind stress over ice and ocean and apply averaged wind
!     stress and internal ice stresses to surface layer of ocean
!YY: make sure stressDivergenceX/Y defined in U/V-grid
!YY: here fu/fv are replaced by utau/vtau
!YY: It is very important that here utau include seaice internal stress
!YY: while in MaCOM utau is just wind stress over grid.
!YY: so, for ocean surface circulation, utau include effect of ice internal stress
!YY: this should be confirmed by ZY
          areaW = 0.5_wp * (AREA(i) + AREA(w)) * SEAICEstressFactor
          utau_ice(i)=(ONE-areaW)*utau(i) + areaW*windTauX(i)     &
              + stressDivergenceX(i) * SEAICEstressFactor
          areaS = 0.5_wp * (AREA(i) + AREA(s)) * SEAICEstressFactor
          vtau_ice(i)=(ONE-areaS)*vtau(i) + areaS*windTauY(i)     &
              + stressDivergenceY(i) * SEAICEstressFactor
        enddo

      ELSE
!     else: useHB87StressCoupling=F

!--   Compute ice-affected wind stress (interpolate to U/V-points)
!     by averaging wind stress and ice-ocean stress according to
!     ice cover
!$acc loop 
!YY: pay attention
        do i = 1,nlpb
          vecFldLoc(i, iu) = uFld(i,kSrf)
          vecFldLoc(i, iv) = vFld(i,kSrf)
          vecIceLoc(i, iu) = uIce(i)
          vecIceLoc(i, iv) = vIce(i)
        end do
!$acc loop private(w,j1, j2, j3, j4,k1, k2, k3, k4,p1, p2, p3, p4,  &
!$acc              fuIceLoc,areaW)

        do i = 1,nlpb
          w = tw(i)
          j1 = au1(i, 1)
          j2 = au2(i, 1)
          j3 = au3(i, 1)
          j4 = au4(i, 1)
          k1 = au1(i, 2)
          k2 = au2(i, 2)
          k3 = au3(i, 2)
          k4 = au4(i, 2)
          p1 = au1(i, 3)
          p2 = au2(i, 3)
          p3 = au3(i, 3)
          p4 = au4(i, 3)

          fuIceLoc=HALF*( DWATN(i)+DWATN(w) )*                        &
              COSWAT *                                                &
              ( vecIceLoc(i,iu)-vecFldLoc(i,iu) )                     &
              - SIGN(SINWAT, fCori(i)) * 0.5_wp *                     &
              ( DWATN(i) *                                            &
              0.5_wp*(vecIceLoc(j4,k4)*p4-vecFldLoc(j4,k4)*p4         &
                       +vecIceLoc(j2,k2)*p2-vecFldLoc(j2,k2)*p2)      &
              + DWATN(w) *                                            &
              0.5_wp*(vecIceLoc(j3,k3)*p3-vecFldLoc(j3,k3)*p3         &
                       +vecIceLoc(j1,k1)*p1-vecFldLoc(j1,k1)*p1) )  
          areaW = 0.5_wp * (AREA(i) + AREA(w)) * SEAICEstressFactor
          !YY: here fu/fv are replaced by u/vtau
          utau_ice(i)=(ONE-areaW)*utau(i)+areaW*fuIceLoc
        enddo

!$acc loop private(s,j1, j2, j3, j4,k1, k2, k3, k4,p1, p2, p3, p4,  &
!$acc              fvIceLoc,areaS)
        do i = 1,nlpb
          s = ts(i)
          j1 = av1(i, 1)
          j2 = av2(i, 1)
          j3 = av3(i, 1)
          j4 = av4(i, 1)
          k1 = av1(i, 2)
          k2 = av2(i, 2)
          k3 = av3(i, 2)
          k4 = av4(i, 2)
          p1 = av1(i, 3)
          p2 = av2(i, 3)
          p3 = av3(i, 3)
          p4 = av4(i, 3)
          
          fvIceLoc=HALF*( DWATN(i)+DWATN(s) )*                       &
              COSWAT *                                               &
              ( vecIceLoc(i,iv)-vecFldLoc(i,iv) )                    &
              + SIGN(SINWAT, fCori(i)) * 0.5_wp *                    &
              ( DWATN(i) *                                           &
              0.5_wp*(vecIceLoc(j1,k1)*p1-vecFldLoc(j1,k1)*p1        &
                      +vecIceLoc(j3,k3)*p3-vecFldLoc(j3,k3)*p3)      &
              + DWATN(s) *                                           &
              0.5_wp*(vecIceLoc(j2,k2)*p2-vecFldLoc(j2,k2)*p2        &
                       +vecIceLoc(j4,k4)*p4-vecFldLoc(j4,k4)*p4) )
          areaS = 0.5_wp * (AREA(i) + AREA(s)) * SEAICEstressFactor
          vtau_ice(i)=(ONE-areaS)*vtau(i)+areaS*fvIceLoc
        enddo

      ENDIF
!$acc end kernels

!YY,ZY: exchange necessary ? only used in mod_csp_dyn_force later
      !call mpi_data_exchange(mpi_sendbuf_1d,utau)
      !call mpi_data_exchange(mpi_sendbuf_1d,vtau)

!$acc end data

   end subroutine seaice_ocean_stress


#ifdef SEAICE_ALLOW_BOTTOMDRAG
   subroutine seaice_bottomdrag_coeffs(uIceLoc, vIceLoc, HEFFMLoc,    &
              HEFFLoc, AREALoc, CbotLoc)
!     *==========================================================*
!     | o Compute the non-linear drag coefficients for ice-bottom
!     |   drag, as a parameterization for grounding fastice
!     |   following
!     |   Lemieux et al. (2015), doi:10.1002/2014JC010678
!     *==========================================================*
!     | written by Martin Losch, Apr 2016
!     *==========================================================*

      implicit none

!     u/vIceLoc :: local copies of the current ice velocity
!     HEFFMLoc  :: local copy of land-sea masks
!     CbotLoc   :: drag coefficients
      real(wp) ::  uIceLoc   (nlpb)
      real(wp) ::  vIceLoc   (nlpb)
      real(wp) ::  HEFFMLoc  (nlpb)
      real(wp) ::  HEFFLoc   (nlpb)
      real(wp) ::  AREALoc   (nlpb)
      real(wp) ::  CbotLoc   (nlpb)

      INTEGER i
      INTEGER kSrf
      real(wp) ::  tmpFld(nlpb)
      real(wp) ::  tmp, htmp, hCrit, recip_k1, u0sq, fac, rFac
      integer :: ibot
      real(wp) :: vecIceLoc(nlpb,2)
      integer(i2) :: maskWS(nlpb,2)
      integer :: uel,uei,uep,vnl,vni,vnp

!$acc data present(uIceLoc, vIceLoc, HEFFMLoc,HEFFLoc, AREALoc, CbotLoc,   &
!$acc              maskW,maskS,ue,vn),   &
!$acc      create(tmpFld,vecIceLoc,maskWS)

      IF (SEAICEbasalDragK2.GT.0.0_wp) THEN
!     avoid this computation for a non-zero coefficient
      kSrf = nk
!     some abbreviations
      u0sq     = SEAICEbasalDragU0*SEAICEbasalDragU0
      recip_k1 = 0.0_wp
      IF ( SEAICEbasalDragK1 .GT. 0.0_wp ) recip_k1 = 1.0_wp/SEAICEbasalDragK1
!     fac scales the soft maximum for more accuracy
      fac = 10.0_wp
      rFac = 1.0_wp/fac

!$acc kernels copyin(kSrf,u0sq,recip_k1,fac,rFac)
!$acc loop
      do i = 1,nlpb
        CbotLoc(i) = 0.0_wp 
        tmpFld (i) = 0.0_wp
      enddo
!$acc loop
      do i = 1,nlpb
        vecIceLoc(i,iu) = uIceLoc(i)
        vecIceLoc(i,iv) = vIceLoc(i)
        maskWS(i,iu) = maskW(i,kSrf)
        maskWS(i,iv) = maskS(i,kSrf)
      enddo
!$acc loop private(tmp,uel,uei,uep,vnl,vni,vnp)
      do i = 1,nlpb
        IF ( AREALoc(i) .GT. 0.01_wp ) THEN
          uel = ue(i, 1)
          uei = ue(i, 2)
          uep = ue(i, 3)
          vnl = vn(i, 1)
          vni = vn(i, 2)
          vnp = vn(i, 3)

    !      tmp = 0.25_wp*(
    ! &         ( vecIceLoc(i,iu)*maskInW( i ,j,bi,bj)
    ! &         + uIceLoc(i+1,j,bi,bj)*maskInW(i+1,j,bi,bj) )**2
    ! &       + ( vIceLoc(i,j  ,bi,bj)*maskInS(i,j  ,bi,bj)
    ! &         + vIceLoc(i,j+1,bi,bj)*maskInS(i,j+1,bi,bj) )**2 )
          tmp = 0.25_wp*(    &
              ( vecIceLoc(i,iu)*maskWS( i,iu) + vecIceLoc(uel,uei)*uep*maskWS(uel,uei) )**2  &
            + ( vecIceLoc(i,iv)*maskWS(i,iv) + vecIceLoc(vnl,vni)*vnp*maskWS(vnl,vni) )**2 )
           tmpFld(i) = SEAICEbasalDragK2 / SQRT(tmp + u0sq)

        ENDIF
      enddo

!$acc loop private(ibot,htmp,hCrit)
      do i = 1,nlpb
        ibot = min(nk, kSurfC(i))    
        IF ( AREALoc(i) .GT. 0.01_wp ) THEN
          htmp = HEFFLoc(i)
          hCrit   = ABS(rF(ibot)/gravityRau0)*AREALoc(i)*recip_k1
          CbotLoc(i) = CbotLoc(i) + tmpFld(i)                    &
              * LOG(EXP( fac*(htmp - hCrit) ) + 1.0_wp)*rFac  &
              * EXP( - SEAICE_cBasalStar                         &
                     *(SEAICE_area_max - AREALoc(i)) )           &
              * HEFFMLoc(i)
        ENDIF
      enddo
!$acc end kernels
      ENDIF    !endif SEAICEbasalDragK2.GT.0.
!$acc end data

   end subroutine seaice_bottomdrag_coeffs
#endif /* SEAICE_ALLOW_BOTTOMDRAG */

end module mitice_dynamics
