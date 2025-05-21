module mitice_itd
use mitice_parameters
use mitice_vars
use mod_csp_basic
#ifdef SeaiceDebug
!use mod_mpi_interfaces
#endif
        
   implicit none

!YY: these vars are memo-consuming, should be dynamically (de)allocated?
#ifdef SEAICE_ITD
!     variables related to ridging schemes
!     ridgingModeNorm :: norm to ensure convervation (N in Lipscomb et al 2007)
!     partFunc   :: participation function (a_n in Lipscomb et al 2007)
!     ridgeRatio :: mean ridge thickness/ thickness of ridging ice
!     hrMin      :: min ridge thickness
!     hrMax      :: max ridge thickness   (SEAICEredistFunc = 0)
!     hrExp      :: ridge e-folding scale (SEAICEredistFunc = 1)
!     hActual    :: HEFFITD/AREAITD
      real(wp), dimension(:), allocatable  ::  ridgingModeNorm 
      real(wp), dimension(:,:), allocatable  ::  partFunc 
      real(wp), dimension(:,:), allocatable  ::  hrMin
      real(wp), dimension(:,:), allocatable  ::  hrMax
      real(wp), dimension(:,:), allocatable  ::  hrExp
      real(wp), dimension(:,:), allocatable  ::  ridgeRatio
      real(wp), dimension(:,:), allocatable  ::  hActual    !YY: note hActual is a local var in seaice_itd_remap
      !$acc declare create(ridgingModeNorm,partFunc,hrMin,hrMax,hrExp,ridgeRatio,hActual)
      
#endif /* SEAICE_ITD */


contains
#ifdef SEAICE_ITD
   subroutine seaice_itd_allocate
     implicit none
     allocate(ridgingModeNorm(nlpb))
     allocate(partFunc(nlpb,0:nITD))
     allocate(hrMin(nlpb,1:nITD))
     allocate(ridgeRatio(nlpb,1:nITD), hActual(nlpb,1:nITD))
     IF ( SEAICEredistFunc .EQ. 0 ) THEN
       allocate(hrMax(nlpb, 1:nITD))
     ELSEIF( SEAICEredistFunc .EQ. 1 ) THEN
       allocate(hrExp(nlpb, 1:nITD))
     ENDIF


   end subroutine seaice_itd_allocate

   subroutine seaice_itd_release
     implicit none
     deallocate(ridgingModeNorm)
     deallocate(partFunc)
     deallocate(hrMin, hrMax,hrExp)
     deallocate(ridgeRatio, hActual)
     IF ( SEAICEredistFunc .EQ. 0 ) THEN
       deallocate(hrMax)
     ELSEIF( SEAICEredistFunc .EQ. 1 ) THEN
       deallocate(hrExp)
     ENDIF
   end subroutine seaice_itd_release
#endif
         

   subroutine seaice_reg_ridge
!     *=================================================================*
!     | o this routine has two purposes:
!     |   (1) clean up after advection (undershoots etc.);
!     |       after advection, the sea ice variables may have unphysical
!     |       values, e.g. < 0 or very thin ice, that are regularized
!     |       here.
!     |   (2) driver for ice ridging;
!     |       concentration as a special case may be > 1 in convergent
!     |       motion and a ridging algorithm redistributes the ice to
!     |       limit the concentration to 1.
!     *=================================================================*

      implicit none
      integer i
      integer it
      real(wp)  recip_deltaTtherm       !reciprocal of time step
      real(wp)  tmpscal1, tmpscal2

!$acc data present(d_HEFFbyNEG,d_HSNWbyNEG,HEFF,AREA,HSNOW,    &
#ifdef SEAICE_ITD
!$acc              HEFFITD,AREAITD, HSNOWITD,     &
#endif
#ifdef SEAICE_VARIABLE_SALINITY
!$acc              saltFluxAdjust,HSALT,             &
#endif
!$acc              TICES)

      recip_deltaTtherm = ONE / SEAICE_deltaTtherm

!$acc kernels loop
      do i = 1,nlpb
        d_HEFFbyNEG(i)    = 0.0_wp
        d_HSNWbyNEG(i)    = 0.0_wp
#ifdef SEAICE_VARIABLE_SALINITY
        saltFluxAdjust(i) = 0.0_wp
#endif /* SEAICE_VARIABLE_SALINITY */
      enddo
!$acc end kernels

! =====================================================================
! ========== PART 1: treat pathological cases (post advdiff) ==========
! =====================================================================

!YY: ifdef EXF_SEAICE_FRACTION, then areamask (seaice fraction) will be
!YY: read in EXF_getffield as an initial observational field, and
!YY: relaxation of AREA is done here towards observations
!YY: TO DO LATER NOW DELETE. REFER TO ORIGINAL MIT CODE

!YY: TREAT NEGATIVE VALUE BY DIVERGENCE
!--   (1) treat the case of negative values:

#ifdef SEAICE_ITD
!$acc kernels
!$acc loop private(tmpscal1,tmpscal2)     !YY: is tmpscal private here, test?
      do i = 1,nlpb
!$acc loop seq
        do it=1,SEAICE_multDim
          tmpscal1=0.0_wp
          tmpscal2=0.0_wp
          tmpscal1=MAX(-HEFFITD(i,it), 0.0_wp)
          !YY: neutralize to zero
          HEFFITD(i,it)=HEFFITD(i,it)+tmpscal1
          !YY: deltaHEFF cumulation by categories
          d_HEFFbyNEG(i)=d_HEFFbyNEG(i)+tmpscal1
          tmpscal2=MAX(-HSNOWITD(i,it), 0.0_wp)
          HSNOWITD(i,it)=HSNOWITD(i,it)+tmpscal2
          d_HSNWbyNEG(i)=d_HSNWbyNEG(i)+tmpscal2
          AREAITD(i,it)=MAX(AREAITD(i,it),0.0_wp)
!     AREA, HEFF, and HSNOW will be updated at end of PART 1
!     by calling SEAICE_ITD_SUM
        enddo
      enddo
!$acc end kernels

!     update mean thicknesses HEFF and HSNOW and total ice
!     concentration AREA to match single category values
      call seaice_itd_sum

#else /* ndef SEAICE_ITD */
!$acc kernels loop
      do i = 1,nlpb
        d_HEFFbyNEG(i)=MAX(-HEFF(i),0.0_wp)
        HEFF(i)=HEFF(i)+d_HEFFbyNEG(i)
        d_HSNWbyNEG(i)=MAX(-HSNOW(i),0.0_wp)
        HSNOW(i)=HSNOW(i)+d_HSNWbyNEG(i)
        AREA(i)=MAX(AREA(i),0.0_wp)
      enddo
!$acc end kernels
#endif /* SEAICE_ITD */

!--   (2) treat the case of very thin ice:

#ifdef SEAICE_ITD
!     Here we risk that even though HEFF may be larger than siEps (=1e-5)
!     HEFFITD can have classes with very small (< siEps) non-zero ice volume.
!     We avoid applying the correction to each class because that leads to
!     funny structures in the net heat and freshwater flux into the ocean.
!     Let us keep our fingers crossed, that the model will be benign!
!$acc kernels loop collapse(2)
      DO IT=1,SEAICE_multDim
        do i = 1,nlpb
          IF (HEFF(i).LE.siEps) THEN
            HEFFITD(i,IT) = 0.0_wp
            HSNOWITD(i,IT) = 0.0_wp
          ENDIF
        enddo
      ENDDO
!$acc end kernels
#endif /* SEAICE_ITD */
!$acc kernels 
!$acc loop independent private(tmpscal1,tmpscal2)
      do i = 1,nlpb
        tmpscal1=0.0_wp
        tmpscal2=0.0_wp
        IF (HEFF(I).LE.siEps) THEN
         tmpscal1=-HEFF(i)
         tmpscal2=-HSNOW(i)
!$acc loop seq
         DO IT=1,SEAICE_multDim
          !YY: surface ice temp set to zero
          TICES(I,IT)=celsius2K
         ENDDO
        ENDIF
        HEFF(I)=HEFF(I)+tmpscal1
        HSNOW(I)=HSNOW(I)+tmpscal2
        d_HEFFbyNEG(I)=d_HEFFbyNEG(I)+tmpscal1
        d_HSNWbyNEG(I)=d_HSNWbyNEG(I)+tmpscal2
      enddo
!$acc end kernels

!--   (3) treat the case of area but no ice/snow:

#ifdef SEAICE_ITD
!$acc kernels loop collapse(2)
      DO IT=1,SEAICE_multDim
        do i = 1,nlpb
          IF ( (HEFFITD(I,IT) .EQ.0.0_wp).AND.(HSNOWITD(I,IT).EQ.0.0_wp))  &
            AREAITD(I,IT)=0.0_wp
        enddo
      ENDDO
!$acc end kernels
#else /* ndef SEAICE_ITD */
!$acc kernels loop
      do i = 1,nlpb
         IF ((HEFF(i).EQ.0.0_wp).AND.(HSNOW(i).EQ.0.0_wp)) AREA(I)=0.0_wp
      enddo
!$acc end kernels
#endif /* SEAICE_ITD */

!--   (4) treat the case of very small area:

#ifdef SEAICE_ITD
!$acc kernels loop collapse(2)
      DO IT=1,SEAICE_multDim
        do i = 1,nlpb
          IF ((HEFFITD(I,IT).GT.0.0_wp).OR.(HSNOWITD(I,IT).GT.0.0_wp)) THEN
!     SEAICE_area_floor*SEAICE_multDim cannot be allowed to exceed 1
!     hence use SEAICE_area_floor devided by SEAICE_multDim
!     (or install a warning in e.g. seaice_readparms.F)
            AREAITD(I,IT)= MAX(AREAITD(I,IT),SEAICE_area_floor*( 1.0_wp / float(SEAICE_multDim) ))
          ENDIF
        enddo
      enddo
!$acc end kernels
#else /* ndef SEAICE_ITD */
!$acc kernels loop
      do i = 1,nlpb
        IF ((HEFF(i).GT.0.0_wp).OR.(HSNOW(i).GT.0.0_wp)) THEN
          AREA(I)=MAX(AREA(I),SEAICE_area_floor)
        ENDIF
      enddo
!$acc end kernels
#endif /* SEAICE_ITD */

!     (5) treat sea ice salinity pathological cases
#ifdef SEAICE_VARIABLE_SALINITY
!$acc kernels copyin(recip_deltaTtherm)
!$acc loop
      do i = 1,nlpb
        IF ( (HSALT(I) .LT. 0.0_wp).OR.(HEFF(I) .EQ. 0.0_wp)  ) THEN
           saltFluxAdjust(I) = - HEFFM(I) * HSALT(I) * recip_deltaTtherm
           HSALT(I) = 0.0_wp
        ENDIF
      enddo
!$acc end kernels
#endif /* SEAICE_VARIABLE_SALINITY */

! =====================================================================
! ========== PART 2: ridging algorithm  ===============================
! =====================================================================

!!YY: TREAT EXCESSIVE >1 VALUE BY CONVERGENCE
!     treat case of excessive ice cover, e.g., due to ridging:

#ifdef SEAICE_ITD

!     ridge ice according to Lipscomb et al. (2007), Bitz et al. (2001)
!     Thorndyke et al. (1975), Hibler (1980)
      call seaice_do_ridging
!     check that all ice thickness categories meet their limits
!     (includes Hibler-type ridging)
      call seaice_itd_redist
!     update mean thicknesses HEFF and HSNOW and total ice
!     concentration AREA to match single category values
      call seaice_itd_sum

#else /* ifndef SEAICE_ITD */

!$acc kernels loop
      do i = 1,nlpb
!     this is the simple Hibler (1979)-type ridging (capping of
!     concentrations > 1) for the non-ITD sea ice model
        AREA(i)=MIN(AREA(i),SEAICE_area_max)
      enddo
!$acc end kernels

#endif /* SEAICE_ITD */

!$acc end data

   end subroutine seaice_reg_ridge



#ifdef SEAICE_ITD
   subroutine seaice_itd_redist
!     *===========================================================*
!     | o checks if absolute ice thickness in any category
!     |   exceeds its category limits
!     | o redistributes sea ice area and volume
!     |   and associated ice properties in thickness space
!     |
!     | Torge Martin, Feb. 2012, torge@mit.edu
!     *===========================================================*

      implicit none

      INTEGER i,  k
!     openwater   :: open water area fraction
      real(wp), dimension(:), allocatable :: openWater


      allocate( openWater(nlpb) )

!$acc data create(openWater) present(AREAITD,HEFFITD,HSNOWITD,Hlimit)

!       calculate area of open water
!$acc kernels
!$acc loop
      do i = 1,nlpb
        openWater(i) = ONE
      enddo
!$acc loop
      do i = 1,nlpb
!$acc loop seq
       do k=1,nITD
         openWater(i) = openWater(i) - AREAITD(i,k)
       enddo
      enddo

!     ----------------------------------------------------
!     | redistribute/"advect" sea ice in thickness space |
!     | as described in Bitz et al. (2001)               |
!     ----------------------------------------------------


!--   Hibler-type "ridging", i.e. cut back excessive ice area fraction ---
!     in case ice concentration exceeds 100% assume that
!     convergence of floe field has eliminated all open water
!     and eventual rafting occured in thinnest category:
!$acc loop
      do i = 1,nlpb
        IF (openWater(i) .lt. 0.0_wp)    &
        !YY: excessive ice area redist to the 1st category
           AREAITD(i,1) = openWater(i) + AREAITD(i,1)
      enddo
!$acc end kernels

!     the following steps only make sense if there are actually
!     multi-categories
      IF (nITD .gt. 1) THEN

!--   check if more thicker ice needs to be rafted to accomodate area excess:
!$acc kernels
!$acc loop
       do i = 1,nlpb
!$acc loop seq
        do k=1,nITD-1
          IF (AREAITD(i,k) .lt. 0.0) THEN
!--   pass concentration deficit up to next thicker category
!--   since all quantities are extensive, we add instead of average
           AREAITD (i,k+1) = AREAITD (i,k+1) + AREAITD (i,k)
           AREAITD (i,k  ) = ZERO
           HEFFITD (i,k+1) = HEFFITD (i,k+1) + HEFFITD (i,k)
           HEFFITD (i,k  ) = ZERO
           HSNOWITD(i,k+1) = HSNOWITD(i,k+1) + HSNOWITD(i,k)
           HSNOWITD(i,k  ) = ZERO
          ENDIF
        enddo
       enddo

!     --- ice thickness redistribution ---
!         now check that ice thickness stays within category limits
!$acc loop
       do i = 1,nlpb
!$acc loop seq
        do k=1,nITD-1
          IF (HEFFITD(i,k) .gt. Hlimit(k)*AREAITD(i,k)) THEN
!--   the upper thickness limit is exceeded: move ice up to next
!     thicker category
           AREAITD (i,k+1) = AREAITD (i,k+1) + AREAITD (i,k)
           AREAITD (i,k  ) = ZERO
           HEFFITD (i,k+1) = HEFFITD (i,k+1) + HEFFITD (i,k)
           HEFFITD (i,k  ) = ZERO
           HSNOWITD(i,k+1) = HSNOWITD(i,k+1) + HSNOWITD(i,k)
           HSNOWITD(i,k  ) = ZERO
          ENDIF
        enddo
       enddo

!$acc loop
       do i = 1,nlpb
!$acc loop seq
        do k=nITD,2,-1
          IF (HEFFITD(i,k) .lt. Hlimit(k-1)*AREAITD(i,k)) THEN
!--   the lower thickness limit is exceeded: move ice down to next thinner
!     category
           AREAITD (i,k-1) = AREAITD (i,k-1) + AREAITD (i,k  )
           AREAITD (i,k  ) = ZERO
           HEFFITD (i,k-1) = HEFFITD (i,k-1) + HEFFITD (i,k)
           HEFFITD (i,k  ) = ZERO
           HSNOWITD(i,k-1) = HSNOWITD(i,k-1) + HSNOWITD(i,k)
           HSNOWITD(i,k  ) = ZERO
          ENDIF
        enddo 
       enddo
!$acc end kernels

!     end nITD>1 constraint
      ENDIF
!$acc end data

      deallocate( openWater )

   end subroutine seaice_itd_redist
#endif /* SEAICE_ITD */

#ifdef SEAICE_ITD
   subroutine seaice_itd_sum
!     *===========================================================*
!     | o sum ice area and volume over all ITD categories
!     |   and write into AREA and HEFF
!     *===========================================================*

      implicit none

      INTEGER i, k

!$acc kernels present(AREA,HEFF,HSNOW,AREAITD,HEFFITD,HSNOWITD)
!$acc loop
      do i = 1,nlpb
         AREA (i) = AREAITD (i,1)
         HEFF (i) = HEFFITD (i,1)
         HSNOW(i) = HSNOWITD(i,1)
      enddo

!$acc loop
      do i = 1,nlpb
!$acc loop seq
        do k=2,nITD
          AREA (i) = AREA (i) + AREAITD (i,k)
          HEFF (i) = HEFF (i) + HEFFITD (i,k)
          HSNOW(i) = HSNOW(i) + HSNOWITD(i,k)
        enddo
      ENDDO
!$acc end kernels

   end subroutine seaice_itd_sum
#endif /* SEAICE_ITD */

   subroutine seaice_do_ridging

!     *===========================================================*
!     | o compute mechanical redistribution of thin (undeformed) into
!     |   thick (deformed, i.e. ridged) ice categories
!     |   according to Thorndyke et al (1975) and Hibler (1980)
!     |   or Bitz et al (2001) and Lipscomb et al (2007)
!     |
!     | Martin Losch, Apr. 2014, Martin.Losch@awi.de
!     *===========================================================*

      implicit none
      integer i
#ifdef SEAICE_ITD
      integer k, l, n

      real(wp) ::  openWater       (nlpb) !open water area fraction
      real(wp) ::  netArea         (nlpb)
!     variables related to ridging schemes
      real(wp) ::  openingRate     (nlpb)
      real(wp) ::  closingRate     (nlpb)
      real(wp) ::  grossClosing    (nlpb)
!     amount of ice that participates in ridging (according to partFunc)
      real(wp) ::  ridgingArea     (nlpb)
      real(wp) ::  ridgingHeff     (nlpb)
      real(wp) ::  ridgingHsnw     (nlpb)
!     fractions of deformed/ridged ice
      real(wp) ::  areaFraction    (nlpb)
      real(wp) ::  volFraction     (nlpb)
!     absolute area/concentration of deformed/ridged ice
      real(wp) ::  ridgedArea      (nlpb)
      LOGICAL  ::  doRidging       (nlpb)
      LOGICAL doRidgeAgain, areaTooLarge

      real(wp) ::  convergence, divergence, shear, divAdv
      real(wp) ::  tmp, tmpFac, hL, hR, expL, expR
      real(wp) ::  areaPR(nlpb,1:nITD)
      real(wp) ::  heffPR(nlpb,1:nITD)
      real(wp) ::  hsnwPR(nlpb,1:nITD)

#endif /* SEAICE_ITD */


#ifndef SEAICE_ITD
!     Hiblers "ridging function" for single category ice
!$acc kernels loop
      do i = 1,nlpb
        AREA(i) = MIN(AREA(i),SEAICE_area_max)
      enddo
!$acc end kernels
#else

!$acc data present(AREAITD,e11,e22,deltaC,opnWtrFrac,fw2ObyRidge,    &
!$acc              ridgingModeNorm,partFunc,ridgeRatio,hrMin,hLimit),    &
!$acc      create(openWater,ridgingArea,ridgingHeff,ridgingHsnw,  &
!$acc            areaFraction,volFraction,ridgedArea, &
!$acc            openingRate,closingRate,grossClosing,netArea,    &
!$acc            areaPR,heffPR,hsnwPR,doRidging,ridgedArea)

!     calculate area of open water
!YY: rafting and ridging only occur when openwater <0 
!$acc kernels
!$acc loop
      do i = 1,nlpb
        openWater(i) = ONE
      enddo
!$acc loop
      do i = 1,nlpb
!$acc loop seq
        do k=1,nITD
          openWater(i) = openWater(i) - AREAITD(i,k)
        enddo
      enddo
!$acc end kernels
      IF ( SEAICEsimpleRidging ) THEN
!--   Hibler-type "ridging", i.e. cut back excessive ice area fraction ---
!     in case ice concentration exceeds 100% assume that
!     convergence of floe field has eliminated all open water
!     and eventual rafting occured in thinnest category:
!$acc kernels loop
        do i = 1,nlpb
          IF (openWater(i) .lt. 0.0_wp) AREAITD(i,1) = openWater(i)+AREAITD(i,1)
        enddo
!$acc end kernels

      ELSE
!     initialisation
!$acc kernels
!$acc loop
      do i = 1,nlpb
        ridgingArea(i)    = 0.0_wp
        ridgingHeff(i)    = 0.0_wp
        ridgingHsnw(i)    = 0.0_wp
        areaFraction(i)   = 0.0_wp
        volFraction(i)    = 0.0_wp
        fw2ObyRidge(i)    = 0.0_wp
        ridgedArea(i)     = 0.0_wp
      enddo
!$acc end kernels

!YY: update hActual,ridgingModeNorm,partFunc,hrMin,hrMax,hrExp,ridgeRatio
      call seaice_prepare_ridging

!     Compute the first strain rate invariant epsilonI (divergence)
!     energy dissipation by convergence = -min (divergence, 0)
!     energy dissipation by shearing    = (1/2) * (Delta - abs(divergence))
!$acc kernels
!$acc loop private(divergence,shear,convergence)
      do i = 1,nlpb
        divergence  = e11(i) + e22(i)
        shear       = 0.5_wp * ( deltaC(i) - ABS(divergence) )
        convergence = - MIN(divergence, 0.0_wp)
        closingRate(i) = SEAICEshearParm*shear + convergence
      enddo
!     we need a new estimate of the total AREA (including the open water
!     fraction, but for computational reason it is not included here)
!$acc loop
      do i = 1,nlpb
        netArea(i) = 0.0_wp
      enddo
!$acc loop
      do i = 1,nlpb
!$acc loop seq
        do k=1,nITD
          netArea(i) = netArea(i) + AREAITD(i,k)
        enddo
      enddo
!$acc loop private(divAdv)
      do i = 1,nlpb
!     divergence rate due to advection; this term need not be zero due
!     to numerical effects
!     (this is copied from CICE but I am not sure about that)
        divAdv = (1.0_wp-netArea(i)-opnWtrFrac(i)) * (1.0_wp/SEAICE_deltaTtherm)
        IF (divAdv .LT. 0.0_wp) closingRate(i) = MAX(closingRate(i), -divAdv)
!     finally compute a non-negative opening rate that will lead to
!     a net area of 1
        openingRate(i) = closingRate(i) + divAdv
      enddo
!$acc end kernels

!     start of the ridging loop

      doRidgeAgain = .TRUE.
      n = 1
      DO WHILE (doRidgeAgain)
!     save pre-ridging ice concentration and ridged ice volume
!$acc kernels
!$acc loop collapse(2)
      DO k=1,nITD
        do i = 1,nlpb
          areaPR(i,k) = AREAITD (i,k)
          heffPR(i,k) = HEFFITD (i,k)
          hsnwPR(i,k) = HSNOWITD(i,k)
        enddo
      ENDDO

!$acc loop private(tmp,tmpFac)
      do i = 1,nlpb
!     Based on the ITD of ridging and ridged ice, convert the net
!     closing rate to a gross closing rate times deltaT.
!     NOTE: 0 < ridgingModeNorm <= 1
        grossClosing(i) = closingRate(i)*SEAICE_deltaTtherm/ridgingModeNorm(i)
!     reduce rates in case more than 100% of open water would be removed
        IF ( partFunc(i,0) .GT. 0.0_wp ) THEN
          tmp = partFunc(i,0)*grossClosing(i)
          IF ( tmp .GT. opnWtrFrac(i) ) THEN
           tmpFac = opnWtrFrac(i)/tmp
           grossClosing(i) = grossClosing(i) * tmpFac
           openingRate(i)  =  openingRate(i) * tmpFac
          ENDIF
        ENDIF
      enddo
!$acc loop private(tmp,tmpFac)
      do i = 1,nlpb
!$acc loop seq
        do k=1,nITD
!     reduce rates in case more than 100% of any ice categroy would be removed
          IF ( areaPR(i,k) .GT. SEAICE_area_reg             &
              .AND. partFunc(i,k) .GT. 0.0_wp ) THEN
            tmp = partFunc(i,k)*grossClosing(i)
            IF ( tmp .GT. AREAITD(i,k) ) THEN
             tmpFac = AREAITD(i,k)/tmp
             grossClosing(i) = grossClosing(i) * tmpFac
             openingRate(i)  =  openingRate(i) * tmpFac
            ENDIF
          ENDIF
        enddo 
      enddo

!     start redistribution
!$acc loop 
      do i = 1,nlpb
!     open water first
        opnWtrFrac(i) = opnWtrFrac(i)         &
            - partFunc(i,0)*grossClosing(i) + openingRate(i)*SEAICE_deltaTtherm
        opnWtrFrac(i) = MAX( 0.0_wp, opnWtrFrac(i) )
      enddo
!$acc end kernels

      DO k=1,nITD
!$acc kernels copyin(k)
!$acc loop independent
        do i = 1,nlpb
          doRidging(i) = areaPR(i,k) .GT. SEAICE_area_reg   &
              .AND. partFunc(i,k) .GT. 0.0_wp               &
              .AND. grossClosing(i) .GT. 0.0_wp             &
              .AND. HEFFM(i) .GT. 0.0_wp
!     this would be safety catch only
!     &        .AND. netArea(i,j) .GT. 1.0_wp
         IF ( doRidging(i) ) THEN
           ridgingArea(i) = partFunc(i,k)*grossClosing(i)
           IF ( ridgingArea(i) .GT. areaPR(i,k) ) THEN
             ridgingArea(i) = areaPR(i,k)
           ENDIF
           areaFraction(i) = ridgingArea(i)/areaPR(i,k)
           !YY: ridgeratio is a >1 value, ice which participates will shrink its area by ratio
           ridgedArea(i)   = ridgingArea(i)/ridgeRatio(i,k)
!     compute ice volume (HEFF) and snow volume to be removed from this
!     ridging category
           ridgingHEFF(i) = heffPR(i,k) * areaFraction(i)
           ridgingHsnw(i) = hsnwPR(i,k) * areaFraction(i)
!     part of the snow mass is pushed into the ocean during ridging;
!     this freshwater flux will be added to the net feshwater flux into
!     the ocean in seaice_growth
           fw2ObyRidge(i) = fw2ObyRidge(i)      &
               + SEAICE_rhoSnow*ridgingHsnw(i)  &
               *(1.0_wp - SEAICEsnowFracRidge)
!     reduce the snow volume that is left for redistribution
           ridgingHsnw(i) = ridgingHsnw(i) * SEAICEsnowFracRidge
!     remove ice concentration, volume (HEFF), and snow volume from
!     this ridging category
           AREAITD(i,k) = AREAITD(i,k) - ridgingArea(i)
           HEFFITD(i,k) = HEFFITD(i,k) - ridgingHeff(i)
           HSNOWITD(i,k)=HSNOWITD(i,k) - ridgingHsnw(i)
         ENDIF
       enddo

!     inner loop over categories: distribute what has been removed from the
!     kth category to all categories according to area/volFraction
       IF ( SEAICEredistFunc .EQ. 0 ) THEN
!     Assume ridged ice is uniformly distributed between hrmin and hrmax
!     (Hibler, 1980), see also s/r seaice_prepare_ridging.
!$acc loop private(hL,hR) independent    !k has been copyin
        do i = 1,nlpb       
!$acc loop seq
         DO l=1,nITD
!     initialising these is essential, because here the ridging-mask doRidging
!     to area/volFraction, and applied via these fields
          areaFraction(i) = 0.0_wp
          volFraction (i) = 0.0_wp

          IF ( doRidging(i) ) THEN
           IF ( hrMin(i,k) .GE. hLimit(l) .OR. hrMax(i,k) .LE. hLimit(l-1) ) THEN
             areaFraction(i) = 0.0_wp
             volFraction (i) = 0.0_wp
           ELSE
             hL = MAX(hrMin(i,k), hLimit(l-1))
             hR = MIN(hrMax(i,k), hLimit(l))
             areaFraction(i) = ( hR - hL ) / ( hrMax(i,k) - hrMin(i,k) )
             volFraction (i) = areaFraction(i)*( hR + hL )      &
               / ( hrMax(i,k) + hrMin(i,k) )
           ENDIF
          ENDIF
!         after computing the fraction ridged ice for this category, apply it
         !YY: here ridgedarea is for outer category k, while areafraction is updated in inner loop
          AREAITD(i,l) = AREAITD(i,l) +areaFraction(i)*ridgedArea(i)
          HEFFITD(i,l) = HEFFITD(i,l) +volFraction(i)*ridgingHeff(i)
          HSNOWITD(i,l) = HSNOWITD(i,l)+volFraction(i)*ridgingHsnw(i)*SEAICEsnowFracRidge
         ENDDO   !enddo for inner l-loop
        enddo    !enddo for nlpb

       ELSEIF ( SEAICEredistFunc .EQ. 1 ) THEN
!     Follow Lipscomb et al. (2007) and model ridge ITD as an exponentially
!     decaying function, see also s/r seaice_prepare_ridging.
!$acc loop private(hL,hR,expL,expR) independent    !k has been copyin
        do i = 1,nlpb       
!$acc loop seq
         DO l=1,nITD
!     initialising these is essential, because here the ridging-mask doRidging
!     to area/volFraction, and applied via these fields
          areaFraction(i) = 0.0_wp
          volFraction (i) = 0.0_wp

          IF ( l.LT.nITD ) THEN
            IF ( doRidging(i) .AND. hrMin(i,k) .LT. hLimit(l)    &
                        .AND. hrExp(i,k) .NE. 0.0_wp ) THEN
             hL   = MAX( hrMin(i,k), hLimit(l-1) )
             hR   = hLimit(l)
             expL = EXP(-( hL - hrMin(i,k) )/hrExp(i,k) )
             expR = EXP(-( hR - hrMin(i,k) )/hrExp(i,k) )
             areaFraction(i) = expL - expR
             volFraction (i) =  &
                 ( ( hL + hrExp(i,k) ) * expL     &
                 - ( hR + hrExp(i,k) ) * expR )   &
                 / ( hrMin(i,k) + hrExp(i,k) )
            ENDIF
          ELSE
            IF ( doRidging(i) .AND. hrExp(i,k) .NE. 0.0_wp ) THEN
              hL   = MAX( hrMin(i,k), hLimit(l-1) )
              expL = EXP(-( hL - hrMin(i,k) )/hrExp(i,k) )
              areaFraction(i) = expL
              volFraction (i) = ( hL + hrExp(i,k) ) * expL     &
                  / ( hrMin(i,k) + hrExp(i,k) )
            ENDIF
          ENDIF
!         after computing the fraction ridged ice for this category, apply it
          !YY: here ridgedarea is for outer category k, while areafraction is updated in inner loop
          AREAITD(i,l) = AREAITD(i,l) +areaFraction(i)*ridgedArea(i)
          HEFFITD(i,l) = HEFFITD(i,l) +volFraction(i)*ridgingHeff(i)
          HSNOWITD(i,l) = HSNOWITD(i,l)+volFraction(i)*ridgingHsnw(i)*SEAICEsnowFracRidge
         ENDDO   !enddo for inner l-loop
        enddo    !enddo for nlpb
       ENDIF   ! endif for SEAICEredistFunc


!!     inner loop over categories: distribute what has been removed from the
!!     kth category to all categories according to area/volFraction
!!$acc loop seq
!       DO l=1,nITD
!!     initialising these is essential, because here the ridging-mask doRidging
!!     to area/volFraction, and applied via these fields
!         do i = 1,nlpb
!           areaFraction(i) = 0.0_wp
!           volFraction (i) = 0.0_wp
!         enddo
!         IF ( SEAICEredistFunc .EQ. 0 ) THEN
!!     Assume ridged ice is uniformly distributed between hrmin and hrmax
!!     (Hibler, 1980), see also s/r seaice_prepare_ridging.
!           do i = 1,nlpb
!             IF ( doRidging(i) ) THEN
!              IF ( hrMin(i,k) .GE. hLimit(l) .OR. hrMax(i,k) .LE. hLimit(l-1) ) THEN
!                areaFraction(i) = 0.0_wp
!                volFraction (i) = 0.0_wp
!              ELSE
!                hL = MAX(hrMin(i,k), hLimit(l-1))
!                hR = MIN(hrMax(i,k), hLimit(l))
!                areaFraction(i) = ( hR - hL ) / ( hrMax(i,k) - hrMin(i,k) )
!                volFraction (i) = areaFraction(i)*( hR + hL )      &
!                  / ( hrMax(i,k) + hrMin(i,k) )
!              ENDIF
!             ENDIF
!           enddo
!        ELSEIF ( SEAICEredistFunc .EQ. 1 ) THEN
!!     Follow Lipscomb et al. (2007) and model ridge ITD as an exponentially
!!     decaying function, see also s/r seaice_prepare_ridging.
!          IF ( l.LT.nITD ) THEN
!            do i = 1,nlpb
!              IF ( doRidging(i) .AND. hrMin(i,k) .LT. hLimit(l)    &
!                          .AND. hrExp(i,k) .NE. 0.0_wp ) THEN
!               hL   = MAX( hrMin(i,k), hLimit(l-1) )
!               hR   = hLimit(l)
!               expL = EXP(-( hL - hrMin(i,k) )/hrExp(i,k) )
!               expR = EXP(-( hR - hrMin(i,k) )/hrExp(i,k) )
!               areaFraction(i) = expL - expR
!               volFraction (i) =  &
!                   ( ( hL + hrExp(i,k) ) * expL     &
!                   - ( hR + hrExp(i,k) ) * expR )   &
!                   / ( hrMin(i,k) + hrExp(i,k) )
!              ENDIF
!            enddo
!          ELSE
!            do i = 1,nlpb
!              IF ( doRidging(i) .AND. hrExp(i,k) .NE. 0.0_wp ) THEN
!                hL   = MAX( hrMin(i,k), hLimit(l-1) )
!                expL = EXP(-( hL - hrMin(i,k) )/hrExp(i,k) )
!                areaFraction(i) = expL
!                volFraction (i) = ( hL + hrExp(i,k) ) * expL     &
!                    / ( hrMin(i,k) + hrExp(i,k) )
!              ENDIF
!            enddo
!          ENDIF
!        ENDIF   ! endif for SEAICEredistFunc
!!     after computing the fraction ridged ice for this category, apply it
!        do i = 1,nlpb
!         !YY: here ridgedarea is for outer category k, while areafraction is updated in inner loop
!          AREAITD(i,l) = AREAITD(i,l) +areaFraction(i)*ridgedArea(i)
!          HEFFITD(i,l) = HEFFITD(i,l) +volFraction(i)*ridgingHeff(i)
!          HSNOWITD(i,l) = HSNOWITD(i,l)+volFraction(i)*ridgingHsnw(i)*SEAICEsnowFracRidge
!        enddo
!!     category l-loop
!       ENDDO

!$acc end kernels
      ENDDO    !!     enddo k-loop

!     determine if the ridging process needs to be repeated
!     we need a new estimate of the total AREA
!$acc kernels 
!$acc loop
      do i = 1,nlpb
        netArea(i) = 0.0_wp
      enddo 
!$acc loop
      do i = 1,nlpb
!$acc loop seq
        do k=1,nITD
          netArea(i) = netArea(i) + AREAITD(i,k)
        enddo
      enddo
!$acc end kernels
      doRidgeAgain   = .FALSE.
!$acc parallel copy(doRidgeAgain)       !YY: test if doridgeagain updated 
!$acc loop gang,vector private(tmp,areaTooLarge,divAdv) independent
      do i = 1,nlpb
        tmp = netArea(i)+opnWtrFrac(i)
        areaTooLarge = tmp - 1.0_wp .GT. 1.0d-11
        IF ( HEFFM(i) .GT. 0.0_wp .AND. areaTooLarge ) THEN
          doRidging(i) = .TRUE.
          doRidgeAgain   = .TRUE.
          divAdv = (1.0_wp-tmp)*(1.0_wp/SEAICE_deltaTtherm)
          closingRate(i) = MAX( 0.0_wp, -divAdv)
          openingRate(i) = MAX( 0.0_wp,  divAdv)
        ELSE
!     set to zero avoid going through this grid point again
          closingRate(i) = 0.0_wp
          openingRate(i) = 0.0_wp
          doRidging(i)   = .FALSE.
        ENDIF
      enddo
!$acc end parallel   !make sure doridgeagin is copyout

      IF ( doRidgeAgain .AND. n.GE.SEAICEridgingIterMax ) THEN
        write(mitice_runmsg_unit,'(A)') 'SEAICE_DO_RIDGING: *** WARNING ***'
        write(mitice_runmsg_unit,*)'SEAICE_DO_RIDGING: DONOT CONVERGE IN : ',SEAICEridgingIterMax,'iterations' 
      ENDIF
      doRidgeAgain = doRidgeAgain .AND. n.LT.SEAICEridgingIterMax
      IF ( doRidgeAgain ) THEN
        WRITE(mitice_runmsg_unit,'(A,I2,A,I10)')                    &
           'SEAICE_DO_RIDGING: Repeat ridging after iteration ',    &
           n, ' in timestep ', myIter
      ENDIF
      IF ( doRidgeAgain ) CALL seaice_prepare_ridging
      n = n + 1
!     ridging iteration
      ENDDO
      ENDIF     ! ENDIF FOR .not. SEAICEsimpleRidging

!$acc end data

#endif /* SEAICE_ITD */

   end subroutine seaice_do_ridging


#ifdef SEAICE_ITD

   subroutine seaice_itd_remap(heffitdpre, areaitdpre)
!     *===========================================================*
!     | o checks if absolute ice thickness in any category
!     |   exceeds its category limits
!     | o remaps sea ice area and volume
!     |   and associated ice properties in thickness space
!     |   following the remapping scheme of Lipscomb (2001), JGR
!     |
!     | Martin Losch, started in May 2014, Martin.Losch@awi.de
!     | with many fixes by Mischa Ungermann (MU)
!     *===========================================================*

      implicit none

!     === Global variables to be checked and remapped ===
!     AREAITD   :: sea ice area      by category
!     HEFFITD   :: sea ice thickness by category
!     === Global variables to be remappped ===
!     HSNOWITD  :: snow thickness    by category

      real(wp)  :: heffitdPre(nlpb,nITD)
      real(wp)  :: areaitdPre(nlpb,nITD)

!     === Local variables ===
      INTEGER i, k
      INTEGER kDonor, kRecvr
      real(wp) slope
      real(wp) etaMin, etaMax, etam, etap, eta2
      real(wp) dh0, da0, daMax
      
      real(wp)  :: dhActual    (nlpb,nITD)
      real(wp)  :: hActual     (nlpb,nITD)
      real(wp)  :: hActualPre  (nlpb,nITD)
      real(wp)  dheff, darea, dhsnw

      real(wp)  :: hLimitNew   (nlpb,0:nITD)
!     coefficients for represent g(h)
!     g0 :: constant coefficient in g(h)
!     g1 :: linear  coefficient in g(h)
!     hL :: left end of range over which g(h) > 0
!     hL :: right end of range over which g(h) > 0
      real(wp)  :: g0 (nlpb,0:nITD)
      real(wp)  :: g1 (nlpb,0:nITD)
      real(wp)  :: hL (nlpb,0:nITD)
      real(wp)  :: hR (nlpb,0:nITD)
!     local copy of AREAITD
      !real(wp) aLoc(nlpb)
      LOGICAL doRemapping (nlpb)
      LOGICAL debug_remapping     !YY: adjust deltaT to make sure displaced boundaries among Hn-1 and Hn+1


      debug_remapping = debuglevel .ge. 1

!$acc data present(heffitdpre, areaitdpre,hlimit),     &
!$acc      create(dhActual,hActual,hActualPre,hLimitNew,g0,g1,hL,hR,doRemapping)

!     initialisation
!$acc kernels
!$acc loop
      do i = 1,nlpb
        doRemapping(i) = .FALSE.
        IF ( HEFFM(i) .NE. 0.0_wp ) doRemapping(i) = .TRUE.
      enddo
!     do not compute regularized hActual as in seaice_growth, because
!     with regularization, hActual deviates too much from the actual
!     category boundaries and the boundary computation fails too often.
!$acc loop collapse(2)
      DO k=1,nITD
        do i = 1,nlpb
          hActualPre (i,k) = 0.0_wp
          hActual    (i,k) = 0.0_wp
          dhActual   (i,k) = 0.0_wp
          IF ( areaitdPre(i,k) .GT. SEAICE_area_reg ) THEN
           hActualPre(i,k) = heffitdPre(i,k)/areaitdPre(i,k)
          ENDIF
          IF ( AREAITD(i,k) .GT. SEAICE_area_reg ) THEN
          !YY: heff is a grid-average value (equivalent to area=1)
           hActual(i,k) = HEFFITD(i,k)/AREAITD(i,k)
          ENDIF
          dhActual(i,k) = hActual(i,k) - hActualPre(i,k)
        enddo
      ENDDO

!     compute new category boundaries

!$acc loop
      do i = 1,nlpb
        hLimitNew(i,0) = hLimit(0)
      enddo
!YY: pay attention. test k dependency 
!$acc loop collapse(2) private(slope)
      DO k=1,nITD-1
        do i = 1,nlpb
          IF ( hActualPre(i,k)  .GT.SEAICE_eps .AND.    &
              hActualPre(i,k+1).GT.SEAICE_eps ) THEN
           slope = ( dhActual(i,k+1) - dhActual(i,k) )  &
               /( hActualPre(i,k+1) - hActualPre(i,k) )
           hLimitNew(i,k) =   hLimit(k) + dhActual(i,k) &
               +     slope * ( hLimit(k) - hActualPre(i,k) )
          ELSEIF ( hActualPre(i,k)  .GT.SEAICE_eps ) THEN
           hLimitNew(i,k) = hLimit(k) + dhActual(i,k)
          ELSEIF ( hActualPre(i,k+1).GT.SEAICE_eps ) THEN
           hLimitNew(i,k) = hLimit(k) + dhActual(i,k+1)
          ELSE
           hLimitNew(i,k) = hLimit(k)
          ENDIF
!     After computing the new boundary, check
!     (1) if it is between two adjacent thicknesses
          IF ( ( AREAITD(i,k).GT.SEAICE_area_reg .AND.    &
                hActual(i,k) .GE. hLimitNew(i,k) ) .OR.   &
              ( AREAITD(i,k+1).GT.SEAICE_area_reg .AND.   &
                hActual(i,k+1) .LE. hLimitNew(i,k) ) )    &
              doRemapping(i) = .FALSE.
!     (2) that it is been the old boudnaries k-1 and k+1
!     (Note from CICE: we could allow this, but would make the code
!     more complicated)
          IF ( ( hLimitNew(i,k) .GT. hLimit(k+1) ) .OR.   &
               ( hLimitNew(i,k) .LT. hLimit(k-1) ) )      &
               doRemapping(i) = .FALSE.
        enddo
      ENDDO
!$acc end kernels

!     Report problems, if there are any. Because this breaks optimization
!     do not do it by default.
!     Where doRemapping is false, the rebinning of seaice_itd_redist
!     (called at the end) will take care of shifting the ice.
      IF ( debug_remapping )  then 
!$acc update self(AREAITD, hActual, hActualPre, hLimitNew, doRemapping)
           CALL SEAICE_ITD_REMAP_CHECK_BOUNDS(                          &
             AREAITD, hActual, hActualPre, hLimitNew, doRemapping)
      ENDIF
!     computing the upper limit of the thickest category does not require
!     any checks and can be computed now
!$acc kernels loop
      do i = 1,nlpb
        hLimitNew(i,nITD) = hLimit(nITD)
        IF ( AREAITD(i,nITD).GT.SEAICE_area_reg )                          &
            hLimitNew(i,nITD) = MAX( 3.0_wp*hActual(i,nITD)                   &
            - 2.0_wp * hLimitNew(i,nITD-1), hLimit(nITD-1) )
      enddo
!$acc end kernels
      
!     end of limit computation, now compute the coefficients of the
!     linear approximations of g(h) => g(eta) = g0 + g1*eta

!     CICE does something specical for the first category.
!     compute coefficients for 1st category
      k = 1
!$acc kernels copyin(k)
!$acc loop
      do i = 1,nlpb
!     initialisation
       ! aLoc(i) = AREAITD(i,k)
!     initialise hL and hR
!     this single line is different from the code that follows below
!     for all categories
        hL(i,k) = hLimitNew(i,k-1)
        hR(i,k) = hLimit(k)
      enddo
!$acc end kernels

      !CALL seaice_itd_remap_linear(g0(:,k), g1(:,k),hL(:,k), hR(:,k), &
      !    hActual(:,k), aLoc,doRemapping)
      CALL seaice_itd_remap_linear(g0(:,k), g1(:,k),hL(:,k), hR(:,k), &
          hActual(:,k), AREAITD(:,k),doRemapping)

!     Find area lost due to melting of thin (category 1) ice

!$acc kernels copyin(k)
!$acc loop private(dh0,etaMax,da0,daMax)
      do i = 1,nlpb
        IF ( doRemapping(i) .AND.AREAITD(i,k) .GT. SEAICE_area_reg ) THEN
!MU if melting of ice in category 1
         IF ( dhActual(i,k) .LT. 0.0_wp ) THEN
!     integrate g(1) from zero to abs(dhActual)
!MU dh0 is max thickness of ice in first category that is melted
          dh0    = MIN(-dhActual(i,k),hLimit(k))
          etaMax = MIN(dh0,hR(i,k)) - hL(i,k)
          IF ( etaMax > 0.0_wp ) THEN
!MU da0 is /int_0^dh0 g dh
           da0 = g0(i,k)*etaMax + g1(i,k)*etaMax*etaMax*0.5_wp
           daMax = AREAITD(i,k) * ( 1.0_wp - hActual(i,k)/hActualPre(i,k))
           da0 = MIN( da0, daMax )
!MU adjust thickness to conserve volume
           IF ( (AREAITD(i,k)-da0) .GT. SEAICE_area_reg ) THEN
             hActual(i,k) = hActual(i,k) * AREAITD(i,k)/( AREAITD(i,k) - da0 )
           ELSE
             hActual(i,k) = ZERO
             da0 = AREAITD(i,k)
           ENDIF
!MU increase open water fraction
           AREAITD(i,k) = AREAITD(i,k) - da0
          ENDIF
         ELSE
!MU H_0* = F_0 * dT
          hLimitNew(i,k-1) = MIN( dhActual(i,k), hLimit(k) )
         ENDIF
        ENDIF
      enddo
!$acc end kernels

!     compute all coefficients
!$acc kernels loop collapse(2)
      DO k=1,nITD
        do i = 1,nlpb
!     initialisation
         !aLoc(i) = AREAITD(i,k)
!     initialise hL and hR
         hL(i,k) = hLimitNew(i,k-1)
         hR(i,k) = hLimitNew(i,k)
        enddo
      enddo
!$acc end kernels
      DO k=1,nITD
        CALL SEAICE_ITD_REMAP_LINEAR(g0(:,k), g1(:,k),hL(:,k), hR(:,k), &
           hActual(:,k), AREAITD(:,k), doRemapping)
      ENDDO

!$acc kernels
!$acc loop private(dheff,darea,dhsnw,etaMin,etaMax,kDonor,kRecvr,etam,etap,eta2)
      do i = 1,nlpb
!$acc loop seq
        DO k=1,nITD-1
          dheff = 0.0_wp
          darea = 0.0_wp
          IF ( doRemapping(i) ) THEN
!     compute integration limits in eta space
            IF ( hLimitNew(i,k) .GT. hLimit(k) ) THEN
              etaMin = MAX(       hLimit(k), hL(i,k)) - hL(i,k)
              etaMax = MIN(hLimitNew(i,k), hR(i,k)) - hL(i,k)
              kDonor = k
              kRecvr = k+1
            ELSE
              etaMin = 0.0_wp
              etaMax = MIN(hLimit(k), hR(i,k+1)) - hL(i,k+1)
              kDonor = k+1
              kRecvr = k
            ENDIF
!     compute the area and volume to be moved
            IF ( etaMax .GT. etaMin ) THEN
              etam  = etaMax-etaMin
              etap  = etaMax+etaMin
              eta2  = 0.5*etam*etap
              darea = g0(i,kDonor)*etam + g1(i,kDonor)*eta2
              dheff = g0(i,kDonor)*eta2                         &
                   +  g1(i,kDonor)*(etaMax**3-etaMin**3)*third  &
                   +  darea*hL(i,kDonor)
            ENDIF
!     ... or shift entire category, if nearly all ice is to be shifted.
            IF ( (darea .GT.AREAITD(i,kDonor)-SEAICE_eps).OR.    &
                 (dheff .GT.HEFFITD(i,kDonor)-SEAICE_eps) ) THEN
              darea = AREAITD(i,kDonor)
              dheff = HEFFITD(i,kDonor)
            ENDIF
!     regularize: reset to zero, if there is too little ice to be shifted ...
            IF ( (darea .LT. SEAICE_eps).OR.(dheff .LT. SEAICE_eps) ) THEN
              darea  = 0.0_wp
              dheff  = 0.0_wp
            ENDIF
!     snow scaled by area
            IF ( AREAITD(i,kDonor) .GT. SEAICE_area_reg ) THEN
!     snow scaled by area (why not volume?), CICE also does it in this way
              dhsnw = darea/AREAITD(i,kDonor) * HSNOWITD(i,kDonor)
            ELSE
              dhsnw = HSNOWITD(i,kDonor)
            ENDIF
!     apply increments
            HEFFITD(i,kRecvr) = HEFFITD(i,kRecvr) + dheff
            HEFFITD(i,kDonor) = HEFFITD(i,kDonor) - dheff
            AREAITD(i,kRecvr) = AREAITD(i,kRecvr) + darea
            AREAITD(i,kDonor) = AREAITD(i,kDonor) - darea
            HSNOWITD(i,kRecvr)=HSNOWITD(i,kRecvr) + dhsnw
            HSNOWITD(i,kDonor)=HSNOWITD(i,kDonor) - dhsnw
!     end if doRemapping
          ENDIF
        enddo
      ENDDO
!$acc end kernels
!$acc end data

   end subroutine seaice_itd_remap


   subroutine seaice_itd_remap_linear( g0, g1, hL, hR,      &
     hActual, area, doRemapping)

!     *===========================================================*
!     | o compute coefficients g0, g1 for piece-wise linear fit
!     |    g(h) = g0 + g1*h
!     | o compute range boundaries hL, hR for this linear fit
!     |
!     | Martin Losch, May 2014, Martin.Losch@awi.de
!     *===========================================================*

      implicit none

!     OUTPUT: coefficients for representing g(h)
!     g0 :: constant coefficient in g(h)
!     g1 :: linear  coefficient in g(h)
!     hL :: left end of range over which g(h) > 0
!     hL :: right end of range over which g(h) > 0
      real(wp) :: g0 (nlpb)
      real(wp) :: g1 (nlpb)
      real(wp) :: hL (nlpb)
      real(wp) :: hR (nlpb)
!     INPUT:
!     hActual :: ice thickness of current category
!     area    :: ice concentration of current category
      real(wp) :: hActual (nlpb)
      real(wp) :: area    (nlpb)
!     regularization constants
!      real(wp) ::  SEAICE_area_reg
!      real(wp) ::  SEAICE_eps
!     doRemapping :: mask where can be done, excludes points where
!                    new category limits are outside certain bounds
      LOGICAL  ::  doRemapping (nlpb)

      INTEGER i
!     recip_etaR :: reciprocal of range interval in eta space
!     etaNoR   :: ratio of distance to lower limit over etaR
      real(wp) ::  auxCoeff, recip_etaR, etaNoR
      real(wp), parameter ::  sixth = 0.666666666666666666666666666_wp

!$acc kernels present(g0, g1, hL, hR,hActual, area,doRemapping)
!$acc loop private(recip_etaR,etaNoR,auxCoeff)
      do i = 1,nlpb
        g0(i) = 0.0_wp
        g1(i) = 0.0_wp
        IF ( doRemapping(i) .AND. area(i) .GT. SEAICE_area_reg .AND.    &
            hR(i) - hL(i) .GT. SEAICE_eps ) THEN
!     change hL and hR if hActual falls outside the central third of the range
          IF ( hActual(i) .LT. (2.0_wp*hL(i) + hR(i))*third ) THEN
           hR(i) = 3.0_wp * hActual(i) - 2.0_wp * hL(i)
          ELSEIF ( hActual(i).GT.(hL(i)+2.0_wp*hR(i))*third ) THEN
           hL(i) = 3.0_wp * hActual(i) - 2.0_wp * hR(i)
          ENDIF
!     calculate new etaR = hR - hL;
!     catch the case of hR=hL, which can happen when hActual=hR or hL
!     before entering this routine; in this case g0=g1=0.
         recip_etaR = 0.0_wp
!MU         IF ( hR(i,j) .GT. hL(i,j) ) ! crucial change; lets the model explode
         IF ( hR(i) - hL(i) .GT. SEAICE_eps )           &
             recip_etaR = 1.0_wp / (hR(i) - hL(i))
!     some abbreviations to avoid computing the same thing multiple times
         etaNoR     = (hActual(i) - hL(i))*recip_etaR
         auxCoeff   = 6.0_wp * area(i)*recip_etaR
!     equations (14) of Lipscomb (2001), JGR
         g0(i) = auxCoeff*( sixth - etaNoR )
         g1(i) = 2.0_wp * auxCoeff*recip_etaR*( etaNoR - 0.5_wp )
        ELSE
!     not doRemapping
!     reset hL and hR
         hL(i) = 0.0_wp
         hR(i) = 0.0_wp
        ENDIF
      enddo
!$acc end kernels

   end subroutine seaice_itd_remap_linear


   subroutine seaice_itd_remap_check_bounds(                  &
         AREAITD, hActual, hActualPre, hLimitNew, doRemapping)

!     *===========================================================*
!     | o where doRemapping = .FALSE. print a warning
!     |
!     | Martin Losch, May 2014, Martin.Losch@awi.de
!     *===========================================================*

      implicit none

!INPUT PARAMETERS: ===================================================
!     hActual :: ice thickness of current category
      real(wp) :: hActual   (nlpb,1:nITD)
      real(wp) :: hActualPre(nlpb,1:nITD)
!     hLimitNew :: new "advected" category boundaries after seaice_growth
      real(wp) :: hLimitNew (nlpb,0:nITD)
!     AREAITD :: ice concentration of current category
      real(wp) :: AREAITD   (nlpb,nITD)
!     doRemapping :: mask where can be done, excludes points where
!                    new category limits are outside certain bounds
      LOGICAL ::  doRemapping (nlpb)

      INTEGER i, k
      CHARACTER(39) tmpBuf

      do i = 1,nlpb
        IF (.NOT.doRemapping(i) ) THEN
          DO k=1,nITD-1
            WRITE(tmpBuf,'(A,I5,A,I10)') ' at (', i, ') in timestep ', myIter
            IF ( AREAITD(i,k).GT.SEAICE_area_reg .AND.     &
                 hActual(i,k) .GE. hLimitNew(i,k) ) THEN
              WRITE(mitice_runmsg_unit,'(A,I3,A)')         &
              'SEAICE_ITD_REMAP: hActual(k) >= hLimitNew(k) for category ', k, tmpBuf
!             PRINT *, hActual(i,j,k),hLimitNew(i,j,k), hLimit(k)
            ENDIF
            IF ( AREAITD(i,k+1).GT.SEAICE_area_reg .AND.   &
                 hActual(i,k+1) .LE. hLimitNew(i,k) ) THEN
              WRITE(mitice_runmsg_unit,'(A,I3,A)')         &
              'SEAICE_ITD_REMAP: hActual(k+1) <= hLimitNew(k) for category ', k, tmpBuf
!              PRINT '(8(1X,E10.4))',
!     &            AREAITD(i,j,k+1,bi,bj), hActual(i,j,k+1),
!     &            hActualPre(i,j,k+1),
!     &            AREAITD(i,j,k,bi,bj), hActual(i,j,k),
!     &            hActualPre(i,j,k),
!     &            hLimitNew(i,j,k), hLimit(k)
            ENDIF
            IF ( hLimitNew(i,k) .GT. hLimit(k+1) ) THEN
              WRITE(mitice_runmsg_unit,'(A,I3,A)')         &
              'SEAICE_ITD_REMAP: hLimitNew(k) > hLimitNew(k+1) for category ', k, tmpBuf
            ENDIF
            IF ( hLimitNew(i,k) .LT. hLimit(k-1) ) THEN
              WRITE(mitice_runmsg_unit,'(A,I3,A)')         &
              'SEAICE_ITD_REMAP: hLimitNew(k) < hLimitNew(k-1) for category ', k, tmpBuf
            ENDIF
          ENDDO
        ENDIF
      enddo

   end subroutine seaice_itd_remap_check_bounds


   subroutine seaice_prepare_ridging

!  *===========================================================*
!   o compute ridging parameters according to Thorndyke et al
!     (1975), Hibler (1980), Bitz et al (2001) or
!     Lipscomb et al (2007)
!  
!   Martin Losch, Apr. 2014, Martin.Losch@awi.de
!  *===========================================================*

      implicit none

!     === Local variables ===
      INTEGER i, k
!     variables related to ridging schemes
!     gSum        :: cumulative distribution function G
      real(wp) :: gSum(nlpb,-1:nITD)
      real(wp)    recip_gStar, recip_aStar, tmp
!     Regularization values squared
      real(wp) area_reg_sq, hice_reg_sq

!$acc data present(hActual,HEFFITD,AREAITD,opnWtrFrac,partFunc,      &
!$acc             hrMin,ridgeRatio,HEFFM,ridgingModeNorm),     &
!$acc      create(gSum)


!     regularization constants
      area_reg_sq = SEAICE_area_reg * SEAICE_area_reg
      hice_reg_sq = SEAICE_hice_reg * SEAICE_hice_reg
!$acc kernels copyin(area_reg_sq,hice_reg_sq)
!$acc loop collapse(2) private(tmp)
      DO k=1,nITD
       do i = 1,nlpb
         hActual(i,k) = 0.0_wp
         IF ( HEFFITD(i,k) .GT. 0.0_wp ) THEN
!     regularize as in seaice_growth: compute hActual with regularized
!     AREA and regularize from below with a minimum thickness
          tmp = HEFFITD(i,k) / SQRT( AREAITD(i,k)**2 + area_reg_sq )
          hActual(i,k) = SQRT(tmp * tmp + hice_reg_sq)
         ENDIF
       enddo
      ENDDO

!     compute the cumulative thickness distribution function gSum
!$acc loop
      do i = 1,nlpb
        gSum(i,-1) = 0.0_wp
        gSum(i, 0)  = 0.0_wp
        IF (opnWtrFrac(i).GT.SEAICE_area_floor) gSum(i,0) = opnWtrFrac(i)
      enddo
!$acc loop
      do i = 1,nlpb
!$acc loop seq
       do k = 1, nITD
         gSum(i,k) = gSum(i,k-1)
         IF ( AREAITD(i,k) .GT. SEAICE_area_floor ) gSum(i,k) = gSum(i,k) + AREAITD(i,k)
       enddo
      enddo
!     normalize
!$acc loop
      do i = 1,nlpb
!$acc loop seq
       do k = 0, nITD
         IF ( gSum(i,nITD).NE.0.0_wp ) gSum(i,k) = gSum(i,k) / gSum(i,nITD)
       enddo 
      enddo
!$acc end kernels

!     Compute the participation function

      IF ( SEAICEpartFunc .EQ. 0 ) THEN
!     Thorndike et al. (1975) discretize b(h) = (2/Gstar) * (1 - G(h)/Gstar)
!     The expressions for the partition function partFunc are found by
!     integrating b(h)g(h) between the category boundaries.
       recip_gStar = 1.0_wp / SEAICEgStar
!$acc kernels copyin(recip_gStar)
!$acc loop collapse(2) independent
       DO k = 0, nITD
        do i = 1,nlpb
          partFunc(i,k) = 0.0_wp
          IF ( gSum(i,k) .LT. SEAICEgStar ) THEN
           partFunc(i,k) = (gSum(i,k)-gSum(i,k-1)) * recip_gStar     &
               *( 2.0_wp - (gSum(i,k-1)+gSum(i,k))*recip_gStar)
          ELSEIF (  gSum(i,k-1) .LT. SEAICEgStar .AND. gSum(i,k) .GE. SEAICEgStar ) THEN
           partFunc(i,k) = (SEAICEgStar-gSum(i,k-1)) * recip_gStar   &
               *( 2.0_wp - (gSum(i,k-1)+SEAICEgStar)*recip_gStar)
          ENDIF
        enddo
       ENDDO
!$acc end kernels
      ELSEIF  ( SEAICEpartFunc .EQ. 1 ) THEN
!     Lipscomb et al. (2007) discretize b(h) = exp(-G(h)/astar) into
!     partFunc(n) = [exp(-G(n-1)/astar - exp(-G(n)/astar] / [1-exp(-1/astar)].
!     The expression is found by integrating b(h)g(h) between the category
!     boundaries.
       recip_astar = 1.0_wp / SEAICEaStar
       tmp = 1.0_wp / ( 1.0_wp - EXP( -recip_astar ) )
!$acc kernels copyin(recip_astar,tmp)
!     abuse gSum as a work array
!$acc loop collapse(2)
!!YY: potential error .check later
       do k = -1, nITD    !YY:not k is not independent
        do i = 1,nlpb
          gSum(i,k)     = EXP(-gSum(i,k)*recip_astar) * tmp
        enddo
       enddo
!$acc loop collapse(2) independent
       do k = 0, nITD    !YY:not k is not independent
        do i = 1,nlpb
          partFunc(i,k) = gSum(i,k-1) - gSum(i,k)
        enddo
       enddo
!$acc end kernels
      ELSE
       STOP 'WARNING:  SEAICEpartFunc > 1 not implemented'
      ENDIF

!     Compute variables of ITD of ridged ice
!     ridgeRatio :: mean ridge thickness/ thickness of ridging ice
!     hrMin      :: min ridge thickness
!     hrMax      :: max ridge thickness   (SEAICEredistFunc = 0)
!     hrExp      :: ridge e-folding scale (SEAICEredistFunc = 1)
      IF ( SEAICEredistFunc .EQ. 0 ) THEN
!$acc kernels loop collapse(2)
        DO k = 1, nITD
         do i = 1,nlpb
           hrMin(i,k)  = 0.0_wp
           hrMax(i,k)  = 0.0_wp
!       avoid divisions by zero
           ridgeRatio(i,k) = 1.0_wp
         enddo
        ENDDO
!$acc end kernels
      ELSEIF ( SEAICEredistFunc .EQ. 1 ) THEN
!$acc kernels loop collapse(2)
        DO k = 1, nITD
         do i = 1,nlpb
           hrMin(i,k)  = 0.0_wp
           hrExp(i,k)  = 0.0_wp
!       avoid divisions by zero
           ridgeRatio(i,k) = 1.0_wp
         enddo
        ENDDO
!$acc end kernels
      ENDIF
      IF ( SEAICEredistFunc .EQ. 0 ) THEN
!     Assume ridged ice is uniformly distributed between hrmin and hrmax.
!     (Hibler, 1980)
!$acc kernels loop collapse(2)
       DO k = 1, nITD
        do i = 1,nlpb
          IF ( hActual(i,k) .GT. 0.0_wp ) THEN
!     This is the original Hibler (1980) scheme:
           hrMin(i,k) = 2.0_wp * hActual(i,k)
           hrMax(i,k) = 2.0_wp * SQRT(hActual(i,k)*SEAICEhStar)
!     CICE does this in addition, so that thick ridging ice is not required
!     to raft:
           hrMin(i,k) = MIN(hrMin(i,k),hActual(i,k)+SEAICEmaxRaft)
           hrMax(i,k) = MAX(hrMax(i,k),  hrMin(i,k)+SEAICE_hice_reg)
!
           ridgeRatio(i,k) = 0.5_wp * (hrMax(i,k)+hrMin(i,k))/hActual(i,k)
          ENDIF
        enddo
       ENDDO
!$acc end kernels
      ELSEIF ( SEAICEredistFunc .EQ. 1 ) THEN
!     Follow Lipscomb et al. (2007) and model ridge ITD as an exponentially
!     decaying function
!$acc kernels loop collapse(2) private(tmp)
       DO k = 1, nITD
        do i = 1,nlpb
          IF ( hActual(i,k) .GT. 0.0_wp ) THEN
           tmp = hActual(i,k)
           hrMin(i,k) = MIN(2.0_wp * tmp, tmp+SEAICEmaxRaft)
           hrExp(i,k) = SEAICEmuRidging*SQRT(tmp)
!     arent we missing a factor 0.5 here?
           ridgeRatio(i,k)=(hrMin(i,k)+hrExp(i,k))/tmp
          ENDIF
        enddo
       ENDDO
!$acc end kernels
      ELSE
       STOP 'WARNING: SEAICEredistFunc > 1 not implemented'
      ENDIF

!     Compute the norm of the ridging mode N (in Lipscomp et al 2007)
!     or omega (in Bitz et al 2001):
!     rigdingModeNorm = net ice area removed / total area participating.
!     For instance, if a unit area of ice with thickness = 1 participates in
!     ridging to form a ridge with a = 1/3 and thickness = 3, then
!     rigdingModeNorm = 1 - 1/3 = 2/3.
!YY: new-formed ridge with area 1/3 and thickness 3m. then removed or redistributed
!YY: ice area from the thickness 1m is  1-1/3.
!$acc kernels
!$acc loop
      do i = 1,nlpb
        ridgingModeNorm(i) = partFunc(i,0)
      ENDDO
!$acc loop collapse(2)
      DO k = 1, nITD
       do i = 1,nlpb
         partFunc(i,k) = partFunc(i,k) * HEFFM(i)
       enddo
      ENDDO
!$acc loop
      do i = 1,nlpb
!$acc loop seq
       DO k = 1, nITD
         ridgingModeNorm(i) = ridgingModeNorm(i)                 &
             + partFunc(i,k)*( 1.0_wp - 1.0_wp/ridgeRatio(i,k) )
       enddo
      ENDDO
!     avoid division by zero
!$acc loop
      do i = 1,nlpb
        IF (ridgingModeNorm(i).LE.0.0_wp) ridgingModeNorm(i) = 1.0_wp
      enddo
!$acc end kernels
!$acc end data

   end subroutine seaice_prepare_ridging
#endif /* SEAICE_ITD */


   subroutine  seaice_calc_ice_strength

!     *==========================================================*
!     | o compute ice strengh PRESS0
!     *==========================================================*

      implicit none

      INTEGER i,  k
      real(wp) tmpscal1, tmpscal2

#ifdef SEAICE_ITD
!$acc data present(opnWtrFrac,AREA,HEFF,HEFFM,PRESS0,ZMAX,ZMIN,     &
!$acc         partFunc,hActual,hrMin,ridgeRatio,ridgingModeNorm)
#else
!$acc data present(AREA,HEFF,HEFFM,PRESS0,ZMAX,ZMIN)
#endif

#ifdef SEAICE_ITD
!$acc kernels loop
      do i = 1,nlpb
        opnWtrFrac(i) = 1.0_wp - AREA(i)
      enddo
!$acc end kernels

      IF ( useHibler79IceStrength ) THEN
#else
      IF ( .TRUE. ) THEN
#endif /* SEAICE_ITD */
!$acc kernels loop private(tmpscal1,tmpscal2)
       do i = 1,nlpb
!--   now set up ice pressure and viscosities
         IF ( (HEFF(i).LE.SEAICEpresH0).AND.(SEAICEpresPow0.NE.1) ) THEN
           tmpscal1=MAX(HEFF(i)/SEAICEpresH0,ZERO)
           tmpscal2=SEAICEpresH0*(tmpscal1**SEAICEpresPow0)
         ELSEIF ( (HEFF(i).GT.SEAICEpresH0).AND.(SEAICEpresPow1.NE.1) ) THEN
           tmpscal1=MAX(HEFF(i)/SEAICEpresH0,ZERO)
           tmpscal2=SEAICEpresH0*(tmpscal1**SEAICEpresPow1)
         ELSE
           tmpscal2=HEFF(i)
         ENDIF
         PRESS0(i)=SEAICE_strength*tmpscal2                 &
             *EXP(-SEAICE_cStar*(SEAICE_area_max-AREA(i)))
         ZMAX(i)  = SEAICE_zetaMaxFac*PRESS0(I)
         ZMIN(i)  = SEAICE_zetaMin
         PRESS0(i)= PRESS0(i)*HEFFM(i)
       enddo
!$acc end kernels
#ifdef SEAICE_ITD
      ELSE
!     not useHiber79IceStrength
!$acc kernels loop
       do i = 1,nlpb
         PRESS0(i) = 0.0_wp
       enddo
!$acc end kernels

       call seaice_prepare_ridging

       IF ( SEAICEredistFunc .EQ. 0 ) THEN
         tmpscal1 = 1.0_wp / 3.0_wp
!$acc kernels copyin(tmpscal1)
!$acc loop
         do i = 1,nlpb
!$acc loop seq
          do k = 1, nITD
!     replace (hrMax**3-hrMin**3)/(hrMax-hrMin) by identical
!     hrMax**2+hrMin**2 + hrMax*hrMin to avoid division by potentially
!     small number
           IF ( partFunc(i,k) .GT. 0.0_wp ) then                   
             PRESS0(i) = PRESS0(i) + partFunc(i,k) *               &
             ( - hActual(i,k)**2 + ( hrMax(i,k)**2 + hrMin(i,k)**2 &
             + hrMax(i,k)*hrMin(i,k) )*tmpscal1 / ridgeRatio(i,k) )
           ENDIF
          enddo 
         enddo
!$acc end kernels
       ELSEIF ( SEAICEredistFunc .EQ. 1 ) THEN
!$acc kernels
!$acc loop
         do i = 1,nlpb
!$acc loop seq
          do k = 1, nitd
           PRESS0(i) = PRESS0(i)     &
               + partFunc(i,k) * ( - hActual(i,k)**2 +   &
               (           hrMin(i,k)*hrMin(i,k)         &
               + 2.0_wp *  hrMin(i,k)*hrExp(i,k)         &
               + 2.0_wp *  hrExp(i,k)*hrExp(i,k)         &
               )/ridgeRatio(i,k) )
          enddo 
         enddo
!$acc end kernels
       ENDIF

       tmpscal1 = SEAICE_cf*0.5_wp*gravity*(rau0-SEAICE_rhoIce)*SEAICE_rhoIce/rau0
!$acc kernels loop copyin(tmpscal1)
       do i = 1,nlpb
         PRESS0(i) = PRESS0(i)/ridgingModeNorm(i) * tmpscal1
         ZMAX(i)  = SEAICE_zetaMaxFac*PRESS0(i)
         ZMIN(i)  = SEAICE_zetaMin
         PRESS0(i)= PRESS0(i)*HEFFM(i)
       enddo
!$acc end kernels
#endif /* SEAICE_ITD */
      ENDIF
!$acc end data

   end subroutine  seaice_calc_ice_strength


   subroutine seaice_calc_strainrates(uFld, vFld, e11Loc, e22Loc, e12Loc, iStep)

!     *==========================================================*
!     | o compute strain rates from ice velocities
!     *==========================================================*
!     | written by Martin Losch, Apr 2007
!     *==========================================================*

      implicit none

!     !INPUT/OUTPUT PARAMETERS:
!     e11Loc :: strain rate tensor, component 1,1
!     e22Loc :: strain rate tensor, component 2,2
!     e12Loc :: strain rate tensor, component 1,2
!     iStep  :: Sub-time-step number
      real(wp), intent(in)  ::  uFld(nlpb)
      real(wp), intent(in)  ::  vFld(nlpb)
      real(wp)  ::  e11Loc(nlpb)
      real(wp)  ::  e22Loc(nlpb)
      real(wp)  ::  e12Loc(nlpbz)
      INTEGER iStep

      INTEGER i
!     hFacU, hFacV :: determine the no-slip boundary condition
      real(wp) hFacU, hFacV, noSlipFac
!     auxillary variables that help writing code that
!     vectorizes even after TAFization
      real(wp)  ::  dudx (nlpb)
      real(wp)  ::  dvdy (nlpb)
      !real(wp)  ::  dudy (nlpb)
      !real(wp)  ::  dvdx (nlpb)
      !real(wp)  ::  uave (nlpb)
      !real(wp)  ::  vave (nlpb)

      real(wp)  ::  vecFld(nlpb,2),SIMaskUV(nlpb,2) 
      real(wp)  ::  dlc(nlpb,2)
      integer   ::  uel,uei,uep,vnl,vni,vnp
      integer :: i1, i2, i3, i4
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4

!$acc data present(uFld, vFld, e11Loc, e22Loc, e12Loc,dxC,dyC,SIMaskU,SIMaskV,  &
!$acc              ue,vn,recip_dxW,recip_dyS,recip_rAz,HEFFM,recip_dyW,recip_dxS),    &
!$acc      create(dudx,dvdy,vecFld,SIMaskUV,dlc)

      noSlipFac = 0.0_wp
      IF ( SEAICE_no_slip ) noSlipFac = 1.0_wp

!YY,ZY: confirm ?? 
!$acc kernels copyin(noSlipFac)
!$acc loop
      do i = 1,nlpb
        vecFld(i, iu) = uFld(i)
        vecFld(i, iv) = vFld(i)
        dlc(i, iu) = dxC(i)
        dlc(i, iv) = dyC(i)
        SIMaskUV(i,iu) = SIMaskU(i)
        SIMaskUV(i,iv) = SIMaskV(i)
      end do
!$acc loop private(uel,uei,uep,vnl,vni,vnp)
      do i = 1,nlpb
        uel = ue(i, 1)
        uei = ue(i, 2)
        uep = ue(i, 3)
        vnl = vn(i, 1)
        vni = vn(i, 2)
        vnp = vn(i, 3)
        dudx(i) = recip_dxW(i) * (vecFld(uel,uei)*uep-vecFld(i,iu))
        dvdy(i) = recip_dyS(i) * (vecFld(vnl,vni)*vnp-vecFld(i,iv))
      enddo
!     evaluate strain rates at C-points
!$acc loop
      do i = 1,nlpb
        !YY: here metric terms are removed. 
        e11Loc(i) = dudx(i)
        e22Loc(i) = dvdy(i)
      enddo
!YY: comment for now
!!--     for OBCS: assume no gradient beyong OB
!      do i = 1,nlpb
!        e11Loc(i) = e11Loc(i)*maskInC(i)
!        e22Loc(i) = e22Loc(i)*maskInC(i)
!      enddo

!     abbreviations at Z-points
!YY,ZY:  e12 is defined in Z-point
!YY,ZY: hfacu, hfacV definition may be revisited 
!$acc loop private(i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,k4,hFacU,hFacV)
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
!     evaluate strain rates at Z-points
!YY, ZY : to confirm later
        hFacU = SIMaskUV(i,iu) - SIMaskUV(i3,j3)
        hFacV = SIMaskUV(i,iv) - SIMaskUV(i4,j4)
        e12Loc(i) = 0.5_wp* (VecFld(i1, j1)*dlc(i1, j1)*k1 + VecFld(i2, j2)*dlc(i2, j2)*k2 &
               - VecFld(i3, j3)*dlc(i3, j3)*k3 - VecFld(i4, j4)*dlc(i4, j4)*k4)*recip_rAz(i) &
           *HEFFM(i1)*HEFFM(i3)*HEFFM(i4)*HEFFM( tw(i3) ) &
           + noSlipFac * (                            &
                (VecFld(i1,j1)*k1+VecFld(i3,j3)*k3) *recip_dyW(i1) * hFacU  &
              + (VecFld(i2,j2)*k2+VecFld(i4,j4)*k4) *recip_dxS(i2) * hFacV )
      end do
!$acc end kernels
!$acc end data

   end subroutine seaice_calc_strainrates





end module mitice_itd
