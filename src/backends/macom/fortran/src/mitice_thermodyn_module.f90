module mitice_thermodyn
use mitice_parameters
use mitice_vars
use mitice_itd
use mod_csp_basic
use mod_misc_basic
use mod_mpi_variables
use mod_mpi_interfaces

   implicit none

contains

   subroutine seaice_growth
!     *==========================================================*
!     | o Updata ice thickness and snow depth
!     *==========================================================*

      implicit none

! unit/sign convention:
!    Within the thermodynamic computation all stocks, except HSNOW,
!      are in 'effective ice meters' units, and >0 implies more ice.
!    This holds for stocks due to ocean and atmosphere heat,
!      at the outset of 'PART 2: determine heat fluxes/stocks'
!      and until 'PART 7: determine ocean model forcing'
!    This strategy minimizes the need for multiplications/divisions
!      by ice fraction, heat capacity, etc. The only conversions that
!      occurs are for the HSNOW (in effective snow meters) and
!      PRECIP (fresh water m/s).
!
! HEFF is effective Hice thickness (m3/m2)
! HSNOW is Heffective snow thickness (m3/m2)
! HSALT is Heffective salt content (g/m2)
! AREA is the seaice cover fraction (0<=AREA<=1)
! Q denotes heat stocks -- converted to ice stocks (m3/m2) early on
!
! For all other stocks/increments, such as d_HEFFbyATMonOCN
! or a_QbyATM_cover, the naming convention is as follows:
!    The prefix 'a_' means available, the prefix 'd_' means delta
!       (i.e. increment), and the prefix 'r_' means residual.
!    The suffix '_cover' denotes a value for the ice covered fraction
!       of the grid cell, whereas '_open' is for the open water fraction.
!    The main part of the name states what ice/snow stock is concerned
!       (e.g. QbyATM or HEFF), and how it is affected (e.g. d_HEFFbyATMonOCN
!       is the increment of HEFF due to the ATMosphere extracting heat from the
!       OCeaN surface, or providing heat to the OCeaN surface).

      INTEGER i
      INTEGER kSurface
      INTEGER IT
      CHARACTER(lc) msgBuf
!     constants
      real(wp)  ::  tempFrz, ICE2SNOW, SNOW2ICE
      real(wp)  ::  QI, QS, recip_QI
      real(wp)  ::  lhSublim

! conversion factors to go from Q (W/m2) to HEFF (ice meters)
      real(wp) ::  convertQ2HI, convertHI2Q
! conversion factors to go from precip (m/s) unit to HEFF (ice meters)
      real(wp) :: convertPRECIP2HI, convertHI2PRECIP
!     Factor by which we increase the upper ocean friction velocity (u*) when
!     ice is absent in a grid cell  (dimensionless)
      real(wp) :: MixedLayerTurbulenceFactor
!     wind speed square
      real(wp) :: SPEED_SQ
!     Regularization values squared
      real(wp) :: area_reg_sq, hice_reg_sq
!     pathological cases thresholds
      real(wp) :: heffTooHeavy
!     local copy of surface layer thickness in meters
      real(wp) :: dzSurf
      real(wp) :: recip_multDim, recip_deltaTtherm
      real(wp) :: recip_HO
!     local value (=1/ice thickness)
      real(wp) :: recip_HH
#ifndef SEAICE_ITD
!     facilitate multi-category snow implementation
      real(wp) :: pFac, pFacSnow
#endif
!     additional factors accounting for a non-uniform sea-ice PDF
      real(wp) :: denominator, recip_denominator, areaPDFfac
      real(wp) :: tmpscal0, tmpscal1, tmpscal2, tmpscal3, tmpscal4


!==   local arrays ==
!--   TmixLoc        :: ocean surface/mixed-layer temperature (in K)
      real(wp) ::  TmixLoc       (nlpb)

#ifndef SEAICE_ITD
!     actual ice thickness (with upper and lower limit)
      real(wp) :: heffActual          (nlpb)
!     actual snow thickness
      real(wp) :: hsnowActual         (nlpb)
#endif
!     actual ice thickness (with lower limit only) Reciprocal
      real(wp) :: recip_heffActual    (nlpb)

!     AREA_PRE :: hold sea-ice fraction field before any seaice-thermo update
      real(wp) :: AREApreTH           (nlpb)
      real(wp) :: HEFFpreTH           (nlpb)
      real(wp) :: HSNWpreTH           (nlpb)
#ifdef SEAICE_ITD
      real(wp) :: AREAITDpreTH        (nlpb,1:nITD)
      real(wp) :: HEFFITDpreTH        (nlpb,1:nITD)
      real(wp) :: HSNWITDpreTH        (nlpb,1:nITD)
      real(wp) :: areaFracFactor      (nlpb,1:nITD)
#endif

!     wind speed
      real(wp) :: UG                  (nlpb)


      real(wp) :: ticeInMult          (nlpb,nITD)
      real(wp) :: ticeOutMult         (nlpb,nITD)
      real(wp) :: heffActualMult      (nlpb,nITD)
      real(wp) :: hsnowActualMult     (nlpb,nITD)
#ifdef SEAICE_ITD
      real(wp) :: recip_heffActualMult(nlpb,nITD)
#endif
      real(wp) :: a_QbyATMmult_cover  (nlpb,nITD)
      real(wp) :: a_QSWbyATMmult_cover(nlpb,nITD)
      real(wp) :: a_FWbySublimMult    (nlpb,nITD)
#ifdef SEAICE_ITD
      real(wp) :: r_QbyATMmult_cover  (nlpb,nITD)
      real(wp) :: r_FWbySublimMult    (nlpb,nITD)
! for lateral melt parameterization:
      real(wp) :: latMeltFrac         (nlpb,nITD)
      real(wp) :: latMeltRate         (nlpb,nITD)
      real(wp) :: floeAlpha, floeDiameter, floeDiameterMin, floeDiameterMax
#endif

!     a_QbyATM_cover :: available heat (in W/m^2) due to the interaction of
!             the atmosphere and the ocean surface - for ice covered water
!     a_QbyATM_open  :: same but for open water
!     r_QbyATM_cover :: residual of a_QbyATM_cover after freezing/melting processes
!     r_QbyATM_open  :: same but for open water
      real(wp) ::  a_QbyATM_cover      (nlpb)
      real(wp) ::  a_QbyATM_open       (nlpb)
      real(wp) ::  r_QbyATM_cover      (nlpb)
      real(wp) ::  r_QbyATM_open       (nlpb)
!     a_QSWbyATM_open   - short wave heat flux over ocean in W/m^2
!     a_QSWbyATM_cover  - short wave heat flux under ice in W/m^2
      real(wp) ::  a_QSWbyATM_open     (nlpb)
      real(wp) ::  a_QSWbyATM_cover    (nlpb)
!     a_QbyOCN :: available heat (in W/m^2) due to the
!             interaction of the ice pack and the ocean surface
!     r_QbyOCN :: residual of a_QbyOCN after freezing/melting
!             processes have been accounted for
      real(wp) ::  a_QbyOCN            (nlpb)
      real(wp) ::  r_QbyOCN            (nlpb)

!     The change of mean ice thickness due to turbulent ocean-sea ice heat fluxes
      real(wp) ::  d_HEFFbyOCNonICE    (nlpb)

!     The sum of mean ice thickness increments due to atmospheric fluxes over
!     the open water fraction and ice-covered fractions of the grid cell
      real(wp) ::  d_HEFFbyATMonOCN    (nlpb)
!     The change of mean ice thickness due to flooding by snow
      real(wp) ::  d_HEFFbyFLOODING    (nlpb)

!     The mean ice thickness increments due to atmospheric fluxes over the open
!     water fraction and ice-covered fractions of the grid cell, respectively
      real(wp) ::  d_HEFFbyATMonOCN_open(nlpb)
      real(wp) ::  d_HEFFbyATMonOCN_cover(nlpb)

      real(wp) ::  d_HSNWbyATMonSNW    (nlpb)
      real(wp) ::  d_HSNWbyOCNonSNW    (nlpb)
      real(wp) ::  d_HSNWbyRAIN        (nlpb)

      real(wp) ::  d_HFRWbyRAIN        (nlpb)

!     a_FWbySublim :: fresh water flux implied by latent heat of
!                     sublimation to atmosphere, same sign convention
!                     as EVAP (positive upward)
      real(wp) ::  a_FWbySublim        (nlpb)
      real(wp) ::  r_FWbySublim        (nlpb)
      real(wp) ::  d_HEFFbySublim      (nlpb)
      real(wp) ::  d_HSNWbySublim      (nlpb)

#ifdef SEAICE_CAP_SUBLIM
!     The latent heat flux which will sublimate all snow and ice
!     over one time step
      real(wp) ::  latentHeatFluxMax   (nlpb)
      real(wp) ::  latentHeatFluxMaxMult(nlpb,nITD)
#endif

#ifdef SEAICE_ITD
      real(wp) ::  d_HEFFbySublim_ITD         (nlpb,1:nITD)
      real(wp) ::  d_HSNWbySublim_ITD         (nlpb,1:nITD)
      real(wp) ::  d_HEFFbyOCNonICE_ITD       (nlpb,1:nITD)
      real(wp) ::  d_HSNWbyATMonSNW_ITD       (nlpb,1:nITD)
      real(wp) ::  d_HEFFbyATMonOCN_ITD       (nlpb,1:nITD)
      real(wp) ::  d_HEFFbyATMonOCN_cover_ITD (nlpb,1:nITD)
      real(wp) ::  d_HEFFbyATMonOCN_open_ITD  (nlpb,1:nITD)
      real(wp) ::  d_HSNWbyRAIN_ITD           (nlpb,1:nITD)
      real(wp) ::  d_HSNWbyOCNonSNW_ITD       (nlpb,1:nITD)
      real(wp) ::  d_HEFFbyFLOODING_ITD       (nlpb,1:nITD)
#endif

      real(wp) ::  SItflux     (nlpb)
      real(wp) ::  SIatmQnt    (nlpb)
      real(wp) ::  SIatmFW     (nlpb)

!#ifdef ALLOW_BALANCE_FLUXES
      logical allow_balance_fluxes
      logical balanceQnet, balanceEmPmR
      real(wp) ::  FWFsiTile
      real(wp) ::  HFsiTile
      real(wp) ::  FWF2HFsiTile
!#endif

      real(wp) :: vecFldLoc(nlpb,2), vecIceLoc(nlpb,2)
      integer  :: uel,uei,uep,vnl,vni,vnp


!YY: may delete later. this is for removing global mean of Qnet and freshwater flux
      allow_balance_fluxes = .false.
      balanceQnet = .false.
      balanceEmPmR= .false.

      kSurface = nk
      !YY: donot consider hFacC now
      !dzSurf   = drF(kSurface)*hFacC(i,kSurface)/gravityRau0
      dzSurf   = drF(kSurface)/gravityRau0

      recip_multDim        = ONE / dble(SEAICE_multDim)
      recip_deltaTtherm = ONE / SEAICE_deltaTtherm

!     Cutoff for iceload
      heffTooHeavy = dzSurf * 0.2_wp
!     RATIO OF SEA ICE DENSITY to SNOW DENSITY
      ICE2SNOW     = SEAICE_rhoIce/SEAICE_rhoSnow
      SNOW2ICE     = ONE / ICE2SNOW

!     HEAT OF FUSION OF ICE (J/m^3)
      QI           = SEAICE_rhoIce*SEAICE_lhFusion
      recip_QI     = ONE / QI
!     HEAT OF FUSION OF SNOW (J/m^3)
      QS           = SEAICE_rhoSnow*SEAICE_lhFusion

!     ICE LATENT HEAT CONSTANT: solid-liquid, liquid-vap
      lhSublim = SEAICE_lhEvap + SEAICE_lhFusion

!     regularization constants
      area_reg_sq = SEAICE_area_reg * SEAICE_area_reg
      hice_reg_sq = SEAICE_hice_reg * SEAICE_hice_reg

! conversion factors to go from Q (W/m2) to HEFF (ice meters)
      convertQ2HI=SEAICE_deltaTtherm/QI
      convertHI2Q = ONE/convertQ2HI
! conversion factors to go from precip (m/s) unit to HEFF (ice meters)
      convertPRECIP2HI=SEAICE_deltaTtherm*rhoConstFresh/SEAICE_rhoIce
      convertHI2PRECIP = ONE/convertPRECIP2HI

!     compute parameters for thickness pdf
      denominator = 0.0_wp
      DO IT=1,SEAICE_multDim
       denominator = denominator + IT * SEAICE_pdf(IT)
      ENDDO
      denominator = (2.0_wp * denominator) - 1.0_wp
      recip_denominator = 1.0_wp / denominator
#ifdef SEAICE_ITD
      areaPDFfac  = 1.0_wp
#else
      areaPDFfac  = denominator * recip_multDim
#endif /* SEAICE_ITD */

#ifdef SEAICE_ITD
! constants for lateral melt parameterization:
! following Steele (1992), Equ. 2
      floeAlpha                  = 0.66_wp
! typical mean diameter used in CICE 4.1:
! (this is currently computed as a function of ice concentration
!  following a suggestion by Luepkes at al. (2012))
!      floeDiameter               = 300. _d 0
! parameters needed for variable floe diameter following Luepkes et al. (2012):
      floeDiameterMin            = 8.0_wp
      floeDiameterMax            = 300.0_wp
#endif


!$acc data create(    &
!$acc        a_QbyATM_cover,a_QbyATM_open,r_QbyATM_cover,r_QbyATM_open,           &
!$acc        a_QSWbyATM_open,a_QSWbyATM_cover,a_QbyOCN,r_QbyOCN,                  &
!$acc        d_HEFFbyOCNonICE,d_HEFFbyATMonOCN,d_HEFFbyFLOODING,                  &
!$acc        d_HEFFbyATMonOCN_open,d_HEFFbyATMonOCN_cover,d_HSNWbyATMonSNW,       &
!$acc        d_HSNWbyOCNonSNW,d_HSNWbyRAIN,a_FWbySublim,r_FWbySublim,             &
!$acc        d_HEFFbySublim,d_HSNWbySublim,                                       &  
!$acc        HEFFpreTH,HSNWpreTH,AREApreTH, heffActualMult,hsnowActualMult,       &
!$acc        recip_heffActual,vecFldLoc,vecIceLoc,UG,                             &

!YY: comment for now
!!!$acc        SItflux,SIatmQnt,SIatmFW ,    &

#ifdef SEAICE_CAP_SUBLIM
!$acc        latentHeatFluxMax, latentHeatFluxMaxMult,   & 
#endif
#ifdef SEAICE_ITD
!$acc        d_HEFFbySublim_ITD,d_HSNWbySublim_ITD,d_HEFFbyOCNonICE_ITD,         &
!$acc        d_HSNWbyATMonSNW_ITD,d_HEFFbyATMonOCN_ITD,                          &
!$acc        d_HEFFbyATMonOCN_cover_ITD,d_HEFFbyATMonOCN_open_ITD,               &
!$acc        d_HSNWbyRAIN_ITD,d_HSNWbyOCNonSNW_ITD,                              &
!$acc        d_HEFFbyFLOODING_ITD,r_QbyATMmult_cover,r_FWbySublimMult,           &
!$acc        latMeltFrac,latMeltRate,HEFFITDpreTH,HSNWITDpreTH,AREAITDpreTH,     &
!$acc        areaFracFactor,heffActualMult,hsnowActualMult,recip_heffActualMult, &
#endif
!$acc        d_HFRWbyRAIN,ticeInMult,ticeOutMult,               &
!$acc        a_QbyATMmult_cover,a_QSWbyATMmult_cover,a_FWbySublimMult),  &  
!$acc        present(HEFF,AREA,Qnet,Qsw,qns_ice,qsr_ice,uice,vice,       &
!$acc        u10,v10,saltFlux,zevap,EmPmR_ice),    &
!$acc        copyin(kSurface,ICE2SNOW,SNOW2ICE,recip_QI,dzSurf,     &
!$acc               convertPRECIP2HI,convertHI2PRECIP,recip_deltaTtherm, &
!$acc               convertQ2HI,convertHI2Q,recip_denominator,areaPDFfac, &
!$acc               area_reg_sq,hice_reg_sq)

! ===================================================================
! =======================PART 0: initializations=====================
! ===================================================================



! array initializations
!$acc kernels
!$acc loop
      do i = 1,nlpb
        a_QbyATM_cover (i)       = 0.0_wp
        a_QbyATM_open(i)         = 0.0_wp
        r_QbyATM_cover (i)       = 0.0_wp
        r_QbyATM_open (i)        = 0.0_wp

        a_QSWbyATM_open (i)      = 0.0_wp
        a_QSWbyATM_cover (i)     = 0.0_wp

        a_QbyOCN (i)             = 0.0_wp
        r_QbyOCN (i)             = 0.0_wp

        d_HEFFbyOCNonICE(i)      = 0.0_wp
        d_HEFFbyATMonOCN(i)      = 0.0_wp
        d_HEFFbyFLOODING(i)      = 0.0_wp

        d_HEFFbyATMonOCN_open(i) = 0.0_wp
        d_HEFFbyATMonOCN_cover(i)= 0.0_wp

        d_HSNWbyATMonSNW(i)      = 0.0_wp
        d_HSNWbyOCNonSNW(i)      = 0.0_wp
        d_HSNWbyRAIN(i)          = 0.0_wp
        a_FWbySublim(i)          = 0.0_wp
        r_FWbySublim(i)          = 0.0_wp
        d_HEFFbySublim(i)        = 0.0_wp
        d_HSNWbySublim(i)        = 0.0_wp
#ifdef SEAICE_CAP_SUBLIM
        latentHeatFluxMax(i)     = 0.0_wp
#endif
        d_HFRWbyRAIN(i)          = 0.0_wp
      enddo
!$acc loop collapse(2)
      do IT=1,SEAICE_multDim
        do i = 1,nlpb
          ticeInMult(i,IT)            = 0.0_wp
          ticeOutMult(i,IT)           = 0.0_wp
          a_QbyATMmult_cover(i,IT)    = 0.0_wp
          a_QSWbyATMmult_cover(i,IT)  = 0.0_wp
          a_FWbySublimMult(i,IT)      = 0.0_wp
#ifdef SEAICE_CAP_SUBLIM
          latentHeatFluxMaxMult(i,IT) = 0.0_wp
#endif
#ifdef SEAICE_ITD
          d_HEFFbySublim_ITD(i,IT)         = 0.0_wp
          d_HSNWbySublim_ITD(i,IT)         = 0.0_wp
          d_HEFFbyOCNonICE_ITD(i,IT)       = 0.0_wp
          d_HSNWbyATMonSNW_ITD(i,IT)       = 0.0_wp
          d_HEFFbyATMonOCN_ITD(i,IT)       = 0.0_wp
          d_HEFFbyATMonOCN_cover_ITD(i,IT) = 0.0_wp
          d_HEFFbyATMonOCN_open_ITD(i,IT)  = 0.0_wp
          d_HSNWbyRAIN_ITD(i,IT)           = 0.0_wp
          d_HSNWbyOCNonSNW_ITD(i,IT)       = 0.0_wp
          d_HEFFbyFLOODING_ITD(i,IT)       = 0.0_wp
          r_QbyATMmult_cover(i,IT)         = 0.0_wp
          r_FWbySublimMult(i,IT)           = 0.0_wp
! for lateral melt parameterization:
          latMeltFrac(i,IT)                = 0.0_wp
          latMeltRate(i,IT)                = 0.0_wp
#endif
        enddo
      enddo

! =====================================================================
! ===========PART 1: treat pathological cases (post advdiff)===========
! =====================================================================
!     This part has been mostly moved to S/R seaice_reg_ridge, which is
!     called before S/R seaice_growth

!     store regularized values of heff, hsnow, area at the onset of thermo.
!$acc loop
      do i = 1,nlpb
        HEFFpreTH(I)=HEFF(I)
        HSNWpreTH(I)=HSNOW(I)
        AREApreTH(I)=AREA(I)
      enddo
#ifdef SEAICE_ITD
!$acc loop collapse(2)
      DO IT=1,SEAICE_multDim
        do i = 1,nlpb
          HEFFITDpreTH(I,IT)=HEFFITD(I,IT)
          HSNWITDpreTH(I,IT)=HSNOWITD(I,IT)
          AREAITDpreTH(I,IT)=AREAITD(I,IT)

!     keep track of areal and volume fraction of each ITD category
          IF (AREA(I) .GT. ZERO) THEN
            areaFracFactor(I,IT)=AREAITD(I,IT)/AREA(I)
          ELSE
!          if there is no ice, potential growth starts in 1st category
            IF (IT .EQ. 1) THEN
              areaFracFactor(I,IT)=ONE
            ELSE
              areaFracFactor(I,IT)=ZERO
            ENDIF
          ENDIF
        enddo
      ENDDO
#endif /* SEAICE_ITD */
!$acc end kernels

!     COMPUTE ACTUAL ICE/SNOW THICKNESS; USE MIN/MAX VALUES
!     TO REGULARIZE SEAICE_SOLVE4TEMP/d_AREA COMPUTATIONS

!$acc kernels 
#ifdef SEAICE_ITD
!$acc loop collapse(2) private(tmpscal1,tmpscal2)
        DO IT=1,SEAICE_multDim
          do i = 1,nlpb
            IF (HEFFITDpreTH(I,IT) .GT. ZERO) THEN
!     regularize AREA with SEAICE_area_reg
              tmpscal1 = SQRT(AREAITDpreTH(I,IT) * AREAITDpreTH(I,IT) + area_reg_sq)
!     heffActual calculated with the regularized AREA
              tmpscal2 = HEFFITDpreTH(I,IT) / tmpscal1
!     regularize heffActual with SEAICE_hice_reg (add lower bound)
              heffActualMult(I,IT) = SQRT(tmpscal2 * tmpscal2 + hice_reg_sq)
!     hsnowActual calculated with the regularized AREA
              hsnowActualMult(I,IT) = HSNWITDpreTH(I,IT) / tmpscal1
!     regularize the inverse of heffActual by hice_reg
              recip_heffActualMult(I,IT)  = AREAITDpreTH(I,IT) /     &
                sqrt(HEFFITDpreTH(I,IT) * HEFFITDpreTH(I,IT)         &
                + hice_reg_sq)
!     Do not regularize when HEFFpreTH = 0
            ELSE
              heffActualMult(I,IT) = ZERO
              hsnowActualMult(I,IT) = ZERO
              recip_heffActualMult(I,IT)  = ZERO
            ENDIF
          enddo
        ENDDO
#else /* ndef SEAICE_ITD */
!$acc loop private(tmpscal1,tmpscal2)
        do i = 1,nlpb
          IF (HEFFpreTH(I) .GT. ZERO) THEN
!if        regularize AREA with SEAICE_area_reg
           tmpscal1 = SQRT(AREApreTH(I)* AREApreTH(I) + area_reg_sq)
!if        heffActual calculated with the regularized AREA
           tmpscal2 = HEFFpreTH(I) / tmpscal1
!if        regularize heffActual with SEAICE_hice_reg (add lower bound)
           heffActual(I) = SQRT(tmpscal2 * tmpscal2 + hice_reg_sq)
!if        hsnowActual calculated with the regularized AREA
           hsnowActual(I) = HSNWpreTH(I) / tmpscal1
!if        regularize the inverse of heffActual by hice_reg
           recip_heffActual(I)  = AREApreTH(I) /                    &
                      sqrt(HEFFpreTH(I)*HEFFpreTH(I) + hice_reg_sq)
!if       Do not regularize when HEFFpreTH = 0
          ELSE
           heffActual(I) = ZERO
           hsnowActual(I) = ZERO
           recip_heffActual(I)  = ZERO
          ENDIF
        enddo
#endif /* SEAICE_ITD */
!$acc end kernels

#ifdef SEAICE_CAP_SUBLIM
!     COMPUTE MAXIMUM LATENT HEAT FLUXES FOR THE CURRENT ICE
!     AND SNOW THICKNESS
!     The latent heat flux over the sea ice which
!     will sublimate all of the snow and ice over one time
!     step (W/m^2)
!$acc kernels copyin(lhSublim)
#ifdef SEAICE_ITD
!$acc loop collapse(2)
        DO IT=1,SEAICE_multDim
          do i = 1,nlpb
            IF (HEFFITDpreTH(I,IT) .GT. ZERO) THEN
             latentHeatFluxMaxMult(I,IT) = lhSublim*recip_deltaTtherm *  &
                 (HEFFITDpreTH(I,IT)*SEAICE_rhoIce +                     &
                  HSNWITDpreTH(I,IT)*SEAICE_rhoSnow)                     &
                 /AREAITDpreTH(I,IT)
            ELSE
             latentHeatFluxMaxMult(I,IT) = ZERO
            ENDIF
          enddo
        ENDDO
#else /* ndef SEAICE_ITD */
!$acc loop
        do i = 1,nlpb
          IF (HEFFpreTH(I) .GT. ZERO) THEN
           latentHeatFluxMax(I) = lhSublim * recip_deltaTtherm *         &
               (HEFFpreTH(I) * SEAICE_rhoIce +                           &
                HSNWpreTH(I) * SEAICE_rhoSnow)/AREApreTH(I)
          ELSE
           latentHeatFluxMax(I) = ZERO
          ENDIF
        enddo
#endif /* SEAICE_ITD */
!$acc end kernels
#endif /* SEAICE_CAP_SUBLIM */

! ===================================================================
! ================PART 2: determine heat fluxes/stocks===============
! ===================================================================

! determine available heat due to the atmosphere -- for open water
! ================================================================
!YY: the compuation of qns/qsr refer to mod_csp_forcing
!$acc kernels
!$acc loop
        do i = 1,nlpb
          !a_QbyATM_open  (i) = -qns(i) -qsr(i)
          !a_QSWbyATM_open(i) = -qsr(i)
          a_QbyATM_open  (i) = Qnet(i)
          a_QSWbyATM_open(i) = Qsw (i)
        enddo

! determine available heat due to the atmosphere -- for ice covered water
! =======================================================================
!     Compute relative wind speed over sea ice.
!$acc loop
        do i = 1,nlpb
          vecFldLoc(i, iu) = u10(i,4)
          vecFldLoc(i, iv) = v10(i,4)
          vecIceLoc(i, iu) = uice(i)
          vecIceLoc(i, iv) = vice(i)
        end do
!$acc loop private(uel,uei,uep,vnl,vni,vnp,SPEED_SQ)
        do i = 1,nlpb
          uel = ue(i, 1)
          uei = ue(i, 2)
          uep = ue(i, 3)
          vnl = vn(i, 1)
          vni = vn(i, 2)
          vnp = vn(i, 3)

          SPEED_SQ=(0.5_wp*(vecFldLoc(i,iu)-vecIceLoc(i,iu)+vecFldLoc(uel,uei)*uep-vecIceLoc(uel,uei)*uep))**2  &
                  +(0.5_wp*(vecFldLoc(i,iv)-vecIceLoc(i,iv)+vecFldLoc(vnl,vni)*vnp-vecIceLoc(vnl,vni)*vnp))**2
          IF ( SPEED_SQ .LE. SEAICE_EPS_SQ ) THEN
            UG(i)=SEAICE_EPS
          ELSE
            UG(i)=SQRT(SPEED_SQ)
          ENDIF
        enddo

!--   Start loop over multi-categories
!$acc loop collapse(2)
!!$acc loop collapse(2) private(pFac,pFacSnow)
        DO IT=1,SEAICE_multDim
          do i = 1,nlpb
            ticeInMult(I,IT)  = TICES(I,IT)
            ticeOutMult(I,IT) = TICES(I,IT)
            TICES(I,IT) = ZERO
#ifndef SEAICE_ITD
!     for SEAICE_ITD heffActualMult and latentHeatFluxMaxMult have been
!     calculated above (instead of heffActual and latentHeatFluxMax)
!--   assume homogeneous distribution between 0 and 2 x heffActual
            pFac = (2.0_wp*IT - 1.0_wp)*recip_denominator
            pFacSnow = 1.0_wp
            IF ( SEAICE_useMultDimSnow ) pFacSnow=pFac
            heffActualMult(I,IT)        = heffActual(I)*pFac
            hsnowActualMult(I,IT)       = hsnowActual(I)*pFacSnow
#ifdef SEAICE_CAP_SUBLIM
            latentHeatFluxMaxMult(I,IT) = latentHeatFluxMax(I)*pFac
#endif
#endif /* ndef SEAICE_ITD */
         enddo
        ENDDO
!$acc end kernels

        DO IT=1,SEAICE_multDim
         call seaice_solve4temp(UG, heffActualMult(:,IT), hsnowActualMult(:,IT), &
#ifdef SEAICE_CAP_SUBLIM
              latentHeatFluxMaxMult(:,IT),                                      &
#endif
              ticeInMult(:,IT),ticeOutMult(:,IT), a_QbyATMmult_cover(:,IT),      &
              a_QSWbyATMmult_cover(:,IT),a_FWbySublimMult(:,IT) )
        ENDDO

!$acc kernels
!$acc loop
        do i = 1,nlpb
!$acc loop seq
          DO IT=1,SEAICE_multDim
!     update TICES
            TICES(I,IT) = ticeOutMult(I,IT)
!     average over categories
#ifdef SEAICE_ITD
!     calculate area weighted mean
!     (fluxes are per unit (ice surface) area and are thus area weighted)
            a_QbyATM_cover   (I) = a_QbyATM_cover(I)              &
                + a_QbyATMmult_cover(I,IT)*areaFracFactor(I,IT)
            a_QSWbyATM_cover (I) = a_QSWbyATM_cover(I)            &
                + a_QSWbyATMmult_cover(I,IT)*areaFracFactor(I,IT)
            a_FWbySublim     (I) = a_FWbySublim(I)                &
                + a_FWbySublimMult(I,IT)*areaFracFactor(I,IT)
#else
            a_QbyATM_cover   (I) = a_QbyATM_cover(I)              &
                + a_QbyATMmult_cover(I,IT)*SEAICE_PDF(IT)
            a_QSWbyATM_cover (I) = a_QSWbyATM_cover(I)            &
                + a_QSWbyATMmult_cover(I,IT)*SEAICE_PDF(IT)
            a_FWbySublim     (I) = a_FWbySublim(I)                &
                + a_FWbySublimMult(I,IT)*SEAICE_PDF(IT)
#endif
          enddo
        ENDDO

! switch heat fluxes from W/m2 to 'effective' ice meters
#ifdef SEAICE_ITD
!$acc loop collapse(2)
        DO IT=1,SEAICE_multDim
          do i = 1,nlpb
           a_QbyATMmult_cover(I,IT) = a_QbyATMmult_cover(I,IT)     &
               * convertQ2HI * AREAITDpreTH(I,IT)   !YY: convert to grid-average,assume area as 1
           a_QSWbyATMmult_cover(I,IT) = a_QSWbyATMmult_cover(I,IT) &
               * convertQ2HI * AREAITDpreTH(I,IT)
! and initialize r_QbyATMmult_cover
           r_QbyATMmult_cover(I,IT)=a_QbyATMmult_cover(I,IT)
!     Convert fresh water flux by sublimation to 'effective' ice meters.
!     Negative sublimation is resublimation and will be added as snow.
           a_FWbySublimMult(I,IT) = SEAICE_deltaTtherm*(1.0_wp/SEAICE_rhoIce) &
                 * a_FWbySublimMult(I,IT)*AREAITDpreTH(I,IT)
           r_FWbySublimMult(I,IT)=a_FWbySublimMult(I,IT)
          enddo
        ENDDO
!$acc loop
        do i = 1,nlpb
          a_QbyATM_open(I)    = a_QbyATM_open(I) * convertQ2HI * ( ONE - AREApreTH(I) )
          a_QSWbyATM_open(I)  = a_QSWbyATM_open(I) * convertQ2HI * ( ONE - AREApreTH(I) )
! and initialize r_QbyATM_open
          r_QbyATM_open(I)=a_QbyATM_open(I)
        enddo
#else /* SEAICE_ITD */
!$acc loop
        do i = 1,nlpb
          a_QbyATM_cover(I) = a_QbyATM_cover(I) * convertQ2HI * AREApreTH(I)
          a_QSWbyATM_cover(I) = a_QSWbyATM_cover(I) * convertQ2HI * AREApreTH(I)
          a_QbyATM_open(I)    = a_QbyATM_open(I) * convertQ2HI * ( ONE - AREApreTH(I) )
          a_QSWbyATM_open(I)  = a_QSWbyATM_open(I) * convertQ2HI * ( ONE - AREApreTH(I) )
! and initialize r_QbyATM_cover/r_QbyATM_open
          r_QbyATM_cover(I)=a_QbyATM_cover(I)
          r_QbyATM_open(I)=a_QbyATM_open(I)
!     Convert fresh water flux by sublimation to 'effective' ice meters.
!     Negative sublimation is resublimation and will be added as snow.
          a_FWbySublim(I) = SEAICE_deltaTtherm*(1.0_wp/SEAICE_rhoIce)    &
                * a_FWbySublim(I)*AREApreTH(I)
          r_FWbySublim(I)=a_FWbySublim(I)
        enddo
#endif /* SEAICE_ITD */

! determine available heat due to the ice pack tying the
! underlying surface water temperature to freezing point
! ======================================================

!$acc loop private(tempFrz,tmpscal1,MixedLayerTurbulenceFactor,tmpscal2)
        do i = 1,nlpb
!         FREEZING TEMP. OF SEA WATER (deg C)
          tempFrz = SEAICE_tempFrz0 + SEAICE_dTempFrz_dS *tFld(I,kSurface,2)
! efficiency of turbulent fluxes : dependency to sign of THETA-TBC
          IF ( tFld(I,kSurface,1) .GE. tempFrz ) THEN
           tmpscal1 = SEAICE_mcPheePiston
          ELSE
           tmpscal1 =SEAICE_frazilFrac*dzSurf/SEAICE_deltaTtherm
          ENDIF
! efficiency of turbulent fluxes : dependency to AREA (McPhee cases)
          IF ( (AREApreTH(I) .GT. 0.0_wp).AND.(.NOT.SEAICE_mcPheeStepFunc) ) THEN
           MixedLayerTurbulenceFactor = ONE - SEAICE_mcPheeTaper * AREApreTH(I)
          ELSEIF ( (AREApreTH(I) .GT. 0.0_wp).AND.(SEAICE_mcPheeStepFunc) ) THEN
           MixedLayerTurbulenceFactor = ONE - SEAICE_mcPheeTaper
          ELSE
           MixedLayerTurbulenceFactor = ONE
          ENDIF
! maximum turbulent flux, in ice meters
!YY: note original HeatCapacity_Cp is equal to rcp in MaCOM
          !tmpscal2= - (rcp*rhoConst * recip_QI)
          tmpscal2= - (rcp*rau0 * recip_QI)           &
              * (tFld(I,kSurface,1)-tempFrz)          &
              * SEAICE_deltaTtherm * HEFFM(i)
! available turbulent flux
          a_QbyOCN(i) = tmpscal1 * tmpscal2 * MixedLayerTurbulenceFactor
          r_QbyOCN(i) = a_QbyOCN(i)
        enddo
!$acc end kernels

#ifdef SEAICE_ITD
! determine lateral melt rate at floe edges based on an
! average floe diameter or a floe size distribution
! following Steele (1992, Tab. 2)
! ======================================================
!$acc kernels copyin(floeAlpha,floeDiameterMin,floeDiameterMax)
!$acc loop collapse(2) private(floeDiameter,tempFrz,tmpscal1,tmpscal2)
        DO IT=1,SEAICE_multDim
          do i = 1,nlpb
           tempFrz = SEAICE_tempFrz0 +                 &
               SEAICE_dTempFrz_dS *tFld(I,kSurface,2)
           tmpscal1=(tFld(I,kSurface,1)-tempFrz)
           tmpscal2=sqrt(0.87_wp + 0.067_wp*UG(i)) * UG(i)

!     variable floe diameter following Luepkes et al. (2012, JGR, Equ. 26)
!     with beta=1
           floeDiameter = floeDiameterMin * floeDiameterMax    &
               / ( floeDiameterMax*( 1.0_wp - AREApreTH(I) )   &
                 + floeDiameterMin*AREApreTH(I) )
!     following the idea of SEAICE_areaLossFormula == 1:
           IF (a_QbyATMmult_cover(i,it).LT.ZERO .OR.           &
              a_QbyATM_open(i) .LT.ZERO .OR.                   &
              a_QbyOCN(i)      .LT.ZERO) THEN
!     lateral melt rate as suggested by Perovich, 1983 (PhD thesis)
!           latMeltRate(i,j,it) = 1.6 _d -6 * tmpscal1**1.36
!     The following for does the same, but is faster
            latMeltRate(i,it) = ZERO
            IF (tmpscal1 .GT. ZERO) &
               latMeltRate(i,it) = 1.6d-6 * exp(1.36_wp*log(tmpscal1))
!     lateral melt rate as suggested by Maykut and Perovich, 1987
!     (JGR 92(C7)), Equ. 24
!           latMeltRate(i,j,it) = 13.5 _d -6 * tmpscal2 * tmpscal1**1.3
!     further suggestion by Maykut and Perovich to avoid
!     latMeltRate -> 0 for UG -> 0
!           latMeltRate(i,j,it) = (1.6 _d -6 + 13.5 _d -6 * tmpscal2)
!    &                          * tmpscal1**1.3
!     factor determining fraction of area and ice volume reduction
!     due to lateral melt
            latMeltFrac(i,it) =                           &
            latMeltRate(i,it)*SEAICE_deltaTtherm*rpi /    &
            (floeAlpha * floeDiameter)
            latMeltFrac(i,it)=max(ZERO, min(latMeltFrac(i,it),ONE))
           ELSE
            latMeltRate(i,it)=0.0_wp
            latMeltFrac(i,it)=0.0_wp
           ENDIF
          enddo
        ENDDO
!$acc end kernels
#endif /* SEAICE_ITD */

! ===================================================================
! =========PART 3: determine effective thicknesses increments========
! ===================================================================

! compute snow/ice tendency due to sublimation
#ifdef SEAICE_ITD
!$acc kernels
!$acc loop private(tmpscal2)
        do i = 1,nlpb
!$acc loop seq
         DO IT=1,SEAICE_multDim
!     First sublimate/deposite snow
          tmpscal2 =                                      &
          MAX(MIN(r_FWbySublimMult(I,IT),HSNOWITD(I,IT) *SNOW2ICE),ZERO)
          d_HSNWbySublim_ITD(I,IT) = - tmpscal2 * ICE2SNOW
!         accumulate change over ITD categories
          d_HSNWbySublim(I)     = d_HSNWbySublim(I) - tmpscal2*ICE2SNOW
          r_FWbySublimMult(I,IT)= r_FWbySublimMult(I,IT) - tmpscal2

!     If anything is left, sublimate ice
          tmpscal2 =                                     &
          MAX(MIN(r_FWbySublimMult(I,IT),HEFFITD(I,IT)),ZERO)
          d_HEFFbySublim_ITD(I,IT) = - tmpscal2
!         accumulate change over ITD categories
          d_HEFFbySublim(I)      = d_HEFFbySublim(I)      - tmpscal2
          r_FWbySublimMult(I,IT) = r_FWbySublimMult(I,IT) - tmpscal2
         enddo
        enddo
!$acc loop collapse(2)
        do IT=1,SEAICE_multDim
         do i = 1,nlpb
!     If anything is left, it will be evaporated from the ocean rather than sublimated.
!     Since a_QbyATM_cover was computed for sublimation, not simple evaporation, we need to
!     remove the fusion part for the residual (that happens to be precisely r_FWbySublim).
          a_QbyATMmult_cover(I,IT) = a_QbyATMmult_cover(I,IT)  &
                                     - r_FWbySublimMult(I,IT)
          r_QbyATMmult_cover(I,IT) = r_QbyATMmult_cover(I,IT)  &
                                     - r_FWbySublimMult(I,IT)
         enddo
        enddo     !endIT loop
!$acc end kernels

! else if for SEAICE_ITD
#else

!$acc kernels
!$acc loop private(tmpscal2)
         do i = 1,nlpb
!     First sublimate/deposite snow
          tmpscal2 =                                      &
          MAX(MIN(r_FWbySublim(I),HSNOW(I)*SNOW2ICE),ZERO)
          d_HSNWbySublim(I) = - tmpscal2 * ICE2SNOW
          HSNOW(I)    = HSNOW(I)  - tmpscal2*ICE2SNOW
          r_FWbySublim(I)   = r_FWbySublim(I) - tmpscal2
         enddo
!$acc loop private(tmpscal2)
         do i = 1,nlpb
!     If anything is left, sublimate ice
          tmpscal2 = MAX(MIN(r_FWbySublim(I),HEFF(I)),ZERO)
          d_HEFFbySublim(I) = - tmpscal2
          HEFF(I)     = HEFF(I)   - tmpscal2
          r_FWbySublim(I)   = r_FWbySublim(I) - tmpscal2
         enddo
!$acc loop
         do i = 1,nlpb
!     If anything is left, it will be evaporated from the ocean rather than sublimated.
!     Since a_QbyATM_cover was computed for sublimation, not simple evaporation, we need to
!     remove the fusion part for the residual (that happens to be precisely r_FWbySublim).
          a_QbyATM_cover(I) = a_QbyATM_cover(I)-r_FWbySublim(I)
          r_QbyATM_cover(I) = r_QbyATM_cover(I)-r_FWbySublim(I)
         enddo
!$acc end kernels

#endif

! compute ice thickness tendency due to ice-ocean interaction

       IF (.NOT.SEAICE_growMeltByConv) THEN

#ifdef SEAICE_ITD
!$acc kernels 
!$acc loop private(tmpscal1)
        do i = 1,nlpb
!$acc loop seq
          do IT=1,SEAICE_multDim
!          ice growth/melt due to ocean heat r_QbyOCN (W/m^2) is
!          equally distributed under the ice and hence weighted by
!          fractional area of each thickness category
           tmpscal1=MAX(r_QbyOCN(i)*areaFracFactor(I,IT), -HEFFITD(I,IT))
           d_HEFFbyOCNonICE_ITD(I,IT)=tmpscal1
           d_HEFFbyOCNonICE(I) = d_HEFFbyOCNonICE(I) + tmpscal1
          enddo
        enddo
!$acc loop
        do i = 1,nlpb
          r_QbyOCN(I)=r_QbyOCN(I)-d_HEFFbyOCNonICE(I)
        enddo
!$acc end kernels
#else /* SEAICE_ITD */
!$acc kernels loop
        do i = 1,nlpb
          d_HEFFbyOCNonICE(I)=MAX(r_QbyOCN(i), -HEFF(I))
          r_QbyOCN(I)=r_QbyOCN(I)-d_HEFFbyOCNonICE(I)
          HEFF(I)=HEFF(I) + d_HEFFbyOCNonICE(I)
        enddo
!$acc end kernels
#endif /* SEAICE_ITD */

      ENDIF !SEAICE_growMeltByConv

! compute snow melt tendency due to snow-atmosphere interaction

#ifdef SEAICE_ITD
!$acc kernels 
!$acc loop private(tmpscal1,tmpscal2)
        do i = 1,nlpb
!$acc loop seq
          do IT=1,SEAICE_multDim
!     Convert to standard units (meters of ice) rather than to meters
!     of snow. This appears to be more robust.
           tmpscal1=MAX(r_QbyATMmult_cover(I,IT), -HSNOWITD(I,IT)*SNOW2ICE)
           tmpscal2=MIN(tmpscal1,0.0_wp)
           d_HSNWbyATMonSNW_ITD(I,IT) = tmpscal2*ICE2SNOW
           d_HSNWbyATMonSNW(I) = d_HSNWbyATMonSNW(I) + tmpscal2*ICE2SNOW
           r_QbyATMmult_cover(I,IT)=r_QbyATMmult_cover(I,IT) - tmpscal2
          enddo
        enddo
!$acc end kernels
#else /* SEAICE_ITD */
!$acc kernels 
!$acc loop private(tmpscal1,tmpscal2)
        do i = 1,nlpb
!     Convert to standard units (meters of ice) rather than to meters
!     of snow. This appears to be more robust.
          tmpscal1=MAX(r_QbyATM_cover(I),-HSNOW(I)*SNOW2ICE)
          tmpscal2=MIN(tmpscal1,0.0_wp)
          d_HSNWbyATMonSNW(I)= tmpscal2*ICE2SNOW
          !HSNOW(I) = HSNOW(I) + tmpscal2*ICE2SNOW
          !!IMPORTANT,YUAN: avoid extremely small negative HSNOW values by snow2ice/ice2snow conv.
          HSNOW(I) = MAX(0.0_wp,HSNOW(I) + tmpscal2*ICE2SNOW)
          r_QbyATM_cover(I)=r_QbyATM_cover(I) - tmpscal2
        enddo
!$acc end kernels
#endif /* SEAICE_ITD */

! compute ice thickness tendency due to the atmosphere

#ifdef SEAICE_ITD
!$acc kernels 
!$acc loop private(tmpscal1,tmpscal2)
        do i = 1,nlpb
!$acc loop seq
          DO IT=1,SEAICE_multDim
           tmpscal1 = HEFFITDpreTH(I,IT)          &
               + d_HEFFbySublim_ITD(I,IT) + d_HEFFbyOCNonICE_ITD(I,IT)
           tmpscal2 = MAX(-tmpscal1, r_QbyATMmult_cover(I,IT)   &
!         Limit ice growth by potential melt by ocean
                   + AREAITDpreTH(I,IT) * r_QbyOCN(I))
           d_HEFFbyATMonOCN_cover_ITD(I,IT) = tmpscal2
           d_HEFFbyATMonOCN_cover(I) = d_HEFFbyATMonOCN_cover(I) + tmpscal2
           d_HEFFbyATMonOCN_ITD(I,IT) = d_HEFFbyATMonOCN_ITD(I,IT) + tmpscal2
           d_HEFFbyATMonOCN(I)       = d_HEFFbyATMonOCN(I) + tmpscal2
           r_QbyATMmult_cover(I,IT)  = r_QbyATMmult_cover(I,IT) - tmpscal2
          enddo
        ENDDO
!$acc end kernels
#else /* ndef SEAICE_ITD */
!$acc kernels 
!$acc loop private(tmpscal2)
        do i = 1,nlpb
          tmpscal2 = MAX(-HEFF(I),r_QbyATM_cover(I)+        &
!         Limit ice growth by potential melt by ocean
             AREApreTH(I) * r_QbyOCN(I))

          d_HEFFbyATMonOCN_cover(I)=tmpscal2
          d_HEFFbyATMonOCN(I)=d_HEFFbyATMonOCN(I)+tmpscal2
          r_QbyATM_cover(I)=r_QbyATM_cover(I)-tmpscal2
          HEFF(I) = HEFF(I) + tmpscal2
        enddo
!$acc end kernels
#endif /* SEAICE_ITD */

! add snow precipitation to HSNOW.
!$acc kernels loop
        do i = 1,nlpb
          d_HSNWbyRAIN(I) = convertPRECIP2HI * ICE2SNOW *      &
              snow(i,4) /rhoConstFresh * AREApreTH(I) 
          d_HFRWbyRAIN(I) = -convertPRECIP2HI *                &
              (prec(i, 4) - snow(i, 4))/rhoConstFresh * AREApreTH(I)
          HSNOW(I) = HSNOW(I) + d_HSNWbyRAIN(I)
        enddo
!$acc end kernels

! compute snow melt due to heat available from ocean.

      IF (.NOT.SEAICE_growMeltByConv) THEN

#ifdef SEAICE_ITD
!$acc kernels 
!$acc loop private(tmpscal1,tmpscal2,tmpscal4)
        do i = 1,nlpb
!$acc loop seq
          DO IT=1,SEAICE_multDim
           tmpscal4 = HSNWITDpreTH(I,IT)            &
                   + d_HSNWbySublim_ITD(I,IT)       &
                   + d_HSNWbyATMonSNW_ITD(I,IT)     &
                   + d_HSNWbyRAIN_ITD(I,IT)
           tmpscal1=MAX(r_QbyOCN(i)*ICE2SNOW*areaFracFactor(I,IT), -tmpscal4)
           tmpscal2=MIN(tmpscal1,0.0_wp)
           d_HSNWbyOCNonSNW_ITD(I,IT) = tmpscal2
           d_HSNWbyOCNonSNW(I) = d_HSNWbyOCNonSNW(I) + tmpscal2
           r_QbyOCN(I)=r_QbyOCN(I) - tmpscal2*SNOW2ICE
          enddo
        ENDDO
!$acc end kernels
#else /* ndef SEAICE_ITD */
!$acc kernels 
!$acc loop private(tmpscal1,tmpscal2)
        do i = 1,nlpb
          tmpscal1=MAX(r_QbyOCN(i)*ICE2SNOW, -HSNOW(I))
          tmpscal2=MIN(tmpscal1,0.0_wp)
          d_HSNWbyOCNonSNW(I) = tmpscal2
          r_QbyOCN(I)=r_QbyOCN(I)-d_HSNWbyOCNonSNW(I)*SNOW2ICE
          HSNOW(I) = HSNOW(I)+d_HSNWbyOCNonSNW(I)
        enddo
!$acc end kernels
#endif /* SEAICE_ITD */

      ENDIF !SEAICE_growMeltByConv

! gain of new ice over open water
!$acc kernels loop private(tmpscal4,tmpscal1,tmpscal2,tmpscal3)
      do i = 1,nlpb
#ifdef SEAICE_ITD
!       HEFF will be updated at the end of PART 3,
!       hence sum of tendencies so far is needed
        tmpscal4 = HEFFpreTH(I)           &
                + d_HEFFbySublim(I)       &
                + d_HEFFbyOCNonICE(I)     &
                + d_HEFFbyATMonOCN(I)
#else /* ndef SEAICE_ITD */
!      HEFF is updated step by step throughout seaice_growth
       tmpscal4 = HEFF(I)
#endif /* SEAICE_ITD */
!      Initial ice growth is triggered by open water
!      heat flux overcoming potential melt by ocean
       tmpscal1=r_QbyATM_open(I)+r_QbyOCN(i) * (1.0_wp - AREApreTH(I))
!      Penetrative shortwave flux beyond first layer
!      that is therefore not available to ice growth/melt
       tmpscal2=SWFracB * a_QSWbyATM_open(I)
!      impose -HEFF as the maxmum melting if SEAICE_doOpenWaterMelt
!      or 0. otherwise (no melting if not SEAICE_doOpenWaterMelt)
       tmpscal3=facOpenGrow*MAX(tmpscal1-tmpscal2,   &
               -tmpscal4*facOpenMelt)*HEFFM(I)
#ifdef SEAICE_ITD
!      ice growth in open water adds to first category
       d_HEFFbyATMonOCN_open_ITD(I,1)=tmpscal3
       d_HEFFbyATMonOCN_ITD(I,1)=d_HEFFbyATMonOCN_ITD(I,1)+tmpscal3
#endif /* SEAICE_ITD */
       d_HEFFbyATMonOCN_open(I)=tmpscal3
       d_HEFFbyATMonOCN(I)=d_HEFFbyATMonOCN(I)+tmpscal3
       r_QbyATM_open(I)=r_QbyATM_open(I)-tmpscal3
       HEFF(I) = HEFF(I) + tmpscal3
      enddo
!$acc end kernels


! convert snow to ice if submerged.

      IF ( SEAICEuseFlooding ) THEN
#ifdef SEAICE_ITD
!$acc kernels
!$acc loop collapse(2) private(tmpscal0,tmpscal1,tmpscal3,tmpscal4)
        DO IT=1,SEAICE_multDim
          do i = 1,nlpb
            tmpscal3 = HEFFITDpreTH(I,IT)           &
                    + d_HEFFbySublim_ITD(I,IT)      &
                    + d_HEFFbyOCNonICE_ITD(I,IT)    &
                    + d_HEFFbyATMonOCN_ITD(I,IT)
            tmpscal4 = HSNWITDpreTH(I,IT)           &
                    + d_HSNWbySublim_ITD(I,IT)      &
                    + d_HSNWbyATMonSNW_ITD(I,IT)    &
                    + d_HSNWbyRAIN_ITD(I,IT)
            tmpscal0 = (tmpscal4*SEAICE_rhoSnow     &
                    +  tmpscal3*SEAICE_rhoIce)      &
                    / rhoConstFresh
            tmpscal1 = MAX( 0.0_wp, tmpscal0 - tmpscal3)
            d_HEFFbyFLOODING_ITD(I,IT) = tmpscal1
          enddo
        ENDDO
!$acc loop
        do i = 1,nlpb
!$acc loop seq
          DO IT=1,SEAICE_multDim
            d_HEFFbyFLOODING(I) = d_HEFFbyFLOODING(I)  + d_HEFFbyFLOODING_ITD(I,IT)
          enddo
        enddo
!$acc end kernels
#else /* ndef SEAICE_ITD */
!$acc kernels
!$acc loop private(tmpscal0,tmpscal1)
        do i = 1,nlpb
          tmpscal0 = (HSNOW(I)*SEAICE_rhoSnow+HEFF(I)*SEAICE_rhoIce) &
              / rhoConstFresh  !YY:divided by rhofresh?
          tmpscal1 = MAX( 0.0_wp, tmpscal0 - HEFF(I))
          d_HEFFbyFLOODING(I)=tmpscal1
          HEFF(I) = HEFF(I)+d_HEFFbyFLOODING(I)
          !YY: avoid negative hsnow values
          HSNOW(I) = MAX(0.0_wp, HSNOW(I)- d_HEFFbyFLOODING(I)*ICE2SNOW)
        enddo
!$acc end kernels
#endif /* SEAICE_ITD */
      ENDIF

#ifdef SEAICE_ITD
! apply ice and snow thickness changes

!$acc kernels loop collapse(2)
      DO IT=1,SEAICE_multDim
        do i = 1,nlpb
         HEFFITD(I,IT) = HEFFITD(I,IT)                      &
                              + d_HEFFbySublim_ITD(I,IT)    &
                              + d_HEFFbyOCNonICE_ITD(I,IT)  &
                              + d_HEFFbyATMonOCN_ITD(I,IT)  &
                              + d_HEFFbyFLOODING_ITD(I,IT)
        !avoid negative hsnowitd values. here my understanding is itd redistribution only take effect in HEFFITD
        !, rather than HSNOWITD, so HSNOWITD dont permit negativel values
        ! HSNOWITD(I,IT) = HSNOWITD(I,IT)                    &
        !                      + d_HSNWbySublim_ITD(I,IT)    &
        !                      + d_HSNWbyATMonSNW_ITD(I,IT)  &
        !                      + d_HSNWbyRAIN_ITD(I,IT)      &
        !                      + d_HSNWbyOCNonSNW_ITD(I,IT)  &
        !                      - d_HEFFbyFLOODING_ITD(I,IT) * ICE2SNOW
         HSNOWITD(I,IT) = MAX( 0.0_wp, HSNOWITD(I,IT)                    &
                              + d_HSNWbySublim_ITD(I,IT)    &
                              + d_HSNWbyATMonSNW_ITD(I,IT)  &
                              + d_HSNWbyRAIN_ITD(I,IT)      &
                              + d_HSNWbyOCNonSNW_ITD(I,IT)  &
                              - d_HEFFbyFLOODING_ITD(I,IT) * ICE2SNOW )
        enddo
      ENDDO
!$acc end kernels
#endif /* SEAICE_ITD */

! ===================================================================
! ==========PART 4: determine ice cover fraction increments=========-
! ===================================================================

!$acc kernels
#ifdef SEAICE_ITD
!--   in thinnest category account for lateral ice growth and melt the
!--   "non-ITD" way, so that the ITD simulation with SEAICE_multDim=1
!--   is identical with the non-ITD simulation;
!--   use HEFF, ARE, HSNOW, etc. as temporal storage for 1st category
!$acc loop
      do i = 1,nlpb
        HEFF(I)=HEFFITD(I,1)
        AREA(I)=AREAITD(I,1)
        HSNOW(I)=HSNOWITD(I,1)
        HEFFpreTH(I)=HEFFITDpreTH(I,1)
        AREApreTH(I)=AREAITDpreTH(I,1)
        recip_heffActual(I)=recip_heffActualMult(I,1)
      enddo
#endif /* SEAICE_ITD */
!$acc loop private(recip_HO,recip_HH,tmpscal4,tmpscal3,tmpscal0,tmpscal1,tmpscal2)
      do i = 1,nlpb

        IF ( latC(I) .LT. ZERO ) THEN
         recip_HO=1.0_wp / HO_south
        ELSE
         recip_HO=1.0_wp / HO
        ENDIF

        recip_HH = recip_heffActual(I)

! gain of ice over open water : computed from
!   (SEAICE_areaGainFormula.EQ.1) from growth by ATM
!   (SEAICE_areaGainFormula.EQ.2) from predicted growth by ATM
        IF (SEAICE_areaGainFormula.EQ.1) THEN
          tmpscal4 = MAX(ZERO,d_HEFFbyATMonOCN_open(I))
        ELSE
          tmpscal4=MAX(ZERO,a_QbyATM_open(I))
        ENDIF

! loss of ice cover by melting : computed from
!   (SEAICE_areaLossFormula.EQ.1) from all but only melt conributions by ATM and OCN
!   (SEAICE_areaLossFormula.EQ.2) from net melt-growth>0 by ATM and OCN
!   (SEAICE_areaLossFormula.EQ.3) from predicted melt by ATM
        IF (SEAICE_areaLossFormula.EQ.1) THEN
          tmpscal3 = MIN( 0.0_wp , d_HEFFbyATMonOCN_cover(I) )    &
           + MIN( 0.0_wp , d_HEFFbyATMonOCN_open(I) )             &
           + MIN( 0.0_wp , d_HEFFbyOCNonICE(I) )
        ELSEIF (SEAICE_areaLossFormula.EQ.2) THEN
          tmpscal3 = MIN( 0.0_wp , d_HEFFbyATMonOCN_cover(I)      &
            + d_HEFFbyATMonOCN_open(I) + d_HEFFbyOCNonICE(I) )
        ELSE
!         compute heff after ice melt by ocn:
          tmpscal0=HEFF(I) - d_HEFFbyATMonOCN(I)
!         compute available heat left after snow melt by atm:
          tmpscal1= a_QbyATM_open(I)+a_QbyATM_cover(I)            &
                    - d_HSNWbyATMonSNW(I)*SNOW2ICE
!         could not melt more than all the ice
          tmpscal2 = MAX(-tmpscal0,tmpscal1)
          tmpscal3 = MIN(ZERO,tmpscal2)
        ENDIF

! apply tendency
        IF ( (HEFF(i).GT.0.0_wp).OR.(HSNOW(i).GT.0.0_wp) ) THEN
         AREA(I)=MAX(0.0_wp, MIN( SEAICE_area_max, AREA(I)  &
           + recip_HO*tmpscal4+HALF*recip_HH*tmpscal3 * areaPDFfac ))
        ELSE
         AREA(I)=0.0_wp
        ENDIF
      enddo
#ifdef SEAICE_ITD
!     transfer 1st category values back into ITD variables
!$acc loop
      do i = 1,nlpb
        HEFFITD(I,1)=HEFF(I)
        AREAITD(I,1)=AREA(I)
        HSNOWITD(I,1)=HSNOW(I)
      enddo
!     now melt ice laterally in all other thickness categories
!     (areal growth, i.e. new ice formation, only occurrs in 1st category)
      IF (SEAICE_multDim .gt. 1) THEN
!$acc loop collapse(2) private(tmpscal1)
       DO IT=2,SEAICE_multDim
         do i = 1,nlpb
          IF (HEFFITD(I,IT).LE.ZERO) THEN
!     when thickness is zero, area should be zero, too:
           AREAITD(I,IT)=ZERO
          ELSE
!     tmpscal1 is the minimal ice concentration after lateral melt that will
!     not lead to an unphysical increase of ice thickness by lateral melt;
!     estimated as the concentration before thermodynamics scaled by the
!     ratio of new ice thickness and ice thickness before thermodynamics
           IF ( HEFFITDpreTH(I,IT).LE.ZERO ) THEN
            tmpscal1=0.0_wp
           ELSE
            tmpscal1=AREAITDpreTH(I,IT)* HEFFITD(I,IT)/HEFFITDpreTH(I,IT)
           ENDIF
!     melt ice laterally based on an average floe sice
!     following Steele (1992)
           AREAITD(I,IT) = AREAITD(I,IT) * (ONE - latMeltFrac(I,IT))
!     limit area reduction so that actual ice thickness does not increase
           AREAITD(I,IT) = max(AREAITD(I,IT), tmpscal1)
          ENDIF
         enddo
       ENDDO
      ENDIF
#endif /* SEAICE_ITD */

!$acc end kernels

#ifdef SEAICE_ITD
! check categories for consistency with limits after growth/melt ...
      IF ( SEAICEuseLinRemapITD ) THEN
        CALL seaice_itd_remap(heffitdPreTH, areaitdPreTH)
      ENDIF
      call seaice_itd_redist
! ... and update total AREA, HEFF, HSNOW
!     (the updated HEFF is used below for ice salinity increments)
      call seaice_itd_sum
#endif /* SEAICE_ITD */


! ===================================================================
! =============PART 5: determine ice salinity increments=============
! ===================================================================

#ifndef SEAICE_VARIABLE_SALINITY
!$acc kernels
!$acc loop private(tmpscal1,tmpscal2,tmpscal3)
      do i = 1,nlpb
        tmpscal1 = d_HEFFbyNEG(I) + d_HEFFbyOCNonICE(I) +     &
                  d_HEFFbyATMonOCN(I) + d_HEFFbyFLOODING(I)   &
                + d_HEFFbySublim(I)
!#ifdef EXF_SEAICE_FRACTION
!!YY: in MaCOM, there is no icefraction input file. is it necessary?
!!YY: in reg_ridge, there is a relaxation of heff by obs.
!     &             + d_HEFFbyRLX(I,J,bi,bj)
!#endif
!atn: can not take out more that surface salinity when SSS<SEAICE_salt0
        tmpscal3 = max( 0.0_wp, min(SEAICE_salt0,tFld(I,kSurface,2)) )
        tmpscal2 = tmpscal1 * tmpscal3 * HEFFM(I)      &
               * recip_deltaTtherm * SEAICE_rhoIce
        saltFlux(I) = tmpscal2
      enddo
!$acc end kernels
#endif /* ndef SEAICE_VARIABLE_SALINITY */

#ifdef SEAICE_VARIABLE_SALINITY
!$acc kernels loop private(tmpscal1,tmpscal2)
      do i = 1,nlpb
! sum up the terms that affect the salt content of the ice pack
        tmpscal1=d_HEFFbyOCNonICE(I)+d_HEFFbyATMonOCN(I)

! recompute HEFF before thermodynamic updates (which is not AREApreTH in legacy code)
        tmpscal2=HEFF(I)-tmpscal1-d_HEFFbyFLOODING(I)
! tmpscal1 > 0 : m of sea ice that is created
!YY: saltfrac:percentage of newly formed ice salinity
        IF ( tmpscal1 .GE. 0.0_wp ) THEN
           saltFlux(I) = HEFFM(I)*recip_deltaTtherm      &
               *SEAICE_saltFrac*tFld(I,kSurface,2)       & 
               *tmpscal1*SEAICE_rhoIce

! tmpscal1 < 0 : m of sea ice that is melted
        ELSE
           saltFlux(I) = HEFFM(I)*recip_deltaTtherm*HSALT(I)*tmpscal1/tmpscal2
        ENDIF

! update HSALT based on surface saltFlux
        HSALT(I) = HSALT(I) + saltFlux(I) * SEAICE_deltaTtherm
        saltFlux(I) = saltFlux(I) + saltFluxAdjust(I)
      enddo
!$acc end kernels
#endif /* SEAICE_VARIABLE_SALINITY */


! ===================================================================
! ==============PART 7: determine ocean model forcing================
! ===================================================================

! compute net heat flux leaving/entering the ocean,
! accounting for the part used in melt/freeze processes

!$acc kernels
#ifdef SEAICE_ITD
! compute total of "mult" fluxes for ocean forcing
!$acc loop
      do i = 1,nlpb
        a_QbyATM_cover(I)   = 0.0_wp
        r_QbyATM_cover(I)   = 0.0_wp
        a_QSWbyATM_cover(I) = 0.0_wp
        r_FWbySublim(I)     = 0.0_wp
      enddo
!$acc loop
      do i = 1,nlpb
!$acc loop seq
        DO IT=1,SEAICE_multDim
! if fluxes in effective ice meters, i.e. ice volume per area, then use:
         a_QbyATM_cover(I)=a_QbyATM_cover(I) + a_QbyATMmult_cover(I,IT)
         r_QbyATM_cover(I)=r_QbyATM_cover(I) + r_QbyATMmult_cover(I,IT)
         a_QSWbyATM_cover(I)=a_QSWbyATM_cover(I) + a_QSWbyATMmult_cover(I,IT)
         r_FWbySublim(I)=r_FWbySublim(I) + r_FWbySublimMult(I,IT)
        enddo
      ENDDO
#endif /* SEAICE_ITD */
!$acc loop
      do i = 1,nlpb
        QNET(I) = r_QbyATM_cover(I) + r_QbyATM_open(I) + a_QSWbyATM_cover(I)  &
            - ( d_HEFFbyOCNonICE(I) + d_HSNWbyOCNonSNW(I)*SNOW2ICE + d_HEFFbyNEG(I)  &
!!YY: relaxation of seaice fraction TO DO LATER
!ifdef EXF_SEAICE_FRACTION
!             + d_HEFFbyRLX(I)   &
!endif
              + d_HSNWbyNEG(I)*SNOW2ICE       &
              - convertPRECIP2HI * snow(i,4) /rhoConstFresh * (ONE-AREApreTH(I))   &
              ) * HEFFM(I)
      enddo
!$acc loop
      do i = 1,nlpb
        QSW(I)  = a_QSWbyATM_cover(I) + a_QSWbyATM_open(I)
      enddo

! switch heat fluxes from 'effective' ice meters to W/m2
!$acc loop
      do i = 1,nlpb
        QNET(I) = QNET(I)*convertHI2Q
        QSW(I)  = QSW(I)*convertHI2Q
      enddo
!$acc end kernels

#ifndef SEAICE_DISABLE_HEATCONSFIX
! treat advective heat flux by ocean-to-ice water exchange (at 0decC)
! ===================================================================
!gf Unlike for evap and precip, the temperature of gained/lost
! ocean liquid water due to melt/freeze of solid water cannot be chosen
! arbitrarily to be e.g. the ocean SST. Indeed the present seaice model
! implies a constant ice temperature of 0degC. If melt/freeze water is exchanged
! at a different temperature, it leads to a loss of conservation in the
! ocean+ice system. While this is mostly a serious issue in the
! real fresh water + non linear free surface framework, a mismatch
! between ice and ocean boundary condition can result in all cases.
! Below we therefore anticipate on external_forcing_surf.F
! to diagnoze and/or apply the correction to QNET.

!$acc kernels
!$acc loop private(tmpscal1,tmpscal3)
        do i = 1,nlpb
!atn: initialize tmpscal1
           tmpscal1 = ZERO
! ocean water going to ice/snow, in precip units
           tmpscal3=rhoConstFresh*HEFFM(I)*(              &
            ( d_HSNWbyATMonSNW(I)*SNOW2ICE                &
            + d_HSNWbyOCNonSNW(I)*SNOW2ICE                &
            + d_HEFFbyOCNonICE(I) + d_HEFFbyATMonOCN(I)   &
            + d_HEFFbyNEG(I)+ d_HSNWbyNEG(I)*SNOW2ICE )   &
            * convertHI2PRECIP                            &
            - snow(i,4)/rhoConstFresh * (ONE-AREApreTH(I)) )
!YY: snow falling on open water (1-area), and snow falling on ice-covered fraction hsnow 
! factor in the heat content as done in external_forcing_surf.F
           IF ( useRealFreshWaterFlux ) THEN     !nonlinear freesurface
             tmpscal1 = - tmpscal3* rcp * tFld(I,kSurface,1)
           ELSE
             tmpscal1 = ZERO
           ENDIF
! remove the mismatch when real fresh water is exchanged (at 0degC here)
           IF ( useRealFreshWaterFlux.AND.SEAICEheatConsFix ) &
             QNET(I)=QNET(I)+tmpscal1
        enddo
!$acc end kernels
#endif /* ndef SEAICE_DISABLE_HEATCONSFIX */

        if(allow_balance_fluxes) then
! compute the net heat flux, including adv. by water, entering ocean+ice

!$acc kernels 
!$acc loop private(tmpscal1,tmpscal2)
          do i = 1,nlpb
!gf 1) SIatmQnt (analogous to qnet; excl. adv. by water exch.)
!ML If I consider the atmosphere above the ice, the surface flux
!ML which is relevant for the air temperature dT/dt Eq
!ML accounts for sensible and radiation (with different treatment
!ML according to wave-length) fluxes but not for "latent heat flux",
!ML since it does not contribute to heating the air.
!ML So this diagnostic is only good for heat budget calculations within
!ML the ice-ocean system.
            SIatmQnt(I) = HEFFM(I)*convertHI2Q*(              &
              a_QSWbyATM_cover(I) + a_QbyATM_cover(I) + a_QbyATM_open(I) )
!gf 2) SItflux (analogous to tflux; includes advection by water
!             exchanged between atmosphere and ocean+ice)
! solid water going to atm, in precip units
            tmpscal1 = rhoConstFresh*HEFFM(I)                   &
             * convertHI2PRECIP * ( - d_HSNWbyRAIN(I)*SNOW2ICE  &
             + a_FWbySublim(I) - r_FWbySublim(I) )
! liquid water going to atm, in precip units
            tmpscal2=HEFFM(I)* (               &
             ( ( zevap(I)-prec(I,4) )*( ONE-AREApreTH(I) ) -runoff(i,4) )   &
               + ( d_HFRWbyRAIN(I) + r_FWbySublim(I) )*convertHI2PRECIP*rhoConstFresh )
! In real fresh water flux + nonlinFS, we factor in the advected specific
! energy (referenced to 0 for 0deC liquid water). In virtual salt flux or
! linFS, rain/evap get a special treatment (see external_forcing_surf.F).
            tmpscal1= - tmpscal1*( -SEAICE_lhFusion + rcp * ZERO )
            IF ( useRealFreshWaterFlux ) THEN     !nonlinear freesurface
              tmpscal2= - tmpscal2*( ZERO + rcp * tFld(I,kSurface,1) )
            ELSE
              tmpscal2= ZERO
            ENDIF
            SItflux(I)= SIatmQnt(I)-tmpscal1-tmpscal2
          enddo
!$acc end kernels
        endif !endif for allow_balance_flux

! compute net fresh water flux leaving/entering
! the ocean, accounting for fresh/salt water stocks.

!$acc kernels
!$acc loop private(tmpscal1)
        do i = 1,nlpb
          tmpscal1= d_HSNWbyATMonSNW(I)*SNOW2ICE     &
                  +d_HFRWbyRAIN(I)                   &
                  +d_HSNWbyOCNonSNW(I)*SNOW2ICE      &
                  +d_HEFFbyOCNonICE(I)               &
                  +d_HEFFbyATMonOCN(I)               &
                  +d_HEFFbyNEG(I)                    &
!ifdef EXF_SEAICE_FRACTION
!                 +d_HEFFbyRLX(I,J,bi,bj)             &
!endif
                  +d_HSNWbyNEG(I)*SNOW2ICE           &
!     If r_FWbySublim>0, then it is evaporated from ocean.
                  +r_FWbySublim(I)
          EmPmR_ice(I)  = HEFFM(I)*(                     &
              (zevap(i)-prec(i,4))* ( ONE - AREApreTH(I) ) -runoff(i,4)  &
              + tmpscal1*convertHI2PRECIP*rhoConstFresh)       !kg/m2/s 
#ifdef SEAICE_ITD
!     beware of the sign: fw2ObyRidge is snow mass moved into the ocean
!     by ridging, so requires a minus sign
          EmPmR_ice(I) = EmPmR_ice(I) - fw2ObyRidge(I)*recip_deltaTtherm * HEFFM(I)
#endif /* SEAICE_ITD */
          if(allow_balance_fluxes) then
! and the flux leaving/entering the ocean+ice
             SIatmFW(I) = HEFFM(I)*(                 &
                 zevap(I)*( ONE - AREApreTH(I) ) - prec(I,4) -runoff(i,4) )  &
                 + a_FWbySublim(I) * SEAICE_rhoIce * recip_deltaTtherm
          endif

        enddo
!$acc end kernels

! Sea Ice Load on the sea surface.
! =================================

        IF ( useRealFreshWaterFlux ) THEN
!$acc kernels loop private(tmpscal1,tmpscal2)
          do i = 1,nlpb
            tmpscal1 = HEFF(I)*SEAICE_rhoIce + HSNOW(I)*SEAICE_rhoSnow
            tmpscal2 = MIN(tmpscal1,heffTooHeavy*rau0)
!           tmpscal2 = HEFF(I)*SEAICE_rhoIce + HSNOW(I)*SEAICE_rhoSnow
            sIceLoad(i) = tmpscal2
          enddo
!$acc end kernels
        ENDIF


! ===================================================================
! =========PART 9: HF/FWF global integrals and balancing=============
! ===================================================================

      if(allow_balance_fluxes) then
! Compute tile integrals of heat/fresh water fluxes to/from atm.
! ==============================================================
        FWFsiTile = 0.0_wp
        IF ( balanceEmPmR ) THEN
!$acc kernels
!$acc loop reduction(+:FWFsiTile)
          do i = 1,nlpb
            FWFsiTile = FWFsiTile + SIatmFW(i) * rAc(i) * maskC(i,kSurface)
          enddo
!$acc end kernels
          CALL mpi_comp_dp_allsum(FWFsiTile)
        ENDIF
! to translate global mean FWF adjustements (see below) we may need :
        FWF2HFsiTile = 0.0_wp
        IF ( balanceEmPmR) THEN
!$acc kernels
!$acc loop reduction(+:FWF2HFsiTile)
          do i = 1,nlpb
          !YY: maskinC is replaced with maskC(ksurface)
            FWF2HFsiTile = FWF2HFsiTile + rcp * tFld(I,kSurface,1) * rAc(i) * maskC(i,kSurface)
          enddo
!$acc end kernels
          CALL mpi_comp_dp_allsum(FWF2HFsiTile)
        ENDIF
        HFsiTile = 0.0_wp
        IF ( balanceQnet ) THEN
!$acc kernels
!$acc loop reduction(+:HFsiTile)
          do i = 1,nlpb
           HFsiTile = HFsiTile + SItflux(i) * rAc(i) * maskC(i,kSurface)
          enddo
!$acc end kernels
          CALL mpi_comp_dp_allsum(HFsiTile)
        ENDIF

! mean SIatmFW
        tmpscal0=FWFsiTile / globalArea
! corresponding mean advection by atm to ocean+ice water exchange
!        (if mean SIatmFW was removed uniformely from ocean)
        tmpscal1=FWFsiTile / globalArea * FWF2HFsiTile / globalArea
! mean SItflux (before potential adjustement due to SIatmFW)
        tmpscal2=HFsiTile / globalArea
! mean SItflux (after potential adjustement due to SIatmFW)
        IF ( balanceEmPmR ) tmpscal2=tmpscal2-tmpscal1

! 3) balancing adjustments
        IF ( balanceEmPmR ) THEN
!$acc kernels loop copyin(tmpscal0,tmpscal1)
          do i = 1,nlpb
            EmPmR_ice(i) = EmPmR_ice(i) - tmpscal0
            SIatmFW(i) = SIatmFW(i) - tmpscal0
!           adjust SItflux consistently
            IF ( useRealFreshWaterFlux ) THEN      !nonlinear freesurface
              tmpscal1=( ZERO + rcp * tFld(I,kSurface,1) )
            ELSE
              tmpscal1=ZERO
            ENDIF
            SItflux(i) = SItflux(i) - tmpscal0*tmpscal1
!           no qnet or tflux adjustement is needed
          enddo
!$acc end kernels
        
          IF (mpi_rank == 0) THEN
             print*, 'rm Global mean of',' SIatmFW = ', tmpscal0, '  @ it=', myIter
          ENDIF
        ENDIF
        IF ( balanceQnet ) THEN
!$acc kernels loop copyin(tmpscal2)
          do i = 1,nlpb
              SItflux(i) = SItflux(i) - tmpscal2
              Qnet(i) = Qnet(i) - tmpscal2
              SIatmQnt(i) = SIatmQnt(i) - tmpscal2
          enddo
!$acc end kernels
          IF (mpi_rank == 0) THEN
             print*,'rm Global mean of', ' SItflux = ', tmpscal2, '  @ it=', myIter
          ENDIF
        ENDIF
      endif ! endif for allow_balance_fluxes

      !YY: update qns/qsr
!$acc kernels loop
      do i = 1,nlpb
        qns_ice(i) = Qsw(i) - Qnet(i)
        qsr_ice(i) = -Qsw(i)
      enddo
!$acc end kernels
!$acc end data

   end subroutine seaice_growth


   subroutine seaice_solve4temp(UG, HICE_ACTUAL, HSNOW_ACTUAL,  &
#ifdef SEAICE_CAP_SUBLIM
     F_lh_max,                                &
#endif
     TSURFin,TSURFout,F_ia, IcePenetSW,FWsublim)

!     *==========================================================*
!     | o Calculate ice growth rate, surface fluxes and
!     |   temperature of ice surface.
!     |   see Hibler, MWR, 108, 1943-1973, 1980
!     *==========================================================*

      implicit none

!     !INPUT PARAMETERS:
!     UG           :: atmospheric wind speed (m/s)
!     HICE_ACTUAL  :: actual ice thickness
!     HSNOW_ACTUAL :: actual snow thickness
!     TSURFin    :: surface temperature of ice/snow in Kelvin
!     !OUTPUT PARAMETERS:
!     TSURFout   :: updated surface temperature of ice/snow in Kelvin
!     F_ia       :: upward seaice/snow surface heat flux to atmosphere (W/m^2)
!     IcePenetSW :: short wave heat flux transmitted through ice (+=upward)
!     FWsublim   :: fresh water (mass) flux due to sublimation (+=up)(kg/m^2/s)
!---- Notes:
!     1) should add IcePenetSW to F_ia to get the net surface heat flux
!        from the atmosphere (IcePenetSW not currently included in F_ia)
!     2) since zero ice/snow heat capacity is assumed, all the absorbed Short
!        -Wave is used to warm the ice/snow surface (heating profile ignored).
!----------
      real(wp) ::  UG (nlpb)
      real(wp) ::  HICE_ACTUAL (nlpb)
      real(wp) ::  HSNOW_ACTUAL(nlpb)
#ifdef SEAICE_CAP_SUBLIM
      real(wp) ::  F_lh_max    (nlpb)
#endif
      real(wp) ::  TSURFin     (nlpb)
      real(wp) ::  TSURFout    (nlpb)
      real(wp) ::  F_ia        (nlpb)
      real(wp) ::  IcePenetSW  (nlpb)
      real(wp) ::  FWsublim    (nlpb)

!     === Local variables ===
!     kSurface  :: vertical index of surface layer
      INTEGER i
      INTEGER kSurface
      INTEGER ITER
      real(wp) ::  D1, D1I
      real(wp) ::  D3(nlpb)
      real(wp) ::  TMELT, XKI, XKS, HCUT, recip_HCUT, XIO
!     SurfMeltTemp :: Temp (K) above which wet-albedo values are used
      real(wp) ::  SurfMeltTemp
!     effConduct :: effective conductivity of combined ice and snow
      real(wp) ::  effConduct(nlpb)
!     lhSublim :: latent heat of sublimation (SEAICE_lhEvap + SEAICE_lhFusion)
      real(wp) ::  lhSublim
!     t1,t2,t3,t4 :: powers of temperature
      real(wp) ::   t1, t2, t3, t4

!-    Constants to calculate Saturation Vapor Pressure
!     Maykut Polynomial Coeff. for Sat. Vapor Press
      !real(wp) ::  C1, C2, C3, C4, C5, QS1
!     Extended temp-range expon. relation Coeff. for Sat. Vapor Press
      real(wp) ::  lnTEN
      real(wp) ::  aa1,aa2,bb1,bb2,Ppascals,cc0,cc1,cc2,cc3t
!     specific humidity at ice surface variables
      real(wp) ::  mm_pi,mm_log10pi

!     tempFrz  :: ocean temperature in contact with ice
!                 (=seawater freezing point) (K)
!     F_c      :: conductive heat flux through seaice+snow (+=upward)
!     F_lwu    :: upward long-wave surface heat flux (+=upward)
!     F_sens   :: sensible surface heat flux         (+=upward)
!     F_lh     :: latent heat flux (sublimation) (+=upward)
!     qhice    :: saturation vapor pressure of snow/ice surface
!     dqh_dTs  :: derivative of qhice w.r.t snow/ice surf. temp
!     dFia_dTs :: derivative of surf heat flux (F_ia) w.r.t surf. temp
      real(wp) ::  tempFrz    (nlpb)
      real(wp) ::  F_c        (nlpb)
      real(wp) ::  F_lwu      (nlpb)
      real(wp) ::  F_sens     (nlpb)
      real(wp) ::  F_lh       (nlpb)
      real(wp) ::  qhice      (nlpb)
      real(wp) ::  dqh_dTs    (nlpb)
      real(wp) ::  dFia_dTs   (nlpb)
      real(wp) ::  absorbedSW (nlpb)
      real(wp) ::  penetSWFrac
      real(wp) ::  delTsurf

!     local copies of global variables
      real(wp) ::  tsurfLoc   (nlpb)
      real(wp) ::  tsurfPrev  (nlpb)
      real(wp) ::  atempLoc   (nlpb)
      real(wp) ::  lwdownLoc  (nlpb)
      real(wp) ::  ALB        (nlpb)
      real(wp) ::  ALB_ICE    (nlpb)
      real(wp) ::  ALB_SNOW   (nlpb)
!     iceOrNot :: this is HICE_ACTUAL.GT.0.
      LOGICAL  ::  iceOrNot(nlpb)


!$acc data present(UG, HICE_ACTUAL, HSNOW_ACTUAL,  &
#ifdef SEAICE_CAP_SUBLIM
!$acc             F_lh_max,                        &
#endif
!$acc             TSURFin,TSURFout,F_ia, IcePenetSW,FWsublim),           &
!$acc      create(D3,effConduct,tempFrz,F_c,F_lwu,F_sens,F_lh,qhice,     &
!$acc            dqh_dTs,dFia_dTs,absorbedSW,penetSWFrac,delTsurf,       &
!$acc            tsurfLoc,tsurfPrev,atempLoc,lwdownLoc,ALB,ALB_ICE,      &
!$acc            ALB_SNOW,iceOrNot)
!-    MAYKUT CONSTANTS FOR SAT. VAP. PRESSURE TEMP. POLYNOMIAL
      !C1=    2.7798202d-06
      !C2=   -2.6913393d-03
      !C3=    0.97920849_wp
      !C4= -158.63779_wp
      !C5= 9653.1925_wp
      !QS1=0.622_wp/1013.0_wp
!-    Extended temp-range expon. relation Coeff. for Sat. Vapor Press
      lnTEN = LOG(10.0_wp)
      aa1 = 2663.5_wp
      aa2 = 12.537_wp
      bb1 = 0.622_wp
      bb2 = 1.0_wp - bb1
      Ppascals = 100000.0_wp
      cc0 = EXP(aa2*lnTEN)
      cc1 = cc0*aa1*bb1*Ppascals*lnTEN
      cc2 = cc0*bb2

      kSurface = nk

!     SENSIBLE HEAT CONSTANT
      D1=SEAICE_dalton*SEAICE_cpAir*SEAICE_rhoAir

!     ICE LATENT HEAT CONSTANT
      lhSublim = SEAICE_lhEvap + SEAICE_lhFusion
      D1I=SEAICE_dalton*lhSublim*SEAICE_rhoAir

!     MELTING TEMPERATURE OF ICE
      TMELT        = celsius2K

!     ICE CONDUCTIVITY
      XKI=SEAICE_iceConduct

!     SNOW CONDUCTIVITY
      XKS=SEAICE_snowConduct

!     CUTOFF SNOW THICKNESS
!     Snow-Thickness above HCUT: SW optically thick snow (=> snow-albedo).
!     Snow-Thickness below HCUT: linear transition to ice-albedo
      HCUT = SEAICE_snowThick
      recip_HCUT = 0.0_wp
      IF ( HCUT.GT.0.0_wp ) recip_HCUT = 1.0_wp / HCUT

!     PENETRATION SHORTWAVE RADIATION FACTOR
      XIO=SEAICE_shortwave

!     Temperature Threshold for wet-albedo:
      SurfMeltTemp = TMELT + SEAICE_wetAlbTemp

!     Initialize variables
!$acc kernels copyin(kSurface,SurfMeltTemp,HCUT,recip_HCUT,XIO,XKI,XKS)
!$acc loop
      do i = 1,nlpb
!     initialise output arrays:
        TSURFout (I) = TSURFin(I)
        F_ia     (I) = 0.0_wp
        IcePenetSW(I)= 0.0_wp
        FWsublim (I) = 0.0_wp
        iceOrNot  (I) = HICE_ACTUAL(I) .GT. 0.0_wp
        absorbedSW(I) = 0.0_wp
        qhice     (I) = 0.0_wp
        dqh_dTs   (I) = 0.0_wp
        F_lh      (I) = 0.0_wp
        F_lwu     (I) = 0.0_wp
        F_sens    (I) = 0.0_wp
!     Make a local copy of LW, surface & atmospheric temperatures
        tsurfLoc (I) = TSURFin(I)
        !lwdownLoc(I,J) = MAX( MIN_LWDOWN, LWDOWN(I,J,bi,bj) )
        !YY: LWDOWN replaced by lwdn(i,4), ATEMP replaced by t10
        lwdownLoc(I) = MAX( MIN_LWDOWN, lwdn(I,4) )
        atempLoc (I) = MAX( celsius2K+MIN_ATEMP, t10(I,4) )

!     FREEZING TEMP. OF SEA WATER (K)
        tempFrz(I) = SEAICE_dTempFrz_dS *tFld(I,kSurface,2)    &
          + SEAICE_tempFrz0 + celsius2K

!     Now determine fixed (relative to tsurf) forcing term in heat budget

        IF(HSNOW_ACTUAL(I).GT.0.0_wp) THEN
!     Stefan-Boltzmann constant times emissivity
!YY: note that snow_emiss is used here, while in MaCOM no emiss thing anywhere
         D3(I)=SEAICE_snow_emiss*SEAICE_boltzmann
!     This is now [(1-emiss)*lwdown - lwdown]
         lwdownLoc(I) = SEAICE_snow_emiss*lwdownLoc(I)
         !lwdownLoc(I) = 0.97_wp*lwdownLoc(I)
        ELSE
!     Stefan-Boltzmann constant times emissivity
         D3(I)=SEAICE_ice_emiss*SEAICE_boltzmann
!     This is now [(1-emiss)*lwdown - lwdown]
         lwdownLoc(I) = SEAICE_ice_emiss*lwdownLoc(I)
         !lwdownLoc(I) = 0.97_wp*lwdownLoc(I)
        ENDIF
      enddo
!$acc loop private(penetSWFrac)
      do i = 1,nlpb
!     DECIDE ON ALBEDO
        IF ( iceOrNot(i) ) THEN

         IF ( latC(I) .LT. 0.0_wp ) THEN
           IF (tsurfLoc(I) .GE. SurfMeltTemp) THEN
             ALB_ICE (I)   = SEAICE_wetIceAlb_south
             ALB_SNOW(I)   = SEAICE_wetSnowAlb_south
           ELSE                  ! no surface melting
             ALB_ICE (I)   = SEAICE_dryIceAlb_south
             ALB_SNOW(I)   = SEAICE_drySnowAlb_south
           ENDIF
         ELSE                   !/ Northern Hemisphere
           IF (tsurfLoc(I) .GE. SurfMeltTemp) THEN
             ALB_ICE (I)   = SEAICE_wetIceAlb
             ALB_SNOW(I)   = SEAICE_wetSnowAlb
           ELSE                  ! no surface melting
             ALB_ICE (I)   = SEAICE_dryIceAlb
             ALB_SNOW(I)   = SEAICE_drySnowAlb
           ENDIF
         ENDIF                  !/ Albedo for snow and ice

!     If actual snow thickness exceeds the cutoff thickness, use snow albedo
         IF (HSNOW_ACTUAL(I) .GT. HCUT) THEN
           ALB(I) = ALB_SNOW(I)
         ELSEIF ( HCUT.LE.ZERO ) THEN
           ALB(I) = ALB_ICE(I)
         ELSE
!     otherwise, use linear transition between ice and snow albedo
           ALB(I) = MIN( ALB_ICE(I) + HSNOW_ACTUAL(I)*recip_HCUT   &
                      *(ALB_SNOW(I) -ALB_ICE(I)), ALB_SNOW(I) )
         ENDIF

!     Determine the fraction of shortwave radiative flux remaining
!     at ocean interface after scattering through the snow and ice.
!     If snow is present, no radiation penetrates through snow+ice
         IF (HSNOW_ACTUAL(I) .GT. 0.0_wp) THEN
           penetSWFrac = 0.0_wp
         ELSE
           penetSWFrac = XIO*EXP(-1.5_wp * HICE_ACTUAL(I))
         ENDIF
!     The shortwave radiative flux leaving ocean beneath ice (+=up).
         IcePenetSW(I) = -(1.0_wp - ALB(I))*penetSWFrac * swdn(I,4)
!     The shortwave radiative flux convergence in the seaice.
         absorbedSW(I) =  (1.0_wp - ALB(I))*(1.0_wp - penetSWFrac)* swdn(I,4)

!     The effective conductivity of the two-layer snow/ice system.
!     Set a minimum sea ice thickness of 5 cm to bound
!     the magnitude of conducive heat fluxes.
!if   * now taken care of by SEAICE_hice_reg in seaice_growth
!        hice_tmp = max(HICE_ACTUAL(I),0.05_wp)
         effConduct(I) = XKI * XKS /           &
             (XKS * HICE_ACTUAL(I) + XKI * HSNOW_ACTUAL(I))

        ENDIF                   !/* iceOrNot */
      enddo
!$acc end kernels

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      DO ITER=1,IMAX_TICE   !YY:  number of iterations for ice surface temp

!$acc kernels copyin(aa1,aa2,lnTEN,bb1,cc1,cc2,Ppascals,D1I,D1,TMELT)
!$acc loop private(t1,t2,t3,t4,mm_log10pi,mm_pi,cc3t)
        do i = 1,nlpb

!-    save tsurf from previous iter
         tsurfPrev(I) = tsurfLoc(I)
         IF ( iceOrNot(I) ) THEN

          t1 = tsurfLoc(I)
          t2 = t1*t1
          t3 = t2*t1
          t4 = t2*t2

!--   Calculate the specific humidity in the BL above the snow/ice
         ! IF ( useMaykutSatVapPoly ) THEN     !YY: default is false
!-    Use the Maykut polynomial
         !  qhice(I)=QS1*(C1*t4+C2*t3 +C3*t2+C4*t1+C5)
         !  dqh_dTs(I) = 0.0_wp
         ! ELSE
!-    Use exponential relation approx., more accurate at low temperatures
!     log 10 of the sat vap pressure
           mm_log10pi = -aa1 / t1 + aa2
!     The saturation vapor pressure (SVP) in the surface
!     boundary layer (BL) above the snow/ice.
!          mm_pi = TEN **(mm_log10pi)
!     The following form does the same, but is faster
           mm_pi = EXP(mm_log10pi*lnTEN)
           qhice(I) = bb1*mm_pi/( Ppascals -(1.0_wp - bb1)*mm_pi )
!     A constant for SVP derivative w.r.t TICE
!          cc3t = TEN **(aa1 / t1)
!     The following form does the same, but is faster
           cc3t = EXP(aa1 / t1 * lnTEN)
!     d(qh)/d(TICE)
           dqh_dTs(I) = cc1*cc3t/((cc2-cc3t*Ppascals)**2 *t2)
         ! ENDIF

!     Calculate the flux terms based on the updated tsurfLoc
!YY: in MaCOM zq_zu is specific humidity at windheight 10m
!YY: in MITgcm, aqh       :: Surface (2m) specific humidity in kg/kg
!              Typical range: 0 < aqh < 0.02
!in MaCOM, q10 is 10m-height value, right? here aqh is 2m-height

          F_c(I)    = effConduct(I)*(tempFrz(I)-t1)
          !F_lh(I)   = D1I*UG(I)*(qhice(I)-AQH(I,J,bi,bj))
          !YY: upward + here, in MaCOM it is reverse
          !qhice here is equal to zqsatw at sea surface?? so should use q10??
          !remember zq_zu is not a global var, which one to use?? ZY
          F_lh(I)   = D1I*UG(I)*(qhice(I)-q10(I,4))
          !F_lh(I)   = D1I*UG(I)*(qhice(I)-zq_zu(I))
#ifdef SEAICE_CAP_SUBLIM
!     if the latent heat flux implied by tsurfLoc exceeds
!     F_lh_max, cap F_lh and decouple the flux magnitude from tIce (tsurfLoc)
          IF (F_lh(I) .GT. F_lh_max(I)) THEN
             F_lh(I)  = F_lh_max(I)
             dqh_dTs(I) = ZERO
          ENDIF
#endif /* SEAICE_CAP_SUBLIM */

          F_lwu(I) = t4 * D3(I)
          F_sens(I)= D1 * UG(I) * (t1 - atempLoc(I))
          F_ia(I) = -lwdownLoc(I) - absorbedSW(I) + F_lwu(I) + F_sens(I) + F_lh(I)
!     d(F_ia)/d(Tsurf)
          dFia_dTs(I) = 4.0_wp*D3(I)*t3 + D1*UG(I) + D1I*UG(I)*dqh_dTs(I)

!-    Update tsurf as solution of : Fc = Fia + d/dT(Fia - Fc) *delta.tsurf
          tsurfLoc(I) = tsurfLoc(I) + ( F_c(I)-F_ia(I) ) / ( effConduct(I)+dFia_dTs(I) )

         ENDIF
       enddo
      ! IF ( useMaykutSatVapPoly ) THEN     !YY: default is false
      !   do i = 1,nlpb
      !     tsurfLoc(I) = MAX( celsius2K+MIN_TICE, tsurfLoc(I) )
      !   enddo
      ! ENDIF
!$acc loop
       do i = 1,nlpb
         tsurfLoc(I) = MIN( tsurfLoc(I), TMELT )
       enddo
!
!$acc end kernels
      ENDDO                     !/* Iterations */
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( postSolvTempIter.EQ.2 ) THEN
!$acc kernels copyin(aa1,aa2,lnTEN,bb1,Ppascals,D1I,D1)
!$acc loop private(t1,t2,t4,mm_log10pi,mm_pi)
        do i = 1,nlpb
         IF ( iceOrNot(I) ) THEN
!     Recalculate the fluxes based on the (possibly) adjusted TSURF
          t1 = tsurfLoc(I)
          t2 = t1*t1
          !t3 = t2*t1
          t4 = t2*t2

         ! IF ( useMaykutSatVapPoly ) THEN
         !  qhice(I)=QS1*(C1*t4+C2*t3 +C3*t2+C4*t1+C5)
         ! ELSE
!     log 10 of the sat vap pressure
           mm_log10pi = -aa1 / t1 + aa2
!     saturation vapor pressure
           mm_pi = EXP(mm_log10pi*lnTEN)
!     over ice specific humidity
           qhice(I) = bb1*mm_pi/( Ppascals -(1.0_wp - bb1)*mm_pi )
         ! ENDIF
          F_c(I)  = effConduct(I) * (tempFrz(I) - t1)
          !YY: ??
          !F_lh(I) = D1I * UG(I)*(qhice(I)-AQH(I))
          !F_lh(I) = D1I * UG(I)*(qhice(I)-zq_zu(I))
          F_lh(I) = D1I * UG(I)*(qhice(I)-q10(I,4))
#ifdef SEAICE_CAP_SUBLIM
          IF (F_lh(I) .GT. F_lh_max(I)) THEN
             F_lh(I)  = F_lh_max(I)
          ENDIF
#endif /* SEAICE_CAP_SUBLIM */
          F_lwu (I) = t4 * D3(I)
          F_sens(I) = D1 * UG(I) * (t1 - atempLoc(I))
!     The flux between the ice/snow surface and the atmosphere.
          F_ia(I) = -lwdownLoc(I) - absorbedSW(I) + F_lwu(I) + F_sens(I) + F_lh(I)

         ENDIF
        enddo
!$acc end kernels
      ELSEIF ( postSolvTempIter.EQ.1 ) THEN
!$acc kernels copyin(D1I)
!$acc loop private(delTsurf)
        do i = 1,nlpb 
         IF ( iceOrNot(I) ) THEN
!     Update fluxes (consistent with the linearized formulation)
          delTsurf  = tsurfLoc(I)-tsurfPrev(I)
          F_c (I)  = effConduct(I)*(tempFrz(I)-tsurfLoc(I))
          F_ia(I) = F_ia(I) + dFia_dTs(I)*delTsurf
          F_lh(I) = F_lh(I) + D1I*UG(I)*dqh_dTs(I)*delTsurf
         ENDIF
       enddo
!$acc end kernels
      ENDIF
!$acc kernels copyin(lhSublim)
!$acc loop
      do i = 1,nlpb
        IF ( iceOrNot(I) ) THEN

!     Save updated tsurf and finalize the flux terms
         TSURFout(I) = tsurfLoc(I)
!     Fresh water flux (kg/m^2/s) from latent heat of sublimation.
!     F_lh is positive upward (sea ice looses heat) and FWsublim
!     is also positive upward (atmosphere gains freshwater)
         FWsublim(I) = F_lh(I)/lhSublim


        ENDIF                   !/* iceOrNot */
      enddo
!$acc end kernels
!$acc end data

  end subroutine seaice_solve4temp

end module mitice_thermodyn
