module mitice_init
   use mod_misc_basic
   use mod_csp_basic
   use mitice_parameters
   use mitice_vars
   use mitice_utility
   use mitice_debug
   use mitice_io
   use mitice_ave
   use mitice_itd

   implicit none

contains
    subroutine mitice_read_params
!   *==========================================================*
!    Called from MaCOM main program
!    read sea ice model parameters from a seperate namelist: namelist.mitice
!   *==========================================================*

      implicit none

!LOCAL VARIABLES:
      CHARACTER(lc) msg   !Informational/error message buffer
      LOGICAL chkFlag, exst
      INTEGER l
      integer nError
      integer istat
      integer, parameter ::  mitice_nml_unit=36

      real(wp) ::  tmp

!-    Old parameters (to be retired one day):
      real(wp) ::  SEAICE_availHeatTaper
      real(wp) ::  SEAICE_gamma_t, SEAICE_gamma_t_frz, SEAICE_availHeatFracFrz

!--   SEAICE parameters
      NAMELIST /SEAICE_PARAMS/      &
      SEAICEuseDYNAMICS, SEAICEuseFREEDRIFT,           &
      SEAICEuseTEM, SEAICEuseTilt,      &
      useHB87stressCoupling,     &
      usePW79thermodynamics, SEAICEuseFlooding,      &
      SEAICE_growMeltByConv,SEAICEadvHeff, SEAICEadvArea, SEAICEadvSnow,      &
      SEAICEadvSalt, SEAICEaddSnowMass,SEAICEmomAdvection,SEAICE_clipVelocities,      &
      SEAICE_no_slip, SEAICEetaZmethod, IMAX_TICE, postSolvTempIter,      &
      SEAICEdiffKhHeff, SEAICEdiffKhSnow, SEAICEdiffKhArea,      &
      SEAICEdiffKhSalt,SEAICE_deltaTtherm, SEAICE_deltaTdyn,      &
      SEAICE_deltaTevp, SEAICE_elasticParm, SEAICE_evpTauRelax,      &
      SEAICE_evpDampC, SEAICEnEVPstarSteps,SEAICE_evpAlpha, SEAICE_evpBeta,      &
      SEAICEaEVPcoeff, SEAICEaEVPcStar, SEAICEaEVPalphaMin,      &
      SEAICE_zetaMin, SEAICE_zetaMaxFac,SEAICEuseLinRemapITD,      &
      useHibler79IceStrength, SEAICEpartFunc, SEAICEredistFunc,      &
      SEAICEridgingIterMax, SEAICEsimpleRidging, SEAICEsnowFracRidge,      &
      SEAICEgStar, SEAICEhStar, SEAICEaStar, SEAICEshearParm,      &
      SEAICEmuRidging, SEAICEmaxRaft, SEAICE_cf,      &
      SEAICEpresH0, SEAICEpresPow0, SEAICEpresPow1,      &
      SEAICE_initialHEFF, SEAICE_areaGainFormula, SEAICE_areaLossFormula,      &
      SEAICE_doOpenWaterGrowth, SEAICE_doOpenWaterMelt,      &
      SEAICE_rhoAir, SEAICE_rhoIce, SEAICE_rhoSnow, ICE2WATR,      &
      SEAICE_cpAir, SEAICEscaleSurfStress,      &
      SEAICE_drag, SEAICE_waterDrag, SEAICEdWatMin, SEAICE_dryIceAlb,      &
      SEAICE_wetIceAlb, SEAICE_drySnowAlb, SEAICE_wetSnowAlb, HO,      &
      SEAICE_drag_south, SEAICE_waterDrag_south,      &
      SEAICE_dryIceAlb_south, SEAICE_wetIceAlb_south,      &
      SEAICE_drySnowAlb_south, SEAICE_wetSnowAlb_south, HO_south,      &
      SEAICE_cBasalStar, SEAICEbasalDragU0, SEAICEbasalDragK1,      &
      SEAICEbasalDragK2, SEAICE_wetAlbTemp,      &
      SEAICE_strength, SEAICE_cStar, SEAICE_eccen,      &
      SEAICEpressReplFac, SEAICE_tensilFac, SEAICE_tensilDepth,      &
      SEAICE_lhFusion, SEAICE_lhEvap, SEAICE_dalton,      &
      SEAICE_iceConduct, SEAICE_snowConduct,      &
      SEAICE_emissivity, SEAICE_ice_emiss, SEAICE_snow_emiss,      &
      SEAICE_snowThick, SEAICE_shortwave,  OCEAN_drag,      &
      SEAICE_tempFrz0, SEAICE_dTempFrz_dS, SEAICE_salt0,      &
      SEAICE_saltFrac, SEAICEstressFactor, SEAICE_availHeatTaper,      &
      SEAICE_mcPheePiston, SEAICE_frazilFrac, SEAICE_mcPheeTaper,      &
      SEAICE_mcPheeStepFunc, SEAICE_gamma_t, SEAICE_gamma_t_frz,      &
      SEAICE_availHeatFrac, SEAICE_availHeatFracFrz, SEAICE_PDF,      &
      AreaFile, HeffFile, uIceFile, vIceFile, HsnowFile, HsaltFile,      &
      SEAICEheatConsFix, SEAICE_multDim, SEAICE_useMultDimSnow,      &
      SEAICE_deltaMin, SEAICE_area_reg, SEAICE_hice_reg,      &
      SEAICE_area_floor, SEAICE_area_max, SEAICE_tauAreaObsRelax,      &
      SEAICE_airTurnAngle, SEAICE_waterTurnAngle,      &
      MIN_ATEMP, MIN_LWDOWN,  MIN_TICE,SEAICE_EPS, SEAICE_EPS_SQ,      &
      SEAICEwriteState, SEAICEuseEVPpickup,      &
      SEAICEuseEVPstar, SEAICEuseEVPrev,      &
#ifdef SEAICE_ITD
      Hlimit_c1, Hlimit_c2, Hlimit_c3, Hlimit,  &
#endif
      SEAICE_monFreq, SEAICE_dumpFreq, SEAICE_taveFreq, dumpFileIntv

#ifdef SEAICE_ITD
      allocate(Hlimit(0:nITD))
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      if(mpi_rank==0) then

!--   set default sea ice parameters
       SEAICEuseDYNAMICS  = .TRUE.
       SEAICEuseFREEDRIFT = .FALSE.
       SEAICEuseTilt      = .TRUE.
       SEAICEheatConsFix  = .FALSE.
       SEAICEuseTEM       = .FALSE.
       SEAICEuseEVPpickup = .TRUE.    !using seaice_sigma field for hotstart; default T for EVPsolver
       SEAICEuseEVPstar   = .TRUE.
       SEAICEuseEVPrev    = .TRUE.
       SEAICE_growMeltByConv = .FALSE.
       useHB87stressCoupling = .TRUE.
       SEAICEscaleSurfStress = .TRUE.
       SEAICEaddSnowMass     = .TRUE.
       usePW79thermodynamics = .TRUE.    !ZERO-layer thermodyn
!      start of ridging parameters
       useHibler79IceStrength= .TRUE.
       SEAICEsimpleRidging   = .TRUE.
!      The range of this proportionality constant is 2 to 17
       SEAICE_cf             = 17.0_wp
!      By default use partition function of Thorndyke et al. (1975) ...
       SEAICEpartFunc        = 0
!      and redistribution function of Hibler (1980)
       SEAICEredistFunc      = 0
       SEAICEridgingIterMax  = 10
!      This parameter is not well constraint (Lipscomb et al. 2007)
       SEAICEshearParm       = 0.5_wp
!      Thorndyke et al. (1975)
       SEAICEgStar           = 0.15_wp
!      suggested by Hibler (1980), Flato and Hibler (1995)
!      SEAICEhStar           = 100. _d 0
!      but this value of 25 seems to give thinner ridges in better agreement
!      with observations (according to Lipscomb et al 2007)
       SEAICEhStar           =  25.0_wp
!      according to Lipscomb et al. (2007) these values for aStar and mu
!      are approximately equivalent to gStar=0.15 (aStar = gStar/3) for
!      SEAICEpartFunc = 1 ...
       SEAICEaStar           = 0.05_wp
!      ... and hStar=25 for SEAICEredistFunc = 1
!      Libscomb et al. (2007): mu =  3,  4,  5,   6
!      correspond to        hStar = 25, 50, 75, 100
       SEAICEmuRidging       = 3.0_wp
       SEAICEmaxRaft         = 1.0_wp
       SEAICEsnowFracRidge   = 0.5_wp
       SEAICEuseLinRemapITD  = .TRUE.
!      end ridging parampeters
       !useMaykutSatVapPoly = .FALSE.
       SEAICEadvHeff      = .TRUE.
       SEAICEadvArea      = .TRUE.
       SEAICEadvSnow      = .TRUE.
#ifdef SEAICE_VARIABLE_SALINITY
       SEAICEadvSalt      = .TRUE.
#else
      SEAICEadvSalt      = .FALSE.
#endif
      SEAICEmomAdvection       = .FALSE.      !default: false
      SEAICEuseFlooding  = .TRUE.
      SEAICE_no_slip     = .FALSE.    !no slip boundary
      SEAICE_clipVelocities = .FALSE.
      SEAICEetaZmethod   = 3
      SEAICEdiffKhArea   = UNSET_RL
      SEAICEdiffKhHeff   = UNSET_RL
      SEAICEdiffKhSnow   = UNSET_RL
      SEAICEdiffKhSalt   = UNSET_RL
      !YY: default using dTtracer in MaCOM
      SEAICE_deltaTtherm = dTtracer
      SEAICE_deltaTdyn   = dTtracer
      SEAICE_deltaTevp   = UNSET_RL

!     Hunke, JCP, 2001 use 615 kg/m^2 for this, but does not recommend using it
      SEAICE_evpDampC    = -1.0_wp
      SEAICE_zetaMin     = 0.0_wp
      SEAICE_zetaMaxFac  = 2.5d8
      SEAICEpresH0       = 1.0_wp
      SEAICEpresPow0     = 1
      SEAICEpresPow1     = 1
      SEAICE_evpTauRelax = -1.0_wp
      SEAICE_elasticParm = 0.33333333333333333333333333_wp
      SEAICE_evpAlpha    = UNSET_RL
      SEAICE_evpBeta     = UNSET_RL
      SEAICEnEVPstarSteps = UNSET_I
      SEAICEaEVPcoeff    = UNSET_RL
      SEAICEaEVPcStar    = UNSET_RL
      SEAICEaEVPalphaMin = UNSET_RL
      SEAICE_initialHEFF = ZERO
#ifdef SEAICE_ITD
!     Coefficients used to calculate sea ice thickness category limits
!     after Lipscomb et al. (2001, JGR), Equ. 22
!     choose between
!      - original parameters of Lipscomb et al. (2001):
!        c1=3.0/N, c2=15*c1, c3=3.0
!      - and a higher resolution of thin end of ITD:
!        c1=1.5/N, c2=42*c1, c3=3.3
      DO l = 0, nITD
       Hlimit(l) = UNSET_RL
      ENDDO
      Hlimit_c1          = 3.0_wp
      Hlimit_c2          = 15.0_wp
      Hlimit_c3          = 3.0_wp
#endif
      SEAICE_rhoIce      = 0.91d03
      !SEAICE_rhoIce      = rhoic
      !SEAICE_rhoSnow     = 330.0_wp
      SEAICE_rhoSnow     = rhosn
      ICE2WATR           = UNSET_RL
      SEAICE_drag        = 0.001_wp    !default 0.001 
      OCEAN_drag         = 0.001_wp  
      SEAICE_waterDrag   = 0.0055_wp 
      SEAICEdWatMin      = 0.25_wp   
      SEAICE_dryIceAlb   = 0.75_wp   
      SEAICE_wetIceAlb   = 0.66_wp   
      SEAICE_drySnowAlb  = 0.84_wp   
      SEAICE_wetSnowAlb  = 0.7_wp    
      HO                 = 0.5_wp    
      SEAICE_drag_south       = UNSET_RL
      SEAICE_waterDrag_south  = UNSET_RL
      SEAICE_dryIceAlb_south  = UNSET_RL
      SEAICE_wetIceAlb_south  = UNSET_RL
      SEAICE_drySnowAlb_south = UNSET_RL
      SEAICE_wetSnowAlb_south = UNSET_RL
      HO_south                = UNSET_RL
!     basal drag parameters following Lemieux et al. (2015)
      SEAICE_cBasalStar = UNSET_RL
      SEAICEbasalDragU0 =  0.0001   ! YY: here I increase it to 10-4 to test 
      !SEAICEbasalDragU0 =  5.0d-05   ! YY: here I increase it to 10-3 to test 
      SEAICEbasalDragK1 =  8.0_wp
      SEAICEbasalDragK2 =  0.0_wp    !default 0.0
!     Lemieux et al. (2015) recommend: SEAICEbasalDragK2 = 15.0
!
      SEAICE_wetAlbTemp  = -1.0d-3
      SEAICE_strength    = 2.75d04
      SEAICE_cStar       = 20.0_wp
      SEAICEpressReplFac = 1.0_wp 
      SEAICE_eccen       = 2.0_wp 
      SEAICE_tensilFac   = 0.0_wp
      SEAICE_tensilDepth = 0.0_wp
!     coefficients for flux computations/bulk formulae
      SEAICE_dalton      = 1.75d-03
!YY: some parameters using values from MaCOM : mod_misc_basic
      SEAICE_rhoAir     = rhoa
      SEAICE_cpAir      = cpa
      SEAICE_lhEvap     = Lv
      SEAICE_lhFusion   = 3.34d05      !YY, MaCOM lfus
      !YY: in exf package, ice and snow emis is 0.95
      !YY: in MaCOM. emissivity is not used. ask ZY
      SEAICE_emissivity = 5.5d-08 / 5.67d-08
      SEAICE_ice_emiss  = SEAICE_emissivity
      SEAICE_snow_emiss = SEAICE_emissivity

      SEAICE_iceConduct  = 2.1656_wp
      SEAICE_snowConduct = 3.1d-01
      SEAICE_snowThick   = 0.15_wp
      SEAICE_shortwave   = 0.30_wp
      SEAICE_salt0       = 0.0_wp 
      SEAICE_saltFrac    = 0.0_wp 

#ifdef SEAICE_ITD
      SEAICE_multDim     = nITD
      SEAICE_PDF(1)     = 1.0_wp
      DO l=2,nITD
       SEAICE_PDF(l)     = 0.0_wp
      ENDDO
#else
      SEAICE_multDim     = 1
      DO l=1,nITD
       SEAICE_PDF(l)     = UNSET_RL
      ENDDO
#endif
      SEAICE_useMultDimSnow = .TRUE.

!     default to be set later (ocean-seaice turbulent flux coeff):
      SEAICE_mcPheeStepFunc = .FALSE.
      SEAICE_mcPheeTaper    = UNSET_RL
      SEAICE_availHeatTaper = UNSET_RL
      SEAICE_mcPheePiston   = UNSET_RL
      SEAICE_frazilFrac     = UNSET_RL
      SEAICE_gamma_t        = UNSET_RL
      SEAICE_gamma_t_frz    = UNSET_RL
      SEAICE_availHeatFrac  = UNSET_RL
      SEAICE_availHeatFracFrz = UNSET_RL
      SEAICE_doOpenWaterGrowth=.TRUE.
      SEAICE_doOpenWaterMelt=.False.
      SEAICE_areaLossFormula=1
      SEAICE_areaGainFormula=1
      SEAICE_tempFrz0    = 0.0901_wp
      SEAICE_dTempFrz_dS = -0.0575_wp
      SEAICEstressFactor = 1.0_wp
      SEAICE_tauAreaObsRelax = -999.0_wp
      AreaFile   = ' '
      HsnowFile  = ' '
      HsaltFile  = ' '
      HeffFile   = ' '
      uIceFile   = ' '
      vIceFile   = ' '
      IMAX_TICE  = 10
      postSolvTempIter = 2

      SEAICE_area_floor = siEPS
      SEAICE_area_reg   = siEPS
      SEAICE_hice_reg   = 0.05_wp
      SEAICE_area_max   = 1.0_wp

      SEAICE_airTurnAngle   = 0.0_wp
      SEAICE_waterTurnAngle = 0.0_wp
      MIN_ATEMP         = -50.0_wp
      MIN_LWDOWN        = 60.0_wp
      MIN_TICE          = -50.0_wp
      SEAICE_deltaMin   = UNSET_RL
      SEAICE_EPS        = 1.0d-10
      SEAICE_EPS_SQ     = SEAICE_EPS*SEAICE_EPS

      !YY: model output or not
      SEAICEwriteState  = .FALSE.
!YY: MaCOM do not have monitor, TO DO LATER
!      SEAICE_monFreq    = monitorFreq
!YY: the best is to specify in setup file 
      SEAICE_dumpFreq   = 1000.0*SEAICE_deltaTtherm
      dumpFileIntv = 1      !1, new dump file with different month; 2, yearly, when dumpFreq is larger, i.e., daily
      !SEAICE_taveFreq   = 10.0*SEAICE_deltaTtherm
      SEAICE_taveFreq   = -1.0        ! if gt 0, time-averaging output enabled


!     Open and read the seaice setup file
      write(mitice_runmsg_unit,'(A)') ' '
      write(mitice_runmsg_unit,'(A)') ' READ SEAICE NAMELIST: opening NAMELIST.MITICE'
! only mpi_rank=0 read nml, then pass to other thread 
!YY: setup file must using following name
!===================================================
      INQUIRE( FILE='namelist.mitice', EXIST=exst )
!===================================================
      IF (.not. exst) then
         WRITE(mitice_runmsg_unit,'(A)') 'SETUP FILE: namelist.mitice does not exist!'
         STOP 'ABNORMAL END: namelist.mitice missing'
      ENDIF
      open (file='namelist.mitice', unit=mitice_nml_unit, action='read', &
           iostat=istat, iomsg=msg, status='old')
      if (istat /= 0) then
         write (mitice_runmsg_unit, *) 'namelist.mitice open failed, error message = ', msg
         write (mitice_runmsg_unit, *) 'we use default seaice parameters'
      else
         read (unit=mitice_nml_unit, nml=SEAICE_PARAMS, iostat=istat, iomsg=msg)
         if (istat /= 0) then
            write (mitice_runmsg_unit, *) 'namelist.mitice: SEAICE_PARAMS read failed, error message = ', msg
            write (mitice_runmsg_unit, *) 'we use default parameters'
         else
             write (mitice_runmsg_unit, *) 'read parameters for MITICE successful'
          endif
      endif
      !close mitice nml file
      close (unit=mitice_nml_unit, status='keep', iostat=istat, iomsg=msg)
      if (istat /= 0) then
         write (*, *) 'namelist.mitice close failed, error message = ', msg
      end if


#ifdef SEAICE_ITD
! SEAICE_multDim has become a runtime parameter but if SEAICE_ITD is defined
!  it needs to equal nITD because of shared code (mostly in seaice_growth.F).
! nITD is set in SEAICE_SIZE.h
      SEAICE_multDim     = nITD
#endif

!--   Set default values (if not specified in seaice namelist)

      IF ( SEAICE_deltaMin .EQ. UNSET_RL ) SEAICE_deltaMin = SEAICE_EPS

!--   If no PDF was prescribed use the default uniform pdf
      tmp = SEAICE_multDim
      DO l = 1, SEAICE_multDim
       IF (SEAICE_PDF(l).EQ.UNSET_RL) SEAICE_PDF(l) = ONE/tmp
      ENDDO
      DO l = SEAICE_multDim+1, nITD
       IF (SEAICE_PDF(l).EQ.UNSET_RL) SEAICE_PDF(l) = 0.0_wp
      ENDDO

      IF (ICE2WATR.EQ.UNSET_RL) ICE2WATR = SEAICE_rhoIce*r1_rau0
      IF (SEAICE_drag_south       .EQ. UNSET_RL) SEAICE_drag_south       = SEAICE_drag
      IF (SEAICE_waterDrag_south  .EQ. UNSET_RL) SEAICE_waterDrag_south  = SEAICE_waterDrag
      IF (SEAICE_dryIceAlb_south  .EQ. UNSET_RL) SEAICE_dryIceAlb_south  = SEAICE_dryIceAlb
      IF (SEAICE_wetIceAlb_south  .EQ. UNSET_RL) SEAICE_wetIceAlb_south  = SEAICE_wetIceAlb
      IF (SEAICE_drySnowAlb_south .EQ. UNSET_RL) SEAICE_drySnowAlb_south = SEAICE_drySnowAlb
      IF (SEAICE_wetSnowAlb_south .EQ. UNSET_RL) SEAICE_wetSnowAlb_south = SEAICE_wetSnowAlb
      IF (HO_south                .EQ. UNSET_RL) HO_south                = HO

!     Basal drag parameter
      IF (SEAICE_cBasalStar .EQ. UNSET_RL) SEAICE_cBasalStar = SEAICE_cStar

      nError = 0

      IF ( SEAICE_deltaTtherm .NE. dTtracer     .OR.          &
          SEAICE_deltaTdyn   .LT. SEAICE_deltaTtherm .OR.     &
          (SEAICE_deltaTdyn/SEAICE_deltaTtherm) .NE.          &
          INT(SEAICE_deltaTdyn/SEAICE_deltaTtherm) ) THEN
         WRITE(mitice_runmsg_unit,'(A)')  &
             'Unsupported combination of SEAICE_deltaTtherm,'
         WRITE(mitice_runmsg_unit,'(A)')  &
             ' SEAICE_deltaTdyn, and dTtracer'
         nError = nError + 1
      ENDIF

!=======================================
!YY: currently only EVP solver is supported. 
!YY: LSR cannot implemented due to MaCOM's grid system
!=======================================
      SEAICEuseEVP = .True.
!     if EVP is turned on, a couple of parameters need to be computed
      IF ( SEAICEuseEVP ) THEN
        IF ( (SEAICE_deltaTdyn/SEAICE_deltaTevp) .NE.               &
          INT(SEAICE_deltaTdyn/SEAICE_deltaTevp) .AND.              &
          .NOT. (SEAICEuseEVPstar.OR.SEAICEuseEVPrev) ) THEN
          WRITE(mitice_runmsg_unit,'(A)')                           &
             'SEAICE_deltaTevp must be a factor of SEAICE_deltaTdyn.'
          nError = nError + 1
        ENDIF
        IF ( SEAICE_elasticParm .LE. 0.0_wp ) THEN
          WRITE(mitice_runmsg_unit,'(A)')                           &
           'SEAICE_elasticParm must greater than 0.'
          nError = nError + 1
        ENDIF
        IF ( SEAICE_evpTauRelax .LE. 0.0_wp )                       &
           SEAICE_evpTauRelax = SEAICE_deltaTdyn*SEAICE_elasticParm
!     determine number of internal steps
        IF ( SEAICEnEVPstarSteps.EQ.UNSET_I ) THEN
          IF ( SEAICE_deltaTevp.EQ.UNSET_RL ) THEN
            WRITE(mitice_runmsg_unit,'(A,A)')                       &
             'S/R SEAICE_readparms: Either ',                       &
             'SEAICEnEVPstarSteps or SEAICE_deltaTevp need to be set.'
            nError = nError + 1
          ELSE
            SEAICEnEVPstarSteps = INT(SEAICE_deltaTdyn/SEAICE_deltaTevp)
          ENDIF
        ENDIF
!     default: evpAlpha = evpBeta
        IF ( SEAICE_evpAlpha .NE. UNSET_RL .AND.                     &
           SEAICE_evpBeta .EQ. UNSET_RL )                            &
             SEAICE_evpBeta = SEAICE_evpAlpha
        IF ( SEAICE_evpBeta .NE. UNSET_RL .AND.                      &
           SEAICE_evpAlpha .EQ. UNSET_RL )                           &
             SEAICE_evpAlpha = SEAICE_evpBeta
!     derive other parameters
        IF ( SEAICE_evpBeta .EQ. UNSET_RL ) THEN
          SEAICE_evpBeta   = SEAICE_deltaTdyn/SEAICE_deltaTevp
        ELSE
          SEAICE_deltaTevp = SEAICE_deltaTdyn/SEAICE_evpBeta
        ENDIF
        IF ( SEAICE_evpAlpha .EQ. UNSET_RL ) THEN
         SEAICE_evpAlpha = 2.0_wp * SEAICE_evpTauRelax/SEAICE_deltaTevp
        ELSE
         SEAICE_evpTauRelax = 0.5_wp *SEAICE_evpAlpha*SEAICE_deltaTevp
        ENDIF
!     this turns on adaptive EVP
        IF ( SEAICEaEVPcoeff .NE. UNSET_RL ) THEN
         IF ( SEAICEaEVPcStar  .EQ.UNSET_RL) SEAICEaEVPcStar   =4.0_wp
         IF (SEAICEaEVPalphaMin.EQ.UNSET_RL) SEAICEaEVPalphaMin=5.0_wp
!     For adaptive EVP we do not need constant parameters alpha and
!     beta, because they are computed dynamically. Reset them to
!     undefined here, so that we know if something funny is going on.
         SEAICE_evpAlpha     = UNSET_RL
         SEAICE_evpBeta      = UNSET_RL
        ENDIF
      ENDIF     !endif for SEAICE_ALLOW_EVP

      IF ( SEAICEuseFREEDRIFT ) SEAICEuseEVP = .FALSE.
      IF ( SEAICEuseFREEDRIFT ) THEN
        WRITE(mitice_runmsg_unit,'(A,A)')                                  &
            'WARNING FROM S/R MITICE_READPARMS:',                          &
            ' switch seaice from EVP to "free drift"'
      ENDIF

#ifndef SEAICE_ITD
      IF ( .NOT.useHibler79IceStrength ) THEN
       useHibler79IceStrength = .TRUE.
       WRITE(mitice_runmsg_unit,'(A)')   &
           'WARNING FROM S/R SEAICE_READPARMS: resetting useHibler79IceStrength = .TRUE., because'
       WRITE(mitice_runmsg_unit,'(A)')     &
           'WARNING FROM S/R SEAICE_READPARMS: SEAICE_ITD is not defined.'
      ENDIF
#endif /* SEAICE_ITD */

!-    The old ways of specifying mcPheeTaper, mcPheePiston & frazilFrac:
!     a) prevent multiple specification of the same coeff;
!     b) if specified, then try to recover old way of setting & default.
      IF ( SEAICE_mcPheeTaper .EQ. UNSET_RL ) THEN
       IF ( SEAICE_availHeatTaper.EQ.UNSET_RL ) THEN
         SEAICE_mcPheeTaper = 0.0_wp
       ELSE
         SEAICE_mcPheeTaper = SEAICE_availHeatTaper
       ENDIF
      ELSEIF ( SEAICE_availHeatTaper.NE.UNSET_RL ) THEN
         WRITE(mitice_runmsg_unit,'(2A)') 'S/R SEAICE_READPARMS: Cannot specify ',  &
         'both SEAICE_mcPheeTaper & SEAICE_availHeatTaper'
         nError = nError + 1
      ENDIF

!-    set SEAICE_frazilFrac if not yet done
      IF ( SEAICE_gamma_t_frz .NE. UNSET_RL ) THEN
       IF ( SEAICE_frazilFrac .EQ. UNSET_RL ) THEN
         SEAICE_frazilFrac = SEAICE_deltaTtherm/SEAICE_gamma_t_frz
       ELSE
         WRITE(mitice_runmsg_unit,'(2A)') 'S/R SEAICE_READPARMS: Cannot specify ',  &
          'both SEAICE_frazilFrac & SEAICE_gamma_t_frz'
         nError = nError + 1
       ENDIF
      ENDIF
      IF ( SEAICE_availHeatFracFrz.NE.UNSET_RL ) THEN
       IF ( SEAICE_frazilFrac .EQ. UNSET_RL ) THEN
         SEAICE_frazilFrac = SEAICE_availHeatFracFrz
       ELSE
        IF ( SEAICE_gamma_t_frz .EQ. UNSET_RL ) THEN
         WRITE(mitice_runmsg_unit,'(2A)') 'S/R SEAICE_READPARMS: Cannot specify ',  &
           'both SEAICE_frazilFrac  & SEAICE_availHeatFracFrz'
        ELSE
         WRITE(mitice_runmsg_unit,'(2A)') 'S/R SEAICE_READPARMS: Cannot specify ',  &
           'both SEAICE_gamma_t_frz & SEAICE_availHeatFracFrz'
        ENDIF
       nError = nError + 1
       ENDIF
      ENDIF

!     the default for SEAICE_gamma_t_frz use to be SEAICE_gamma_t:
      IF ( SEAICE_gamma_t .NE. UNSET_RL .AND.       &
           SEAICE_frazilFrac .EQ. UNSET_RL ) THEN
         SEAICE_frazilFrac = SEAICE_deltaTtherm/SEAICE_gamma_t
      ENDIF
!     the default for SEAICE_availHeatFracFrz use to be SEAICE_availHeatFrac:
      IF ( SEAICE_availHeatFrac.NE.UNSET_RL .AND.   &
           SEAICE_frazilFrac .EQ. UNSET_RL ) THEN
         SEAICE_frazilFrac = SEAICE_availHeatFrac
      ENDIF
      IF ( SEAICE_frazilFrac .EQ. UNSET_RL ) THEN
         SEAICE_frazilFrac = 1.0_wp
      ENDIF

!-    start by setting SEAICE_availHeatFrac (used in seaice_init_fixed.F
!     to set SEAICE_mcPheePiston once drF is known)
      IF ( SEAICE_gamma_t .NE. UNSET_RL ) THEN
       IF ( SEAICE_availHeatFrac.EQ.UNSET_RL ) THEN
         SEAICE_availHeatFrac = SEAICE_deltaTtherm/SEAICE_gamma_t
       ELSE
        WRITE(mitice_runmsg_unit,'(2A)') 'S/R SEAICE_READPARMS: Cannot specify ',  &
         'both SEAICE_gamma_t & SEAICE_availHeatFrac'
         nError = nError + 1
       ENDIF
      ENDIF
      IF ( SEAICE_mcPheePiston .NE. UNSET_RL .AND.     &
           SEAICE_availHeatFrac.NE. UNSET_RL ) THEN
        IF ( SEAICE_gamma_t .EQ. UNSET_RL ) THEN
         WRITE(mitice_runmsg_unit,'(2A)') 'S/R SEAICE_READPARMS: Cannot specify ',  &
           'both SEAICE_mcPheePiston & SEAICE_availHeatFrac'
        ELSE
         WRITE(mitice_runmsg_unit,'(2A)') 'S/R SEAICE_READPARMS: Cannot specify ',  &
           'both SEAICE_mcPheePiston & SEAICE_gamma_t'
        ENDIF
        nError = nError + 1
      ENDIF

!     Set diffusivity if not done in namelist
      IF ( SEAICEdiffKhArea .EQ. UNSET_RL )SEAICEdiffKhArea = SEAICEdiffKhHeff
      IF ( SEAICEdiffKhArea .EQ. UNSET_RL )SEAICEdiffKhArea = 0.0_wp
      IF ( SEAICEdiffKhHeff .EQ. UNSET_RL )SEAICEdiffKhHeff = SEAICEdiffKhArea
      IF ( SEAICEdiffKhSnow .EQ. UNSET_RL )SEAICEdiffKhSnow = SEAICEdiffKhHeff
      IF ( SEAICEdiffKhSalt .EQ. UNSET_RL )SEAICEdiffKhSalt = SEAICEdiffKhHeff


      IF ( nError .GT. 0 ) THEN
       WRITE(mitice_runmsg_unit,'(2A)') 'S/R SEAICE_READPARMS: ',    &
         'Error reading parameter file namelist.mitice'
       WRITE(mitice_runmsg_unit,'(A,I3,A)') 'S/R SEAICE_READPARMS: ', nError,   &
         ' parameters values are inconsistent or incomplete'
      ENDIF

!--   Now set-up any remaining parameters that result from other params


!     Check the consitency of a few parameters
      IF ( SEAICE_emissivity .LT. 1.0d-04 ) THEN
       WRITE(mitice_runmsg_unit,'(2A)')                              &
           'SEAICE_emissivity is no longer emissivity*(boltzmann ',  &
           'constant) but really an emissivity.'
       WRITE(mitice_runmsg_unit,'(2A)')                              &
           'Typical values are near 1 ',                             &
           '(default is 5.5/5.67=0.9700176...).'
       WRITE(mitice_runmsg_unit,'(A,E13.6,A)')                       &
           'Please change SEAICE_emissivity in data.seaice to ',     &
           SEAICE_emissivity, '/5.67e-8.'
       STOP 'ABNORMAL END: S/R SEAICE_READPARMS'
      ENDIF

!--   Since the default of SEAICE_waterDrag has changed, issue a warning
!     in case of large values
      chkFlag = .FALSE.
      IF ( SEAICE_waterDrag .GT. 1.0_wp ) THEN
       chkFlag = .TRUE.
      ENDIF
      IF ( SEAICE_waterDrag_South .GT. 1.0_wp ) THEN
       chkFlag = .TRUE.
      ENDIF
      IF ( chkFlag ) THEN
       WRITE(mitice_runmsg_unit,'(3A)') '** WARNING ** SEAICE_READPARMS: ', &
            'That is 3 orders of magnitude larger',                         &
            ' than the default of 5.5e-3.'
      ENDIF

!YY: Pack seaice parameters and broadcast to other threads
      CALL mitice_mpi_namelist_pack_and_bcast

     else     ! elseif for mpirank==0
     
      CALL mitice_mpi_namelist_unpack_and_bcast
     
     endif   !endif for mpi_rank==0 

    end subroutine mitice_read_params


    subroutine mitice_init_fixed
!     *==========================================================*
! Called from MaCOM main program
! Initialize the sea ice model of MITgcm
! TO DO : mpi pack&bcast and unpack
!     *==========================================================*
    implicit none

#ifdef SEAICE_ITD
      character(lc) :: msgBuf
!     k - loop counter for ITD categories
      integer k
      real(wp)  ::  tmpVar
      logical computeHlimit
#endif
      integer i
      integer kSrf


!   determine the fraction of sw radiation penetrating the model shallowest layer
      real(wp) swfracba(2)
      real(wp) tmpFac
!     surface layer thickness in meters
      real(wp) dzSurf

!YY: rC_z(kSrf) is surface height in meter
      kSrf = nk

!YY: only master thread set parameters, then pack and broadcast by MPI
      if (mpi_rank == 0) then

!SHORTWAVE_HEATING
       tmpFac    = -1.0
       swfracba(1) = ABS(rF_z(nkp1))    !YY, ZY: using rF_z, which is not available in MaCOM now
       swfracba(2) = ABS(rF_z(nk))      !YY: rf_z is negative?? see SWFRAC
       CALL SWFRAC( 2, tmpFac, swfracba )
       SWFracB = swfracba(2)

!--   Set mcPheePiston coeff (if still unset)
       dzSurf = drF(kSrf)
       ! YY,  ZY:check, divided by gravityRau0
       IF ( .not.ln_bous )dzSurf = drF(kSrf) / gravityRau0
      IF ( SEAICE_mcPheePiston.EQ.UNSET_RL ) THEN
        IF ( SEAICE_availHeatFrac.NE.UNSET_RL ) THEN
          SEAICE_mcPheePiston = SEAICE_availHeatFrac * dzSurf/SEAICE_deltaTtherm
        ELSE
          SEAICE_mcPheePiston = MCPHEE_TAPER_FAC * STANTON_NUMBER * USTAR_BASE
          SEAICE_mcPheePiston = MIN( SEAICE_mcPheePiston, dzSurf/SEAICE_deltaTtherm )
        ENDIF
      ENDIF


#ifdef SEAICE_ITD
!     zeroth category needs to be zero
      Hlimit(0)    = 0.0_wp
!     thickest category is unlimited
      Hlimit(nITD) = 999.9_wp
      computeHlimit=.FALSE.
!     check if Hlimit contains useful values
      DO k = 1, nITD
       IF ( Hlimit(k).EQ.UNSET_RL )    computeHlimit=.TRUE.
       IF ( Hlimit(k).LE.Hlimit(k-1) ) computeHlimit=.TRUE.
      ENDDO
      IF ( computeHlimit ) THEN
       WRITE(mitice_runmsg_unit,'(A,I2,A)') 'SEAICE_INIT_FIXED: Computing ',  &
            nITD, ' thickness category limits.'
!     use Equ. 22 of Lipscomb et al. (2001, JGR) to generate ice
!     thickness category limits:
!     - dependends on given number of categories nITD
!     - choose between original parameters of Lipscomb et al. (2001):
!       c1=3.0/N, c2=15*c1, c3=3.0
!       or emphasize thin end of ITD (in order to enhance ice growth):
!       c1=1.5/N, c2=42*c1, c3=3.3
!       -> HINT: set parameters c1, c2 and c3 in seaice_readparms.F
       IF ( nITD.GT.1 ) THEN
        tmpVar = nITD
        tmpVar = ONE / tmpVar
        Hlimit_c1 = Hlimit_c1*tmpVar
        Hlimit_c2 = Hlimit_c2*Hlimit_c1
        DO k=1,nITD-1
         Hlimit(k) = Hlimit(k-1) + Hlimit_c1 + Hlimit_c2           &
          *( ONE + TANH( Hlimit_c3 *( FLOAT(k-1)*tmpVar - ONE ) ) )
        ENDDO
       ENDIF
      ENDIF
!     thickest category is unlimited
      Hlimit(nITD) = 999.9_wp

#endif /* SEAICE_ITD */

! MPI pack and broadcast
      CALL mitice_mpi_initfixed_pack_and_bcast
    else     ! elseif for mpirank==0
      CALL mitice_mpi_initfixed_unpack_and_bcast
    endif    ! endif for mpirank==0

!--   set the first seaice computation time
      NxtIceTime = int(SEAICE_deltaTdyn, 8)  !first dyn computation time
      NxtTHDTime = int(SEAICE_deltaTtherm, 8)  !first thermdyn computation time

!--   Summarise seaice configuration
      if (mpi_rank == 0) call seaice_summary
      

!-    convert SEAICE_doOpenWaterGrowth/Melt logical switch to numerical
!     facOpenGrow/facOpenMelt
      facOpenGrow = 0.0_wp
      facOpenMelt = 0.0_wp
      IF (SEAICE_doOpenWaterGrowth) facOpenGrow = 1.0_wp
      IF (SEAICE_doOpenWaterMelt)   facOpenMelt = 1.0_wp

    end subroutine mitice_init_fixed


    subroutine mitice_init_allocate

    implicit none

       !Global vars -- masks
       !YUAN, ZY: for vars defined at Z-point, their dimension should be nlpbz
       !??
       allocate(HEFFM(nlpb), SIMaskU(nlpb),SIMaskV(nlpb)) 
       allocate(seaiceMaskU(nlpb),seaiceMaskV(nlpb))
       !Global vars -- dynamic
       allocate(uice(nlpb), vice(nlpb))
       allocate(area(nlpb), heff(nlpb), hsnow(nlpb))
#if defined SEAICE_ITD
       allocate(areaitd(nlpb,nITD), heffitd(nlpb,nITD), hsnowitd(nlpb,nITD))
       allocate(opnWtrFrac(nlpb),fw2ObyRidge(nlpb) )
       !allocate(Hlimit(0:nITD))
       call seaice_itd_allocate
#endif
       allocate(TICES(nlpb,SEAICE_multDim))    !seaice_multDim is 1 if SEAICE_ITD not defined

       allocate(PRESS0(nlpb), ZMAX(nlpb), ZMIN(nlpb))
       allocate(uIceNm1(nlpb), vIceNm1(nlpb))
       allocate(e11(nlpb),e22(nlpb),e12(nlpbz))
       !allocate(etaZ(nlpb),zetaZ(nlpb))
       allocate(etaZ(nlpbz),zetaZ(nlpbz))
       !allocate(eta(nlpb),zeta(nlpb))
       allocate(deltaC(nlpb))
       allocate(tensileStrFac(nlpb))
       allocate(FORCEX0(nlpb), FORCEY0(nlpb))
       allocate(FORCEX(nlpb), FORCEY(nlpb))
       allocate(DWATN(nlpb))
       !YY:seaiceMass only used in dynamic module
       allocate(seaiceMassC(nlpb), seaiceMassU(nlpb),seaiceMassV(nlpb))

       allocate(stressDivergenceX(nlpb), stressDivergenceY(nlpb))

       allocate(d_HEFFbyNEG(nlpb), d_HSNWbyNEG(nlpb))

       if(SEAICEuseEVP) then    !YY: useEVP must be read in as namelist, because init of var is hereafter
         allocate(seaice_sigma1(nlpb), seaice_sigma2(nlpb), seaice_sigma12(nlpbz))
       endif
#ifdef SeaiceDebug
      allocate(tmp_outvar(nlpb,nITD))
      tmp_outvar = 0.0
#endif

#ifdef SEAICE_VARIABLE_SALINITY
       allocate(HSALT(nlpb), saltFluxAdjust(nlpb))
#endif

#if defined SEAICE_ALLOW_FREEDRIFT
       allocate(uice_fd(nlpb), vice_fd(nlpb))
#endif

       allocate( CbotC(nlpb) )
       allocate(sIceLoad(nlpb))
       allocate(utau_ice(nlpb),vtau_ice(nlpb))
       allocate(qns_ice(nlpb),qsr_ice(nlpb),EmPmR_ice(nlpb))
!
    end subroutine mitice_init_allocate


    subroutine mitice_init_vars
!   *==========================================================*
!    Called from MaCOM main program
!    Set initial conditions for dynamic variables and time-dependent arrays
!    TO DO : 
!   *==========================================================*
      implicit none

      INTEGER i
      integer w, s
      INTEGER kSrf
      real(sp)  mask_uice
      INTEGER k,ibot   !ibot is bottom index
      real(wp) recip_tensilDepth

      !YY:in MaCOM, kTopC=nkp1, rf(nkp1) is sea surface interface, rC(nk) is depth in the middle
       kSrf = nk

!$acc update device(maskC,maskW,maskS,tFld,gravityRau0,rF,kSurfC,nk)

!$ACC DATA PRESENT(area,heff,hsnow,uice,vice,seaiceMaskU,seaiceMaskV, &
#ifdef SEAICE_ITD
      !$ACC  areaitd,heffitd,hsnowitd,opnWtrFrac,fw2ObyRidge,         &
#endif
#ifdef SEAICE_ALLOW_FREEDRIFT
      !$ACC  uice_fd,vice_fd,                                         &
#endif
      !$ACC  uIceNm1,vIceNm1,e11,e22,e12,deltaC,DWATN,CbotC,          &
      !$ACC  etaZ,zetaZ,tensileStrFac,ZMAX,ZMIN,PRESS0,               &
      !$ACC  stressDivergenceX,stressDivergenceY,                     &
      !$ACC  seaice_sigma1,seaice_sigma2,seaice_sigma12,              &
      !$ACC  FORCEX,FORCEY,FORCEX0,FORCEY0,                           &
      !$ACC  seaiceMassC,seaiceMassU,seaiceMassV,TICES,sIceLoad,      &
#ifdef SEAICE_VARIABLE_SALINITY
      !$ACC  HSALT,       &
#endif
      !$ACC  HEFFM, SIMaskU,SIMaskV),                                 &
      !$ACC  COPYIN(kSrf) 


!--   Initialise grid parameters that do no change during the integration
!--   Initialise grid info
!$acc kernels
!$acc loop
      do i = 1,nlpb
        !YY, ZY: maskC in MaCOM is type integer. Is type conversion necessary ?
        !HEFFM  (i) = maskC(i,kSrf)
        !SIMaskU(i) = maskW(i,kSrf)
        !SIMaskV(i) = maskS(i,kSrf)
        !YY: or
        HEFFM  (i) = real( maskC(i,kSrf),8 )
        SIMaskU(i) = real( maskW(i,kSrf),8 )
        SIMaskV(i) = real( maskS(i,kSrf),8 )
      enddo
!--   Initialise all variables
!$acc loop
      DO i = 1,nlpb
        HEFF(i) = 0.0_wp
        AREA(i) = 0.0_wp
        HSNOW(i)  = 0.0_wp
        UICE(i)= 0.0_wp
        VICE(i)= 0.0_wp
#ifdef SEAICE_ALLOW_FREEDRIFT
        uice_fd(i)= 0.0_wp
        vice_fd(i)= 0.0_wp
#endif
!
        uIceNm1(i) = 0.0_wp
        vIceNm1(i) = 0.0_wp
        e11    (i) = 0.0_wp
        e22    (i) = 0.0_wp
        deltaC (i) = 0.0_wp
        !ETA    (i) = 0.0_wp
        !ZETA   (i) = 0.0_wp
        FORCEX (i) = 0.0_wp
        FORCEY (i) = 0.0_wp
        seaiceMassC(i)= 0.0_wp
        seaiceMassU(i)= 0.0_wp
        seaiceMassV(i)= 0.0_wp
        stressDivergenceX(i) = 0.0_wp
        stressDivergenceY(i) = 0.0_wp
        if(SEAICEuseEVP) then
          seaice_sigma1 (i) = 0.0_wp
          seaice_sigma2 (i) = 0.0_wp
        endif
        DWATN(i)   = 0.0_wp
        CbotC(i)   = 0.0_wp
        PRESS0(i)  = 0.0_wp
        !PRESS (i)  = 0.0_wp
        FORCEX0(i) = 0.0_wp
        FORCEY0(i) = 0.0_wp
        ZMAX(i)    = 0.0_wp
        ZMIN(i)    = 0.0_wp
        tensileStrFac(i) = 0.0_wp
#ifdef SEAICE_VARIABLE_SALINITY
        HSALT(i)  = 0.0_wp
#endif
        qsr_ice(i)   = 0.0_wp    ! temporary array
        qns_ice(i)   = 0.0_wp
        EmPmR_ice(i) = 0.0_wp
        utau_ice(i)  = 0.0_wp
        vtau_ice(i)  = 0.0_wp
      enddo

!$acc loop collapse(2)
      do k=1,nITD
        do i = 1,nlpb
         TICES(i,k)=0.0_wp
#ifdef SEAICE_ITD
         AREAITD(i,k) = 0.0_wp
         HEFFITD(i,k) = 0.0_wp
         HSNOWITD(i,k)= 0.0_wp
#endif
        enddo
      enddo

! for vars at Z-point, their dimension are nlpbz
!$acc loop
       do i = 1,nlpbz
        if(SEAICEuseEVP) then
          seaice_sigma12(i) = 0.0_wp
        endif
        etaZ   (i) = 0.0_wp
        zetaZ  (i) = 0.0_wp
        e12    (i) = 0.0_wp
       enddo
!$acc end kernels

      
      IF ( SEAICE_taveFreq.GT.0.0_wp ) THEN
!     Initialize averages to zero
        call mitice_ave_allocate
        call mitice_ave_reset
        NxtAveTime = int(SEAICE_taveFreq, 8)
      ENDIF

!YY: seaiceMask is initialized as land-sea mask. however, it is dynamic mask with regards to AREA
!$acc kernels
!$acc loop private(w,s,mask_uice)
      DO i=1,nlpb
        w = tw(i)
        s = ts(i)
        seaiceMaskU(i)=   0.0_wp
        seaiceMaskV(i)=   0.0_wp
        mask_uice=HEFFM(i)+HEFFM(w)
        IF(mask_uice.GT.1.5_wp) seaiceMaskU(i)=1.0_wp
        mask_uice=HEFFM(i)+HEFFM(s)
        IF(mask_uice.GT.1.5_wp) seaiceMaskV(i)=1.0_wp
      ENDDO


!!YY: keep it here in case of regional modeling capability
!!YY: TO DO LATER
!#ifdef OBCS_UVICE_OLD
!#ifdef SEAICE_CGRID
!        IF (useOBCS) THEN
!C--   If OBCS is turned on, close southern and western boundaries
!         DO i=1-OLx,sNx+OLx
!C Southern boundary
!          J_obc = OB_Js(i,bi,bj)
!          IF ( J_obc.NE.OB_indexNone ) THEN
!           seaiceMaskU(i,J_obc,bi,bj)=   0.0 _d 0
!           seaiceMaskV(i,J_obc,bi,bj)=   0.0 _d 0
!          ENDIF
!         ENDDO
!         DO j=1-OLy,sNy+OLy
!C Western boundary
!          I_obc = OB_Iw(j,bi,bj)
!          IF ( I_obc.NE.OB_indexNone ) THEN
!           seaiceMaskU(I_obc,j,bi,bj)=   0.0 _d 0
!           seaiceMaskV(I_obc,j,bi,bj)=   0.0 _d 0
!          ENDIF
!         ENDDO
!        ENDIF
!#endif /* SEAICE_CGRID */
!#endif /* OBCS_UVICE_OLD */

!$acc loop collapse(2)
      DO k=1,nITD
        DO i=1,nlpb
         TICES(i,k)=273.0_wp
        ENDDO
      ENDDO
!$acc loop
      DO i=1,nlpb
        seaiceMassC(i)=1000.0_wp
        seaiceMassU(i)=1000.0_wp
        seaiceMassV(i)=1000.0_wp
      ENDDO
!$acc end kernels

!YY: not necessary
!--   Update overlap regions
!      CALL mpi_data_exchange(mpi_sendbuf_1d, seaiceMaskU)
!      CALL mpi_data_exchange(mpi_sendbuf_1d, seaiceMaskV)

      if (debuglevel.ge.1)then
        !call mpi_debug_output_1d('HEFFM', HEFFM)
        !call mpi_debug_output_1d('seaiceMaskU', seaiceMaskU)
        !call mpi_debug_output_1d('seaiceMaskV', seaiceMaskV)
      endif

! read initial fields from files
!NOTE: hotstart fields are incorporated in MaCOM' restart_in code
      if (.not. restart_in) then
        write(mitice_runmsg_unit, *) 'READ INITIAL FIELD'
!$acc kernels loop present(SEAICE_initialHEFF)
        do i = 1,nlpb
          HEFF(i)=SEAICE_initialHEFF*HEFFM(i)
          UICE(i)=ZERO
          VICE(i)=ZERO
        enddo
!$acc end kernels

!--   Read initial sea-ice velocity from file (if available)
!YY: only support NC files
!YY, ZY: make initial field using MaCOM's method. TO DO BY ZY
        IF ( uIceFile .NE. ' ' )                                                &
          CALL mpi_netcdf_read_exchange(trim(adjustl(uIceFile)), 'uIce', uIce, 1,optional_ts = 1)
        IF ( vIceFile .NE. ' ' )                                                &
          CALL mpi_netcdf_read_exchange(trim(adjustl(vIceFile)), 'vIce', vIce, 1,optional_ts = 1)
        CALL mpi_data_exchange(mpi_sendbuf_1d, uIce)
        CALL mpi_data_exchange(mpi_sendbuf_1d, vIce)
!$acc update device(uIce,vIce)
        IF ( uIceFile .NE. ' ' .AND. vIceFile .NE. ' ' ) THEN
!$acc kernels loop
          do i = 1,nlpb
            uIce(i) = uIce(i)*seaiceMaskU(i)
            vIce(i) = vIce(i)*seaiceMaskV(i)
          enddo
!$acc end kernels
          CALL mpi_data_exchange(mpi_sendbuf_1d, uIce)
          CALL mpi_data_exchange(mpi_sendbuf_1d, vIce)
        ENDIF

!--   Read initial sea-ice thickness from file if available.
       IF ( HeffFile .NE. ' ' ) THEN
        CALL mpi_netcdf_read_exchange(trim(adjustl(HeffFile)), 'HEFF', HEFF, 1,optional_ts = 1)
        !CALL mpi_data_exchange(mpi_sendbuf_1d, HEFF)
!$acc update device(HEFF)
!!$acc kernels loop
!        do i = 1,nlpb
!          HEFF(i) = MAX(HEFF(i),ZERO)
!        enddo
!!$acc end kernels
       ENDIF
!$acc kernels loop
       do i = 1,nlpb
         IF(HEFF(i).GT.ZERO) AREA(i)=ONE
       enddo
!$acc end kernels

!--   Read initial sea-ice area from file if available.
       IF ( AreaFile .NE. ' ' ) THEN
        CALL mpi_netcdf_read_exchange(trim(adjustl(AreaFile)), 'AREA', AREA, 1,optional_ts = 1)
        !CALL mpi_data_exchange(mpi_sendbuf_1d, AREA)
!$acc update device(AREA)
!$acc kernels loop
        do i = 1,nlpb
          IF ( AREA(i) .LE. ZERO ) HEFF(i) = ZERO
          AREA(i) = MAX(AREA(i),ZERO)
          AREA(i) = MIN(AREA(i),ONE)
          IF ( HEFF(i) .LE. ZERO ) AREA(i) = ZERO
          IF ( HEFF(i) .LE. ZERO ) HEFF(i) = ZERO
        enddo 
!$acc end kernels
       ENDIF

!$acc kernels loop
!YY: is this reasonable??
       do i = 1,nlpb
         HSNOW(i) = 0.2_wp * AREA(i)
       enddo
!$acc end kernels

!--   Read initial snow thickness from file if available.
       IF ( HsnowFile .NE. ' ' ) THEN
        CALL mpi_netcdf_read_exchange(trim(adjustl(HsnowFile)), 'HSNOW', HSNOW, 1,optional_ts = 1)
        CALL mpi_data_exchange(mpi_sendbuf_1d, HSNOW)
!$acc kernels loop
        do i = 1,nlpb
          HSNOW(i) = MAX(HSNOW(i),ZERO)
        enddo
!$acc end kernels
       ENDIF

#ifdef SEAICE_ITD
!$acc kernels loop
       do i = 1,nlpb
         AREAITD(i,1)   = AREA(i)
         HEFFITD(i,1)   = HEFF(i)
         HSNOWITD(i,1)  = HSNOW(i)
         opnWtrFrac(i)  = 1.0_wp - AREA(i)
         fw2ObyRidge(i) = 0.0_wp
       enddo
!$acc end kernels
       CALL seaice_itd_redist
#endif

#ifdef SEAICE_VARIABLE_SALINITY
!$acc kernels loop
       do i = 1,nlpb
         HSALT(i)=HEFF(i)* tFld(i,kSrf,2) * SEAICE_rhoIce*SEAICE_saltFrac
       enddo
!$acc end kernels
!--   Read initial sea ice salinity from file if available.
       IF ( HsaltFile .NE. ' ' ) THEN
        CALL mpi_netcdf_read_exchange(trim(adjustl(HsaltFile)), 'HSALT', HSALT, 1,optional_ts = 1)
        CALL mpi_data_exchange(mpi_sendbuf_1d, HSALT)
       ENDIF
#endif /* SEAICE_VARIABLE_SALINITY */

      ENDIF   !YY: endif for initial field 

!---  Complete initialization
!$acc kernels loop
      do i = 1,nlpb
        !ZETA(i)   = HEFF(i)*(1.0d11)
        !ETA (i)   = ZETA(i)/SEAICE_eccen**2
        PRESS0(i) = SEAICE_strength*HEFF(i)*EXP(-SEAICE_cStar*(ONE-AREA(i)))
        ZMAX(i)   = SEAICE_zetaMaxFac*PRESS0(i)
        ZMIN(i)   = SEAICE_zetaMin
        PRESS0(i) = PRESS0(i)*HEFFM(i)
      enddo
!$acc end kernels
      !YY,ZY: if keep useRealFreshWaterFlux switch, then need include it in MaCOM.
      !SO NEED TO remove this switch later. default is admitting real fw flux
      !sIceLoad is used in mod_csp_dyn_force
      IF ( useRealFreshWaterFlux ) THEN
!$acc kernels loop
        do i = 1,nlpb
          sIceLoad(i) = HEFF(i)*SEAICE_rhoIce + HSNOW(i)*SEAICE_rhoSnow
        enddo
!$acc end kernels
      ENDIF

!     compute tensile strength factor k: tensileStrength = k*PRESS
!     can be done in initialisation phase as long as it depends only
!     on depth
      IF ( SEAICE_tensilFac .NE. 0.0_wp ) THEN
       recip_tensilDepth = 0.0_wp
       IF ( SEAICE_tensilDepth .GT. 0.0_wp )              &
         recip_tensilDepth = 1.0_wp / SEAICE_tensilDepth
!$acc kernels copyin(recip_tensilDepth)
!$acc loop private(ibot)
       do i = 1,nlpb
!!YY: here ocean bottom r_Low is replaced by rF_z(kSurfC), rF_z is not present in MaCOM
!!YY: rF_z is better , right?
!ZY to check
         ibot = min(nk, kSurfC(i))    !YY: refer to mod_csp_vdiff_coef, divided by rho*g??
         tensileStrFac(i) = SEAICE_tensilFac*HEFFM(i)     &
               *exp(-ABS( rF(ibot)/gravityRau0 )*recip_tensilDepth)
       enddo
!$acc end kernels
      ENDIF

!$acc end data

!! this is to dump initial field as NC format, ts_dump = 1
!YY: if hot start, then the first frame should be in restart_in read subroutine
      IF (SEAICEwriteState .and. .not.restart_in) THEN     !default is false. set in init.
        write(prefix, '(A,i4.4,i2.2,i2.2,i2.2,i2.2,i2.2)')   &
          trim(adjustl(foldername))//'dump_',nowyear,nowmon,nowday,nowhour,nowmin,nowsec
        call seaice_io_dump(prefix,1)
        call seaice_io_init
        NxtDumpTime = int(SEAICE_dumpFreq, 8)   ! put in the end
      ENDIF

    end subroutine mitice_init_vars

!==============================================================================
    subroutine gpu_seaice_parameters_update
!==============================================================================
!VECTOR VARS
!      !$ACC UPDATE DEVICE(                                                       &
!      !$ACC  HEFFM, SIMaskU,SIMaskV,seaiceMaskU,seaiceMaskV,                     &
!      !$ACC  uice,vice,area,heff,hsnow,                                          &
!#ifdef SEAICE_ITD
!      !$ACC  areaitd,heffitd,hsnowitd,opnWtrFrac,fw2ObyRidge,                    &
!#endif
!      !$ACC  TICES,PRESSO,ZMAX,ZMIN,uIceNm1,vIceNm1,e11,e22,e12,                 &
!      !$ACC  etaZ,zetaZ,deltaC,tensileStrFac,FORCEX0,FORCEY0,FORCEX,FORCEY,      &
!      !$ACC  DWATN,seaiceMassC,seaiceMassU,seaiceMassV,                          &
!      !$ACC  stressDivergenceX,stressDivergenceY,d_HEFFbyNEG,d_HSNWbyNEG,        &
!      !$ACC  seaice_sigma1,seaice_sigma2,seaice_sigma12,                         &
!#ifdef SEAICE_VARIABLE_SALINITY
!      !$ACC  HSALT,saltFluxAdjust,                                               &
!#endif
!#if defined SEAICE_ALLOW_FREEDRIFT
!      !$ACC  uice_fd,vice_fd,                                                    &
!#endif
!      !$ACC CbotC,sIceLoad,fCori )


       !$ACC UPDATE DEVICE(fCori)

!SCALAR VARS AND PARAMETERS

      !these parameters are initialized in subroutine init_fixed

      !$ACC UPDATE DEVICE(                                                       &
#ifdef SEAICE_ITD
      !$ACC  Hlimit,      &
#endif
      !$ACC  SWFracB,SEAICE_mcPheePiston,facOpenGrow,facOpenMelt)

      !Integer or real parameters read from namelist
      !$ACC UPDATE DEVICE(   &
      !$ACC     SEAICE_mcPheeStepFunc,SEAICEheatConsFix,SEAICE_useMultDimSnow,           &
      !$ACC     SEAICEdiffKhHeff,SEAICEdiffKhSnow,SEAICEdiffKhArea,SEAICEdiffKhSalt,     &
      !$ACC     SEAICE_deltaTtherm,SEAICE_deltaTdyn,SEAICE_deltaTevp,                    &
      !$ACC     SEAICE_evpTauRelax,SEAICE_evpAlpha,SEAICE_evpBeta,                       &
      !$ACC     SEAICEaEVPalphaMin,SEAICE_zetaMin,SEAICE_zetaMaxFac,                     &
      !$ACC     SEAICEsnowFracRidge,SEAICEgStar,SEAICEhStar,                             &
      !$ACC     SEAICEaStar,SEAICEshearParm,SEAICEmuRidging,SEAICEmaxRaft,               &
      !$ACC     SEAICEpresH0,SEAICEpresPow0,SEAICEpresPow1,SEAICE_initialHEFF,           &
      !$ACC     SEAICE_areaGainFormula,SEAICE_areaLossFormula,                           &
      !$ACC     SEAICE_rhoAir,SEAICE_rhoIce,SEAICE_rhoSnow,ICE2WATR,                     &
      !$ACC     SEAICE_drag,SEAICE_waterDrag,SEAICEdWatMin,                              &
      !$ACC     SEAICE_dryIceAlb,SEAICE_wetIceAlb,SEAICE_drySnowAlb,SEAICE_wetSnowAlb,   &
      !$ACC     HO,SEAICE_drag_south,SEAICE_waterDrag_south,                             &
      !$ACC     SEAICE_dryIceAlb_south,SEAICE_wetIceAlb_south,SEAICE_drySnowAlb_south,   &
      !$ACC     SEAICE_wetSnowAlb_south,HO_south,SEAICE_cBasalStar,                      &
      !$ACC     SEAICEbasalDragU0,SEAICEbasalDragK1,SEAICEbasalDragK2,SEAICE_wetAlbTemp, &
      !$ACC     SEAICE_strength,SEAICE_cStar,SEAICEpressReplFac,SEAICE_lhFusion,         &
      !$ACC     SEAICE_ice_emiss,SEAICE_snow_emiss,SEAICE_tempFrz0,SEAICE_dTempFrz_dS,   &
      !$ACC     SEAICE_salt0,SEAICE_saltFrac,SEAICEstressFactor,SEAICE_frazilFrac,       &
      !$ACC     SEAICE_PDF,SEAICE_multDim,SEAICE_deltaMin,SEAICE_tensilFac,              &
      !$ACC     SEAICE_area_reg,SEAICE_hice_reg,SEAICE_area_floor, SEAICE_area_max,      &
      !$ACC     MIN_ATEMP,MIN_LWDOWN,MIN_TICE,SEAICE_EPS,SEAICE_EPS_SQ)

      !these pars may not necessary upload to device
      !$ACC UPDATE DEVICE(SEAICEnEVPstarSteps,SEAICEridgingIterMax,SEAICEaEVPcoeff, &
      !$ACC           SEAICEaEVPcStar,SEAICE_cf,SEAICE_cpAir, ICE2WATR,             &
      !$ACC           SEAICE_lhEvap,SEAICE_dalton,SEAICE_iceConduct,                &
      !$ACC           SEAICE_snowConduct,SEAICE_emissivity,SEAICE_snowThick,        &
      !$ACC           SEAICE_shortwave,OCEAN_drag,SEAICE_mcPheeTaper,               &
      !$ACC           SEAICE_availHeatFrac,SEAICE_tauAreaObsRelax,                  &
      !$ACC           SEAICE_airTurnAngle,SEAICE_waterTurnAngle,                    &
      !$ACC           SEAICE_monFreq,SEAICE_dumpFreq,SEAICE_taveFreq,               &
      !$ACC           SEAICE_tensilDepth,SEAICE_eccen)

    end subroutine gpu_seaice_parameters_update

end module mitice_init
