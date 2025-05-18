module mitice_parameters
!  *==========================================================*
!  | o Basic parameter header for sea ice model.
!  *==========================================================*

! modules from MaCOM
use mod_misc_basic
use mod_csp_basic

implicit none

!-- Time
   integer(i8), save, public :: NxtIceTime    !time in second to next seaice dyn calling
   integer(i8), save, public :: NxtTHDTime    !time in second to next thermdyn calling
   integer(i8), save, public :: NxtDumpTime    !time in second to next dump calling
   integer(i8), save, public :: NxtAveTime    !time in second to next average calling
   logical,     public, save :: NxtDumpFile = .false.         !.true. when seaice output time enter a new month, see csp_init_step
!-- seaice load 
   real(wp), parameter  ::  sIceLoadFac = 1.0_wp !factor to scale (and turn off) sIceLoad (sea-ice loading)
!!YY: these are logical switches for MaCOM. consider move them to MaCOM, or turn off
   logical, public, parameter :: useRealFreshWaterFlux=.True.

!--   Sea ice array-size definition, 7 or 5
      integer, public, parameter ::  nITD=7   !number of seaice thickness categories

! - dynamics:
!     SEAICEuseDYNAMICS :: If false, do not use dynamics;
!                          default is to use dynamics.
!     SEAICEuseFREEDRIFT :: If True use free drift velocity instead of EVP
!                           or LSR
!     SEAICEuseEVP      :: If true use elastic viscous plastic solver
!     SEAICEuseEVPstar  :: If true use modified elastic viscous plastic
!                          solver (EVP*) by Lemieux et al (2012)
!     SEAICEuseEVPrev   :: If true use "revisited" elastic viscous plastic
!                          solver following Bouillon et al. (2013), very similar
!                          to EVP*, but uses fewer implicit terms and drops
!                          one 1/e^2 in equations for sigma2 and sigma12
!     SEAICEuseEVPpickup :: Set to false in order to start EVP solver with
!                          non-EVP pickup files.  Default is true.
!                          Applied only if SEAICEuseEVP=.TRUE.
!     SEAICEuseTEM      :: to use the truncated ellipse method (see Geiger et al.
!                          1998) set this parameter to true, default is false
!     SEAICEuseTilt     :: If true then include surface tilt term in dynamics
!     SEAICE_no_slip    :: apply no slip boundary conditions to seaice velocity
!     SEAICE_clipVelocities :: clip velocities to +/- 40cm/s
!     SEAICEaddSnowMass :: in computing seaiceMass, add snow contribution
!                          default is .TRUE.
!     useHB87stressCoupling :: use an intergral over ice and ocean surface
!                          layer to define surface stresses on ocean
!                          following Hibler and Bryan (1987, JPO)
!     useHibler79IceStrength :: if true original ice strength parameterization
!                          other use Rothrock (1975) parameterization based
!                          on energetics and an ice thickness distribution
!                          (default = .true.)
!     SEAICEscaleSurfStress :: if TRUE, scale ice-ocean and ice-atmosphere
!                          stress on ice by concenration (AREA) following
!                          Connolley et al. (2004), JPO. (default = .TRUE.)
!     SEAICEsimpleRidging :: use Hibler(1979) ridging (default=.true.)
!     SEAICEuseLinRemapITD :: use linear remapping (Lipscomb et al. 2001)
!                             .TRUE. by default
! - advection:
!     SEAICEadvHeff     :: turn on advection of effective thickness
!                          (default = .true.)
!     SEAICEadvArea     :: turn on advection of fraction area
!                          (default = .true.)
!     SEAICEadvSnow     :: turn on advection of snow (does not work with
!                          non-default Leap-frog scheme for advection)
!     SEAICEadvSalt     :: turn on advection of salt (does not work with
!                          non-default Leap-frog scheme for advection)
!     SEAICEmomAdvection:: turn on advection of momentum (default = .false.)
!     SEAICEuseAbsVorticity    :: useAbsVorticity, useJamartMomAdv for vector
!     SEAICEuseJamartMomAdv    :: invariant momentum in the ocean
! - thermodynamics:
!     usePW79thermodynamics :: use "0-layer" thermodynamics as described in
!                           Parkinson and Washington (1979) and Hibler (1979)
!     SEAICE_useMultDimSnow :: use same fixed pdf for snow as for
!                              multi-thickness-category ice (default=.TRUE.)
!     SEAICEuseFlooding  :: turn on scheme to convert submerged snow into ice
!     SEAICEheatConsFix  :: If true then fix ocn<->seaice advective heat flux.
!     useMaykutSatVapPoly :: use Maykut Polynomial for saturation vapor pressure
!                         instead of extended temp-range exponential law; def=F.
!     SEAICE_mcPheeStepFunc :: use step function (not linear tapering) in
!                           ocean-ice turbulent flux
!     SEAICE_doOpenWaterGrowth :: use open water heat flux directly to grow ice
!                           (when false cool ocean, and grow later if needed)
!     SEAICE_doOpenWaterMelt   :: use open water heat flux directly to melt ice
!                           (when false warm ocean, and melt later if needed)
!     SEAICE_growMeltByConv :: grow/melt according to convergence of turbulence
!                              and conduction, rather than in two steps (default)
! - other (I/O, ...):
!     SEAICEwriteState  :: If true, dump sea ice state to file;
      logical, public, save ::                                         &
          SEAICEuseDYNAMICS, SEAICEuseFREEDRIFT,                       & 
          SEAICEuseEVP, SEAICEuseEVPstar, SEAICEuseEVPrev,             &
          SEAICEuseEVPpickup,                                          &
          useHibler79IceStrength, SEAICEsimpleRidging,                 &
          SEAICEuseLinRemapITD,SEAICEuseTEM, SEAICEuseTilt,            &
          SEAICE_no_slip,SEAICEscaleSurfStress,                        &
          SEAICE_clipVelocities, SEAICEaddSnowMass,                    &
          useHB87stressCoupling,                &
          SEAICEadvHeff, SEAICEadvArea,                                &
          SEAICEadvSnow, SEAICEadvSalt, SEAICEmomAdvection,            &
          usePW79thermodynamics,                                       &
          SEAICE_useMultDimSnow, SEAICEuseFlooding, SEAICEheatConsFix, &
          !useMaykutSatVapPoly,         & 
          SEAICE_mcPheeStepFunc,                  &
          SEAICE_doOpenWaterGrowth, SEAICE_doOpenWaterMelt,            &
          SEAICE_growMeltByConv,SEAICEwriteState


!     IMAX_TICE         :: number of iterations for ice surface temp
!                          (default=10)
!     postSolvTempIter :: select flux calculation after surf. temp solver
!                         iteration
!                         0 = none, i.e., from last iter
!                         1 = use linearized approx (consistent with tsurf
!                             finding)
!                         2 = full non-lin form
!     SEAICEnEVPstarSteps :: number of evp*-steps
!     SEAICEpresPow0    :: HEFF exponent for ice strength below SEAICEpresH0
!     SEAICEpresPow1    :: HEFF exponent for ice strength above SEAICEpresH0
!     rigding parameters (only active when SEAICE_ITD is defined)
!     SEAICEpartFunc    :: =0 use Thorndyke et al (1975) participation function
!                       :: =1 use Lipscomb et al (2007) participation function
!     SEAICEredistFunc  :: =0 assume ridged ice is uniformly distributed
!                             (Hibler, 1980)
!                          =1 Following Lipscomb et al. (2007), ridged ice is
!                             distributed following an exponentially
!                             decaying function
!     SEAICEridgingIterMax :: maximum number of ridging iterations
!     end ridging parameters
!     SEAICE_areaLossFormula :: selects formula for ice cover loss from melt
!                        :: 1=from all but only melt conributions by ATM and OCN
!                        :: 2=from net melt-growth>0 by ATM and OCN
!                        :: 3=from predicted melt by ATM
!     SEAICE_areaGainFormula :: selects formula for ice cover gain from open
!                               water growth
!                        :: 1=from growth by ATM
!                        :: 2=from predicted growth by ATM
!     SEAICEetaZmethod   :: determines how shear-viscosity eta is computed at
!                           Z-points
!                           0=simple averaging from C-points (default and old)
!                           3=weighted averaging of squares of strain rates
!                             (recommended for energy conservation)
!     SEAICE_multDim     :: number of ice categories

      INTEGER IMAX_TICE, postSolvTempIter
      INTEGER SEAICEnEVPstarSteps
      INTEGER SEAICE_areaLossFormula
      INTEGER SEAICE_areaGainFormula
      INTEGER SEAICEetaZmethod
      INTEGER SEAICE_multDim
      INTEGER SEAICEpresPow0, SEAICEpresPow1
      INTEGER SEAICEpartFunc, SEAICEredistFunc
      INTEGER SEAICEridgingIterMax

!     AreaFile          :: File containing initial sea-ice concentration
!     HsnowFile         :: File containing initial snow thickness
!     HsaltFile         :: File containing initial sea ice salt content
!     HeffFile          :: File containing initial sea-ice thickness
!     uIceFile          :: File containing initial sea-ice U comp. velocity
!     vIceFile          :: File containing initial sea-ice V comp. velocity
!        !!! NOTE !!! Initial sea-ice thickness can also be set using
!        SEAICE_initialHEFF below.  But a constant initial condition
!        can mean large artificial fluxes of heat and freshwater in
!        the surface layer during the first model time step.

      CHARACTER(lc) AreaFile
      CHARACTER(lc) HsnowFile
      CHARACTER(lc) HsaltFile
      CHARACTER(lc) HeffFile
      CHARACTER(lc) uIceFile
      CHARACTER(lc) vIceFile


!     SEAICE_deltaTtherm :: Seaice timestep for thermodynamic equations (s)
!     SEAICE_deltaTdyn   :: Seaice timestep for dynamic solver          (s)
!     SEAICE_deltaTevp   :: Seaice timestep for EVP solver              (s)
!     SEAICE_elasticParm :: parameter that sets relaxation timescale
!                           tau = SEAICE_elasticParm * SEAICE_deltaTdyn
!     SEAICE_evpTauRelax :: relaxation timescale tau                    (s)
!     SEAICE_evpDampC    :: evp damping constant (Hunke,JCP,2001)       (kg/m^2)
!     SEAICE_evpAlpha    :: dimensionless parameter 2*evpTauRelax/deltaTevp
!     SEAICE_evpBeta     :: dimensionless parameter deltaTdyn/deltaTevp
!     SEAICEaEVPcoeff    :: main coefficent for adaptive EVP (largest
!                           stabilized frequency)
!     SEAICEaEVPcStar    :: multiple of stabilty factor: alpha*beta=cstar*gamma
!     SEAICEaEVPalphaMin :: lower limit of alpha and beta, regularisation
!                           to prevent singularities of system matrix,
!                           e.g. when ice concentration is too low.
!     SEAICE_zetaMaxFac  :: factor determining the maximum viscosity    (s)
!                          (default = 5.e+12/2.e4 = 2.5e8)
!     SEAICE_zetaMin     :: lower bound for viscosity (default = 0)    (N s/m^2)
!     SEAICEpresH0       :: HEFF threshold for ice strength            (m)
!     SEAICE_monFreq     :: SEAICE monitor frequency.                   (s)
!     SEAICE_dumpFreq    :: SEAICE dump frequency.                      (s)
!     SEAICE_taveFreq    :: SEAICE time-averaging frequency.            (s)
!     SEAICE_initialHEFF :: initial sea-ice thickness                   (m)
!     SEAICE_rhoAir      :: density of air                              (kg/m^3)
!     SEAICE_rhoIce      :: density of sea ice                          (kg/m^3)
!     SEAICE_rhoSnow     :: density of snow                             (kg/m^3)
!     ICE2WATR           :: ratio of sea ice density to water density
!     SEAICE_cpAir       :: specific heat of air                        (J/kg/K)
!
!     OCEAN_drag         :: unitless air-ocean drag coefficient (default 0.001)
!     SEAICE_drag        :: unitless air-ice drag coefficient   (default 0.001)
!     SEAICE_waterDrag   :: unitless water-ice drag coefficient (default 0.0055)
!     SEAICEdWatMin      :: minimum linear water-ice drag applied to DWATN
!                           (default 0.25 m/s)
!
!     SEAICE_dryIceAlb   :: winter albedo
!     SEAICE_wetIceAlb   :: summer albedo
!     SEAICE_drySnowAlb  :: dry snow albedo
!     SEAICE_wetSnowAlb  :: wet snow albedo
!     HO                 :: AKA "lead closing parameter", demarcation thickness
!                           between thin and thick ice. Alternatively, HO (in
!                           meters) can be interpreted as the thickness of ice
!                           formed in open water.
!                           HO is a key ice-growth parameter that determines
!                           the partition between vertical and lateral growth.
!                           The default is 0.5m, increasing this value leads
!                           slower formation of a closed ice cover and thus to
!                           more ice (and thicker) ice, decreasing to faster
!                           formation of a closed ice cover (leads are closing
!                           faster) and thus less (thinner) ice.
!
!     SEAICE_drag_south       :: Southern Ocean SEAICE_drag
!     SEAICE_waterDrag_south  :: Southern Ocean SEAICE_waterDrag
!     SEAICE_dryIceAlb_south  :: Southern Ocean SEAICE_dryIceAlb
!     SEAICE_wetIceAlb_south  :: Southern Ocean SEAICE_wetIceAlb
!     SEAICE_drySnowAlb_south :: Southern Ocean SEAICE_drySnowAlb
!     SEAICE_wetSnowAlb_south :: Southern Ocean SEAICE_wetSnowAlb
!     HO_south                :: Southern Ocean HO
!
!     Parameters for basal drag of grounded ice following
!     Lemieux et al. (2015), doi:10.1002/2014JC010678
!     SEAICE_cBasalStar (default = SEAICE_cStar)
!     SEAICEbasalDragU0 (default = 5e-5)
!     SEAICEbasalDragK1 (default = 8)
!     SEAICEbasalDragK2  :: if > 0, turns on basal drag
!                           (default = 0, Lemieux suggests 15)
!
!     SEAICE_wetAlbTemp  :: Temp (deg.C) above which wet-albedo values are used
!     SEAICE_strength    :: sea-ice strength Pstar
!     SEAICE_cStar       :: sea-ice strength paramter C* (def: 20)
!     SEAICE_tensilFac   :: sea-ice tensile strength factor, values in [0,1]
!     SEAICE_tensilDepth :: crtical depth for sea-ice tensile strength (def 0.)
!     SEAICEpressReplFac :: interpolator between PRESS0 and regularized PRESS
!                           1. (default): pure pressure replace method (PRESS)
!                           0.          : pure Hibler (1979) method (PRESS0)
!     SEAICE_eccen       :: sea-ice eccentricity of the elliptical yield curve
!     SEAICE_eccfr       :: sea-ice eccentricity of the elliptical flow rule
!     SEAICE_lhFusion    :: latent heat of fusion for ice and snow (J/kg)
!     SEAICE_lhEvap      :: latent heat of evaporation for water (J/kg)
!     SEAICE_dalton      :: Dalton number (= sensible heat transfer coefficient)
!     SEAICE_iceConduct  :: sea-ice conductivity
!     SEAICE_snowConduct :: snow conductivity
!     SEAICE_emissivity  :: longwave ocean-surface emissivity (-)
!     SEAICE_ice_emiss   :: longwave ice-surface emissivity (-)
!     SEAICE_snow_emiss  :: longwave snow-surface emissivity (-)
!     SEAICE_snowThick   :: cutoff snow thickness (for snow-albedo)
!     SEAICE_shortwave   :: ice penetration shortwave radiation factor
!     SEAICE_saltFrac    :: salinity of newly formed seaice defined as a
!                           fraction of the ocean surface salinity at the time
!                           of freezing
!     SEAICE_salt0       :: prescribed salinity of seaice (in g/kg).
!     facOpenGrow        :: 0./1. version of logical SEAICE_doOpenWaterGrowth
!     facOpenMelt        :: 0./1. version of logical SEAICE_doOpenWaterMelt
!     SEAICE_mcPheePiston:: ocean-ice turbulent flux "piston velocity" (m/s)
!                           that sets melt efficiency.
!     SEAICE_mcPheeTaper :: tapering down of turbulent flux term with ice
!                           concentration. The 100% cover turb. flux is
!                           multiplied by 1.-SEAICE_mcPheeTaper
!     SEAICE_frazilFrac  :: Fraction of surface level negative heat content
!                           anomalies (relative to the local freezing point)
!                           may contribute as frazil over one time step.
!     SEAICE_tempFrz0    :: sea water freezing point is
!     SEAICE_dTempFrz_dS :: tempFrz = SEAICE_tempFrz0 + salt*SEAICE_dTempFrz_dS
!     SEAICE_PDF         :: prescribed sea-ice distribution within grid box
!     SEAICEstressFactor :: factor by which ice affects wind stress (default=1)
!     SEAICE_deltaMin    :: small number used to reduce singularities of Delta
!     SEAICE_area_max    :: usually set to 1. Seeting areaMax below 1 specifies
!                           the minimun amount of leads (1-areaMax) in the
!                           ice pack.
!     SEAICE_area_floor  :: usually set to 1x10^-5. Specifies a minimun
!                           ice fraction in the ice pack.
!     SEAICE_area_reg    :: usually set to 1x10^-5. Specifies a minimun
!                           ice fraction for the purposes of regularization
!     SEAICE_hice_reg    :: usually set to 5 cm. Specifies a minimun
!                           ice thickness for the purposes of regularization
!     SEAICEdiffKhArea   :: sets the diffusivity for area (m^2/s)
!     SEAICEdiffKhHeff   :: sets the diffusivity for effective thickness (m^2/s)
!     SEAICEdiffKhSnow   :: sets the diffusivity for snow on sea-ice (m^2/s)
!     SEAICEdiffKhSalt   :: sets the diffusivity for sea ice salinity (m^2/s)
!     SEAICE_airTurnAngle   :: turning angles of air-ice interfacial stress
!     SEAICE_waterTurnAngle :: and ice-water interfacial stress (in degrees)
!     SEAICE_tauAreaObsRelax :: Timescale of relaxation to observed
!                               sea ice concentration (s), default=unset
!     ridging parameters (Lipscomb et al, 2007, Bitz et al. 2001):
!     SEAICE_cf       :: ratio of total energy sinks to gravitational sink
!                        (scales ice strength, suggested values: 2 to 17)
!     SEAICEgStar     :: maximum ice concentration that participates in ridging
!     SEAICEhStar     :: empirical thickness (ridging parameter)
!     SEAICEaStar     :: ice concentration parameter similar to gStar for
!                        exponential distribution (Lipscomb et al 2007)
!     SEAICEshearParm :: <=1 reduces amount of energy lost to ridge building
!     SEAICEmuRidging :: tuning parameter similar to hStar for Lipcomb et al
!                        (2007)-scheme
!     SEAICEmaxRaft   :: regularization parameter (default=1)
!     SEAICEsnowFracRidge :: fraction of snow that remains on ridged
!     dumpFileIntv    :: Generate new snapshot output file of ice state vars:  1, monthly; 2,yearly
!
      real(wp), public, save ::  SEAICE_deltaTtherm, SEAICE_deltaTdyn, SEAICE_deltaTevp
      real(wp), public, save ::  SEAICE_monFreq, SEAICE_dumpFreq, SEAICE_taveFreq
      integer,  public, save ::  dumpFileIntv
      real(wp), public, save ::  SEAICE_initialHEFF
      real(wp), public, save ::  SEAICE_rhoAir, SEAICE_rhoIce, SEAICE_rhoSnow, ICE2WATR
      real(wp), public, save ::  SEAICE_cpAir
      real(wp), public, save ::  SEAICE_drag, SEAICE_waterDrag, SEAICEdWatMin
      real(wp), public, save ::  SEAICE_dryIceAlb, SEAICE_wetIceAlb
      real(wp), public, save ::  SEAICE_drySnowAlb, SEAICE_wetSnowAlb, HO
      real(wp), public, save ::  SEAICE_drag_south, SEAICE_waterDrag_south
      real(wp), public, save ::  SEAICE_dryIceAlb_south, SEAICE_wetIceAlb_south
      real(wp), public, save ::  SEAICE_drySnowAlb_south, SEAICE_wetSnowAlb_south, HO_south
      real(wp), public, save ::  SEAICE_cBasalStar, SEAICEbasalDragU0
      real(wp), public, save ::  SEAICEbasalDragK1, SEAICEbasalDragK2
      real(wp), public, save ::  SEAICE_wetAlbTemp
      real(wp), public, save ::  SEAICE_strength, SEAICE_cStar, SEAICEpressReplFac
      real(wp), public, save ::  SEAICE_tensilFac, SEAICE_tensilDepth
      real(wp), public, save ::  SEAICE_eccen
      real(wp), public, save ::  SEAICE_lhFusion, SEAICE_lhEvap
      real(wp), public, save ::  SEAICE_dalton
      real(wp), public, save ::  SEAICE_iceConduct, SEAICE_snowConduct
      real(wp), public, save ::  SEAICE_emissivity, SEAICE_ice_emiss, SEAICE_snow_emiss
      real(wp), public, save ::  SEAICE_snowThick, SEAICE_shortwave
      real(wp), public, save ::  SEAICE_saltFrac, SEAICE_salt0, SEAICEstressFactor
      real(wp), public, save ::  SEAICE_mcPheeTaper, SEAICE_mcPheePiston
      real(wp), public, save ::  SEAICE_frazilFrac, SEAICE_availHeatFrac
      real(wp), public, save ::  facOpenGrow, facOpenMelt
      real(wp), public, save ::  SEAICE_tempFrz0, SEAICE_dTempFrz_dS
      real(wp), public, save ::  SEAICE_PDF(nITD)
      real(wp), public, save ::  OCEAN_drag
      real(wp), public, save ::  SEAICE_deltaMin
      real(wp), public, save ::  SEAICE_area_reg, SEAICE_hice_reg
      real(wp), public, save ::  SEAICE_area_floor, SEAICE_area_max
      real(wp), public, save ::  SEAICE_airTurnAngle, SEAICE_waterTurnAngle
      real(wp), public, save ::  SEAICE_elasticParm, SEAICE_evpTauRelax
      real(wp), public, save ::  SEAICE_evpAlpha, SEAICE_evpBeta
      real(wp), public, save ::  SEAICE_evpDampC, SEAICE_zetaMin, SEAICE_zetaMaxFac
      real(wp), public, save ::  SEAICEaEVPcoeff, SEAICEaEVPcStar, SEAICEaEVPalphaMin
      real(wp), public, save ::  SEAICEpresH0
      real(wp), public, save ::  SEAICEdiffKhArea, SEAICEdiffKhHeff, SEAICEdiffKhSnow
      real(wp), public, save ::  SEAICEdiffKhSalt
      real(wp), public, save ::  SEAICE_tauAreaObsRelax
      real(wp), public, save ::  SEAICEgStar, SEAICEhStar, SEAICEaStar, SEAICEshearParm
      real(wp), public, save ::  SEAICEmuRidging, SEAICEmaxRaft, SEAICE_cf
      real(wp), public, save ::  SEAICEsnowFracRidge


!     MIN_ATEMP         :: minimum air temperature   (deg C)
!     MIN_LWDOWN        :: minimum downward longwave (W/m^2)
!     MIN_TICE          :: minimum ice temperature   (deg C)
!     SEAICE_EPS        :: small number
!     SEAICE_EPS_SQ     :: small number square

      real(wp), public, save ::  MIN_ATEMP, MIN_LWDOWN, MIN_TICE
      real(wp), public, save ::  SEAICE_EPS, SEAICE_EPS_SQ

#ifdef SEAICE_ITD
!     Hlimit            :: ice thickness category limits (m), array of
!                          size nITD+1,  (0:nITD)
!     Hlimit_c1,_c2,_c3 :: coefficients set in seaice_readparams.F to
!                          calculate Hlimit in seaice_init_fixed.F
      real(wp), allocatable, dimension(:), public, save ::  Hlimit
      real(wp), public, save ::  Hlimit_c1, Hlimit_c2, Hlimit_c3
#endif /* SEAICE_ITD */


      !$ACC DECLARE CREATE(  &
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
      !$ACC     MIN_ATEMP,MIN_LWDOWN,MIN_TICE,SEAICE_EPS,SEAICE_EPS_SQ,                  &
      !$ACC     SEAICE_mcPheePiston,facOpenGrow, facOpenMelt)

#ifdef SEAICE_ITD
      !$ACC DECLARE CREATE(Hlimit)
#endif

      !these pars may not necessary upload to device
      !$ACC DECLARE CREATE(SEAICEnEVPstarSteps,SEAICEridgingIterMax,SEAICEaEVPcoeff, &
      !$ACC           SEAICEaEVPcStar,SEAICE_cf,SEAICE_cpAir, ICE2WATR,             &
      !$ACC           SEAICE_lhEvap,SEAICE_dalton,SEAICE_iceConduct,                &
      !$ACC           SEAICE_snowConduct,SEAICE_emissivity,SEAICE_snowThick,        &
      !$ACC           SEAICE_shortwave,OCEAN_drag,SEAICE_mcPheeTaper,               &
      !$ACC           SEAICE_availHeatFrac,SEAICE_tauAreaObsRelax,                  &
      !$ACC           SEAICE_airTurnAngle,SEAICE_waterTurnAngle,                    &
      !$ACC           SEAICE_monFreq,SEAICE_dumpFreq,SEAICE_taveFreq,               &
      !$ACC           SEAICE_tensilDepth,SEAICE_eccen)

!!Y: mcpheepiston need?,  SEAICE_evpDampC,elasticParm


!--   Constants used by sea-ice model
      real(wp), parameter  ::  ZERO = 0.0_wp, ONE = 1.0_wp, TWO = 2.0_wp
      real(wp), parameter  ::  QUART = 0.25_wp, HALF = 0.5_wp
      real(wp), parameter  ::  third =0.333333333333333333333333333_wp 
      real(wp), parameter  ::  siEps = 1.0d-5

!--   Constants needed by McPhee formulas for turbulent ocean fluxes :
!        Stanton number (dimensionless), typical friction velocity
!        beneath sea ice (m/s), and tapering factor (dimensionless)
      real(wp), parameter  :: STANTON_NUMBER = 0.0056_wp, USTAR_BASE=0.0125_wp, MCPHEE_TAPER_FAC = 12.5_wp


!--   Used to indicate variables that have not been given a value
      real(dp), parameter, public  ::  UNSET_FLOAT8 = 1.234567D5 
      real(sp), parameter, public  ::  UNSET_FLOAT4 = 1.234567E5
      real(dp), parameter, public  ::  UNSET_RL     = 1.234567D5
      real(sp), parameter, public  ::  UNSET_RS     = 1.234567D5
      integer, parameter           ::  UNSET_I      = 123456789

!--  precipitation freshwater density
      real(wp), parameter :: rhoConstFresh       = 999.8_wp
!--convert centigrade (Celsius) degree to Kelvin
      real(wp), parameter :: celsius2K = 273.15_wp

!--SEAICE_boltzmann   :: Stefan-Boltzman constant (not a run time parameter)
!-- Step is from mod_misc_basic
      real(wp), parameter ::  SEAICE_boltzmann = Stef 

end module mitice_parameters
