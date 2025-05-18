module mitice_utility
use mod_misc_basic
use mod_mpi_interfaces
use mod_mpi_variables
use mitice_vars, only: mitice_runmsg_unit
use mitice_parameters

   implicit none

contains
!==============================================================================
   subroutine mitice_run_info_open
!==============================================================================
      implicit none

      integer :: istat
      character(lc) :: msg, filename, domid_str

      mitice_runmsg_unit = idom*10+2

      write (domid_str, "(i2.2)") idom

      IF (mpi_rank == 0) THEN
         filename = "mitice_run_information_d"//trim(domid_str)//".output"
         open (file=trim(filename), unit=mitice_runmsg_unit, action='write', iostat=istat, iomsg=msg)
         if (istat /= 0) then
            write (*, *) 'mitice_run_information.output open failed, error message = ', msg
            stop
         end if
         write (unit=mitice_runmsg_unit, fmt=*) '         NMEFC MaCOM-Mitice coupled run information'
      END IF

   end subroutine mitice_run_info_open

!==============================================================================
   subroutine mitice_run_info_close
!==============================================================================
      implicit none
      integer :: istat
      character(lc) :: msg, time, date, zone, timestamp

      IF (mpi_rank == 0) THEN
         call date_and_time(date=date, time=time, zone=zone)
         timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                           time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//zone
         write (mitice_runmsg_unit,"(a)") ' Seaice Integration Successfully Finished !'
         write (mitice_runmsg_unit,"(a)") trim(timestamp)
         close (unit=mitice_runmsg_unit, status='keep', iostat=istat, iomsg=msg)
         if (istat /= 0) then
            write (*, *) 'mitice_run_information.output close failed, error message = ', msg
         end if
      END IF

   end subroutine mitice_run_info_close

!==============================================================================
   subroutine swfrac(imax,fact,swdk)
!YY: here -200 is depth (m), rF is Pa in MaCOM
!==============================================================================
!     | o Compute solar short-wave flux penetration.
!     | Compute fraction of solar short-wave flux penetrating to
!     | specified depth, swdk, due to exponential decay in
!     | Jerlov water type jwtype.
!     | Reference : Two band solar absorption model of Paulson
!     |             and Simpson (1977, JPO, 7, 952-956)
!     | Notes
!     | =====
!     | Parameter jwtype is hardcoded to 2 for time being.
!     | Below 200m the solar penetration gets set to zero,
!     | otherwise the limit for the exponent (+/- 5678) needs to
!     | be taken care of.
!     | Written by   : Jan Morzel
!     | Date         : July 12, 1995
!     *==========================================================*
      implicit none

!     !INPUT/OUTPUT PARAMETERS:
!     input arguments
!     imax    :: number of vertical grid points
!     fact    :: scale  factor to apply to depth array
      INTEGER  imax
      real(wp) fact
!     input/output arguments
!     swdk    :: on input: vertical depth for desired sw fraction
!               (fact*swdk) is negative distance (m) from surface
!     swdk    :: on output: short wave (radiation) fractional decay
      real(wp) ::    swdk(imax)

!     LOCAL VARIABLES:
!     max number of different water types
      INTEGER   jwtype
      INTEGER, parameter ::   nwtype=5
      real(wp) :: facz
      real(wp) :: rfac(nwtype),a1(nwtype),a2(nwtype)
      INTEGER i

!     Jerlov water type :
!                    I       IA        IB       II      III
!     jwtype :       1       2         3        4       5
      rfac  =   (/ 0.58_wp, 0.62_wp, 0.67_wp, 0.77_wp, 0.78_wp /)
        a1  =   (/ 0.35_wp,  0.6_wp,  1.0_wp,  1.5_wp,  1.4_wp /)
        a2  =   (/ 23.0_wp, 20.0_wp, 17.0_wp, 14.0_wp,  7.9_wp /)

      jwtype=2

      DO i = 1,imax
        facz = fact*swdk(i)
        IF ( facz .LT. -200.0_wp ) THEN
          swdk(i) = 0.0_wp
        ELSE
          swdk(i) =       rfac(jwtype)  * exp( facz/a1(jwtype) )   &
             + (1.0_wp - rfac(jwtype)) * exp( facz/a2(jwtype) )
        ENDIF
      ENDDO
   end subroutine swfrac

   subroutine freeze_surface
!     *==========================================================*
!     | o Check water temperature and limit range of temperature
!     | appropriately.
!     *==========================================================*

      implicit none

!     Tfreezing :: Freezing threshold temperature.
      INTEGER i,kSrf

      kSrf = nk

!     Check for water that should have frozen
!$acc kernels loop copyin(kSrf)
      do i = 1,nlpb
        IF (tFld(i,kSrf,1) .LT. -1.9_wp) THEN
          tFld(i,kSrf,1) = -1.9_wp
        ENDIF
      enddo
!$acc end kernels

   end subroutine freeze_surface

   subroutine seaice_summary

!  *==========================================================*
!  | o Summarize seaice parameters.
!  *==========================================================*

      implicit none

      if (mpi_rank == 0) then


      WRITE(mitice_runmsg_unit,'(A)')'// ======================================================='
      WRITE(mitice_runmsg_unit,'(A)')'// Seaice configuration (SEAICE_PARM01) >>> START <<<'
      WRITE(mitice_runmsg_unit,'(A)')'// ======================================================='

!--  Time-stepping related param.

      WRITE(mitice_runmsg_unit,'(A)') ' '
      WRITE(mitice_runmsg_unit,'(A)')'   Seaice time stepping configuration   > START <  '
      WRITE(mitice_runmsg_unit,'(A)')'   ----------------------------------------------'

     write(mitice_runmsg_unit,*)'SEAICE_deltaTtherm=',SEAICE_deltaTtherm, ' /thermodynamic timestep/'
     write(mitice_runmsg_unit,*)'SEAICE_deltaTdyn  =',SEAICE_deltaTdyn, ' /dynamic timestep/'
     write(mitice_runmsg_unit,*)'SEAICE_deltaTevp  =',SEAICE_deltaTevp, ' /EVP timestep/'

     WRITE(mitice_runmsg_unit,'(A)')'   Seaice dynamics configuration   > START <  '
     WRITE(mitice_runmsg_unit,'(A)')'   ------------------------------------------'

!--  Seaice-Dynamics parameters
     write(mitice_runmsg_unit,*)'SEAICEuseDYNAMICS =',SEAICEuseDYNAMICS, ' /use dynamics/'
     IF (SEAICEuseDYNAMICS) THEN
       write(mitice_runmsg_unit,*)'model grid type is C-GRID'

       write(mitice_runmsg_unit,*)'SEAICEuseEVP = ',SEAICEuseEVP,' /* use EVP solver */'
       write(mitice_runmsg_unit,*)'SEAICEuseFREEDRIFT = ', SEAICEuseFREEDRIFT,' /* use free drift solution */'
       write(mitice_runmsg_unit,*)'OCEAN_drag  = ', OCEAN_drag, ' /* air-ocean drag coefficient */'
       write(mitice_runmsg_unit,*)'SEAICE_drag  = ', SEAICE_drag, ' /* air-ice drag coefficient */'
       write(mitice_runmsg_unit,*)'SEAICE_drag_south =', SEAICE_drag_south, ' /* Southern Ocean SEAICE_drag */'
       write(mitice_runmsg_unit,*)'SEAICE_waterDrag  =', SEAICE_waterDrag, ' /* water-ice drag (no units) */'
       write(mitice_runmsg_unit,*)'SEAICE_waterDrag_south =',SEAICE_waterDrag_south,' /* Southern Ocean waterDrag */'
       write(mitice_runmsg_unit,*)'SEAICEdWatMin =', SEAICEdWatMin,' /* minimum linear water-ice drag (in m/s) */'
       write(mitice_runmsg_unit,*)'SEAICEuseTilt =', SEAICEuseTilt,' /* include surface tilt in dyna. */'
       write(mitice_runmsg_unit,*)'SEAICEuseTEM      =', SEAICEuseTEM,' /* use truncated ellipse rheology */'
       write(mitice_runmsg_unit,*)'SEAICE_strength   =',SEAICE_strength, ' /* sea-ice strength Pstar */'
       write(mitice_runmsg_unit,*)'SEAICE_cStar      =', SEAICE_cStar,' /* sea-ice strength parameter cStar */'
       write(mitice_runmsg_unit,*)'SEAICEpressReplFac=', SEAICEpressReplFac,' /* press. replacement method factor */'
       write(mitice_runmsg_unit,*)'SEAICE_tensilFac  =',SEAICE_tensilFac, ' /* sea-ice tensile strength factor */'
       write(mitice_runmsg_unit,*)'SEAICE_tensilDepth=',SEAICE_tensilDepth, ' /* crit. depth for tensile strength */'
       write(mitice_runmsg_unit,*)'SEAICEpresH0   =', SEAICEpresH0,' /* sea-ice strength Heff threshold */'
       write(mitice_runmsg_unit,*)'SEAICEpresPow0 =', SEAICEpresPow0,' /* exponent for Heff<SEAICEpresH0 */'
       write(mitice_runmsg_unit,*)'SEAICEpresPow1 =', SEAICEpresPow1,' /* exponent for Heff>SEAICEpresH0 */'     
       write(mitice_runmsg_unit,*)'SEAICEetaZmethod =',SEAICEetaZmethod, ' /* method computing eta at Z-point */'    
       write(mitice_runmsg_unit,*)'SEAICE_eccen    =', SEAICE_eccen, ' /* elliptical yield curve eccent */'   
       write(mitice_runmsg_unit,*)'SEAICEstressFactor    =',SEAICEstressFactor,' /* wind stress scaling factor */'   
       write(mitice_runmsg_unit,*)'SEAICE_no_slip    =',SEAICE_no_slip, ' /* no slip boundary conditions */'   
       write(mitice_runmsg_unit,*)'SEAICE_clipVeloctities =', SEAICE_clipVelocities,' /* impose max. vels. */'     
       write(mitice_runmsg_unit,*)'useHB87stressCoupling  =',useHB87stressCoupling, ' /* altern. ice-ocean stress */'    
       write(mitice_runmsg_unit,*)'SEAICEscaleSurfStress  =',SEAICEscaleSurfStress, &
                                  ' /* scale atm. and ocean-surface stress with AREA */'  
       write(mitice_runmsg_unit,*)'SEAICEaddSnowMass =',SEAICEaddSnowMass, ' /* add snow mass to seaiceMassC/U/V */'   
       write(mitice_runmsg_unit,*)'SEAICE_elasticParm=',SEAICE_elasticParm, ' /* EVP elastic parameter */'     
       write(mitice_runmsg_unit,*)'SEAICE_evpTauRelax=',SEAICE_evpTauRelax, ' /* EVP relaxation timescale */'     
       write(mitice_runmsg_unit,*)'SEAICE_evpDampC   =',SEAICE_evpDampC, ' /* EVP damping parameter */'    
       write(mitice_runmsg_unit,*)'SEAICEuseEVPstar  =', SEAICEuseEVPstar,' /* use EVP* solver */'     
       write(mitice_runmsg_unit,*)'SEAICEuseEVPrev   =', SEAICEuseEVPrev,' /* use "revisited EVP" solver */'   
       write(mitice_runmsg_unit,*)'SEAICE_evpAlpha   =',SEAICE_evpAlpha, ' /* EVP* parameter*/'   
       write(mitice_runmsg_unit,*)'SEAICE_evpBeta    =',SEAICE_evpBeta, ' /* EVP*  parameter */'     
       write(mitice_runmsg_unit,*)'SEAICEaEVPcoeff   =', SEAICEaEVPcoeff ,' /* adaptive EVP parameter*/'    
       write(mitice_runmsg_unit,*)'SEAICEaEVPcStar   =',SEAICEaEVPcStar , ' /* adaptive EVP parameter*/'     
       write(mitice_runmsg_unit,*)'SEAICEnEVPstarSteps =',SEAICEnEVPstarSteps, ' /* num. of EVP* steps */'   
       write(mitice_runmsg_unit,*)'SEAICEuseEVPpickup=', SEAICEuseEVPpickup, ' /* start EVP solver with EVP pickup*/'  
  
      ENDIF !end if SEAICEuseDYNAMICS

       write(mitice_runmsg_unit,*)   ' '  
       write(mitice_runmsg_unit,*)'   Seaice advection diffusion config,   > START <  '     
       write(mitice_runmsg_unit,*)'   -----------------------------------------------'   
       write(mitice_runmsg_unit,*)'SEAICEmomAdvection =', SEAICEmomAdvection, ' /* advect sea ice momentum */'       
       write(mitice_runmsg_unit,*)'SEAICEadvHeff =',SEAICEadvHeff, ' /* advect effective ice thickness */'    
       write(mitice_runmsg_unit,*)'SEAICEadvArea =',SEAICEadvArea, ' /* advect fractional ice area */'     
       write(mitice_runmsg_unit,*)'SEAICEadvSnow =',SEAICEadvSnow, ' /* advect snow layer together with ice */'   
#ifdef SEAICE_VARIABLE_SALINITY
       write(mitice_runmsg_unit,*)'SEAICEadvSalt =', SEAICEadvSalt,' /* advect salinity together with ice */'
#endif     
       write(mitice_runmsg_unit,*)'SEAICEdiffKhArea   =',SEAICEdiffKhArea, ' /* diffusivity (m^2/s) for area */'     
       write(mitice_runmsg_unit,*)'SEAICEdiffKhHeff   =',SEAICEdiffKhHeff, ' /* diffusivity (m^2/s) for heff */'   
       write(mitice_runmsg_unit,*)'SEAICEdiffKhSnow   =',SEAICEdiffKhSnow, ' /* diffusivity (m^2/s) for snow */'



#ifdef SEAICE_ITD
!--   ITD parameters
       write(mitice_runmsg_unit,*)'   Seaice ice thickness distribution configuration   > START <  '     
       write(mitice_runmsg_unit,*)'   -----------------------------------------------------------'   
       write(mitice_runmsg_unit,*)'nITD              =',SEAICE_multDim, ' /* number of ice thickness categories */'
       write(mitice_runmsg_unit,*)'SEAICEuseLinRemapITD  =',SEAICEuseLinRemapITD,' /* select linear remapping scheme for ITD */'     
       write(mitice_runmsg_unit,*)'useHibler79IceStrength  =',useHibler79IceStrength,' /* select ice strength parameter */'   
       write(mitice_runmsg_unit,*)'SEAICEsimpleRidging  =',SEAICEsimpleRidging,' /* select ridging scheme */'     
       write(mitice_runmsg_unit,*)'SEAICEpartFunc   =', SEAICEpartFunc,' /* select ridging participation function */'     
       write(mitice_runmsg_unit,*)'SEAICEredistFunc =',SEAICEredistFunc, ' /* select ridging redistribution function */'   
       write(mitice_runmsg_unit,*)'SEAICE_cf  =', SEAICE_cf  ,' /* ice strength parameter */'         
       write(mitice_runmsg_unit,*)'SEAICEshearParm  =',SEAICEshearParm, ' /* amount of energy lost to shear */'
#endif

!--   Thermodynamics parameters       
       write(mitice_runmsg_unit,*)' '
       write(mitice_runmsg_unit,*)'   Seaice thermodynamics configuration   > START <  '
       write(mitice_runmsg_unit,*)'   -----------------------------------------------'     
       write(mitice_runmsg_unit,*)'SEAICE_rhoIce     =',SEAICE_rhoIce , ' /* density of sea ice (kg/m3) */'
       write(mitice_runmsg_unit,*)'SEAICE_rhoSnow    =',SEAICE_rhoSnow ,' /* density of snow (kg/m3) */'     
       write(mitice_runmsg_unit,*)'usePW79thermodynamics  =',usePW79thermodynamics, ' /* default 0-layer TD */'

      IF (.NOT.usePW79thermodynamics) THEN
        write(mitice_runmsg_unit,*)  '   seaice thermodynamics is OFF  '
      ELSE
       write(mitice_runmsg_unit,*)'SEAICE_lhEvap     =', SEAICE_lhEvap, ' /* latent heat of evaporation */' 
       write(mitice_runmsg_unit,*)'SEAICE_lhFusion   =', SEAICE_lhFusion,' /* latent heat of fusion */'
       write(mitice_runmsg_unit,*)'SEAICE_mcPheePiston =',SEAICE_mcPheePiston, &
                                 ' /* turbulent flux "piston velocity" (m/s) */'
       write(mitice_runmsg_unit,*)'SEAICE_mcPheeTaper =',SEAICE_mcPheeTaper, &
                                 ' /* tapering of turbulent flux (0.< <1.) for AREA=1. */'
       write(mitice_runmsg_unit,*)'SEAICE_frazilFrac =',SEAICE_frazilFrac,'  &
                                  /* frazil (T<tempFrz) to seaice conversion rate (0.< <1.) */'
       write(mitice_runmsg_unit,*)'SEAICE_tempFrz0   =',SEAICE_tempFrz0,' /* freezing temp. of sea water (intercept) */'
       write(mitice_runmsg_unit,*)'SEAICE_growMeltByConv  =', SEAICE_growMeltByConv ,' /* grow,melt by vert. conv. */'  
       write(mitice_runmsg_unit,*)'SEAICE_doOpenWaterGrowth =',SEAICE_doOpenWaterGrowth, ' /* grow by open water */'
       write(mitice_runmsg_unit,*)'SEAICE_doOpenWaterMelt =',SEAICE_doOpenWaterMelt, ' /* melt by open water */'
       write(mitice_runmsg_unit,*)'HO                =', HO,' /* nominal thickness of new ice */'
       write(mitice_runmsg_unit,*)'SEAICE_area_max        =',SEAICE_area_max,' /* set to les than 1. to mimic open leads */'
#ifdef SEAICE_VARIABLE_SALINITY
       write(mitice_runmsg_unit,*)'SEAICE_saltFrac =',SEAICE_saltFrac,' /* fraction of ocn salinity in new ice */'
#else
       write(mitice_runmsg_unit,*)'SEAICE_salt0   =', SEAICE_salt0,' /* constant sea ice salinity */'
#endif       
       
       write(mitice_runmsg_unit,*)'SEAICEuseFlooding =',SEAICEuseFlooding,  ' /* turn submerged snow into ice */'
       write(mitice_runmsg_unit,*)' '
       write(mitice_runmsg_unit,*)'   Seaice air-sea fluxes configuration,   > START <  '
       write(mitice_runmsg_unit,*)'   -----------------------------------------------'
       write(mitice_runmsg_unit,*)'SEAICEheatConsFix  =', SEAICEheatConsFix,' /* accound for ocn<->seaice advect. heat flux */'
       write(mitice_runmsg_unit,*)'IMAX_TICE         =', IMAX_TICE,' /* iterations for ice surface temp */'
       write(mitice_runmsg_unit,*)'postSolvTempIter=',postSolvTempIter,' /* flux calculation after surf. temp iter */'
       write(mitice_runmsg_unit,*)'SEAICE_dryIceAlb  =', SEAICE_dryIceAlb  ,' /* winter albedo */'
       write(mitice_runmsg_unit,*)'SEAICE_wetIceAlb  =', SEAICE_wetIceAlb  ,' /* summer albedo */'
       write(mitice_runmsg_unit,*)'SEAICE_drySnowAlb =',SEAICE_drySnowAlb , ' /* dry snow albedo */'
       write(mitice_runmsg_unit,*)'SEAICE_wetSnowAlb =',SEAICE_wetSnowAlb , ' /* wet snow albedo */'
       write(mitice_runmsg_unit,*)'SEAICE_wetAlbTemp=',SEAICE_wetAlbTemp , ' /* Temp (o.C) threshold for wet-albedo */'
       write(mitice_runmsg_unit,*)'SEAICE_snow_emiss =',SEAICE_snow_emiss , ' /* snow emissivity */'
       write(mitice_runmsg_unit,*)'SEAICE_ice_emiss =',SEAICE_ice_emiss , ' /* seaice emissivity */'
       write(mitice_runmsg_unit,*)'SEAICE_iceConduct =', SEAICE_iceConduct ,' /* sea-ice conductivity */'
       write(mitice_runmsg_unit,*)'SEAICE_snowConduct=', SEAICE_snowConduct,' /* snow conductivity */'
       write(mitice_runmsg_unit,*)'SEAICE_snowThick  =',SEAICE_snowThick ,' /* cutoff snow thickness (for albedo) */'
       write(mitice_runmsg_unit,*)'SEAICE_shortwave  =', SEAICE_shortwave  ,' /* penetration shortwave radiation */'
       !write(mitice_runmsg_unit,*)'useMaykutSatVapPoly =',useMaykutSatVapPoly,' /* use Maykut Polynomial for Sat.Vap.Pr */'
       write(mitice_runmsg_unit,*)'MIN_ATEMP         =', MIN_ATEMP,' /* minimum air temperature */'
       write(mitice_runmsg_unit,*)'MIN_LWDOWN        =',MIN_LWDOWN, ' /* minimum downward longwave */'
       write(mitice_runmsg_unit,*)'MIN_TICE          =',MIN_TICE, ' /* minimum ice temperature */'

      ENDIF   !end if usePW79thermodynamics

       write(mitice_runmsg_unit,*) ' '
       write(mitice_runmsg_unit,*)'   Seaice initialization and IO config.,   > START <  '
       write(mitice_runmsg_unit,*)'   -------------------------------------------------'
!--  Initial Condition/Input related param.
       write(mitice_runmsg_unit,*)'SEAICE_initialHEFF=',SEAICE_initialHEFF, ' /* initial sea-ice thickness */'
       IF ( AreaFile .NE. ' ' ) THEN
         write(mitice_runmsg_unit,*)'AreaFile =',AreaFile, ' /* Initial ice concentration File */'
       ELSE
         write(mitice_runmsg_unit,*)'AreaFile is not set. /* Initial ice concentration File */'
       ENDIF
    
       IF ( HeffFile .NE. ' ' ) THEN
         write(mitice_runmsg_unit,*)'HeffFile =',HeffFile, ' /* Initial effective ice thickness File */'
       ELSE
         write(mitice_runmsg_unit,*)'HeffFile is not set. /* Initial effective ice thickness File */'
       ENDIF
       IF ( HsnowFile .NE. ' ' ) THEN
         write(mitice_runmsg_unit,*)'HsnowFile =', HsnowFile,' /* Initial snow thickness File */'
       ELSE
         write(mitice_runmsg_unit,*)'HsnowFile is not set. /* Initial snow thickness File */'
       ENDIF
#ifdef SEAICE_VARIABLE_SALINITY
       IF ( HsaltFile .NE. ' ' ) THEN
         write(mitice_runmsg_unit,*)'HsaltFile =', HsaltFile,' /* Initial HSALT File */'
       ELSE
         write(mitice_runmsg_unit,*)'HsaltFile is not set. /* Initial HSALT File */'
       ENDIF
#endif
       IF ( uIceFile .NE. ' ' ) THEN
         write(mitice_runmsg_unit,*)'uIceFile =',uIceFile, ' /* Initial U-ice velocity File */'
       ELSE
         write(mitice_runmsg_unit,*)'uIceFile is not set. /* Initial U-ice velocity File */'
       ENDIF
       IF ( vIceFile .NE. ' ' ) THEN
         write(mitice_runmsg_unit,*)'vIceFile =',vIceFile, ' /* Initial V-ice velocity File */'
       ELSE
         write(mitice_runmsg_unit,*)'vIceFile is not set. /* Initial V-ice velocity File */'
       ENDIF
    
       write(mitice_runmsg_unit,*) ' '
       write(mitice_runmsg_unit,*)'// ======================================================='
       write(mitice_runmsg_unit,*)'// Seaice configuration (SEAICE_PARM01) >>> END <<<'
               
      endif  !endif for mpi_rank==0

   end subroutine seaice_summary

end module mitice_utility
