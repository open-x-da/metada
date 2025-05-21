module mitice_vars

!  *==========================================================*
!  | o Basic header for sea ice model.
!  |   Contains most sea ice field declarations.
!  *==========================================================*

! modules from MaCOM
use mod_misc_basic
use mod_csp_basic
! modules from MITICE
use mitice_parameters
implicit none

!Mitice module miscellaneous vars
      integer, public, save  ::  mitice_runmsg_unit
#ifdef SeaiceDebug
      real(wp), public, save, allocatable, dimension(:,:)  :: tmp_outvar  !
#endif

!Grid variables for seaice
      real(wp), public, save, allocatable, dimension(:)  :: HEFFM  !
      !$acc declare create(HEFFM)
      
!     static masks (depend only on geometry)
      real(wp), public, save, allocatable, dimension(:)  :: SIMaskU  !
      real(wp), public, save, allocatable, dimension(:)  :: SIMaskV  !
      !$acc declare create(SIMaskU,SIMaskV)

!     dynamic masks (depend on area)
      real(wp), public, save, allocatable, dimension(:) :: seaiceMaskU  !
      real(wp), public, save, allocatable, dimension(:) :: seaiceMaskV  !
      !$acc declare create(seaiceMaskU,seaiceMaskV)

!Dynamic variables
    real(wp), public, save, allocatable, dimension(:)  :: uice  !
    real(wp), public, save, allocatable, dimension(:)  :: vice  !
    !$acc declare create(uice,vice)

    real(wp), public, save, allocatable, dimension(:)  :: area  !
    !$acc declare create(area)

    real(wp), public, save, allocatable, dimension(:)  :: heff  !
    !$acc declare create(heff)

    real(wp), public, save, allocatable, dimension(:)  :: hsnow !
    !$acc declare create(hsnow)
# ifdef SEAICE_ITD
    real(wp), public, save, allocatable, dimension(:,:)  :: areaitd !
    !$acc declare create(areaitd)

    real(wp), public, save, allocatable, dimension(:,:)  :: heffitd !
    !$acc declare create(heffitd)
    real(wp), public, save, allocatable, dimension(:,:)  :: hsnowitd !

    !$acc declare create(hsnowitd)
!   fraction of open water (= 1-AREA) needed for ridging parameterization
    real(wp), public, save, allocatable, dimension(:)  :: opnWtrFrac  !
    !$acc declare create(opnWtrFrac)

    real(wp), public, save, allocatable, dimension(:)  :: fw2ObyRidge  !
    !$acc declare create(fw2ObyRidge)
# endif
!   bulk and shear viscosity of ice
!    real(wp), public, save, allocatable, dimension(:)  :: ETA     !
    real(wp), public, save, allocatable, dimension(:)  :: etaZ    !
    !$acc declare create(etaZ)

!    real(wp), public, save, allocatable, dimension(:)  :: ZETA    !
    real(wp), public, save, allocatable, dimension(:)  :: zetaZ   !
    !$acc declare create(zetaZ)

!   ice strength/pressure term
    !real(wp), public, save, allocatable, dimension(:)  :: PRESS
!     strain rate tensor
    real(wp), public, save, allocatable, dimension(:)  :: e11  !
    real(wp), public, save, allocatable, dimension(:)  :: e22  !
    real(wp), public, save, allocatable, dimension(:)  :: e12  !
    !$acc declare create(e11,e22,e12)

!   deformation rate tensor invariant, for viscous plastic sea ice =
!   sqrt[(e11**2+e22**2)*(1+1/e**2)+ 4./e**2*e12C**2 + 2*e11*e22*(1-1/e**2))
    real(wp), public, save, allocatable, dimension(:)  :: deltaC    !
    !$acc declare create(deltaC)
!
    real(wp), public, save, allocatable, dimension(:)  :: FORCEX   !
    real(wp), public, save, allocatable, dimension(:)  :: FORCEY   !
    !$acc declare create(FORCEX,FORCEY)

    real(wp), public, save, allocatable, dimension(:)  :: uIceNm1  !
    real(wp), public, save, allocatable, dimension(:)  :: vIceNm1  !
    !$acc declare create(uIceNm1,vIceNm1)
!
    real(wp), public, save, allocatable, dimension(:)  :: seaiceMassC  !
    real(wp), public, save, allocatable, dimension(:)  :: seaiceMassU  !
    real(wp), public, save, allocatable, dimension(:)  :: seaiceMassV  !
    !$acc declare create(seaiceMassC,seaiceMassU,seaiceMassV)

!    
    real(wp), public, save, allocatable, dimension(:)  :: DWATN    !
    !$acc declare create(DWATN)
    real(wp), public, save, allocatable, dimension(:)  :: PRESS0   !
    !$acc declare create(PRESS0)
    real(wp), public, save, allocatable, dimension(:)  :: FORCEX0  !
    real(wp), public, save, allocatable, dimension(:)  :: FORCEY0  !
    !$acc declare create(FORCEX0,FORCEY0)
    real(wp), public, save, allocatable, dimension(:)  :: ZMAX   !
    real(wp), public, save, allocatable, dimension(:)  :: ZMIN   !
    !$acc declare create(ZMAX,ZMIN)

!   factor k to compute the maximal tensile stress from k*PRESS0,
!   in analogy to the maximal compressive stress PRESS0
    real(wp), public, save, allocatable, dimension(:)  :: tensileStrFac  ! 
    !$acc declare create(tensileStrFac)
!
    real(wp), public, save, allocatable, dimension(:)  :: CbotC     !
    !$acc declare create(CbotC)

!   The change of mean ice thickness due to out-of-bounds values following
!   sea ice dynamics and advection
    real(wp), public, save, allocatable, dimension(:)  :: d_HEFFbyNEG  !
    real(wp), public, save, allocatable, dimension(:)  :: d_HSNWbyNEG  !
    !$acc declare create(d_HEFFbyNEG,d_HSNWbyNEG)

!   TICES :: Seaice/snow surface temperature for each category
    real(wp), public, save, allocatable, dimension(:,:)  :: TICES   !
    !$acc declare create(TICES)
 
#ifdef SEAICE_ALLOW_FREEDRIFT
    real(wp), public, save, allocatable, dimension(:) ::  uice_fd,vice_fd 
    !$acc declare create(uice_fd,vice_fd)
#endif
    
!   additional fields needed by the EVP solver
!   seaice_sigma1  - sigma11+sigma22, defined at C-points
!   seaice_sigma2  - sigma11-sigma22, defined at C-points
!   seaice_sigma12 - off-diagonal term, defined at Z-points
    real(wp), public, save, allocatable, dimension(:)  :: seaice_sigma1    !
    real(wp), public, save, allocatable, dimension(:)  :: seaice_sigma2    !
    real(wp), public, save, allocatable, dimension(:)  :: seaice_sigma12   !
    !$acc declare create(seaice_sigma1,seaice_sigma2,seaice_sigma12)

!   stressDivergenceX/Y - divergence of stress tensor
    real(wp), public, save, allocatable, dimension(:)  :: stressDivergenceX  !
    real(wp), public, save, allocatable, dimension(:)  :: stressDivergenceY  !
    !$acc declare create(stressDivergenceX,stressDivergenceY)

!   SWFracB :: fraction of surface Short-Wave radiation reaching
!              the bottom of ocean surface level
    real(wp), public, save  ::   SWFracB
    !$acc declare create(SWFracB)

#ifdef SEAICE_VARIABLE_SALINITY
    real(wp), public, save, allocatable, dimension(:)  ::  HSALT          !
    !$acc declare create(HSALT)
! if HEFFM is zero, or HSALT<0, fix it
    real(wp), public, save, allocatable, dimension(:)  ::  saltFluxAdjust !
    !$acc declare create(saltFluxAdjust)
#endif /* SEAICE_VARIABLE_SALINITY */

!sea-ice loading, expressed in Mass of ice+snow / unitarea (kg/m^2)
!Note: only used with Sea-Ice & RealFreshWater formulation
    real(wp), public, save, allocatable, dimension(:)  :: sIceLoad
    !$acc declare create(sIceLoad)

! store temporary forcing data
    real(wp), public, save, allocatable, dimension(:)  :: qns_ice,qsr_ice
    !$acc declare create(qns_ice,qsr_ice)
    real(wp), public, save, allocatable, dimension(:)  :: utau_ice,vtau_ice
    !$acc declare create(utau_ice,vtau_ice)
    real(wp), public, save, allocatable, dimension(:)  :: EmPmR_ice
    !$acc declare create(EmPmR_ice)
    


end module mitice_vars
