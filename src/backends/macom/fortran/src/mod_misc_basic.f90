module mod_misc_basic
   implicit none

   integer, public, parameter  :: hp = SELECTED_REAL_KIND(3)    ! half precision (real 2)
   integer, public, parameter  :: sp = kind(1.0)                ! single precision (real 4)
   integer, public, parameter  :: dp = kind(1.0d0)              ! double precision (real 8)
   integer, public, parameter  :: wp = dp
   integer, public, parameter  :: i1 = selected_int_kind(2)    ! integer byte
   integer, public, parameter  :: i2 = selected_int_kind(4)    ! integer short
   integer, public, parameter  :: i4 = selected_int_kind(9)    ! integer long
   integer, public, parameter  :: i8 = selected_int_kind(14)    ! integer larger
   integer, public, parameter  :: lc = 256                      ! lenght of character strings
   integer, public, parameter  :: max_domain = 20
   integer, public, parameter  :: iu = 1
   integer, public, parameter  :: iv = 2
   integer, public, parameter  :: tide_total = 147
   complex(wp), public, parameter :: unitC = (0,1)
   complex(wp), public, parameter :: unitConj = (0,-1)
   real(wp), public, parameter :: rau0 = 1026.0_wp                ! volumic mass of reference     [kg/m3]
   real(wp), public, parameter :: gravity = 9.8_wp
   real(wp), public, parameter :: recip_gravity = 1.0_wp/gravity
   real(wp), public, parameter :: epsmin = 1.0e-12_wp
   real(wp), public, parameter :: thetaMax = 1.0e+20_wp
   real(wp), public, parameter :: SItoBar = 1.0e-05_wp               ! Pascal to bar
   real(wp), public, parameter :: SItodBar = 1.0e-04_wp               ! Pascal to decibar
   real(wp), public, parameter :: vkarmn = 0.4_wp                   ! von Karman constant
   real(wp), public, parameter :: rpi = 3.141592653589793_wp     ! pi
   real(wp), public, parameter :: radPi = rpi/180.0_wp
   real(wp), public, parameter :: rhoa = 1.22_wp                  ! air density                   [kg/m3]
   real(wp), public, parameter :: rt0 = 273.15_wp                ! freezing point of fresh water [Kelvin]
   real(wp), public, parameter :: rt0_snow = 273.15_wp                ! melting point of snow         [Kelvin]
   real(wp), public, parameter :: albo = 0.066_wp                 ! ocean albedo assumed to be constant
   real(wp), public, parameter :: Stef = 5.67e-08_wp              ! Stefan Boltzmann constant
   real(wp), public, parameter :: Lv = 2.5e+06_wp               ! latent heat of vaporization
   real(wp), public, parameter :: rhosn = 330.0_wp                 ! volumic mass of snow          [kg/m3]
   real(wp), public, parameter :: xlsn = 110.121e+06_wp           ! volumetric latent heat fusion of snow  [J/m3]
   real(wp), public, parameter :: lfus = xlsn/rhosn               ! latent heat of fusion of fresh ice [J/Kg]
   real(wp), public, parameter :: rcp = 3991.86795711963_wp      ! ocean specific heat           [J/Kg/K]
   real(wp), public, parameter :: cpa = 1000.5_wp                ! specific heat of air
   real(wp), public, parameter :: rhoic = 900.0_wp                 ! volumic mass of sea ice       [kg/m3]
   real(wp), public, parameter :: rcpic = 1.8837e+06_wp            ! volumetric specific heat for ice   [J/m3/K]
   real(wp), public, parameter :: cpic = rcpic/rhoic              ! specific heat for ice   [J/Kg/K]
   real(wp), public, parameter :: r1_rau0_rcp = 1.0_wp/(rau0*rcp)
   real(wp), public, parameter :: r1_rau0 = 1.0_wp/rau0
   real(wp), public, parameter :: rsmall = 0.5_wp*epsilon(1.0_wp)   ! smallest real computer value
   real(wp), public, parameter :: rho_c = 0.01_wp                  ! density criterion for mixed layer depth
   real(wp), public, parameter :: avt_c = 5.0e-4_wp                ! Kz criterion for the turbocline depth
   real(dp), public, parameter :: dround = 1000000.0d0
   !YY add sphere radius
   real(dp), public, parameter  :: rSphere = 6370.0d03
   real(dp), public, parameter  :: recip_rSphere = 1.0_wp/rSphere
   !YY add debug level
   integer, public, parameter  :: debuglevel=0    !define 2 levels, 1: basic; 2: more debug, add to namelist later

end module mod_misc_basic
