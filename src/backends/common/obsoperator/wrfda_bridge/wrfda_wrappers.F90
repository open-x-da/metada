! Wrapper module that brings WRFDA transform routines into a module we build
! and link against. It includes the WRFDA-provided .inc implementations for
! the specified observation families.

module metada_wrfda_wrappers
  use da_define_structures
  use da_control
  use da_reporting
  implicit none
contains

  ! Sound (sonde upper-air) and sonde surface
  #include "da_sound/da_transform_xtoy_sound.inc"
  #include "da_sound/da_transform_xtoy_sound_adj.inc"
  #include "da_sound/da_transform_xtoy_sonde_sfc.inc"
  #include "da_sound/da_transform_xtoy_sonde_sfc_adj.inc"

  ! SYNOP
  #include "da_synop/da_transform_xtoy_synop.inc"
  #include "da_synop/da_transform_xtoy_synop_adj.inc"

  ! METAR
  #include "da_metar/da_transform_xtoy_metar.inc"
  #include "da_metar/da_transform_xtoy_metar_adj.inc"

  ! BUOY
  #include "da_buoy/da_transform_xtoy_buoy.inc"
  #include "da_buoy/da_transform_xtoy_buoy_adj.inc"

  ! SHIPS
  #include "da_ships/da_transform_xtoy_ships.inc"
  #include "da_ships/da_transform_xtoy_ships_adj.inc"

  ! AIREP
  #include "da_airep/da_transform_xtoy_airep.inc"
  #include "da_airep/da_transform_xtoy_airep_adj.inc"

  ! PILOT
  #include "da_pilot/da_transform_xtoy_pilot.inc"
  #include "da_pilot/da_transform_xtoy_pilot_adj.inc"

end module metada_wrfda_wrappers


