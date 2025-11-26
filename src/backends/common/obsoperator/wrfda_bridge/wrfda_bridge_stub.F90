! Stub source file for wrfda_bridge library when WRFDA support is disabled
! This ensures the library always has at least one source file

! Empty module to ensure the file compiles
! When WRFDA is disabled, this provides the minimal source needed for the library

module wrfda_bridge_stub
  implicit none
contains
  ! Empty module - no functionality when WRFDA is not available
end module wrfda_bridge_stub

