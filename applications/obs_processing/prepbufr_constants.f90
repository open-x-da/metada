module prepbufr_constants
  implicit none

  ! Missing value for BUFR data
  real(kind=8), parameter :: R8BFMS = 10.0E10

  ! Array dimensions and sizes
  integer, parameter :: NHR8PM = 8      ! Actual number of BUFR parameters in header
  integer, parameter :: MXR8PM = 10     ! Maximum number of BUFR parameters
  integer, parameter :: MXR8LV = 400    ! Maximum number of BUFR levels
  integer, parameter :: MXR8VN = 10     ! Maximum number of BUFR event sequences
  integer, parameter :: MXR8VT = 6      ! Maximum number of BUFR variable types
  integer, parameter :: MXSTRL = 80     ! Maximum size of a string
  
  ! Other constants
  integer, parameter :: STRLN = 180
  integer, parameter :: MAXTYPE = 15
  integer, parameter :: MAXPLAT = 15
  integer, parameter :: MAXPARM = 5
  integer, parameter :: MAXSAID = 15
  integer, parameter :: NFILO = 15
  
  ! Virtual temperature constants
  real(kind=8), parameter :: VIRTMP_PROG_CODE = 8.0
  real(kind=8), parameter :: VIRTMP_REASON_CODE = 3.0

end module prepbufr_constants 