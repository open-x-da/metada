module prepbufr_data
  use prepbufr_constants
  implicit none
  
  ! Main data arrays from former COMMON block
  real(kind=8) :: hdr(NHR8PM)
  real(kind=8) :: evns(MXR8PM, MXR8LV, MXR8VN, MXR8VT)
  integer :: nlev
  
  ! Station data
  character(len=150) :: outstg 
  character(len=8) :: subset
  character(len=300) :: inf, outf, config
  character(len=300) :: argv
  character(len=101) :: crec
  character(len=6), dimension(MAXTYPE) :: type
  character(len=1), dimension(MAXPARM) :: parm
  character(len=4) :: id
  character(len=10) :: idatec
  character(len=5), dimension(MAXSAID) :: said
  character(len=5) :: sid
  character(len=5) :: vtmp
  
  ! File handling
  integer :: iunso(NFILO) = (/51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65/)
  character(len=6), dimension(NFILO) :: filo = (/'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR', &
                                               'VADWND', 'SATBOG', 'SATEMP', 'ADPSFC', 'SFCSHP', &
                                               'SFCBOG', 'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND'/)
  
  ! Variable types
  character(len=1), dimension(MXR8VT) :: var = (/'P','Q','T','Z','U','V'/)
  
  ! Control variables
  logical :: skip_vtmp = .false.
  logical :: found
  integer :: plat(MAXPLAT)
  
  ! Other variables
  integer :: io, stat, n, inlength, np, nt, nplat
  integer :: count, k, flag, pflag, p, ns, s, platflag
  real(kind=8) :: lat1, lat2, lon1, lon2, platform
  integer :: tv_ev_idx, tvflag
  integer :: iuno
  integer :: idate, ierrpb
  
end module prepbufr_data 