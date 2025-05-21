module mitice_debug
!*===========================================================*
! for debug purpose. ONLY 1d or 2d vars printed to a NC file
! All debug or diagnosis subroutine better be docked here
!   usage: call mpi_debug_output_1d("HEFF", HEFF)
!          call mpi_debug_output_2d("HEFFITD", HEFFITD,nITD)
!*===========================================================*
use mod_mpi_test
use mod_mpi_test_variables
USE mod_mpi_variables
use mod_mpi_interfaces
use mod_mpi_csp_io_misc
use mod_io_netcdf

   implicit none

contains

    subroutine mpi_debug_output_1d(varname,var)
!   *===========================================================*
!   For temporary output of NC file (fname: myIter_varname_debug.nc)
!   note : ONLY dp type
!   usage: call mpi_debug_output_1d("HEFF", HEFF)
!   *===========================================================*
       implicit none
       real(wp), intent(in) :: var(:)
       character(*),intent(in) :: varname
       character(lc) :: varnametmp

       write(varnametmp,'(A1,A,A9)')  '_', adjustl(trim(varname)),'_debug.nc'
       CALL mpi_data_1d_output_comp(varnametmp,var)
    end subroutine mpi_debug_output_1d

    subroutine mpi_debug_output_2d(varname,var,d2)
!   *===========================================================*
!   For temporary output of NC file (fname: myIter_varname_debug.nc)
!   note : ONLY dp type
!   *===========================================================*
       implicit none
       real(wp), intent(in) :: var(:,:)
       character(*),intent(in) :: varname
       integer,intent(in)  :: d2
       character(lc) :: varnametmp

       write(varnametmp,'(A1,A,A9)')  '_', adjustl(trim(varname)),'_debug.nc'
       CALL mpi_data_2d_output_comp(varnametmp, var,d2)
    end subroutine mpi_debug_output_2d

end module mitice_debug
