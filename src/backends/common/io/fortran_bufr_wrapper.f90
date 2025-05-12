! fortran_bufr_wrapper.f90
! Fortran wrapper for BUFR file operations with C/C++ interoperability

module fortran_bufr_wrapper
  use, intrinsic :: iso_c_binding
  implicit none

contains
  ! Open a BUFR file with a specific unit number
  ! This wraps the Fortran OPEN, OPENBF, and DATELEN calls
  subroutine open_bufr_file(filename, unit_num, status) bind(C, name="open_bufr_file_")
    character(kind=c_char), intent(in) :: filename(*)
    integer(c_int), intent(in), value :: unit_num
    integer(c_int), intent(out) :: status
    
    character(len=1024) :: f_filename
    integer :: i
    character(len=2) :: io_method

    ! Convert C string to Fortran string
    i = 1
    do
      if (filename(i) == c_null_char) exit
      f_filename(i:i) = filename(i)
      i = i + 1
    end do
    f_filename(i:) = ' '
    
    ! Set the I/O method for OPENBF
    io_method = 'IN'
    
    ! Initialize status to success
    status = 0
    
    ! Attempt to open the file
    open(unit=unit_num, file=f_filename(1:i-1), form='UNFORMATTED', status='OLD', err=100)
    
    ! Call BUFR opening routines
    call openbf(unit_num, io_method, unit_num)
    call datelen(10)
    return
    
100 continue
    ! Error occurred during OPEN
    status = -1
    return
  end subroutine open_bufr_file
  
  ! Close a BUFR file for a specific unit number
  subroutine close_bufr_file(unit_num) bind(C, name="close_bufr_file_")
    integer(c_int), intent(in), value :: unit_num
    
    ! Close the BUFR file
    call closbf(unit_num)
    
    ! Close the Fortran file unit
    close(unit_num)
  end subroutine close_bufr_file
  
end module fortran_bufr_wrapper 