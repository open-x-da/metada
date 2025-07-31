module macom_logger
  use iso_c_binding
  implicit none

  private

  ! Log levels
  integer, parameter, public :: LOG_DEBUG = 1
  integer, parameter, public :: LOG_INFO = 2
  integer, parameter, public :: LOG_WARNING = 3
  integer, parameter, public :: LOG_ERROR = 4

  ! Log level names
  character(len=8), parameter :: level_names(4) = [character(len=8) :: &
                                             "DEBUG", "INFO ", "WARN ", "ERROR"]

  ! Default log level
  integer, save :: current_log_level = LOG_INFO

  ! Public interfaces
  public :: log_set_level
  public :: macom_log_debug, macom_log_info, macom_log_warning, macom_log_error

contains

  !-----------------------------------------------------------------------------
  ! Set the minimum log level
  !-----------------------------------------------------------------------------
  subroutine log_set_level(level)
    integer, intent(in) :: level

    if (level >= LOG_DEBUG .and. level <= LOG_ERROR) then
      current_log_level = level
    end if
  end subroutine log_set_level

  !-----------------------------------------------------------------------------
  ! Internal logging function
  !-----------------------------------------------------------------------------
  subroutine log_message(level, component, message)
    integer, intent(in) :: level
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: message
    character(len=len(component)) :: component_clean

    component_clean = trim(adjustl(component))

    if (level >= current_log_level) then
      write (*, '(A,A,A,1X,A,A,A,1X,A)') '[', trim(level_names(level)), ']', &
        '[', trim(component_clean), ']', trim(message)
    end if
  end subroutine log_message

  !-----------------------------------------------------------------------------
  ! Public logging interfaces - Renamed to avoid conflicts
  !-----------------------------------------------------------------------------
  subroutine macom_log_debug(component, message)
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: message

    call log_message(LOG_DEBUG, component, message)
  end subroutine macom_log_debug

  subroutine macom_log_info(component, message)
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: message

    call log_message(LOG_INFO, component, message)
  end subroutine macom_log_info

  subroutine macom_log_warning(component, message)
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: message

    call log_message(LOG_WARNING, component, message)
  end subroutine macom_log_warning

  subroutine macom_log_error(component, message)
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: message

    call log_message(LOG_ERROR, component, message)
  end subroutine macom_log_error

end module macom_logger
