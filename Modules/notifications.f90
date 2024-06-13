module notifications
  !! Adds colored terminal print formatting for program status reporting
  implicit none

  public :: error, success, warn, info
  private


  ! checks if output is to terminal (my_isatty(1)==1) or is redirected
  interface
    function my_isatty(fd) bind(C, name = 'isatty')
      use, intrinsic :: iso_c_binding, only: c_int
      integer(c_int)        :: my_isatty
      integer(c_int), value :: fd
    end function
  end interface

contains

subroutine error(str)
  implicit none
  character(len=*), intent(in) :: str
  character(len=*), parameter  :: esc = char(27)
  character(len=*), parameter  :: &
    fg_bright_red      = esc//"[91m", &
    color_reset        = esc//"[0m"
  character(len=len(str)+5+4+7) :: str_new

  if (my_isatty(1)==1) then
    str_new = fg_bright_red//'ERROR: '//color_reset//str
  else
    str_new = 'ERROR: '//str
  endif

  print '(A)', str_new
  error stop
end subroutine error

subroutine success(str)
  implicit none
  character(len=*), intent(in) :: str
  character(len=*), parameter  :: esc = char(27)
  character(len=*), parameter  :: &
    fg_bright_green    = esc//"[92m", &
    color_reset        = esc//"[0m"
  character(len=len(str)+5+4+9) :: str_new

  if (my_isatty(1)==1) then
  	str_new = fg_bright_green//'SUCCESS: '//color_reset//str
  else
  	str_new = 'SUCCESS: '//str
  endif
  print '(A)', str_new
end subroutine success

subroutine info(str)
  implicit none
  character(len=*), intent(in) :: str
  character(len=*), parameter  :: esc = char(27)
  character(len=*), parameter  :: &
    fg_bright_blue      = esc//"[94m", &
    color_reset    = esc//"[0m"
  character(len=len(str)+5+4+6) :: str_new

  if (my_isatty(1)==1) then
  	str_new = fg_bright_blue//'INFO: '//color_reset//str
  else
  	str_new = 'INFO: '//str
  end if
  print '(A)', str_new
end subroutine info

subroutine warn(str)
  implicit none
  character(len=*), intent(in)  :: str
  character(len=*), parameter   :: esc = char(27)
  character(len = *), parameter :: &
    fg_bright_yellow      = esc//"[93m", &
    color_reset    = esc//"[0m"
  character(len=len(str)+5+4+9) :: str_new
  
  if (my_isatty(1)==1) then
    str_new = fg_bright_yellow//'WARNING: '//color_reset//str
  else
  	str_new = 'WARNING: '//str
  end if
  print '(A)', str_new
end subroutine warn


end module notifications
