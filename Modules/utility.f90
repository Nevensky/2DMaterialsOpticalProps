module utility
  implicit none

  public :: parseCommandLineArgs
  private

contains 
  subroutine parseCommandLineArgs(config_file)
    implicit none
    character(len=200), intent(out) :: config_file

    if(command_argument_count() > 1) then
      print *, 'ERROR: Please provide a single argument corresponding to the config_file path.'
      stop
    else if (command_argument_count() ==0) then
      config_file ='config.in'
    else
      call get_command_argument(1,config_file) 
    endif
    config_file = trim(config_file)
    print *, 'Config file: ',config_file
  end subroutine parseCommandLineArgs

end module utility