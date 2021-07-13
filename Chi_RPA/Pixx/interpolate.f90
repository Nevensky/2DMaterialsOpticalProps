program Pi_pol
  USE ISO_Fortran_env
  implicit none
  

  INTEGER, PARAMETER :: sp = real32
  INTEGER, PARAMETER :: dp = real64

  character (len=100) :: pol,dato1,dato2,dato3,dato4,dato5,dato6
  character (len=200) :: config_file
  integer :: No, Nlf
  integer :: No_interp, No_tot, counter
  integer :: io, jo, iG, jG
  integer :: ios_conf
  complex(kind=dp), dimension(:,:),   allocatable  :: Qeff      ! effective charge carriers matrix.
  complex(kind=dp), dimension(:,:,:), allocatable  :: S0         ! korelacijska matrica
  complex(kind=dp), dimension(:,:,:), allocatable  :: Pi_tot, Pi_tot_interp
  complex(kind=dp), dimension(:,:),   allocatable  :: Pi_dia, tmp
  complex(kind=dp) :: Pi_inter, Pi_intra
  real(kind=dp)    :: domega, domega_tot
  real(kind=dp)    :: omin, omax, oi
  real(kind=dp)    :: Gamma_intra, ReChi0,ImChi0
  real(kind=dp)    :: c0


  real(kind=dp), parameter :: pi      = 4.D0*atan(1.D0)
  real(kind=dp), parameter :: Hartree = 2.0D0*13.6056923D0
  

  namelist /config/ pol, Nlf, No_interp, No, omin, omax,Gamma_intra, c0
  call parseCommandLineArgs(config_file) ! config file is the first argument passed to the command
  open(10,file=config_file)
  read(10,nml=config,iostat=ios_conf)
  close(10)

  ! pol = 'xx'
  ! call parseCommandLineArgs(pol)
  ! Nlf = 33  ! 5 Ha gr
  ! Nlf = 65  ! 20 Ha gr
  ! pol = 'zz'
  ! Nlf = 29 ! 5 Ha hbn
  ! c0 = 32.5411 ! a.u.
  ! Gamma_intra = 0.025 ! eV
  ! No = 2001
  ! omin = 1.0d-5
  ! omax = 50.0
  ! No_interp = 10  
  
  No_tot = (No_interp-1)*(No-2)
  
  Gamma_intra = Gamma_intra/Hartree
  omin = omin/Hartree ! iz eV u Hartree
  omax = (omax/Hartree + omin) 
  
  domega = (omax-omin)/(No-1)
  domega_tot = (omax-omin)/No_tot

  dato1 = 'Corrfun_'//adjustl(trim(pol))
  dato2 = 'Qeff_'//adjustl(trim(pol))
  dato6 = 'Pi_RPA_'//adjustl(trim(pol))//'_GiGj_dense'
  dato3 = 'Pi_RPA_'//adjustl(trim(pol))//'_dense' ! G=G'=0
  dato4 = 'Pi_RPA_'//adjustl(trim(pol))//'_dense_inter' ! G=G'=0
  dato5 = 'Pi_RPA_'//adjustl(trim(pol))//'_dense_intra' ! G=G'=0

  print *, 'STARTING PI_pol current-ccurent tensor calc using KK rel.'

  allocate(S0(-No:No,Nlf,Nlf)) 
  allocate(Qeff(Nlf,Nlf)) 

  allocate(Pi_dia(Nlf,Nlf))
  allocate(Pi_tot(No-1,Nlf,Nlf))

  open(74,FILE = dato1)
  omega_loop_D: do io=-No,No
    read(74,*) ! dummy
    ! read(74,'(10F15.10)')((S0(io,iG,jG), jG = 1,Nlf), iG = 1,Nlf)
    do iG=1,Nlf
      do jG=1,Nlf
        read(74,'(2F15.10)') S0(io,iG,jG)
      enddo
    enddo
  end do omega_loop_D
  close(74)

  
  open(75,FILE = dato2)
  do iG=1,Nlf
    do jG=1,Nlf
      read(75,'(2F15.10)') Qeff(iG,jG)
    end do
  end do
  ! read(75,'(10F15.10)')((Qeff(iG,jG), jG = 1,Nlf), iG = 1,Nlf)
  close(75)

  print *, "COMPLETED: Read S0 and Qeff matrices."
  
  
  print *, "STARTED: current-current response calculation."
  ! open(77,FILE = dato3)
  ! open(99)
  ! calculate current-current response Pi(omegaxNlfxNlf) for all local field vectors
  ! separately calculate interband and intraband contributions for (1,1) components
  omega_loop_B: do io = 1,No-1
    ! print *, io
    oi = (io-1)*domega

    iG_loop: do iG = 1,Nlf
      jG_loop: do jG = 1,Nlf

        call genReChi0(io,No,Nlf,iG,jG,oi,domega,S0,Rechi0)
        ! print *,'ReChi0: ',ReChi0
        call genImChi0(io,No,Nlf,iG,jG,oi,domega,S0,ImChi0)
        ! print *,'ImChi0: ',ImChi0
        
        if (io == 1) then ! omega=0.0 Ha
          Pi_dia(iG,jG) = -cmplx(ReChi0,0.0) ! neven debug: diamagnetski doprinos ??
        end if

        Pi_tot(io,iG,jG) = cmplx(ReChi0,ImChi0) 
        Pi_tot(io,iG,jG) = Pi_tot(io,iG,jG) + Pi_dia(iG,jG) ! Pi_RPA = Pi_paramagnetski + Pi_diamagnetski

        Pi_inter = Pi_tot(io,1,1)
        Pi_intra = Qeff(1,1)*oi/(oi + cmplx(0.0,1.0)*Gamma_intra)
        
        ! Pi_tot(io,iG,jG) = Pi_tot(io,iG,jG) + Qeff(iG,jG)*oi/(oi + cmplx(0.0,1.0)*Gamma_intra) ! dodavanje intraband clana
        
      end do jG_loop
    end do iG_loop
    
    ! WRITTING TOTAL RESPONSE FUNCTION Pi for a given polarization 'pol' to file
    ! write(77,*) oi*Hartree,Pi_tot(1,1)
    
    ! write(99,*) oi*Hartree,Pi_tot(4,5)

    ! vodljivost u jedinicama 2*pi*e^2/h   
    ! if(io > 1) then
    !   write(401,*) oi*Hartree, real(-cmplx(0.0,1)*c0*Pi_inter/oi)
    !   write(402,*) oi*Hartree, real(-cmplx(0.0,1)*c0*Pi_intra/oi)
    ! endif


  end do omega_loop_B
  print *, "COMPLETED: current-current response calculation."
  ! close(77)
  ! close(99)
  
  ! No_tot = (No-2)*No_interp - (No_interp-1)
  ! No_tot = (No_interp-1)*(No-2)
  ! domega_tot = (omax-omin)/No_tot

  print *, "STARTED: interpoaltion of frequencies, No_interp:",No_interp

  ! linear interpolation of current-current response Pi(omegaxNlfxNlf) on a dense frequency grid
  allocate(tmp(Nlf,Nlf))
  allocate(Pi_tot_interp(No_tot,Nlf,Nlf))
  open(77,file = dato3)
  open(78,file = dato4)
  open(79,file = dato5)
  counter = 0
  omega_loop_E: do io=1,No-2
    ! print *,'debug Pitot: ',Pi_tot(io,iG,jG)
    tmp(:,:) =  ( Pi_tot(io+1,:,:)-Pi_tot(io,:,:) ) / No_interp
    omega_interp_loop: do jo=1,No_interp-1
      counter = counter + 1
      oi = counter*domega_tot
        iG_loop2: do iG = 1,Nlf
          jG_loop2: do jG = 1,Nlf 
            ! print *,counter,iG,jG
            Pi_tot_interp(counter,iG,jG) = Pi_tot(io,iG,jG) + tmp(iG,jG)*jo
            
            ! Pi_inter = Pi_tot(1,1)
            ! Pi_intra = Qeff(1,1)*oi/(oi + cmplx(0.0,1.0)*Gamma_intra)

            Pi_tot_interp(counter,iG,jG) = Pi_tot_interp(counter,iG,jG) + Qeff(iG,jG)*oi/(oi + cmplx(0.0,1.0)*Gamma_intra)
          end do jG_loop2
        end do iG_loop2
        Pi_inter = Pi_tot_interp(counter,1,1)
        Pi_intra = Qeff(1,1)*oi/(oi + cmplx(0.0,1.0)*Gamma_intra)
        ! WRITTING INTERPOLATED TOTAL RESPONSE FUNCTION Pi for a given polarization 'pol' to file for G=G'=0
        write(77,*) oi*Hartree, real( Pi_tot_interp(counter,1,1) ), aimag( Pi_tot_interp(counter,1,1) )
        write(78,*) oi*Hartree, real( Pi_inter ), aimag( Pi_inter )
        write(79,*) oi*Hartree, real( Pi_intra ), aimag( Pi_intra )

    end do omega_interp_loop
  end do omega_loop_E
  print *,'interpolation counter: ',counter
  print *, 'No_tot: ',No_tot

  close(77)
  close(78)
  close(79)
  print *, "COMPLETED: Interpoaltion of frequencies."


  print *, "STARTED: Writting interpoalted current-current response to file. Polarization: ", pol
  ! WRITTING INTERPOLATED TOTAL RESPONSE FUNCTION Pi for a given polarization 'pol' to file for all G,G'
  open(80,file = dato6)
  do io = 1, No_tot
    do iG = 1, Nlf
      do jG = 1, Nlf
        oi = io*domega_tot
        write(80,*) oi*Hartree, Pi_tot_interp(io,iG,jG)
      enddo
    enddo
  enddo
  close(80)

  deallocate(tmp)
  
  deallocate(S0)
  deallocate(Qeff)
  deallocate(Pi_dia)
  deallocate(Pi_tot)

  print *, "COMPLETED: Writting interpoalted current-current response to file. Polarization: ", pol

  print *,'PROGRAM EXECUTION ENDED FOR CALC = 2'

contains

  subroutine parseCommandLineArgs(config_file)
    implicit none
    character(len=200), intent(out) :: config_file
  
    if(command_argument_count() > 1) then
    print *, 'ERROR: Please provide a single argument corresponding to the config_file path.'
      stop
    else if (command_argument_count() ==0) then
      config_file ='config.interp.in'
    else
      call get_command_argument(1,config_file) 
    endif
    config_file = trim(config_file)
    print *, 'Config file: ',config_file
  end subroutine parseCommandLineArgs

  ! subroutine parseCommandLineArgs(pol)
  !   implicit none
  !   character(len=100), intent(out) :: pol

  !   if(command_argument_count() > 1 .or. command_argument_count()==0 ) then
  !     print *, 'ERROR: Please provide a single argument corresponding to polarization component.'
  !     stop
  !   else
  !     call get_command_argument(1,pol) 
  !   endif
  !   pol = trim(pol)
  !   print *, 'Polarization: ',pol
  ! end subroutine parseCommandLineArgs

  subroutine genReChi0(io,No,Nlf,iG,jG,oi,domega,S0,ReChi0)
    implicit none
    integer,          intent(in)  :: io, No, Nlf
    integer,          intent(in)  :: iG, jG
    real(kind=dp),    intent(in)  :: oi, domega
    complex(kind=dp), intent(in)  :: S0(-No:No,Nlf,Nlf)
    real(kind=dp),    intent(out) :: ReChi0

    integer :: jo
    real(kind=dp) :: oj, fact

    ! neven debug
    ! print *,'S0(-5,1,1):', S0(-5,1,1)

    ReChi0 = 0.0 ! real part of the response function
    ! static limit
    if (io == 1) then
      do  jo = 2,No
        oj = (jo-1)*domega
        fact = domega/oj
        if (jo == 2) then
          fact = 3.0/2.0
        else if (jo == No) then
          fact = 0.5*domega/oj
        end if
        ReChi0 = ReChi0 + fact*( real(S0(-jo+1,iG,jG)) - real(S0(jo-1,iG,jG)) )
      end do
    else if (io == 2) then
      do  jo = 1,No
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj)
        endif
        if (jo == 1) then
          fact = 1.0
        else if (jo == 2) then
          fact = 0.0
        else if (jo == 3) then
          fact=-3.0/2.0
        else if (jo == No) then
          fact = 0.5*domega/(oi-oj)
        end if
        ReChi0 = ReChi0 + fact*real(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)
        if (jo == 1 .or. jo == No) then
          fact = 0.5*domega/(oi + oj)
        end if
        ReChi0 = ReChi0 + fact*real(S0(-jo + 1,iG,jG))
      end do
    else if (io == (No-1)) then
      do  jo = 1,No
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj)
        end if
        if (jo == 1) then
          fact = 0.5*domega/(oi-oj)
        else if (jo == (No-2)) then
          fact = 3.0/2.0
        else if (jo == (No-1)) then
          fact = 0.0
        else if (jo == No) then
          fact=-1.0
        end if
        ReChi0 = ReChi0 + fact*real(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)
        if (jo == 1 .or. jo == No) then
          fact = 0.5*domega/(oi + oj)
        end if
        ReChi0 = ReChi0 + fact*real(S0(-jo + 1,iG,jG))
      end do
    else
      do  jo = 1,No
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj)
        end if
        if (jo == 1) then
          fact = 0.5*domega/(oi-oj)
        else if (jo == (io-1)) then
          fact = 3.0/2.0
        else if (jo == io) then
          fact = 0.0
        else if (jo == (io + 1)) then
          fact=-3.0/2.0
        else if (jo == No) then
          fact = 0.5*domega/(oi-oj)
        end if
        ReChi0 = ReChi0 + fact*real(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)
        if (jo == 1 .or. jo == No) then
          fact = 0.5*domega/(oi + oj)
        end if
        ReChi0 = ReChi0 + fact*real(S0(-jo + 1,iG,jG))
      end do
    end if
    ReChi0 = ReChi0 + pi*aimag(S0(io-1,iG,jG))
    
  end subroutine genReChi0

  subroutine genImChi0(io,No,Nlf,iG,jG,oi,domega,S0,ImChi0)
    implicit none
    integer,          intent(in)  :: io, No, Nlf
    integer,          intent(in)  :: iG, jG
    real(kind=dp),    intent(in)  :: oi, domega
    complex(kind=dp), intent(in)  :: S0(-No:No,Nlf,Nlf)
    real(kind=dp),    intent(out) :: ImChi0

    integer :: jo
    real(kind=dp) :: oj, fact

    ImChi0 = 0.0 ! Imaginary part of the response function Im(Chi)
    ! static limit
    if (io == 1) then
      do  jo = 2,No
        oj = (jo-1)*domega
        fact = domega/oj
        if (jo == 2) then
          fact = 3.0/2.0
        else if (jo == No) then
          fact = 0.5*domega/oj
        end if
        ImChi0 = ImChi0 + fact*(aimag(S0(-jo + 1,iG,jG)) - aimag(S0(jo-1,iG,jG)))
      end do
    else if (io == 2) then
      do  jo = 1,No
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj)
        end if
        if (jo == 1) then 
          fact = 1.0
        else if (jo == 2) then 
          fact = 0.0
        else if (jo == 3) then 
          fact=-3.0/2.0
        else if (jo == No) then
          fact = 0.5*domega/(oi-oj)
        end if

        ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)

        if (jo == 1 .or. jo == No) then
          fact = 0.5*domega/(oi + oj)
        end if
        ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jG))
      end do
    else if (io == (No-1)) then
      do  jo = 1,No
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj) 
        end if
        if (jo == 1) then
          fact = 0.5*domega/(oi-oj)
        else if (jo == (No-2)) then
          fact = 3.0/2.0
        else if (jo == (No-1)) then
          fact = 0.0
        else if (jo == No) then
          fact=-1.0
        end if

        ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)

        if (jo == 1 .or. jo == No) then
          fact = 0.5*domega/(oi + oj)
        end if

        ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jG))
      end do
    else
      do  jo = 1,No
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj)
        end if

        if (jo == 1) then
          fact = 0.5*domega/(oi-oj)
        else if (jo == (io-1)) then
          fact = 3.0/2.0
        else if (jo == io) then
          fact = 0.0
        else if (jo == (io + 1)) then
          fact=-3.0/2.0
        else if (jo == No) then
          fact = 0.5*domega/(oi-oj)
        end if

        ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)

        if (jo == 1 .or. jo == No) then
          fact = 0.5*domega/(oi + oj)
        end if

        ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jG))
      end do
    end if
    
    ImChi0 = ImChi0 - pi*real(S0(io-1,iG,jG)) ! ovaj dio je razlicit od Sloss, S0 je kompleksno polje
    
  end subroutine genImChi0
end program Pi_pol