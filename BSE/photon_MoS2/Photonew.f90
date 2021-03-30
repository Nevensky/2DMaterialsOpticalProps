program photon
 
!  program for the calculation of the current-current 
!  reponse  function in MoS2
use ISO_Fortran_env

implicit none

! single / double precision
integer, parameter :: sp = real32
integer, parameter :: dp = real64
integer, parameter :: qp = real128

! file i/o debug
integer :: ist,ist2,ist4,ist5,ist6,ist7,ist8,ist9,ist10,ist11,ist12
integer :: lno,lno2,lno9,lno10,lno11,lno12

integer :: Nthreads

! Nlfd je minimalno 2 zbog 2X2 blok matrice za p mod

character(len=2) :: pol
integer          :: no, no_ladd, nq
integer          :: qmin, qmax
integer          :: Nl
real(kind=dp)    :: h         ! distance between layer
real(kind=dp)    :: eta       ! analitcko produljenje
real(kind=dp)    :: a0        ! [a.u.]  unit cell parameter in parallel direction 
real(kind=dp)    :: c0        ! [a.u.]  unit cell parameter in perependicular direction 


! constants
real(kind=dp),    parameter :: pi      = 4.D0*atan(1.D0)
real(kind=dp),    parameter :: eV      = 1.602176487D-19
real(kind=dp),    parameter :: Hartree = 2.0D0*13.6056923D0
real(kind=dp),    parameter :: Planck  = 6.626196D-34
real(kind=dp),    parameter :: kB      = 1.3806503D-23
real(kind=dp),    parameter :: aBohr   = 0.5291772d0
real(kind=dp),    parameter :: gamma =  1.0D0/137.0D0 ! fine structure constant e^2/(\hbar c)
! complex(kind=dp), parameter :: rone    = cmplx(1.0,0.0)
! complex(kind=dp), parameter :: czero   = cmplx(0.0,0.0)
! complex(kind=dp), parameter :: ione    = cmplx(0.0,1.0)

! scalars
integer          :: i, j, k, m, iq, io
real(kind=dp)    :: q, o, dq
real(kind=dp)    :: omin, omax
real(kind=dp)    :: domega
complex(kind=dp) :: Pi_ladder_down,Pi_ladder_up,Pi_RPA
complex(kind=dp) :: d0, dxx, dyy, dyz, dzz
complex(kind=dp) :: oi
complex(kind=dp) :: beta
complex(kind=dp) :: Pip
complex(kind=dp) :: ieta

! arrays..
complex(kind=dp), dimension(:),   allocatable :: Pixx,Pizz
complex(kind=dp), dimension(:,:), allocatable :: Pi0 
complex(kind=dp), dimension(:,:), allocatable :: eps, Imat

character (len=100) :: rundir, rpa_xx_file, rpa_zz_file, ladd_down_x_file, ladd_up_x_file, ladd_down_z_file, ladd_up_z_file
namelist /directories/ rundir, rpa_xx_file, rpa_zz_file, ladd_down_x_file, ladd_up_x_file, ladd_down_z_file, ladd_up_z_file
namelist /system/ omin, dq, qmin, qmax
namelist /config/ no, no_ladd
namelist /parameters/ Nl, h, eta, a0, c0
namelist /parallel/ Nthreads

! load namelist
open(10,file='config.in')
read(10,nml=directories,iostat=ist4)
read(10,nml=system,iostat=ist5)
read(10,nml=config,iostat=ist6)
read(10,nml=parameters,iostat=ist7)
read(10,nml=parallel,iostat=ist8)
close(10)

nq   = qmax-qmin
omin = omin/Hartree ! iz eV u Hartree
omax = (omax/Hartree + omin) 

h    = h/aBohr        ! distance between layers, converted from Angstrom to Bohr
eta  = eta/Hartree
ieta = cmplx(0.0,eta)


domega = (omax-omin)/(no-1)

allocate(Pixx(no))
allocate(Pizz(no))
allocate(Pi0(Nl,Nl))
allocate(eps(Nl,Nl))
allocate(Imat(Nl,Nl))

open(31,file = trim(rpa_xx_file) )

open(32,file = trim(ladd_down_x_file) )

open(33,file = trim(ladd_up_x_file) )

do  io =  1,no_ladd
  read(31,*) o, Pi_RPA
  if (io <= no_ladd) then
    read(32,*) o, Pi_ladder_down
    read(33,*) o, Pi_ladder_up
  end if
  Pixx(io) = Pi_RPA + Pi_ladder_down + Pi_ladder_up
end do
close(31)
close(32)
close(33)

open(34,file = trim(rpa_zz_file) )
open(35,file = trim(ladd_down_z_file) )
open(36,file = trim(ladd_up_z_file) )

do  io =  1,no_ladd
  read(31,*) o, Pi_RPA
  if (io <= no_ladd) then
    read(32,*) o, Pi_ladder_down
    read(33,*) o, Pi_ladder_up
  end if
  Pizz(io)=Pi_RPA+Pi_ladder_down+Pi_ladder_up
end do

close(34)
close(35)
close(36)


open(13,file='photon_plot.dat')

q_loop: do iq =  qmin,qmax
  print*, ' iq: ', iq

  q = (iq-1)*dq

  omega_loop: do io = 3,no_ladd
    o =  omin+(io-1)*domega
    oi =  cmplx(o,eta)
    beta =  cmplx(gamma**2  *oi**2 -q**2)
    beta =  sqrt(beta)
    
    dxx =  2.0*pi * cmplx(0.0,1.0) * c0 * gamma**2/beta
    dyy =  2.0*pi * cmplx(0.0,1.0) * c0 * beta/(oi**2)
    dyz = -2.0*pi * cmplx(0.0,1.0) * q * c0/(oi**2)
    dzz =  2.0*pi * cmplx(0.0,1.0) * q**2 * c0/(beta * oi**2)
    
    
    if (pol == 'xx') then 
      d0 = dxx ! S-MOD
    else if (pol == 'zz') then
      d0 = dyy ! P-MOD
    else
      print *, 'FATAL ERROR Specified polarization component not supported.'
      stop
    end if
    
    eps =  cmplx(0.0,0.0)
    Imat(1:Nl,1:Nl) = cmplx(0.0,0.0)
    do  i =  1,Nl
      do  j =  1,Nl
        ! za multilayere nalazi efektivni epsilon = relativna permitivnost
        eps(i,j)= -Pixx(io) * d0 * exp( cmplx(0.0,1.0) * h * beta*abs(real(i)-real(j)) )
        Pi0(i,j)= cmplx(0.0,0.0)
      end do
      eps(i,i) = cmplx(1.0,0.0)+eps(i,i)
      Imat(i,i)= cmplx(1.0,0.0)
      Pi0(i,i) = Pixx(io)
    end do
    
    call gjel(eps,Nl,Nl,Imat,Nl,Nl) ! invertiranje (dio Dysonove jedn.)
    
    Pip =  cmplx(0.0,0.0)
    Pip = sum( eps(1,1:Nl) * Pi0(1:Nl,1) ) ! Dysonova eqn.
    ! do  i =  1,Nl
      ! Pip =  Pip+eps(1,i)*Pi0(i,1)  ! Dysonova eqn.
    ! end do 

    
    
    ! optical conductivity in units  pi*e^2/2h s-mode
    if(iq == 1) write(100,*) o*Hartree, real(-cmplx(0.0,1.0) * 4.0 * c0 * Pip/oi)
    ! if(iq.eq.1)write(200,*)o*Hartree,real(-cmplx(0.0,1.0)*4.0*c0*Pip/oi)
    
    
    
    
    write(3,*) 10.0 * q/aBohr, o*Hartree, real(-cmplx(0.0,1.0) * 4.0 * c0 * Pip/oi)
    
    
    
  end do omega_loop
end do q_loop

close(13)

deallocate(Pixx)
deallocate(Pizz)
deallocate(Pi0)
deallocate(eps)
deallocate(Imat)


end program photon





