program photon
  !  program for the calculation of the current - current corelation function in P
  !  XY - plane  Imat cell parameter  a0=  8.7521  a.u.  
  !  Z - dir Imat cell parameter c0=  32.33282 a.u.  

  use ISO_Fortran_env
  implicit none

  ! single / double precision
  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
  integer, parameter :: qp = real128


  integer i, j, n, m, io, iq,itheta, iG, jG, kG, Nlf, Nlf2
  integer ist1,ist2,ist3,ist4,ist5,ist6

  integer :: NG     ! total number of G vectors 
  integer :: NGd    ! minimum number of coefficients CG over all evc.n files
  integer :: NkI    ! number of wave vectors in IBZ
  integer :: Nk     ! = 48*NkI, number of wave vectors in FBZ with No symmetry 
  integer :: No     ! number of frequencies
  integer :: Nlfd   ! dimenzija polja local field vektora, x2 zbog 2X2 blok matrice za p mod 
  integer :: Ntheta ! ? angle ... dq?
  namelist /config/  NG, NGd, NkI, No, Nlfd, Ntheta


  ! ... constants ...
  real(kind=dp),  parameter :: pi      = 4.d0*atan(1.d0)
  real(kind=dp),  parameter :: eV      = 1.602176487d-19
  real(kind=dp),  parameter :: kB      = 1.3806503d-23
  real(kind=dp),  parameter :: Hartree = 2.0D0*13.6056923d0
  real(kind=dp),  parameter :: Planck  = 6.626196d-34
  real(kind=dp),  parameter :: aBohr   = 0.5291772d0
  real(kind=dp),  parameter :: gamma   = 1.0/137.0

  ! ...scalars...

  real(kind=dp) :: Eref, Gcar
  real(kind=dp) :: q, o, domega
  real(kind=dp) :: theta, dtheta
  complex(kind = dp) :: oi, beta
  ! complex(kind = dp) :: Aii, Ajj, Aij, Bij, diagonal1, diagonal2
   

  ! ...arrays...
  integer, allocatable       :: parG(:)
  integer, allocatable       :: Gi(:)

  real(kind=dp), allocatable :: G(:,:)
  real(kind=dp), allocatable :: Glf(:,:)
  real(kind=dp), allocatable :: KC(:,:)

  complex(kind=dp), dimension(:,:), allocatable :: TS, TP, TScheck, TPcheck
  ! complex(kind=dp), dimension(:,:), allocatable :: Imat, Imat2
  complex(kind=dp), dimension(:,:), allocatable :: Dxx0, Dyy0, Dzz0 ,Dyz0, Dzy0
  complex(kind=dp), dimension(:,:), allocatable :: Pixx0, Piyy0, Piyz0, Pizy0, Pizz0
  complex(kind=dp), dimension(:,:), allocatable :: Pixx, Piyy, Piyz, Pizy, Pizz

  character(len=200) :: path
  character(len=35)  :: tag,buffer

  character(len=200) :: config_file

  character(len=200) :: rundir, savedir, scf_file, rpa_xx_file, rpa_yy_file, rpa_zz_file
  namelist /directories/ rundir, savedir, scf_file, rpa_xx_file, rpa_yy_file, rpa_zz_file

  character(len=3)   :: lf     ! crystal local field effects included in 'xyz' or 'z' direction
  integer            :: qmin, qmax 
  real(kind=dp)      :: dq, omin, omax, eta
  real(kind=dp)      :: Ecut   ! [Hartree] cutoff energy for crystal local field vectors
  namelist /system/ lf, qmin, qmax, dq, omin, omax, eta, Ecut

  real(kind=dp)      :: a0     ! [a.u.]  unit cell parameter in parallel direction 
  real(kind=dp)      :: c0     ! [a.u.]  unit cell parameter in perependicular direction 
  real(kind=dp)      :: Vcell  ! [a.u.^3] unit-cell volume 

  namelist /parameters/ a0, c0, Vcell

  integer            :: Nthreads
  namelist  /parallel/ Nthreads


  ! load namelist
  call parseCommandLineArgs(config_file) ! config file is the first argument passed to ./Photon <arg1>
  open(10,file=config_file)
  read(10,nml=directories,iostat=ist1)
  read(10,nml=system,iostat=ist2)
  read(10,nml=parameters,iostat=ist3)
  read(10,nml=parallel,iostat=ist4)
  close(10)

  ! constants
  Nk   = 48*NkI           ! number of wave vectors in FBZ with No symmetry ops.
  Gcar = 2.0*pi/a0   
  eta  = eta/Hartree
  omin = omin/Hartree ! iz eV u Hartree
  omax = (omax/Hartree + omin) 

  ! allocation of arrays
  allocate( G(3,NG) )
  allocate( parG(NG) )
  allocate( Glf(3,Nlfd) )
  allocate( KC(3,3) )
  allocate( Gi(3) )

  ! allocate( Imat(Nlfd,Nlfd)  )
  ! allocate( Imat2(Nlfd,Nlfd) )

  allocate( Dxx0(Nlfd,Nlfd) )
  allocate( Dyy0(Nlfd,Nlfd) )
  allocate( Dzz0(Nlfd,Nlfd) )
  allocate( Dyz0(Nlfd,Nlfd) )
  allocate( Dzy0(Nlfd,Nlfd) )

  allocate( Pixx0(Nlfd,Nlfd) )
  allocate( Piyy0(Nlfd,Nlfd) )
  allocate( Piyz0(Nlfd,Nlfd) )
  allocate( Pizy0(Nlfd,Nlfd) )
  allocate( Pizz0(Nlfd,Nlfd) )
  allocate( Pixx(Nlfd,Nlfd) )
  allocate( Piyy(Nlfd,Nlfd) )
  allocate( Piyz(Nlfd,Nlfd) )
  allocate( Pizy(Nlfd,Nlfd) )
  allocate( Pizz(Nlfd,Nlfd) )

  allocate( TS(Nlfd,Nlfd) )
  allocate( TP(Nlfd,Nlfd) )
  ! allocate( TScheck(Nlfd,Nlfd) )
  ! allocate( TPcheck(Nlfd,Nlfd) )

  dtheta = pi/(2.0*dble(Ntheta))
  domega = (omax-omin) / (No-1)

  ! KC transformation matrix from rec.cryst. axes to cart.koord.
  ! if G' is vector in rec.cryst. axes then a = KC*a' is vector in cart. axes
  path = trim(rundir)//"/"//trim(scf_file)
  call loadKC(path,KC)
  print *,"status: KC transformation matrix (rec.cryst.->cart.) loaded."


  ! Reading the reciprocal vectors in crystal coordinates and transformation
  ! in Cartesian cordinates.
  ! call loadG(NG,KC,G)
  call loadG_QE6(savedir, KC, NG, parG, G)
  print *,"status: G vectors loaded. NG=",NG

  ! Reciprocal vectors for crystal local field effects calculations in array Glf(3,Nlf)
  call genGlfandParity(lf,Ecut,NG,Gcar,G, parG, Nlf,Nlfd,Glf)
  print *, "Glf matrix generated."
  print *, 'Nlf: ',Nlf,' Nlfd: ',Nlfd

  ! convert Glf from crystal to Cartesian coords.
  do iG = 1,Nlf 
    Glf(:,iG) = Gcar*Glf(:,iG)
  enddo

  ! Nlf2 = 2*Nlf
  ! print *, 'Nlf2 (p-mode): ',Nlf2

  !*******************************************************************
  !   START ITERATING OVER TRANSFER WAVEVECTORS AND FREQUNECIES
  !******************************************************************
   
  angle_loop: do itheta = 1,Ntheta
    theta = (itheta-1)*dtheta
    print *, 'theta = ',180.0*theta/pi,'°'

    q_loop: do iq =  qmin, qmax
      print*, ' iq: ', iq
    
      q = (iq-1)*dq

      omega_loop: do io = 1,No - 1
        o = omin + (io-1)*domega  
        ! print*,'omega = ', o*Hartree,' [eV]'
        if(itheta /= 1) then 
          q = gamma*o*sin(theta)
        elseif (Ntheta /=1 .and. (qmax-qmin)/=0) then
          print *,'FATAL ERROR: Ntheta is >1 while qmin!=qmax. Choose one or the other.'
          stop
        else
          q = (iq-1)*dq
        endif
        
        oi = cmplx(o,eta) 
        beta = cmplx(gamma**2 * oi**2 - q**2)
        beta = sqrt(beta)

        ! matrix elements of free photon propagator D^0_{\mu\nu}(G,G')
        call genD0(q, oi, beta, c0, Nlf, parG, Glf, Dxx0, Dyy0, Dzz0 ,Dyz0, Dzy0)

        ! Reading unscreened current current response tensor Pi^0_{\mu\nu}(G,G')
        call readPi0(io, No, Nlf, rpa_xx_file, rpa_yy_file, rpa_zz_file, Pixx0, Piyy0, Piyz0, Pizy0, Pizz0)

        ! KONSTRUKCIJA TS MATRICE S - MOD
        call genTS(Nlf, Dxx0, Pixx0, TS)

        ! invertiranje matrice TS
        ! call gjel(TS,Nlf,Nlfd,Imat,Nlf,Nlfd)
       call invert(TS)


        ! KONSTRUKCIJA TP MATRICE P - MOD
        call genTP(Nlf, Dyy0, Dyz0, Dzy0, Dzz0, Piyy0, Piyz0, Pizy0, Pizz0, TP)
        

        ! invertiranje matrice TP
        ! call gjel(TP,Nlf2,Nlfd,Imat2,Nlf2,Nlfd)
        call invert(TP)


        ! SCREENED CURRENT - CURRENT MATRICES S-mod: 'Pixx'; P-mod:'Piyy, Piyz, Pizy, Pizz' 
        call genScreenedPi(Nlf, TS, TP, Pixx0, Piyy0, Piyz0, Pizy0, Pizz0,Pixx, Piyy, Piyz, Pizy, Pizz)
        

        ! DEBUG: mozda prepraviti?
        ! unscreened conductivity [ pi*e^2/2h ]
        write(301,*)o*Hartree, real(-cmplx(0.0,1.0)*4.0*c0*Pixx0(1,1)/o)
        write(302,*)o*Hartree, real(-cmplx(0.0,1.0)*4.0*c0*Piyy0(1,1)/o)
        write(303,*)o*Hartree, real(-cmplx(0.0,1.0)*4.0*c0*Pizz0(1,1)/o)

        !***********************************
        !       MACROSCOPIC CONDUCTIVITY 
        !*********************************** 


        ! KONSTRUKCIJA TS MATRICE S - MOD

        call genTS(Nlf, Dxx0, Pixx, TS)

        ! invertiranje matrice TS
        ! call gjel(TS,Nlf,Nlfd,Imat,Nlf,Nlfd)
        call invert(TS)


        ! KONSTRUKCIJA TP MATRICE P - MOD
        call genTP(Nlf, Dyy0, Dyz0, Dzy0, Dzz0, Piyy, Piyz, Pizy, Pizz, TP)


        ! invertiranje matrice TP
        ! call gjel(TP,Nlf2,Nlfd,Imat2,Nlf2,Nlfd)
        call invert(TP)


        ! MAKROSKOPSKI SIGMA
        call writeSigma_macroscopic(o, c0, Nlf, Pixx, Piyy, Pizz, TS, TP)


      enddo omega_loop 
    enddo q_loop
  enddo angle_loop

  ! deallocate all arrays
  if (allocated(Dxx0))  deallocate(Dxx0)
  if (allocated(Dyy0))  deallocate(Dyy0)
  if (allocated(Dzz0))  deallocate(Dzz0)
  if (allocated(Dyz0))  deallocate(Dyz0)
  if (allocated(Dzy0))  deallocate(Dzy0)
  if (allocated(Pixx0)) deallocate(Pixx0)
  if (allocated(Piyy0)) deallocate(Piyy0)
  if (allocated(Piyz0)) deallocate(Piyz0)
  if (allocated(Pizy0)) deallocate(Pizy0)
  if (allocated(Pizz0)) deallocate(Pizz0)
  if (allocated(TP))    deallocate(TP)
  if (allocated(TS))    deallocate(TS)
  if (allocated(G))     deallocate(G)
  if (allocated(Glf))   deallocate(Glf)
  if (allocated(KC))    deallocate(KC)
  if (allocated(parG))  deallocate(parG)
  if (allocated(Gi))    deallocate(Gi)
  if (allocated(Pixx))  deallocate(Pixx)
  if (allocated(Piyy))  deallocate(Piyy)
  if (allocated(Piyz))  deallocate(Piyz)
  if (allocated(Pizy))  deallocate(Pizy)
  if (allocated(Pizz))  deallocate(Pizz)

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

  subroutine invert(A)
      implicit none

      ! integer,          intent(in)    :: Nlf
      ! integer,          intent(in)    :: Nlfd
      complex(kind=dp), intent(inout) :: A(:,:)

      integer :: Nlf
      complex(kind=dp), allocatable :: Acheck(:,:)


      ! MKL vars
      integer :: info_trf, info_tri
      integer :: lwork
      integer, allocatable :: ipiv(:)
      integer, allocatable :: work(:)

      ! DEBUG: dgemm vars (inversion success check)
      integer :: K
      real(kind=dp)    :: alpha = 1.0
      real(kind=dp)    :: beta = 0.0
      complex(kind=dp) :: checkIdentity

      Nlf = size(A,1)
      lwork = Nlf
      K = Nlf      
      allocate(Acheck(size(A,1),size(A,2)))
      Acheck = A

      allocate(ipiv(max(1,Nlf)))
      call zgetrf( Nlf, Nlf, A, Nlf, ipiv, info_trf)
      if (info_trf/=0) then
        print *, 'FATAL ERROR: LU decomposition failed.'
        print *, 'info_trf: ',info_trf
        stop
      endif

      allocate(work(max(1,lwork)))
      call zgetri( Nlf, A, Nlf, ipiv, work, lwork, info_tri )

      if (info_tri/=0) then
        print *, 'FATAL ERROR: Matrix inversion failed.'
        print *, 'info_trf: ',info_tri
        stop
      endif
      
      deallocate(ipiv)
      deallocate(work)

      ! DEBUG: provjera inverzije 
      call zgemm('N','N', Nlf, Nlf, K, alpha, A, Nlf, Acheck, K, beta, Acheck, Nlf)
      checkIdentity = sum(abs(A)) - sum( (/ ( abs(a(i,i)), i=1, size(a, 1)) /) )
      if (real(checkIdentity)>10d-8 .or. imag(checkIdentity)>10d-8) then
        print *, 'FATAL ERROR: Matrix inversion failed.'
        stop
      endif

      deallocate(Acheck)

  end subroutine invert

  subroutine loadKC(path,KC)
    character(len=100), intent(in)  :: path
    real(kind=dp),      intent(out) :: KC(3,3)
    
    integer           :: j
    integer           :: ios, lno=0
    character(len=35) :: tag, buffer

    tag='     reciprocal axes: (cart. coord.'

    open(300,FILE=path,status='old')
    do  i = 1,100000
      read(300,'(a) ') buffer
      lno = lno+1
      if (buffer == tag) then
        do  j = 1,3
          read(300,'(23X,3F10.3) ',err=10001,iostat=ios,end=20001) KC(1,j), KC(2,j), KC(3,j)
        end do
        EXIT
      end if
    end do
    close(300)
 
    goto 8000
    10001   write(*,*) 'Error reading line ',lno+1,', iostat = ',ios
    20001   write(*,*) 'Number of lines read = ',lno
    8000 continue
  
  end subroutine loadKC

  subroutine loadG_QE6(savedir,KC,NG,parG,G)
    ! Reading the reciprocal vectors in crystal coordinates and transformation
    ! in Cartesian cordinates.
    implicit none
    character(len=*), intent(in)   :: savedir
    real(kind=dp),    intent(in)   :: KC(3,3)
    integer,          intent(out)  :: parG(:) ! paritet svakog valnog vektora G
    integer,          intent(out)  :: NG
    real(kind=dp), allocatable, intent(out)  :: G(:,:)  ! polje valnih vektora G u recp. prost. za wfn.

    integer :: ios0, ios1, ios2
    integer :: iuni
    integer :: n, m
    integer :: Nspin
    logical :: gamma_only
    integer, allocatable :: Gi(:,:)
    character(len=218)   :: fname      

    fname = trim(savedir)//'/charge-density.dat'
    print *,'status: Reading Gvecs from file: ',adjustl(trim(fname))
    
    open(newunit=iuni,file=fname,form = 'unformatted',status='old',iostat=ios0,err=199)

    read(iuni, iostat=ios1) gamma_only, NG, Nspin
    ! print *, 'Number of Gvecs (NG):', NG

    read(iuni, iostat=ios1) ! dummy for b1, b2, b3 rec.latt.vecs. 

    allocate(Gi(3,NG))
    read (iuni, iostat=ios2) Gi(1:3,1:NG)
    close(iuni)

    if (Gi(1,1) /= 0  .or.  Gi(2,1) /= 0  .or.  Gi(3,1) /= 0) then
      print*,'WARNING G vectors input is wrong.'
      print*,'G(1:3,1) is not (0,0,0)!'
      stop
    end if

    ! transformation to cart.coord (now all G components are in 2pi/a0 units)
    allocate(G(3,NG))
    G(1:3,1:NG) = 0.0
    do iG=1,NG
      do n = 1,3
        do m = 1,3
          G(n,iG) = G(n,iG) + KC(n,m)*dble( Gi(m,iG) )
        end do
      end do
      parG(iG) = Gi(3,iG)
    end do
    close(iuni)

    deallocate(Gi)

    goto 5000
    199 write(*,*) 'error cant open file id 199, ist=',ios0
    stop
    200 write(*,*) 'error cant read file id 200, ist=',ios1
    stop
    201   write(*,*) '201 buffer1 read. Error reading miller indices ',', iostat = ',ios2
    stop
    202   write(*,*) '202 buffer1 read.'
    5000 continue 
  end subroutine loadG_QE6

  subroutine genGlfandParity(lf,Ecut,NG,Gcar,G, parG, Nlf,Nlfd,Glf)
    ! Generate Reciprocal vectors for crystal local field 
    ! effects calculations in array Glf(3,Nlf)
    implicit none
    character(len=*), intent(in)    :: lf
    integer,          intent(in)    :: NG
    real(kind=dp),    intent(in)    :: Ecut
    real(kind=dp),    intent(in)    :: Gcar
    real(kind=dp),    intent(in)    :: G(:,:)
    integer,          intent(inout) :: Nlfd
    integer,          intent(inout) :: parG(:)
    integer,          intent(out)   :: Nlf
    real(kind=dp),    intent(out)   :: Glf(:,:)
  
    integer     :: iG
    real(kind=dp) :: Eref
  
    Nlf = 0
    if (lf == 'z') then
      ! local field efekti samo u okomitom smjeru (z)
      do iG = 1, NG
        if (G(1,iG) == 0.0 .and. G(2,iG) == 0.0) then
          Eref = Gcar**2 * G(3,iG)**2 / 2.0
          if (Eref <= Ecut) then
            Nlf = Nlf + 1
            Glf(1:2,Nlf) = 0.0
            Glf(3,Nlf) = G(3,iG)
            if ( (parG(iG)/2)*2 == parG(iG) ) then
              parG(Nlf) = 1
            else
              parG(Nlf) = -1
            endif
          end if
        end if
      end do
    else if (lf == 'xyz') then
      ! local field efekti samo u svim smjerovima (xyz)
      do  iG = 1, NG
        Eref = Gcar**2 * sum(G(1:3,iG)**2) / 2.0
        if (Eref <= Ecut) then
          Nlf = Nlf+1
          Glf(1:3,Nlf) = G(1:3,iG)
          if ( (parG(iG)/2)*2 == parG(iG) ) then
            parG(Nlf) = 1
          else
            parG(Nlf) = -1
          end if
        end if
      end do
    end if
    if (Nlf > Nlfd) then
      print*,'Nlf is bigger than Nlfd'
      stop
    else if(2*Nlf<Nlfd) then
      ! DEBUG this line is different than in other programs where Nlfd=Nlf
      Nlfd = 2*Nlf
    end if
  
  end subroutine genGlfandParity
  
  subroutine genD0(q, oi, beta, c0, Nlf, parG, Glf, Dxx0, Dyy0, Dzz0 ,Dyz0, Dzy0)
    ! matrix elements of free photon propagator D^0_{\mu\nu}(G,G')
    
    implicit none
    integer,            intent(in)   :: Nlf
    real(kind=dp),      intent(in)   :: c0
    real(kind=dp),      intent(in)   :: q
    complex(kind=dp),   intent(in)   :: oi, beta
    integer,            intent(in)   :: parG(:)
    real(kind=dp),      intent(in)   :: Glf(:,:)
    complex(kind=dp),   intent(out), dimension(:,:) :: Dxx0, Dyy0, Dzz0 ,Dyz0, Dzy0
    
    complex(kind=dp) :: Aii, Ajj, Aij, Bij, diagonal1, diagonal2
    
    
    integer :: iG,jG
    real(kind=dp),  parameter :: pi    = 4.d0*atan(1.d0)
    real(kind=dp),  parameter :: gamma = 1.0/137.0
    
    do iG = 1,Nlf 
      Aii = parG(iG) * (cmplx(1.0,0.0) &
          & - exp(cmplx(0.0,1.0)*beta*c0)) / (beta**2 - Glf(3,iG)**2)
      do jG = 1,Nlf
        Ajj = Aii*parG(jG) / ( beta**2 - Glf(3,jG)**2 )
        Aij = 2.0 * ( beta**2 + Glf(3,iG)*Glf(3,jG) )* Ajj / c0
        Bij = 2.0 * beta*( Glf(3,iG) + Glf(3,jG) ) * Ajj / c0
        if(jG == iG) then
          diagonal1 = 2.0*cmplx(0.0,1.0)*beta/(beta**2 - Glf(3,iG)**2)
          diagonal2 = 2.0*cmplx(0.0,1.0)*Glf(3,iG) &
                    & / (beta**2 - Glf(3,iG)**2)       
        else 
          diagonal1 = cmplx(0.0,0.0)
          diagonal2 = cmplx(0.0,0.0) 
        endif
        ! matricni element Dxx0 i Dyy0
        Dxx0(iG,jG) = 2.0*pi*cmplx(0.0,1.0) * gamma**2 *(Aij + diagonal1)/beta 
        Dyy0(iG,jG) = 2.0*pi*cmplx(0.0,1.0) * beta * (Aij + diagonal1)/(oi**2)
        ! matricni elementi Dyz0 i Dzy0
        Dyz0(iG,jG) = -2.0*pi*cmplx(0.0,1.0)*Q*(Bij + diagonal2) / (oi**2)
        Dzy0(iG,jG) = Dyz0(iG,jG)      
        ! matricni element Dzz0
        Dzz0(iG,jG) = 2.0*pi*cmplx(0.0,1.0)* Q**2 *(Aij + diagonal1) / (beta*oi**2)
      enddo
      Dzz0(iG,iG) = Dzz0(iG,iG)-4.0*pi/(oi**2)
    enddo
  end subroutine genD0
  
  subroutine readPi0(io, No, Nlf, file_xx, file_yy, file_zz, Pixx0, Piyy0, Piyz0, Pizy0, Pizz0)
    ! Reading unscreened  current current response tensor Pi^0_{\mu\nu}(G,G')
    ! DEBUG: mozda izbaciti iz omega_loopa i ucitati ih sve u RAM odjednom
    implicit none

    integer,            intent(in) :: io, No, Nlf
    character(len=200), intent(in) :: file_xx, file_yy, file_zz
    complex(kind=dp),   intent(out),  dimension(:,:) :: Pixx0, Piyy0, Piyz0, Pizy0, Pizz0


    integer :: io2, iG, jG
    integer :: iuni1, iuni2, iuni3

    open(newunit=iuni1,file=adjustl(trim(file_xx)))
    open(newunit=iuni2,file=adjustl(trim(file_yy)))
    open(newunit=iuni3,file=adjustl(trim(file_zz)))

    do io2 = 1,No - 1
      if (io2==io) then
        read(iuni1,*) ! skip line
        read(iuni1,'(10F15.10)')((Pixx0(iG,jG),jG = 1,Nlf),iG = 1,Nlf)
        read(iuni2,*) ! skip line
        read(iuni2,'(10F15.10)')((Piyy0(iG,jG),jG = 1,Nlf),iG = 1,Nlf)
        read(iuni3,*) ! skip line
        read(iuni3,'(10F15.10)')((Pizz0(iG,jG),jG = 1,Nlf),iG = 1,Nlf)

        ! approximation current current response tensor components Pi^0_{yz}=Pi^0_{zy} = 0 
        Piyz0(:,:) = cmplx(0.0,0.0)        
        Pizy0(:,:) = cmplx(0.0,0.0)
        exit
      else
        read(iuni1,*) ! skip line
        read(iuni2,*)
        read(iuni3,*)
        do iG=1,Nlf
          do jG=1,Nlf
            read(iuni1,*) ! skip line
            read(iuni2,*) ! skip line
            read(iuni3,*) ! skip line
          enddo
        enddo
      endif

    enddo  
    close(iuni1)
    close(iuni2)
    close(iuni3)
      
  end subroutine readPi0


  subroutine genScreenedPi(Nlf,TS, TP,Pixx0, Piyy0, Piyz0, Pizy0, Pizz0,Pixx, Piyy, Piyz, Pizy, Pizz)
    ! SCREENED CURRENT - CURRENT MATRICES S-mod: 'Pixx'; P-mod:'Piyy, Piyz, Pizy, Pizz' 
    implicit none
  
    integer,          intent(in) :: Nlf
    complex(kind=dp), intent(in),  dimension(:,:) :: TS, TP
    complex(kind=dp), intent(in),  dimension(:,:) :: Pixx0, Piyy0, Piyz0, Pizy0, Pizz0
    complex(kind=dp), intent(out), dimension(:,:) :: Pixx, Piyy, Piyz, Pizy, Pizz
  
    integer :: iG, jG, kG
  
    Pixx = cmplx(0.0,0.0) 
    Piyy = cmplx(0.0,0.0)
    Piyz = cmplx(0.0,0.0)
    Pizy = cmplx(0.0,0.0)
    Pizz = cmplx(0.0,0.0)
    do iG = 1,Nlf
      do jG = 1,Nlf       
        do kG = 1,Nlf
          ! Pi_{xx} (s-mod)
          Pixx(iG,jG) = Pixx(iG,jG) + TS(iG,kG)*Pixx0(kG,jG)
          
          ! Pi_{yy} (p-mod)
          Piyy(iG,jG) = Piyy(iG,jG) + TP(iG,kG)*Piyy0(kG,jG) &
                      & + TP(iG,kG+Nlf)*Pizy0(kG,jG)
  
          ! Pi_{yz} (p-mod)
          Piyz(iG,jG) = Piyz(iG,jG) + TP(iG,kG)*Piyz0(kG,jG) &
                      & + TP(iG,kG+Nlf)*Pizz0(kG,jG)
  
          ! Pi_{zy} (p-mod)
          Pizy(iG,jG) = Pizy(iG,jG) + TP(iG+Nlf,kG)*Piyy0(kG,jG) &
                      & + TP(iG+Nlf,kG+Nlf)*Pizy0(kG,jG)
          ! Pi_{zz} (p-mod)
          Pizz(iG,jG) = Pizz(iG,jG) + TP(iG+Nlf,kG)*Piyz0(kG,jG) &
                      & + TP(iG+Nlf,kG+Nlf)*Pizz0(kG,jG)
        enddo
      enddo
    enddo
  end subroutine genScreenedPi

  subroutine genTS(Nlf, Dxx0, Pixx, TS)
    implicit none

    integer,          intent(in)  :: Nlf
    complex(kind=dp), intent(in)  :: Dxx0(:,:)
    complex(kind=dp), intent(in)  :: Pixx(:,:)
    complex(kind=dp), intent(out) :: TS(:,:)

    integer :: iG, jG !, kG
    complex(kind=dp), allocatable :: Imat(:,:)

    allocate(Imat(Nlf,Nlf))

    Imat = cmplx(0.0,0.0)
    TS = cmplx(0.0,0.0)
    do iG = 1,Nlf
      Imat(iG,iG) = cmplx(1.0,0.0)
      do jG = 1,Nlf    
        ! do kG = 1,Nlf   
        TS(iG,jG) = sum( Pixx(iG,:)*Dxx0(:,jG) )
        ! enddo
      enddo
    enddo
    TS = Imat - TS

    deallocate(Imat)

  end subroutine genTS

  subroutine genTP(Nlf, Dyy0, Dyz0, Dzy0, Dzz0, Piyy, Piyz, Pizy, Pizz, TP)
    implicit none

    integer,          intent(in)  :: Nlf
    complex(kind=dp), intent(in),  dimension(:,:)  :: Dyy0, Dyz0, Dzy0, Dzz0
    complex(kind=dp), intent(in),  dimension(:,:)  :: Piyy, Piyz, Pizy, Pizz
    complex(kind=dp), intent(out), dimension(:,:)  :: TP


    integer :: Nlf2  
    integer :: iG, jG !, kG
    complex(kind=dp), allocatable :: Imat(:,:)

    Nlf2 = 2*Nlf ! TP matrix has a 2x2 block for z and y

    allocate(Imat(Nlf2,Nlf2))

    Imat(1:Nlf2,1:Nlf2) = cmplx(0.0,0.0)
    TP(1:Nlf2,1:Nlf2) = cmplx(0.0,0.0) 
    do iG = 1,Nlf2
      Imat(iG,iG) = cmplx(1.0,0.0)
    enddo
    

    do iG = 1,Nlf2
      do jG = 1,Nlf2 
        if(iG <= Nlf .and. jG <= Nlf) then 
          ! do kG = 1,Nlf   
          TP(iG,jG) = sum(Piyy(iG,:)*Dyy0(:,jG) + Piyz(iG,:)*Dzy0(:,jG) )
          ! enddo
        elseif (iG <= Nlf.and.jG > Nlf) then
          ! do kG = 1,Nlf
          TP(iG,jG) = sum( Piyy(iG,:)*Dyz0(:,jG-Nlf) + Piyz(iG,:)*Dzz0(:,jG-Nlf) )
          ! enddo
        elseif (iG > Nlf .and. jG <= Nlf) then
          ! do kG = 1,Nlf
          TP(iG,jG) = sum( Pizy(iG-Nlf,:)*Dyy0(:,jG) + Pizz(iG-Nlf,:)*Dzy0(:,jG) )
          ! enddo
        else
          ! do  kG = 1,Nlf
          TP(iG,jG) = sum( Pizy(iG-Nlf,:)*Dyz0(:,jG-Nlf) + Pizz(iG-Nlf,:)*Dzz0(:,jG-Nlf) )
          ! enddo
        endif
      enddo
    enddo
    TP = Imat - TP    
      
  end subroutine genTP

  subroutine writeSigma_macroscopic(o, c0, Nlf, Pixx, Piyy, Pizz, TS, TP)
    ! Calculates macroscopic conducitivities from current-current response functions Pi_{\mu\nu}
    ! Write the conducitivities and response functions to files.

    implicit none
    integer,          intent(in) :: Nlf
    real(kind=dp),    intent(in) :: c0
    real(kind=dp),    intent(in) :: o
    complex(kind=dp), intent(in), dimension(:,:) ::  Pixx, Piyy, Pizz, TS, TP

    
    real(kind=dp),  parameter :: Hartree = 2.0d0*13.6056923D0


    integer :: Nlf2
    integer :: iuni1, iuni2, iuni3, iuni4, iuni5, iuni6
    complex(kind = dp) :: Sigma_xx,Sigma_yy,Sigma_zz

    Nlf2 = 2*Nlf

    Sigma_xx = sum( Pixx(1,1:Nlf) * TS(1:Nlf,1) )                                        ! s-mod
    Sigma_yy = sum( Piyy(1,1:Nlf) * TP(1:Nlf,1)     + Piyz(1,1:Nlf)*TP(Nlf:Nlf2,1) )     ! p-mod
    Sigma_zz = sum( Pizy(1,1:Nlf) * TP(1:Nlf,Nlf+1) + Pizz(1,1:Nlf)*TP(Nlf:Nlf2,Nlf+1) ) ! p-mod

    ! write sigma_{\mu\nu}^\text{macro} [ pi*e^2/2h ]
    open(newunit=iuni1, file='sigma_macro_xx')
    write(iuni1,*) o*Hartree, real(-cmplx(0.0,1.0)*4.0*c0*Sigma_xx/o)
    close(iuni1)

    open(newunit=iuni2, file='sigma_macro_yy')
    write(iuni2,*) o*Hartree, real(-cmplx(0.0,1.0)*4.0*c0*Sigma_yy/o)
    close(iuni2)

    open(newunit=iuni3, file='sigma_macro_zz')
    write(iuni3,*) o*Hartree, real(-cmplx(0.0,1.0)*4.0*c0*Sigma_zz/o)
    close(iuni3)

    ! write Pi_{\mu\nu}(1,1)
    open(newunit=iuni4, file='Pi_11_xx')
    write(iuni4,*) o*Hartree, real(-cmplx(0.0,1.0)*4.0*c0*Pixx(1,1)/o)
    close(iuni4)

    open(newunit=iuni5, file='Pi_11_yy')
    write(iuni5,*) o*Hartree, real(-cmplx(0.0,1.0)*4.0*c0*Piyy(1,1)/o)
    close(iuni5)

    open(newunit=iuni6, file='Pi_11_zz')
    write(iuni6,*) o*Hartree, real(-cmplx(0.0,1.0)*4.0*c0*Pizz(1,1)/o)
    close(iuni6)

      
  end subroutine writeSigma_macroscopic

end program photon


