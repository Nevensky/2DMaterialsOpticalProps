program photon
  !  program for the calculation of the current - current corelation function in P
  !  XY - plane  Imat cell parameter  a0=  8.7521  a.u.  
  !  Z - dir Imat cell parameter c0=  32.33282 a.u.  
  use OMP_lib
  use ISO_Fortran_env
  use utility, only: int2str
  !use mkl_service

  implicit none

  ! single / double precision
  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
  integer, parameter :: qp = real128

  integer :: parallelCount = 0
  integer :: i, io, iq,itheta, iG, jG, kG, Nlf, Nlf2
  integer :: ist1, ist2, ist3, ist4, ist5, ist6
  integer :: Ntheta, Nq ! number of angles, number of transfer wavevectors

  integer :: No     ! number of frequencies
  integer :: NG     ! total number of G vectors 
  integer :: NGd    ! minimum number of coefficients CG over all evc.n files
  integer :: NkI    ! number of wave vectors in IBZ
  integer :: Nk     ! = 48*NkI, number of wave vectors in FBZ with No symmetry 
  integer :: Nlfd   ! dimenzija polja local field vektora, x2 zbog 2X2 blok matrice za p mod 
  namelist /config/  No, NG, NGd, NkI, Nlfd


  ! ... constants ...
  real(kind=dp),  parameter :: pi      = 4.0_dp*atan(1.0_dp)
  real(kind=dp),  parameter :: eV      = 1.602176487d-19
  real(kind=dp),  parameter :: kB      = 1.3806503d-23
  real(kind=dp),  parameter :: Hartree = 2.0_dp*13.6056923_dp
  real(kind=dp),  parameter :: Planck  = 6.626196d-34
  real(kind=dp),  parameter :: aBohr   = 0.5291772_dp
  real(kind=dp),  parameter :: gamma   = 1.0_dp/137.0_dp

  ! ...scalars...

  real(kind=dp) :: Gcar
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
  complex(kind=dp), dimension(:,:), allocatable :: Dxx, Dyy, Dzz ,Dyz, Dzy
  complex(kind=dp), dimension(:,:,:), allocatable :: Pixx0, Piyy0, Piyz0, Pizy0, Pizz0
  complex(kind=dp), dimension(:,:), allocatable :: Pixx, Piyy, Piyz, Pizy, Pizz
  complex(kind=dp) :: sPixx, sPiyy, sPiyz, sPizy, sPizz

  character(len=200) :: path

  character(len=200) :: config_file

  character(len=200) :: rundir, savedir, scf_file, rpa_xx_file, rpa_yy_file, rpa_zz_file
  namelist /directories/ rundir, savedir, scf_file, rpa_xx_file, rpa_yy_file, rpa_zz_file

  character(len=3)   :: lf     ! crystal local field effects included in 'xyz' or 'z' direction
  integer            :: qmin, qmax
  integer            :: theta_min, theta_max
  real(kind=dp)      :: dq, omin, omax, eta
  real(kind=dp)      :: Ecut   ! [Hartree] cutoff energy for crystal local field vectors
  namelist /system/ lf, qmin, qmax, theta_min, theta_max, dq, omin, omax, eta, Ecut

  real(kind=dp)      :: a0     ! [a.u.]  unit cell parameter in parallel direction 
  real(kind=dp)      :: c0     ! [a.u.]  unit cell parameter in perependicular direction 
  real(kind=dp)      :: Vcell  ! [a.u.^3] unit-cell volume 

  namelist /parameters/ a0, c0, Vcell

  integer            :: Nthreads, NthreadsMKL
  namelist  /parallel/ Nthreads, NthreadsMKL

  integer, save :: thread_id
  !$omp threadprivate(thread_id)

  print *, "PROGRAM PHOTON RUN STARTED"

  ! load namelist
  call parseCommandLineArgs(config_file) ! config file is the first argument passed to ./Photon <arg1>
  open(10,file=config_file)
  read(10,nml=directories,iostat=ist1)
  read(10,nml=system,iostat=ist2)
  read(10,nml=config,iostat=ist3)
  read(10,nml=parameters,iostat=ist4)
  read(10,nml=parallel,iostat=ist5)
  close(10)

  ! constants
  Nk   = 48*NkI           ! number of wave vectors in FBZ with No symmetry ops.
  Gcar = 2.0_dp*pi/a0   
  eta  = eta/Hartree
  omin = omin/Hartree ! iz eV u Hartree
  omax = (omax/Hartree + omin) 

  Ntheta = theta_max-theta_min
  Nq = qmax-qmin


  ! OpenMP set number of threads
  call omp_set_num_threads(Nthreads)
  ! MKL set number of threads
  call mkl_set_num_threads(NthreadsMKL)

  ! allocation of arrays
  allocate( G(3,NG) )
  allocate( parG(NG) )
  allocate( Glf(3,Nlfd) )
  allocate( KC(3,3) )
  allocate( Gi(3) )

  ! allocate( Imat(Nlfd,Nlfd) ) ! dim NlfxNlf
  ! allocate( Imat2(2*Nlfd,2*Nlfd) ) ! dim Nlf*2 x Nlf*2

  allocate( Dxx0(Nlfd,Nlfd) )
  allocate( Dyy0(Nlfd,Nlfd) )
  allocate( Dzz0(Nlfd,Nlfd) )
  allocate( Dyz0(Nlfd,Nlfd) )
  allocate( Dzy0(Nlfd,Nlfd) )
  allocate( Dxx(Nlfd,Nlfd) )
  allocate( Dyy(Nlfd,Nlfd) )
  allocate( Dzz(Nlfd,Nlfd) )
  allocate( Dyz(Nlfd,Nlfd) )
  allocate( Dzy(Nlfd,Nlfd) )

  allocate( Pixx0(No,Nlfd,Nlfd) )
  allocate( Piyy0(No,Nlfd,Nlfd) )
  allocate( Piyz0(No,Nlfd,Nlfd) )
  allocate( Pizy0(No,Nlfd,Nlfd) )
  allocate( Pizz0(No,Nlfd,Nlfd) )
  allocate( Pixx(Nlfd,Nlfd) )
  allocate( Piyy(Nlfd,Nlfd) )
  allocate( Piyz(Nlfd,Nlfd) )
  allocate( Pizy(Nlfd,Nlfd) )
  allocate( Pizz(Nlfd,Nlfd) )

  allocate( TS(Nlfd,Nlfd) ) ! dim NlfxNlf
  allocate( TP(2*Nlfd,2*Nlfd) )     ! dim Nlf*2 x Nlf*2
  ! allocate( TScheck(Nlfd,Nlfd) )
  ! allocate( TPcheck(2*Nlfd,2*Nlfd) )

  domega = (omax-omin) / (No-1)
  if (Ntheta/=0) then
    dtheta = pi/(2.0_dp*dble(Ntheta))
  else
    dtheta = 1.0_dp
  endif

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
  call genGlfandParity(lf,Ecut,NG,Gcar, G, parG, Nlf, Nlfd, Glf)
  print *, "status: Glf matrix generated."
  print *, 'Nlf: ',Nlf,' Nlfd: ',Nlfd

  ! convert Glf from crystal to Cartesian coords.
  do iG = 1,Nlf 
    Glf(:,iG) = Gcar*Glf(:,iG)
  enddo
  print *, "status: Glf matrix converted to Cartesian coords."

  ! Nlf2 = 2*Nlf
  ! print *, 'Nlf2 (p-mode): ',Nlf2

  ! Reading unscreened current current response tensor Pi^0_{\mu\nu}(G,G')
  call loadPi0(No, Nlf, rpa_xx_file, rpa_yy_file, rpa_zz_file, Pixx0, Piyy0, Piyz0, Pizy0, Pizz0)
  print *, "status: Current-current response tensor Pi0(omega,G,G\') loaded."
  print *, "Pixx0 shape: ",shape(Pixx0)
  print *, "Piyy0 shape: ",shape(Piyy0)
  print *, "Pizz0 shape: ",shape(Pizz0)
  !*******************************************************************
  !   START ITERATING OVER TRANSFER WAVEVECTORS AND FREQUNECIES
  !******************************************************************

  angle_loop: do itheta = theta_min, theta_max
    theta = (itheta-1)*dtheta
    print *, 'theta = ',180.0*theta/pi,'°'

    print *, 'DEBUG: entering parallel region'
    !$omp parallel shared(No,Nq,Ntheta,domega,dq,dtheta,qmin,qmax,omin,omax,eta,c0,Nlf,parG,Glf,Pixx0, Piyy0, Piyz0, Pizy0, Pizz0, parallelCount) firstprivate(itheta) num_threads(Nthreads) default(private)
    !$omp master
    print *, 'Requested threads: ',Nthreads, 'Available threads: ',omp_get_num_threads()
    !$omp end master
    thread_id =  omp_get_thread_num()

    !$omp do    
    do iq =  qmin, qmax ! q_loop: 
      !$omp critical(parallelInfo)
      parallelCount = parallelCount + 1
      if (Nq/=0) then
        print*, 'thread id: ',thread_id, ' iq = ', iq
        write (*,'(A11,I6,A2,I6,A5,F5.1,A4)') 'progress: ',parallelCount, ' /',Nq+1,' (',(real(parallelCount)/real(Nq+1))*100.0,'% )'
      endif
      !$omp end critical(parallelInfo)

      ! write (*,'(A11,I6,A2,I6,A5,F5.1,A4)') 'progress: ',iq+1-qmin, ' /',Nq,' (',(real(iq+1-qmin)/real(Nq))*100.0,'% )'
    
      q = (iq-1)*dq

      omega_loop: do io = 1, No-1
        o = omin + (io-1)*domega  
        ! print*,'omega = ', o*Hartree,' [eV]'
        if (Ntheta /=0 .and. Nq/=0) then
          print *,'FATAL ERROR: Both Ntheta>1 and Nq>1. Choose one or the other.'
          stop
        else if(itheta /= 1 .and. Nq==0) then 
          q = gamma*o*sin(theta)
          ! print *, 'q = γ·ω·sin(𝜗)' ! DEBUG
        else
          q = (iq-1)*dq
          ! print *, 'q = (qᵢ-1)·dq' ! DEBUG
        endif
        
        oi = dcmplx(o,eta) 
        beta = dcmplx(gamma**2 * oi**2 - q**2)
        beta = sqrt(beta)

        ! matrix elements of free photon propagator D^0_{\mu\nu}(G,G')
        call genD0(q, oi, beta, c0, Nlf, parG, Glf, Dxx0, Dyy0, Dzz0 ,Dyz0, Dzy0)

        ! Reading unscreened current current response tensor Pi^0_{\mu\nu}(G,G')
        ! DEBUG: moved outside the loop
        ! call loadPi0_omega(io, No, Nlf, rpa_xx_file, rpa_yy_file, rpa_zz_file, Pixx0, Piyy0, Piyz0, Pizy0, Pizz0)

        ! KONSTRUKCIJA TS MATRICE S - MOD
        call genTS(Nlf, Dxx0, Pixx0(io,:,:), TS)

        ! invertiranje matrice TS
        ! call gjel(TS,Nlf,Nlfd,Imat,Nlf,Nlfd)
        call invert(TS)

        ! KONSTRUKCIJA TP MATRICE P - MOD
        call genTP(Nlf, Dyy0, Dyz0, Dzy0, Dzz0, Piyy0(io,:,:), Piyz0(io,:,:), Pizy0(io,:,:), Pizz0(io,:,:), TP)
        

        ! invertiranje matrice TP
        ! call gjel(TP,Nlf2,2*Nlfd,Imat2,Nlf2,Nlfd)
        call invert(TP)


        ! SCREENED CURRENT - CURRENT MATRICES S-mod: 'Pixx'; P-mod:'Piyy, Piyz, Pizy, Pizz' 
        call genScreenedPi(Nlf, TS, TP, Pixx0(io,:,:), Piyy0(io,:,:), Piyz0(io,:,:), Pizy0(io,:,:), Pizz0(io,:,:), Pixx, Piyy, Piyz, Pizy, Pizz)
        

        ! call sumScreenedPi(c0, Glf, Pixx, Piyy, Piyz, Pizy, Pizz, sPixx, sPiyy, sPiyz, sPizy, sPizz)
        call sumIntScreenedPi(c0, Glf, Pixx, Piyy, Piyz, Pizy, Pizz, sPixx, sPiyy, sPiyz, sPizy, sPizz)

        ! DEBUG: mozda prepraviti?
        ! !$omp critical(writeOutputs)
        ! unscreened current-current response tensor
        call writePi(o, Pixx0(io,1,1),trim("pixx0_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writePi(o, Piyy0(io,1,1),trim("piyy0_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writePi(o, Pizz0(io,1,1),trim("pizz0_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        ! screened current-current response tensor
        call writePi(o, Pixx(1,1),trim("pixx_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writePi(o, Piyy(1,1),trim("piyy_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writePi(o, Pizz(1,1),trim("pizz_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        ! unscreened conductivity (1st component) [ pi*e^2/2h ]
        call writeSigma(o, c0, Pixx0(io,1,1),trim("sigma0_xx_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writeSigma(o, c0, Piyy0(io,1,1),trim("sigma0_yy_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writeSigma(o, c0, Pizz0(io,1,1),trim("sigma0_zz_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        ! screened conductivity (1st component) [ pi*e^2/2h ]
        call writeSigma(o, c0, Pixx(1,1),trim("sigmaSc_xx_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writeSigma(o, c0, Piyy(1,1),trim("sigmaSc_yy_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writeSigma(o, c0, Pizz(1,1),trim("sigmaSc_zz_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))

        ! screened current-current response tensor summed over all Nlf for a symmetric trilayer
        call writePi(o, sPixx,trim("pixx_sumGiGj_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writePi(o, sPiyy,trim("piyy_sumGiGj_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writePi(o, sPizz,trim("pizz_sumGiGj_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))

        ! screened conductivity [ pi*e^2/2h ] summed over all Nlf for a symmetric trilayer
        call writeSigma(o, c0, sPixx,trim("sigmaSc_sumGiGj_xx_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writeSigma(o, c0, sPiyy,trim("sigmaSc_sumGiGj_yy_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        call writeSigma(o, c0, sPizz,trim("sigmaSc_sumGiGj_zz_q#"//int2str(iq))//trim("_theta#"//int2str(itheta)))
        ! !$omp end critical(writeOutputs)

        ! DEBUG doesn't work 
        ! calculation of reflected, transmited and absorbed coefficients 
        call genSpectra(o, oi, beta, iq, Nq, itheta, theta, Ntheta, c0, Nlf, parG, Glf, Pixx, Piyy, Pizz, Piyz, Pizy)


        !***********************************
        !       MACROSCOPIC CONDUCTIVITY 
        !*********************************** 


        ! KONSTRUKCIJA TS MATRICE S - MOD

        ! call genTS(Nlf, Dxx0, Pixx, TS)

        ! invertiranje matrice TS
        ! call gjel(TS,Nlf,Nlfd,Imat,Nlf,Nlfd)
        ! call invert(TS)


        ! KONSTRUKCIJA TP MATRICE P - MOD
        ! call genTP(Nlf, Dyy0, Dyz0, Dzy0, Dzz0, Piyy, Piyz, Pizy, Pizz, TP)


        ! invertiranje matrice TP
        ! call gjel(TP,Nlf2,2*Nlfd,Imat2,Nlf2,Nlfd)
        ! call invert(TP)


        ! MAKROSKOPSKI SIGMA
        ! call writeSigma_macroscopic(o, c0, Nlf, Pixx, Piyy, Pizz, TS, TP)

      enddo omega_loop 
    enddo !q_loop
    !$omp end do
    !$omp end parallel
  enddo angle_loop

  ! deallocate all arrays
  if (allocated(Dxx0))  deallocate(Dxx0)
  if (allocated(Dyy0))  deallocate(Dyy0)
  if (allocated(Dzz0))  deallocate(Dzz0)
  if (allocated(Dyz0))  deallocate(Dyz0)
  if (allocated(Dzy0))  deallocate(Dzy0)
  if (allocated(Dxx))   deallocate(Dxx)
  if (allocated(Dyy))   deallocate(Dyy)
  if (allocated(Dzz))   deallocate(Dzz)
  if (allocated(Dyz))   deallocate(Dyz)
  if (allocated(Dzy))   deallocate(Dzy)
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

  print *, "PROGRAM PHOTON RUN ENDED"

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
      complex(kind=dp), allocatable :: Icheck(:,:)


      ! MKL vars
      integer :: info_trf, info_tri
      integer :: lwork
      integer,          allocatable :: ipiv(:)
      complex(kind=dp), allocatable :: work(:)

      ! DEBUG: dgemm vars (inversion success check)
      integer :: K
      real(kind=dp)    :: alpha = 1.0_dp
      real(kind=dp)    :: beta  = 0.0_dp
      complex(kind=dp) :: checkIdentity, checkIdentity2
      complex(kind=dp) :: traceIdentity

      Nlf = size(A,1)
      lwork = 16*Nlf ! often results in segmentation fault due to lwork being too small
      K = Nlf   
      allocate(Acheck(size(A,1),size(A,2)))
      allocate(Icheck(size(A,1),size(A,2)))
      Acheck = A

      allocate(ipiv(max(1,Nlf)))
      call zgetrf( Nlf, Nlf, A, Nlf, ipiv, info_trf)
      if (info_trf/=0) then
        print *, 'FATAL ERROR: LU decomposition failed.'
        print *, 'info_trf: ', info_trf
        print *, 'Nlf: ', Nlf
        stop
      endif

      allocate(work(max(1,lwork)))
      call zgetri( Nlf, A, Nlf, ipiv, work, lwork, info_tri )

      if (info_tri/=0) then
        print *, 'FATAL ERROR: Matrix inversion failed.'
        print *, 'info_trf: ', info_tri
        print *, 'Nlf: ', Nlf
        stop
      endif
      
      deallocate(ipiv)
      deallocate(work)

      ! DEBUG: check correctness of matrix inversion
      call zgemm('N','N', Nlf, Nlf, K, alpha, A, Nlf, Acheck, K, beta, Icheck, Nlf)
      traceIdentity = sum( (/ ( abs(Icheck(i,i)), i=1, size(Icheck, 1)) /) )
      checkIdentity = sum(abs(Icheck)) - traceIdentity
      checkIdentity2 = dcmplx(Nlf,0.0_dp) - traceIdentity
      if ( dble(checkIdentity)>10d-4 .or. aimag(checkIdentity)>10d-4 .or. &
         & dble(checkIdentity2)>10d-4 .or. aimag(checkIdentity2)>10d-4 ) then
        print *, 'FATAL ERROR: Matrix inversion failed. (A⁻¹ · A ≠ 𝟙)'
        print *, 'Nlf: ', Nlf
        print *, 'size Ainv: ', size(A), 'size A:',size(Acheck)
        print *, 'Re(check): ', dble(checkIdentity)
        print *, 'Im(check): ', aimag(checkIdentity)
        print *, 'Re(check2): ', dble(checkIdentity2)
        print *, 'Im(check2): ', aimag(checkIdentity2)
        stop
      endif

      deallocate(Acheck)
      deallocate(Icheck)

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
        exit
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
          Eref = Gcar**2 * G(3,iG)**2 / 2.0_dp
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
        Eref = Gcar**2 * sum(G(1:3,iG)**2) / 2.0_dp
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
    real(kind=dp),  parameter :: pi    = 4.0_dp*atan(1.0_dp)
    real(kind=dp),  parameter :: gamma = 1.0_dp/137.0_dp
    
    do iG = 1,Nlf 
      Aii = parG(iG) * (dcmplx(1.0_dp,0.0_dp) &
          & - exp(dcmplx(0.0_dp,1.0_dp)*beta*c0)) / (beta**2 - Glf(3,iG)**2)
      do jG = 1,Nlf
        Ajj = Aii*parG(jG) / ( beta**2 - Glf(3,jG)**2 )
        Aij = 2.0_dp * ( beta**2 + Glf(3,iG)*Glf(3,jG) )* Ajj / c0
        Bij = 2.0_dp * beta*( Glf(3,iG) + Glf(3,jG) ) * Ajj / c0
        if(jG == iG) then
          diagonal1 = 2.0_dp*dcmplx(0.0_dp,1.0_dp)*beta/(beta**2 - Glf(3,iG)**2)
          diagonal2 = 2.0_dp*dcmplx(0.0_dp,1.0_dp)*Glf(3,iG) &
                    & / (beta**2 - Glf(3,iG)**2)       
        else 
          diagonal1 = dcmplx(0.0_dp,0.0_dp)
          diagonal2 = dcmplx(0.0_dp,0.0_dp) 
        endif
        ! matricni element Dxx0 i Dyy0
        Dxx0(iG,jG) = 2.0_dp*pi*dcmplx(0.0_dp,1.0_dp) * gamma**2 *(Aij + diagonal1)/beta 
        Dyy0(iG,jG) = 2.0_dp*pi*dcmplx(0.0_dp,1.0_dp) * beta * (Aij + diagonal1)/(oi**2)
        ! matricni elementi Dyz0 i Dzy0
        Dyz0(iG,jG) = -2.0_dp*pi*dcmplx(0.0_dp,1.0_dp)*Q*(Bij + diagonal2) / (oi**2)
        Dzy0(iG,jG) = Dyz0(iG,jG)      
        ! matricni element Dzz0
        Dzz0(iG,jG) = 2.0_dp*pi*dcmplx(0.0_dp,1.0_dp)* Q**2 *(Aij + diagonal1) / (beta*oi**2)
      enddo
      Dzz0(iG,iG) = Dzz0(iG,iG)-4.0_dp*pi/(oi**2)
    enddo
  end subroutine genD0

subroutine loadPi0(No, Nlf, file_xx, file_yy, file_zz, Pixx0, Piyy0, Piyz0, Pizy0, Pizz0)
    ! Reading unscreened  current current response tensor Pi^0_{\mu\nu}(G,G')
    ! DEBUG: mozda izbaciti iz omega_loopa i ucitati ih sve u RAM odjednom
    implicit none

    integer,            intent(in) :: No, Nlf
    character(len=200), intent(in) :: file_xx, file_yy, file_zz
    complex(kind=dp),   intent(out),  dimension(:,:,:) :: Pixx0, Piyy0, Piyz0, Pizy0, Pizz0

    integer       :: io, iG, jG
    integer       :: iuni1, iuni2, iuni3
    real(kind=dp) :: o ! frequency tmp var, not used

    open(newunit=iuni1,file=adjustl(trim(file_xx)))
    open(newunit=iuni2,file=adjustl(trim(file_yy)))
    open(newunit=iuni3,file=adjustl(trim(file_zz)))
    do io = 1,No
      do iG = 1,Nlf
        do jG = 1,Nlf
            read(iuni1,*) o, Pixx0(io,iG,jG)
            read(iuni2,*) o, Piyy0(io,iG,jG)
            read(iuni3,*) o, Pizz0(io,iG,jG)
        enddo
      enddo
    enddo
    close(iuni1)
    close(iuni2)
    close(iuni3)

    Piyz0(:,:,:) = dcmplx(0.0_dp,0.0_dp)
    Pizy0(:,:,:) = dcmplx(0.0_dp,0.0_dp)
  end subroutine loadPi0
  
  ! subroutine loadPi0_omega(io_in, No, Nlf, file_xx, file_yy, file_zz, Pixx0, Piyy0, Piyz0, Pizy0, Pizz0)
  !   ! Reading unscreened  current current response tensor Pi^0_{\mu\nu}(G,G')
  !   ! DEBUG: mozda izbaciti iz omega_loopa i ucitati ih sve u RAM odjednom
  !   implicit none

  !   integer,            intent(in) :: io_in, No, Nlf
  !   character(len=200), intent(in) :: file_xx, file_yy, file_zz
  !   complex(kind=dp),   intent(out),  dimension(:,:) :: Pixx0, Piyy0, Piyz0, Pizy0, Pizz0


  !   integer       :: io, iG, jG
  !   integer       :: iuni1, iuni2, iuni3
  !   real(kind=dp) :: o ! frequency tmp var, not used

  !   open(newunit=iuni1,file=adjustl(trim(file_xx)))
  !   open(newunit=iuni2,file=adjustl(trim(file_yy)))
  !   open(newunit=iuni3,file=adjustl(trim(file_zz)))

  !   do io = 1,No
  !     if (io==io_in) then
  !       do iG = 1,Nlf
  !         do jG = 1,Nlf
  !             read(iuni1,*) o, Pixx0(iG,jG)
  !             read(iuni2,*) o, Piyy0(iG,jG)
  !             read(iuni3,*) o, Pizz0(iG,jG)
  !         enddo
  !       enddo
  !     else
  !       do iG=1,Nlf
  !         do jG=1,Nlf
  !             ! skip lines
  !             read(iuni1,*) 
  !             read(iuni2,*) 
  !             read(iuni3,*) 
  !         enddo
  !       enddo
  !     endif
  !   enddo
  !   Piyz0(:,:) = dcmplx(0.0_dp,0.0_dp)
  !   Pizy0(:,:) = dcmplx(0.0_dp,0.0_dp)
  !   close(iuni1)
  !   close(iuni2)
  !   close(iuni3)
  ! end subroutine loadPi0_omega


  subroutine genScreenedPi(Nlf, TS, TP, Pixx0, Piyy0, Piyz0, Pizy0, Pizz0, Pixx, Piyy, Piyz, Pizy, Pizz)
    ! SCREENED CURRENT - CURRENT MATRICES S-mod: 'Pixx'; P-mod:'Piyy, Piyz, Pizy, Pizz' 
    implicit none
  
    integer,          intent(in) :: Nlf
    complex(kind=dp), intent(in),  dimension(:,:) :: TS, TP
    complex(kind=dp), intent(in),  dimension(:,:) :: Pixx0, Piyy0, Piyz0, Pizy0, Pizz0
    complex(kind=dp), intent(out), dimension(:,:) :: Pixx, Piyy, Piyz, Pizy, Pizz
  

    integer :: iG, jG, kG
    integer :: Nlf2

    Nlf2 = 2*Nlf ! TP matrix has a 2x2 block for z and y
  
    Pixx = dcmplx(0.0_dp,0.0_dp) 
    Piyy = dcmplx(0.0_dp,0.0_dp)
    Piyz = dcmplx(0.0_dp,0.0_dp)
    Pizy = dcmplx(0.0_dp,0.0_dp)
    Pizz = dcmplx(0.0_dp,0.0_dp)
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

  subroutine sumScreenedPi(c0, Glf, Pixx, Piyy, Piyz, Pizy, Pizz, sPixx, sPiyy, sPiyz, sPizy, sPizz)
    ! non general function, works only for a symmetric trilayer
    ! sums Pi over all local field vectors to observe sigma^sc in upper layer
    ! z=d/2, z'=d/2
    ! Pi = 1/L \sum_Gz \sum_Gz' Pi(Gz,Gz') * exp(i Gz * d/2) * exp(-i Gz * d/2) 
    implicit none
    real(kind=dp),    intent(in)  :: c0
    real(kind=dp),    intent(in)  :: Glf(:,:)
    complex(kind=dp), intent(in), dimension(:,:) :: Pixx, Piyy, Piyz, Pizy, Pizz
    complex(kind=dp), intent(out) :: sPixx, sPiyy, sPiyz, sPizy, sPizz

    real(kind=dp) :: dist = 6.0654 ! distance between layers [bohr]

    sPixx = dcmplx(0.0_dp,0.0_dp)
    sPiyy = dcmplx(0.0_dp,0.0_dp)
    sPizz = dcmplx(0.0_dp,0.0_dp)
    sPiyz = dcmplx(0.0_dp,0.0_dp)
    sPizy = dcmplx(0.0_dp,0.0_dp)

    do iG = 1,Nlf
      do jG = 1,Nlf 
        sPixx = sPixx + Pixx(iG,jG) * exp( dcmplx(0.0, (Glf(3,iG) - Glf(3,jG) ) * dist) )
        sPiyy = sPiyy + Piyy(iG,jG) * exp( dcmplx(0.0, (Glf(3,iG) - Glf(3,jG) ) * dist) )
        sPizz = sPizz + Pizz(iG,jG) * exp( dcmplx(0.0, (Glf(3,iG) - Glf(3,jG) ) * dist) )
        sPiyz = sPiyz + Piyz(iG,jG) * exp( dcmplx(0.0, (Glf(3,iG) - Glf(3,jG) ) * dist) )
        sPizy = sPizy + Pizy(iG,jG) * exp( dcmplx(0.0, (Glf(3,iG) - Glf(3,jG) ) * dist) )
      enddo
    enddo 

    sPixx = 1/c0 * sPixx
    sPiyy = 1/c0 * sPiyy
    sPizz = 1/c0 * sPizz
    sPiyz = 1/c0 * sPiyz
    sPizy = 1/c0 * sPizy

  end subroutine sumScreenedPi

  subroutine sumIntScreenedPi(c0, Glf, Pixx, Piyy, Piyz, Pizy, Pizz, sPixx, sPiyy, sPiyz, sPizy, sPizz)
    ! non general function, works only for a symmetric trilayer
    ! sums Pi over all local field vectors to observe sigma^sc in upper layer
    ! integrates each component over dz,dz' from 0 to L/2
    ! Pi = 1/L \sum_Gz \sum_Gz' \int_0^{L/2}\int_0^{L/2} dz dz' Pi(Gz,Gz') * exp(-i Gz * z) * exp(i Gz' * z') 
    implicit none
    real(kind=dp),    intent(in)  :: c0
    real(kind=dp),    intent(in)  :: Glf(:,:)
    complex(kind=dp), intent(in), dimension(:,:) :: Pixx, Piyy, Piyz, Pizy, Pizz
    complex(kind=dp), intent(out) :: sPixx, sPiyy, sPiyz, sPizy, sPizz

    ! real(kind=dp) :: dist = 6.0654 ! distance between layers [bohr]
    ! instead of +/- L/2 we could integrate from +/- dist/2

    sPixx = dcmplx(0.0_dp,0.0_dp)
    sPiyy = dcmplx(0.0_dp,0.0_dp)
    sPizz = dcmplx(0.0_dp,0.0_dp)
    sPiyz = dcmplx(0.0_dp,0.0_dp)
    sPizy = dcmplx(0.0_dp,0.0_dp)

    do iG = 1,Nlf
      do jG = 1,Nlf
        if (iG==1 .and. jG/=1) then
          sPixx = sPixx + Pixx(iG,jG) * ( -dcmplx(0.0_dp,1.0_dp) * ( -1 + exp(dcmplx(0.0,Glf(3,jG)*c0/2)) ) ) * c0/(2*Glf(3,jG))
          sPiyy = sPiyy + Piyy(iG,jG) * ( -dcmplx(0.0_dp,1.0_dp) * ( -1 + exp(dcmplx(0.0,Glf(3,jG)*c0/2)) ) ) * c0/(2*Glf(3,jG))
          sPizz = sPizz + Pizz(iG,jG) * ( -dcmplx(0.0_dp,1.0_dp) * ( -1 + exp(dcmplx(0.0,Glf(3,jG)*c0/2)) ) ) * c0/(2*Glf(3,jG))
          sPiyz = sPiyz + Piyz(iG,jG) * ( -dcmplx(0.0_dp,1.0_dp) * ( -1 + exp(dcmplx(0.0,Glf(3,jG)*c0/2)) ) ) * c0/(2*Glf(3,jG))
          sPizy = sPizy + Pizy(iG,jG) * ( -dcmplx(0.0_dp,1.0_dp) * ( -1 + exp(dcmplx(0.0,Glf(3,jG)*c0/2)) ) ) * c0/(2*Glf(3,jG))
        else if (iG/=1 .and. jG==1) then
          sPixx = sPixx + Pixx(iG,jG) * ( -dcmplx(0.0_dp,1.0_dp) * ( 1 - exp(-dcmplx(0.0,Glf(3,iG)*c0/2)) ) ) * c0/(2*Glf(3,iG))
          sPiyy = sPiyy + Piyy(iG,jG) * ( -dcmplx(0.0_dp,1.0_dp) * ( 1 - exp(-dcmplx(0.0,Glf(3,iG)*c0/2)) ) ) * c0/(2*Glf(3,iG))
          sPizz = sPizz + Pizz(iG,jG) * ( -dcmplx(0.0_dp,1.0_dp) * ( 1 - exp(-dcmplx(0.0,Glf(3,iG)*c0/2)) ) ) * c0/(2*Glf(3,iG))
          sPiyz = sPiyz + Piyz(iG,jG) * ( -dcmplx(0.0_dp,1.0_dp) * ( 1 - exp(-dcmplx(0.0,Glf(3,iG)*c0/2)) ) ) * c0/(2*Glf(3,iG))
          sPizy = sPizy + Pizy(iG,jG) * ( -dcmplx(0.0_dp,1.0_dp) * ( 1 - exp(-dcmplx(0.0,Glf(3,iG)*c0/2)) ) ) * c0/(2*Glf(3,iG))
        else if (iG==1 .and. jG==1) then
          sPixx = sPixx + Pixx(iG,jG) * c0**2/4
          sPiyy = sPiyy + Piyy(iG,jG) * c0**2/4
          sPizz = sPizz + Pizz(iG,jG) * c0**2/4
          sPiyz = sPiyz + Piyz(iG,jG) * c0**2/4
          sPizy = sPizy + Pizy(iG,jG) * c0**2/4
        else
          sPixx = sPixx - Pixx(iG,jG) * ( -1 + exp(dcmplx(0.0,Glf(3,iG)*c0/2)) ) * ( -1 + exp(dcmplx(0.0,Glf(3,jG)*c0/2)) ) &
                & / ( Glf(3,iG) * Glf(3,jG) * exp(dcmplx(0.0,Glf(3,iG)*c0/2)) )
          sPiyy = sPiyy - Piyy(iG,jG) * ( -1 + exp(dcmplx(0.0,Glf(3,iG)*c0/2)) ) * ( -1 + exp(dcmplx(0.0,Glf(3,jG)*c0/2)) ) &
                & / ( Glf(3,iG) * Glf(3,jG) * exp(dcmplx(0.0,Glf(3,iG)*c0/2)) )
          sPizz = sPizz - Pizz(iG,jG) * ( -1 + exp(dcmplx(0.0,Glf(3,iG)*c0/2)) ) * ( -1 + exp(dcmplx(0.0,Glf(3,jG)*c0/2)) ) & 
                & / ( Glf(3,iG) * Glf(3,jG) * exp(dcmplx(0.0,Glf(3,iG)*c0/2)) )
          sPiyz = sPiyz - Piyz(iG,jG) * ( -1 + exp(dcmplx(0.0,Glf(3,iG)*c0/2)) ) * ( -1 + exp(dcmplx(0.0,Glf(3,jG)*c0/2)) ) & 
                & / ( Glf(3,iG) * Glf(3,jG) * exp(dcmplx(0.0,Glf(3,iG)*c0/2)) )
          sPizy = sPizy - Pizy(iG,jG) * ( -1 + exp(dcmplx(0.0,Glf(3,iG)*c0/2)) ) * ( -1 + exp(dcmplx(0.0,Glf(3,jG)*c0/2)) ) & 
                & / ( Glf(3,iG) * Glf(3,jG) * exp(dcmplx(0.0,Glf(3,iG)*c0/2)) )
        end if

      enddo
    enddo 

    sPixx = 1/c0 * sPixx
    sPiyy = 1/c0 * sPiyy
    sPizz = 1/c0 * sPizz
    sPiyz = 1/c0 * sPiyz
    sPizy = 1/c0 * sPizy

  end subroutine sumIntScreenedPi

  subroutine genTS(Nlf, Dxx0, Pixx, TS)
    implicit none

    integer,          intent(in)  :: Nlf
    complex(kind=dp), intent(in)  :: Dxx0(:,:)
    complex(kind=dp), intent(in)  :: Pixx(:,:)
    complex(kind=dp), intent(out) :: TS(:,:)

    integer :: iG, jG !, kG
    complex(kind=dp), allocatable :: Imat(:,:)

    allocate(Imat(Nlf,Nlf))

    Imat = dcmplx(0.0_dp,0.0_dp)
    TS = dcmplx(0.0_dp,0.0_dp)
    do iG = 1,Nlf
      Imat(iG,iG) = dcmplx(1.0,0.0)
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

    Imat(1:Nlf2,1:Nlf2) = dcmplx(0.0_dp,0.0_dp)
    TP(1:Nlf2,1:Nlf2) = dcmplx(0.0_dp,0.0_dp) 
    do iG = 1,Nlf2
      Imat(iG,iG) = dcmplx(1.0,0.0)
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

  subroutine  genSpectra(o, oi, beta, iq, Nq, itheta, theta, Ntheta, c0, Nlf, parG, Glf, Dxx, Dyy, Dzz, Dyz, Dzy)
    ! Calculation of reflected, transmited and absorbed coefficients
    ! DEBUG: should rename D to Pi because this is the current-current respone tensor

    integer,          intent(in) :: iq, itheta, Ntheta, Nq
    integer,          intent(in) :: Nlf
    real(kind=dp),    intent(in) :: o, theta
    real(kind=dp),    intent(in) :: c0
    complex(kind=dp), intent(in) :: oi, beta
    integer,          intent(in) :: parG(:)
    real(kind=dp),    intent(in) :: Glf(:,:)
    complex(kind=dp), intent(in), dimension(:,:) :: Dxx, Dyy, Dzz, Dyz, Dzy

    real(kind=dp),  parameter :: gamma   = 1.0_dp/137.0_dp
    real(kind=dp),  parameter :: pi      = 4.0_dp*atan(1.0_dp)
    real(kind=dp),  parameter :: Hartree = 2.0_dp*13.6056923_dp

    integer          :: iG, jG
    real(kind=dp)    :: Ff1, Ff2, Gf1, Gf2
    real(kind=dp)    :: A_s, Tr_s, R_s, A_p, Tr_p, R_p 
    complex(kind=dp) :: Dxx_r, Dyy_r, Dzz_r, Dyz_r
    complex(kind=dp) :: Dxx_t, Dyy_t, Dzz_t, Dyz_t
    complex(kind=dp) :: Dxx_a, Dyy_a, Dzz_a, Dyz_a
    complex(kind=dp) :: Tran_s, Tran_p,  Ref_s, Ref_p, Abso_s, Abso_p
    complex(kind=dp) :: Dxx0, Dyy0
  
    Dxx0 = 2*pi * gamma**2 * dcmplx(0.0_dp,1.0_dp) / beta 
    Dyy0 = 2*pi * gamma * dcmplx(0.0_dp,1.0_dp) / oi 
    Dxx_r = dcmplx(0.0_dp,0.0_dp)
    Dyy_r = dcmplx(0.0_dp,0.0_dp)
    Dzz_r = dcmplx(0.0_dp,0.0_dp)
    Dyz_r = dcmplx(0.0_dp,0.0_dp)
    Dxx_t = dcmplx(0.0_dp,0.0_dp)
    Dyy_t = dcmplx(0.0_dp,0.0_dp)
    Dzz_t = dcmplx(0.0_dp,0.0_dp)
    Dyz_t = dcmplx(0.0_dp,0.0_dp)
    Dxx_a = dcmplx(0.0_dp,0.0_dp)
    Dyy_a = dcmplx(0.0_dp,0.0_dp)
    Dzz_a = dcmplx(0.0_dp,0.0_dp)
    Dyz_a = dcmplx(0.0_dp,0.0_dp)

    do iG = 1,Nlf 
      Ff1 = (2.0_dp*parG(iG)/sqrt(c0))*sin(beta*c0/2.0_dp)/(beta+Glf(3,iG)) 
      Gf1 = (2.0_dp*parG(iG)/sqrt(c0))*sin(beta*c0/2.0_dp)/(beta-Glf(3,iG)) 
      do jG = 1,Nlf
        Ff2 = (2.0_dp*parG(jG)/sqrt(c0))*sin(beta*c0/2.0_dp)/(beta+Glf(3,jG))
        Gf2 = (2.0_dp*parG(jG)/sqrt(c0))*sin(beta*c0/2.0_dp)/(beta-Glf(3,jG))
 
        ! reflection
        Dxx_r = Dxx_r + Ff1 * Dxx(iG,jG) * Gf2
        Dyy_r = Dyy_r + Ff1 * Dyy(iG,jG) * Gf2 
        Dzz_r = Dzz_r + Ff1 * Dzz(iG,jG) * Gf2 
        Dyz_r = Dyz_r + Ff1 * Dyz(iG,jG) * Gf2
        ! transmission
        Dxx_t = Dxx_t + Gf1 * Dxx(iG,jG) * Ff2
        Dyy_t = Dyy_t + Gf1 * Dyy(iG,jG) * Ff2
        Dzz_t = Dzz_t + Gf1 * Dzz(iG,jG) * Ff2
        Dyz_t = Dyz_t + Gf1 * Dyz(iG,jG) * Ff2
        ! absorption
        Dxx_a = Dxx_a + Ff1 * Dxx(iG,jG) * Ff2
        Dyy_a = Dyy_a + Ff1 * Dyy(iG,jG) * Ff2
        Dzz_a = Dzz_a + Ff1 * Dzz(iG,jG) * Ff2
        Dyz_a = Dyz_a + Ff1 * Dyz(iG,jG) * Ff2
      enddo
    enddo
 
    ! s-mode 
    Ref_s  = Dxx0 * Dxx_r
    Tran_s = Dxx0 * Dxx_t
    Abso_s = Dxx_a
 
    ! p-mode
    Ref_p = Dyy0 * ( Dyy_r * cos(theta) &
          & - Dzz_r * sin(theta)**2 / cos(theta) )

    Tran_p = Dyy0 * ( Dyy_t * cos(theta) - 2.0_dp * Dyz_t * sin(theta) &
           & + Dzz_t * sin(theta)**2 / cos(theta) )
 
    Abso_p = cos(theta) * Dyy_a * cos(theta) &
           & + sin(theta) * Dzz_a * sin(theta) &
           & - 2.0_dp*cos(theta) * Dyz_a * sin(theta)
 
 
    ! Absorbption
    A_s = (4.0_dp*pi*gamma/o) * aimag(Abso_s)
    A_p = (4.0_dp*pi*gamma/o) * aimag(Abso_p)
    ! Reflection 
    R_s = cos(theta) * real( Ref_s*conjg(Ref_s) ) 
    R_p = real( Ref_p * conjg(Ref_p) )  
    !  Tran_smission
    Tr_s = 1.0_dp - cos(theta) * ( real(Tran_s*conjg(Tran_s) ) - 2.0_dp * real(Tran_s) ) 
    Tr_p = 1.0_dp - real( Tran_p*conjg(Tran_p) ) - cos(theta) * 2.0_dp * real(Tran_p) 

    ! write A,T,R spectra to file for each angle itheta
    call writeSpectra(iq, itheta, Ntheta, Nq, o, A_p, R_p)

  end subroutine genSpectra

  subroutine writeSpectra(iq, itheta, Ntheta, Nq, o, A_p, R_p)
    implicit none

    integer,          intent(in) :: iq, itheta, Nq, Ntheta
    real(kind=dp),    intent(in) :: o, A_p, R_p

    real(kind=dp),  parameter :: Hartree = 2.0_dp*13.6056923_dp

    integer :: iuni1, iuni2, iuni3
    logical :: exist_a, exist_t, exist_r
    character(len=20) :: id

    if (Ntheta/=0 .and. Nq==0) then
      id = "_theta=#"//int2str(itheta)
    else if (Ntheta==0 .and. Nq/=0) then
      id = "_q#"//int2str(iq)
    else
      id = ""
    endif

    inquire(file="absorption"//id, exist=exist_a)
    if (exist_a) then
      open(newunit=iuni1, file="absorption"//id, status="old", position="append", action="write")
      write(iuni1,*) o*Hartree, A_p
      close(iuni1)
    else
      open(newunit=iuni1, file="absorption"//id, status="new", action="write")
      write(iuni1,*) o*Hartree, A_p
      close(iuni1)
    end if

    inquire(file="transmission"//id, exist=exist_t)
    if (exist_t) then
      open(newunit=iuni2, file="transmission"//id, status="old", position="append", action="write")
      write(iuni2,*) o*Hartree, 1 - A_p - R_p
      close(iuni2)
    else
      open(newunit=iuni2, file="transmission"//id, status="new", action="write")
      write(iuni2,*) o*Hartree, 1 - A_p - R_p
      close(iuni2)
    end if    

    inquire(file="reflection"//id, exist=exist_r)
    if (exist_r) then
      open(newunit=iuni3, file="reflection"//id, status="old", position="append", action="write")
      write(iuni3,*) o*Hartree, R_p
      close(iuni3)
    else
      open(newunit=iuni3, file="reflection"//id, status="new", action="write")
      write(iuni3,*) o*Hartree, R_p
      close(iuni3)
    end if   

  end subroutine writeSpectra

  subroutine writePi(o, Pi, filename)
    implicit none
    real(kind=dp),     intent(in) :: o
    complex(kind=dp),  intent(in) :: Pi
    character(len=*), intent(in) :: filename


    real(kind=dp),  parameter :: Hartree = 2.0_dp*13.6056923_dp

    integer :: iuni
    logical :: exist

    inquire(file=trim(adjustl(filename)), exist=exist)
    if (exist) then
      open(newunit=iuni, file=trim(adjustl(filename)), status="old", position="append", action="write")
      write(iuni,*) o*Hartree, real(Pi), aimag(Pi)
      close(iuni)
    else
      open(newunit=iuni, file=trim(adjustl(filename)), status="new", action="write")
      write(iuni,*) o*Hartree, real(Pi), aimag(Pi)
      close(iuni)
    end if
      
  end subroutine writePi

  subroutine writeSigma(o, c0, Pi, filename)
    implicit none
    real(kind=dp),     intent(in) :: o, c0
    complex(kind=dp),  intent(in) :: Pi
    character(len=*), intent(in)  :: filename


    real(kind=dp),  parameter :: Hartree = 2.0_dp*13.6056923_dp

    integer :: iuni
    logical :: exist

    inquire(file=trim(adjustl(filename)), exist=exist)
    if (exist) then
      open(newunit=iuni, file=trim(adjustl(filename)), status="old", position="append", action="write")
      write(iuni,*) o*Hartree, real(-dcmplx(0.0_dp,1.0_dp)*4.0_dp*c0*Pi/o), &
                             & aimag(-dcmplx(0.0_dp,1.0_dp)*4.0_dp*c0*Pi/o)
      close(iuni)
    else
      open(newunit=iuni, file=trim(adjustl(filename)), status="new", action="write")
      write(iuni,*) o*Hartree, real(-dcmplx(0.0_dp,1.0_dp)*4.0_dp*c0*Pi/o), &
                             & aimag(-dcmplx(0.0_dp,1.0_dp)*4.0_dp*c0*Pi/o)
      close(iuni)
    end if
      
  end subroutine writeSigma

  subroutine writeSigma_macroscopic(o, c0, Nlf, Pixx, Piyy, Pizz, TS, TP)
    ! Calculates macroscopic conducitivities from current-current response functions Pi_{\mu\nu}
    ! Write the conducitivities and response functions to files.

    implicit none
    integer,          intent(in) :: Nlf
    real(kind=dp),    intent(in) :: c0
    real(kind=dp),    intent(in) :: o
    complex(kind=dp), intent(in), dimension(:,:) ::  Pixx, Piyy, Pizz, TS, TP

    
    real(kind=dp),  parameter :: Hartree = 2.0_dp*13.6056923_dp


    integer :: Nlf2
    integer :: iuni1, iuni2, iuni3, iuni4, iuni5, iuni6
    complex(kind = dp) :: Sigma_xx,Sigma_yy,Sigma_zz

    Nlf2 = 2*Nlf

    Sigma_xx = sum( Pixx(1,1:Nlf) * TS(1:Nlf,1) )                                        ! s-mod
    Sigma_yy = sum( Piyy(1,1:Nlf) * TP(1:Nlf,1)     + Piyz(1,1:Nlf)*TP(Nlf:Nlf2,1) )     ! p-mod
    Sigma_zz = sum( Pizy(1,1:Nlf) * TP(1:Nlf,Nlf+1) + Pizz(1,1:Nlf)*TP(Nlf:Nlf2,Nlf+1) ) ! p-mod

    ! write sigma_{\mu\nu}^\text{macro} [ pi*e^2/2h ]
    open(newunit=iuni1, file='sigma_macro_xx')
    write(iuni1,*) o*Hartree, real(-dcmplx(0.0_dp,1.0_dp)*4.0_dp*c0*Sigma_xx/o)
    close(iuni1)

    open(newunit=iuni2, file='sigma_macro_yy')
    write(iuni2,*) o*Hartree, real(-dcmplx(0.0_dp,1.0_dp)*4.0_dp*c0*Sigma_yy/o)
    close(iuni2)

    open(newunit=iuni3, file='sigma_macro_zz')
    write(iuni3,*) o*Hartree, real(-dcmplx(0.0_dp,1.0_dp)*4.0_dp*c0*Sigma_zz/o)
    close(iuni3)

    ! write Pi_{\mu\nu}(1,1)
    open(newunit=iuni4, file='Pi_11_xx')
    write(iuni4,*) o*Hartree, real(-dcmplx(0.0_dp,1.0_dp)*4.0_dp*c0*Pixx(1,1)/o)
    close(iuni4)

    open(newunit=iuni5, file='Pi_11_yy')
    write(iuni5,*) o*Hartree, real(-dcmplx(0.0_dp,1.0_dp)*4.0_dp*c0*Piyy(1,1)/o)
    close(iuni5)

    open(newunit=iuni6, file='Pi_11_zz')
    write(iuni6,*) o*Hartree, real(-dcmplx(0.0_dp,1.0_dp)*4.0_dp*c0*Pizz(1,1)/o)
    close(iuni6)

      
  end subroutine writeSigma_macroscopic

end program photon


