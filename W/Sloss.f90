program surface_loss
 
!---------------------------------------------------------------------------
!        program FOR ab initio-SURFACE LOSS CALCULATION FOR LAYERED SYSTEMS
!        USING SUPERCELL METHOD
!---------------------------------------------------------------------------
!        Quantum Esspresso:
!        verbosity           = 'high'
!       WARNING: VALID ONLY FOR NORM-CONSERVING PSEUDOPOTENTIALS
!---------------------------------------------------------------------------


!USE ISO_Fortran_env ! precision kind, fortran 2008
use OMP_lib
use iotk_module
use ModPointR

implicit none

character (len=iotk_attlenx) :: attr
logical :: found

! single / double precision
! integer, parameter :: sp = real32
! integer, parameter :: dp = real64
! integer, parameter :: qp = real128

! integer, parameter :: sp = SELECTED_REAL_KIND(p=6)
! integer, parameter :: dp = SELECTED_REAL_KIND(p=15)
! integer, parameter :: qp = SELECTED_REAL_KIND(p=33)

! NOT used, read directly from QE
! integer :: nMPx,nMPy,nMPz
! parameter(nMPx=201,nMPy=201,nMPz=1)


character (len=100) :: rundir, savedir, band_file, scf_file
namelist /directories/ rundir, savedir, scf_file, band_file

integer :: debugCount = 0



integer :: ik,i,j,jk,it,lk,Ntot,iG0,Nsymm,iq, &
           io,n,m,iG,R1,K1,R2,K2, Nlf,NG1,   &
           NG2,iG1,iG2,jG,kG,jo,jump,loss,   &
           iGfast,ikmin,kG1,kG2,nord, lf


integer :: Nk     ! = 48*NkI, number of wave vectors in FBZ with no symmetry 
integer :: NkI    ! number of wave vectors in IBZ
integer :: Nband  ! number of bands
integer :: Nocc   ! Number of occupied bands (unit cell)
integer :: NelQE  ! Number of electrons(unit cell)
integer :: NGd    ! number of coefficients CG shulod be less than minimum number of coefficients all over all evc.n files ... moglo bi se dinamicki alocirati 
integer :: NG     ! total number of G vectors  
integer :: no     ! number of frequencies
integer :: nq     ! broj valnih vektora tu je 2 jer je rucno paralelizirano!
integer :: Nlfd   ! dimenzija polja za local field zasto prozivoljno 50, ne moze se znati unaprijed
namelist /config/ NkI, Nband, Nocc, NelQE, NGd, NG, no, nq, Nlfd


! file i/o debug
integer :: ist,ist2,ist4,ist5,ist6,ist7,ist8,ist9,ist10,ist11,ist12
integer :: lno,lno2,lno9,lno10,lno11,lno12

! constants
real(kind=dp),    parameter :: pi = 4.D0*ATAN(1.D0)
real(kind=dp),    parameter :: eV = 1.602176487D-19
real(kind=dp),    parameter :: Hartree = 2.0D0*13.6056923D0
real(kind=dp),    parameter :: Planck = 6.626196D-34
real(kind=dp),    parameter :: three = 3.0d0 
real(kind=dp),    parameter :: aBohr = 0.5291772d0
complex(kind=dp), parameter :: rone  = cmplx(1.0,0.0)
complex(kind=dp), parameter :: czero = cmplx(0.0,0.0)
complex(kind=dp), parameter :: ione  = cmplx(0.0,1.0)

! scalars
real(kind=dp) :: kx,ky,kz
real(kind=dp) :: KQx,KQy,KQz
real(kind=dp) :: qgx,qgy,qgz
real(kind=dp) :: omin,omax
real(kind=dp) :: qx,qy,qz
real(kind=dp) :: kmin
real(kind=dp) :: domega
real(kind=dp) :: o
real(kind=dp) :: K11,K22,K33
real(kind=dp) :: Lor ! lorentzian
real(kind=dp) :: De 
real(kind=dp) :: Gabs ! 
real(kind=dp) :: kref ! trazi najmanju k-tocku sampliranu u MP meshu u kojem se moze izracunati ILS
real(kind=dp) :: Eref 
real(kind=dp) :: Gxx1,Gyy1,Gzz1
real(kind=dp) :: Gxx2,Gyy2,Gzz2
real(kind=dp) :: fact
real(kind=dp) :: oi,oj
real(kind=dp) :: ImChi0, ReChi0
! real(kind=dp) :: q ! ne koristi se
! real(kind=dp) :: qmax ! ne koristi se
real(kind=dp) :: Nel ! Number of electrons(1BZ integration)
real(kind=dp) :: absq
real(kind=dp) :: error
real(kind=dp) :: W1,W2
real(kind=dp) :: ImW
real(kind=dp) :: Wind
real(kind=dp) :: W2KK
real(kind=dp) :: KKS = 0.0
real(kind=dp) :: SKK = 0.0
real(kind=dp) :: WindKK
real(kind=dp) :: krefM
complex(kind=dp) :: G0


! parameters
integer :: qmin,qmax
real(kind=dp) :: Gcar   ! unit cell norm.
real(kind=dp) :: Efermi ! [eV] Fermi en. 
real(kind=dp) :: a0     ! [a.u.]  unit cell parameter in parallel direction 
real(kind=dp) :: c0     ! [a.u.]  unit cell parameter in perependicular direction 
real(kind=dp) :: eps    ! 1.0D-4 threshold
real(kind=dp) :: T      ! [eV] temperature 
real(kind=dp) :: eta    ! damping i\eta
real(kind=dp) :: Ecut   ! [Hartree] cutoff energy for crystal local field calculations , for Ecut=0 S matrix is a scalar ?
real(kind=dp) :: Vcell  ! [a.u.^3] unit-cell volume 
namelist /parameters/ Efermi, a0, c0, eps, T, eta, Ecut, Vcell
namelist /system/ lf, loss, jump, omin, omax, qmin, qmax

! scalar arrays
integer,       dimension(3) :: Gi                           ! pomocna funkcija
integer,       dimension(:),      allocatable :: Gfast      ! pomocna funkcija
integer,       dimension(:),      allocatable :: parG       ! paritet svakog valnog vektora
real(kind=dp), dimension(:),      allocatable :: factMatrix

! multidim arrays
real(kind=dp), dimension(48,3,3)  :: R                      ! matr. simetrijskih op.
real(kind=dp), dimension(48,3,3)  :: RI                     ! inverz od R
real(kind=dp), dimension(3,3)     :: KC                     ! pomocna funkcija
real(kind=dp), dimension(:,:),    allocatable    :: kI
real(kind=dp), dimension(:,:),    allocatable    :: E       ! vl. vr. danog k-i i band-i
! real(kind=dp), dimension(:,:),    allocatable    :: k
real(kind=dp), dimension(:,:),    allocatable    :: ktot    ! polje jedinstvenih k-tocaka u FBZ
real(kind=dp), dimension(:,:),    allocatable    :: G       ! polje valnih vektora G u recp. prost. za wfn.
real(kind=dp), dimension(:,:),    allocatable    :: V       ! matr. gole coulomb. int.
real(kind=dp), dimension(:,:,:),  allocatable    :: S0      ! korelacijska matrica
real(kind=dp), dimension(:,:,:),  allocatable    :: S0_partial ! pomocna var. za S0 redukciju
real(kind=dp), dimension(:,:),    allocatable    :: Glf     ! local field effect polje valnih vekt. u rec. prost.
real(kind=dp), dimension(:,:),    allocatable    :: GlfV    ! local field effect za nezasjenjenu (golu) interakciju V


! Nlfd dimenziju ne znamo a priori zapravo, trebalo bi staviti sve te matrice allocatable i 
! naknadno ih alocirati

complex(kind=dp), dimension(:,:),   allocatable  :: Imat ! jedinicna matr.
complex(kind=dp), dimension(:,:),   allocatable  :: diel_epsilon ! Epsilon (GG')  = I - V(GG')Chi0
complex(kind=dp), dimension(:,:),   allocatable  :: Chi ! (eq. 2.88 nakon invertiranja) ;oprez bio je double precision

complex(kind=dp), dimension(:)    , allocatable  :: MnmK1K2 ! nabojni vrhovi
complex(kind=dp), dimension(:,:)  , allocatable  :: Chi0 ! (eq. 2.89)
complex(kind=dp), dimension(:,:,:), allocatable  :: WT ! time ordered RPA screened coulomb int. (eq. 2.93)
complex(kind=dp), dimension(:,:)  , allocatable  :: Gammap ! omega>0 ,eq....(skripta 5) \sum_{q,m} \int \dd omega' S(\omega')/{(\omega-\omega'-e_{k+q,m} +i\eta}) za GW se koristi se za ovaj dio 
complex(kind=dp), dimension(:,:)  , allocatable  :: Gammam ! omega<0, allocatable  

character (len=100) :: bandn,bandm,dummy,pathk1,pathk2,dato, path
character (len=35)  :: tag,buffer


! complex(kind=dp), dimension(:), pointer,save :: C1,C2 ! Fourierovi koef. u razvoju wfn. iz QE 
complex(kind=dp), dimension(:), allocatable :: C1,C2

! OpenMP vars
integer      :: Nthreads
integer,save :: thread_id
!$omp threadprivate(thread_id)
namelist  /parallel/ Nthreads


! MKL matrix inversion vars
integer :: info_trf, info_tri
integer,allocatable :: ipiv(:)
integer :: lwork
integer,allocatable :: work(:)

! load namelist
open(10,file='config.in')
read(10,nml=directories,iostat=ist4)
read(10,nml=system,iostat=ist5)
read(10,nml=config,iostat=ist6)
read(10,nml=parameters,iostat=ist7)
read(10,nml=parallel,iostat=ist8)
close(10)

! constants
Nk     = 48*NkI                   ! number of wave vectors in FBZ with no symmetry 
T      = T/Hartree                ! convert temperature from eV to Hartree
Efermi = Efermi/Hartree           ! convert Fermi en. from eV to Hartree
eta    = eta/Hartree
Gcar   = 2.0*pi/a0                ! unit cell norm.

! scalar arrays
allocate(parG(NG))                ! paritet svakog valnog vektora
allocate(Gfast(Nlfd*NGd))
allocate(factMatrix(no))
allocate(MnmK1K2(Nlfd))            ! nabojni vrhovi

! multidim arrays
allocate(kI(3,NkI))
allocate(E(NkI,Nband))       ! vl. vr. danog k-i i band-i
! allocate(k(3,Nk))
allocate(ktot(3,Nk))               ! ukupno jedinstvenih k-tocaka u FBZ
allocate(G(3,NG))                  ! polje valnih vektora G u recp. prost. za wfn.
allocate(V(Nlfd,Nlfd))             ! matr. gole coulomb. int.
allocate(S0(no,Nlfd,Nlfd))         ! korelacijska matrica
! allocate(S0_partial(no,Nlfd,Nlfd))         ! korelacijska matrica
allocate(Glf(3,Nlfd))              ! local field effect polje valnih vekt. u rec. prost.
allocate(GlfV(3,Nlfd))             ! local field effect za nezasjenjenu (golu) interakciju V

allocate(Imat(Nlfd,Nlfd))          ! jedinicna matr.
allocate(diel_epsilon(Nlfd,Nlfd) ) ! Epsilon (GG')  = I - V(GG')Chi0
allocate(Chi(Nlfd,Nlfd))           ! (eq. 2.88 nakon invertiranja) ;oprez bio je double precision

allocate(Chi0(Nlfd,Nlfd))          ! (eq. 2.89)
allocate(WT(no,Nlfd,Nlfd))         ! time ordered RPA screened coulomb int. (eq. 2.93)
allocate(Gammap(Nlfd,Nlfd))        ! omega>0 ,eq....(skripta 5) \sum_{q,m} \int \dd omega' S(\omega')/{(\omega-\omega'-e_{k+q,m} +i\eta}) za GW se koristi se za ovaj dio 
allocate(Gammam(Nlfd,Nlfd))        ! omega<0



domega = (omax-omin)/(no-1)



! call FOR POINT GROUP TRANSFORMATIONS
! Point group transformations are in Cartesian coordinate
path = trim(rundir)//trim(scf_file)
call PointR(path,Nsymm,R,RI)
print *,"status: PointR done."


! Upis valnih vektora iz irreducibilne Brillouinove
! zone i pripadnih energijskih nivoa iz filea '****.band'.
! wave vectors are in Cartesian coordinate
path=trim(rundir)//trim(band_file)
call loadkIandE(path, NkI, Nband, Nocc, kI, E)
print *,"status: kI and E loaded."


! generator 1.B.Z.
! Dio programa koji pomocu operacija tockaste grupe i vektora iz
! I.B.Z. generira sve (MEDJUSOBNO RAZLICITE!!!) valne vektore u 1.B.Z.
! Ntot-Tot number of different points ''ktot'' inside 1.B.Z

call genFBZ(Nk,NkI,Nsymm,eps,kI,R,Ntot,ktot)
print *,"status: FBZ generated."

! ! Checking 1BZ integration
! ! provjeri je li broj el. u FBZ (Nel) odgovara stvarnom broju el. u jed. cel. (NelQE)
call checkFBZintegration(Nband,NkI,Nsymm,Ntot,eps,kI,RI,Efermi,E,NelQE,Nel)
print *,"status: FBZ integration correct."


! KC transformation matrix from rec.cryst. axes to cart.koord.
! If G' is vector in rec.cryst. axes then a=KC*a' is vector in cart. axes
path = TRIM(rundir)//TRIM(scf_file)
call loadKC(path,KC)
print *,"status: KC transformation matrix (rec.cryst.->cart.) loaded."


! Reading the reciprocal vectors in crystal coordinates and transformation
! in Cartesian cordinates.
call loadG(KC,parG,G)
print *,"status: G vectors loaded."


! Generate Reciprocal vectors for crystal local field effects calculations in array Glf(3,Nlf)
call genGlfandParity(lf,Ecut,NG,Gcar,G,Nlf,Nlfd,parG,Glf)

! MKL matrix inversion vars
allocate(ipiv(MAX(1,MIN(Nlf, Nlf))))
lwork = Nlf
allocate(work(Nlf))



! IBZ q LOOP STARTS HERE!!!
! iq=0 ne moze biti nula, opticki racun
! iq=2 do iq=...cutoff transfer q vektor!
! ikmin = min. valni vektor u BZ svi veci su visekratnici tog minimalnog
q_loop: do  iq = qmin,qmax ! 42,61
  
  ! searching min. q=(qx,qy,qz) in Gamma - M direction
  call findMinQ(Ntot, ktot, qx, qy, qz)
  
  ! Info file
  call writeInfo(qx, qy, qz, Gcar, Nsymm, Nlf, Ntot, NkI, Nband, eta, T, Nel, NelQE )

  ! intialize correlation matrix
  S0(1:no,1:Nlf,1:Nlf) = cmplx(0.0,0.0)
  
  
! 1.B.Z  LOOP STARTS HERE !!!!

print *, 'DEBUG: entering parallel region'
!$omp parallel shared(S0,iq,kI,ktot,RI,eps,E,G,NkI,Nsymm,NG,Ntot,Nocc,Nband,NGd,Nlf,Nlfd,eta,Vcell) private(ik, S0_partial,MnmK1K2,K11,K22,K33,kx,ky,kz,i,j,it,R1,R2,iG0,KQx,KQy,KQz,iG,jG,jk,K1,K2,n,m,pathk1,pathk2,bandn,bandm,NG1,NG2,io,o,De,Lor,Gxx1,Gxx2,Gyy1,Gyy2,Gzz1,Gzz2,Gfast,iGfast, iG1, iG2,attr,C1,C2) firstprivate(savedir,jump,domega) num_threads(Nthreads) !  default(private) 
thread_id =  omp_get_thread_num()

!$omp do 
do ik = 1, Ntot   ! k_loop_FBZ_2nd: 
  ! neven debug
  ! print *, 'thread id:',thread_id

  ! debug vito (prepraviti za paralelnu izvedbu)
  ! open(122,FILE='status')
  ! write(122,*) 'iq=',iq
  ! write(122,*) 'ik=',ik
  ! close(122)
  
  
  kx = ktot(1,ik)
  ky = ktot(2,ik)
  kz = ktot(3,ik)
  
  ! trazenje (kx,ky,kz) u ireducibilnoj zoni
  call findKinIBZ(ik, NkI, Nsymm, eps, kx, ky, kz, RI, kI, R1, K1)


  KQx = kx + qx
  KQy = ky + qy
  KQz = kz + qz

  !$omp critical(printWaveVector)
  ! thread_id =  omp_get_thread_num()
  print *,'thread id:',thread_id,'ik: ',ik
  print *, 'KQx,KQy,KQz:',KQx,KQy,KQz
  !$omp end critical(printWaveVector)

  ! trazenje (KQx,KQy) prvo u 1.B.Z a onda u I.B.Z.
  call findKQinBZ(KQx, KQy, KQz, eps, Nsymm, NkI, Ntot, NG, ktot, kI, RI, G, iG0, R2, K2)
  
  
  ! R1 - integer, redni broj point operacije R1 u transformaciji K=R1*K1.
  ! K1 - integer, redni broj valnog vektora K1 u transformaciji K=R1*K1.
  ! iG0 i R2-integeri, redni broj vektora reciprocne restke G0 i point operacije R2 u transformaciji K+Q=G0+R2*K2.
  ! K2 - integer, redni broj valnog vektora K2 u transformaciji  K+Q=G0+R2*K2.
  

  allocate(S0_partial(no,Nlf,Nlf))   ! pomocna var. za redukciju S0
  S0_partial(1:no,1:Nlf,1:Nlf) = cmplx(0.0)
  
  bands_n_loop: do  n = 1, Nocc         ! filled bands loop
    bands_m_loop: do  m = Nocc+1, Nband ! empty bands loop
      ! ucitavanje evc.dat binarnih datoteka za fiksni K1,K2, i vrpce n i m
      
      !$omp critical(pathk_read)
      
      ! otvara save/K.000x/evc.dat u atributu <evc band> ispod CnK(G) koef.
      call paths(savedir,K1,K2,n,m,pathk1,pathk2,bandn,bandm) 

      ! print *,'pathk1',pathk1
      call iotk_open_read(10+ik,pathk1)
      call iotk_scan_empty(10+ik,"INFO",attr=attr) ! Otvaranje atribute za INFO
      call iotk_scan_attr(attr,"igwx",NG1)
      
      allocate (C1(NG1))                           ! Alociranje polja C1
      call iotk_scan_dat(10+ik,bandn,C1)           ! Ucitavanje podataka iza evc.n
      call iotk_close_read(10+ik)

      ! print *,'pathk2',pathk2
      call iotk_open_read(10+ik,pathk2)
      call iotk_scan_empty(10+ik,"INFO",attr=attr)
      call iotk_scan_attr(attr,"igwx",NG2)
      allocate (C2(NG2))                           ! Alociranje polja C1
      call iotk_scan_dat(10+ik,bandm,C2)           ! Ucitavanje podataka iza evc.n
      call iotk_close_read(10+ik)

      !$omp end critical(pathk_read)

      if (NGd > NG1) then
        write(*,*) 'NGd is bigger than NG1=',NG1
        STOP
      else if (NGd > NG2) then
        write(*,*) 'NGd is bigger than NG2=',NG2
        STOP
      end if


      ! Konstrukcija stupca matricnih elementa nabojnih vrhova MnmK1K2(G)      
      call genMnmK1K2(jump, eps, Nlf, iG0, NG1, NG2, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2)

      deallocate(C1)
      deallocate(C2)

      do  io = 1,no
        o = (io-1)*domega
        De = o + E(K1,n) - E(K2,m) 
        Lor = -eta/(De**2 + eta**2)     ! ovo bi analticki bila delta funkcija imag. dio od 1/De
        if (abs(Lor) >= 1.0D-3/eta) then
          do  iG = 1,Nlf
            do  jG = 1,Nlf
              S0_partial(io,iG,jG) = S0_partial(io,iG,jG) - 2.0*Lor*MnmK1K2(iG)*conjg(MnmK1K2(jG)) / (pi*Ntot*Vcell)
            end do
          end do
        end if
      end do
              
    end do bands_m_loop ! end of m do loop
  end do bands_n_loop   ! end of n do loop


  !$omp critical(sumS0)      
  ! neven debug
  ! print *, 'K1: ',K1,'K2: ',K2,'R1: ',R1, 'R2: ',R2
  ! print *, 'sum(S0): ', sum(S0(1:no,1:Nlf,1:Nlf)
  S0(1:no,1:Nlf,1:Nlf) = S0(1:no,1:Nlf,1:Nlf) + S0_partial(1:no,1:Nlf,1:Nlf)
  !$omp end critical(sumS0)

  deallocate(S0_partial)
  jump = 1


end do ! k_loop_FBZ_2nd !  end of 1.B.Z do loop
 !$omp end do
 !$omp end parallel
 
 print *, 'DEBUG: exiting parallel region'

 ! STOP

  ! Convert (qx,qy,qz) and Glf from Cartesian coords. to atomic units
  qx = Gcar*qx ! convert qx *2p/a0
  qy = Gcar*qy
  qz = Gcar*qz
  

  GlfV(1:3,1:Nlf) = Gcar*Glf(1:3,1:Nlf)

   ! Kramers-Kroning relacije
  
  omega_loop_A: do  io = 1,no-1
    call genChi0(io,no,Nlf,domega,S0,Chi0)
    
    ! MATRIX V(G,G')
    call genV(eps,qx,qy,Nlf,parG,Glf,GlfV,V)
    
    call genDielectricEpsilon(Nlf,Chi0,V,diel_epsilon)
    
    !  invertiranje matrice ''diel_epsilon = 1-Chi_0*V''
    call gjel(diel_epsilon,Nlf,Nlfd,Imat,Nlf,Nlfd)
    ! call dgetrf( Nlf,Nlfd, diel_epsilon, Nlf, ipiv, info_trf)
    ! call dgetri( Nlf, diel_epsilon, Nlf, ipiv, work, lwork, info_tri )

    ! nezasjenjeni Chi
    call genChi(Nlf,diel_epsilon,Chi0,Chi)

    !  SCREENED COULOMB INTERACTION W^T_GG'(Q,\omega)
    call genWT(io,Nlf,V,Chi,WT)


    ! neven debug WT
    oi = (io-1)*domega ! koristi se u ispisu WT(io,1,1)
    write(20008,*) oi*Hartree,aimag(WT(io,1,1))
    write(10008,*) oi*Hartree,real(WT(io,1,1))

    write(22308,*) oi*Hartree,aimag(WT(io,2,2))
    write(12308,*) oi*Hartree,real(WT(io,2,3))

  end do omega_loop_A
  

  
  ! ispis time ordered zasjenjene kulonske interakcije W_GG'^T(Q,\omega)
  call writeWT_Qi(iq,Nlf,domega,WT)


  S0(1:no-1,1:Nlf,1:Nlf) = -(1.0/pi)*aimag( WT(1:no-1,1:Nlf,1:Nlf) )
  
  ! podatci za GW sve na dalje

  KKS = 0.0
  SKK = 0.0  
  omega_loop_B: do  io = 1,no-1
    ! print*,io
    oi = (io-1)*domega

    do  iG = 1,Nlf
      do  jG = 1,Nlf
        call genW1W2(io,no,domega,S0,W1,W2)
        
        ImW = -pi * S0(io,iG,jG)
        ! stvari vezane uz GW...
        Gammap(iG,jG) = cmplx(W1,ImW)
        Gammam(iG,jG) = cmplx(-W2,0.0)
        if (iG == 1 .and. jG == 1) then
          W2KK = W2
        end if
        if (io == 1) then
          G0 = Gammap(1,1)
        end if
        
      end do
    end do
    
    !  Provjera Kramers-Kroning relacija

    Wind = real(WT(io,1,1)-V(1,1))
    WindKK = real(Gammap(1,1)) - W2KK

    fact = domega
    if (io == 1 .or. io == no-1) then
      fact = 0.5*domega
    end if
    KKS = KKS + fact*(WindKK-Wind)*(WindKK-Wind)
    SKK = SKK + fact*Wind*Wind
    


  end do omega_loop_B

  ! write Kramers-Kroning relations check
  call writeKramKron_Qi(iq,qx, qy, qz, Gcar, KKS, SKK, G0, WT, V)

  
end do q_loop

! MKL matrix inversion vars
deallocate(ipiv)
deallocate(work)

! deallocate NGd related vars
deallocate(Gfast)


! deallocaate scalar arrays
deallocate(parG)         
deallocate(factMatrix)
deallocate(MnmK1K2)    
! deallocaate multidim arrays
deallocate(kI)
deallocate(E)     
! deallocate(k)
deallocate(ktot)       
deallocate(G)          
deallocate(V)     
deallocate(S0) 
deallocate(Glf)      
deallocate(GlfV)     
deallocate(Imat)  
deallocate(diel_epsilon)
deallocate(Chi)   
deallocate(Chi0)  
deallocate(WT) 
deallocate(Gammap)
deallocate(Gammam)


contains
  subroutine findMinQ(Ntot, ktot, qx, qy, qz)
    ! searching min. q=(qx,qy,qz) in Gamma -> M direction
    integer,       intent(in)  :: Ntot
    real(kind=dp), intent(in)  :: ktot(:,:)
    real(kind=dp), intent(out) :: qx, qy, qz

    integer       :: i, ikmin
    real(kind=dp) :: kmin, kref, krefM ! , absq

    kmin = 1.0
    Ntot_loop: do  i = 1, Ntot ! loop over different k-points in FBZ
      kref = sqrt(sum(ktot(1:3,i)**2))
      ! neven debug
      print *,'i=',i,' kref: ',kref
      if (kref == 0.0) then
        CYCLE Ntot_loop
      else if (kref < kmin) then
        kmin = kref
        ikmin = i
        krefM = kmin
      end if
    end do Ntot_loop
    ! neve debug
    ! print *,'ikmin=',ikmin,'kmin=',kmin,'ktot(1:3,ikmin)',ktot
    qx = (iq-1) * ktot(1,ikmin)
    qy = (iq-1) * ktot(2,ikmin)
    qz = (iq-1) * ktot(3,ikmin)
    ! absq = sqrt(qx**2 + qy**2 + qz**2)
  end subroutine findMinQ


  subroutine findKinIBZ(ik, NkI, Nsymm, eps, kx, ky, kz, RI, kI, R1, K1)
    ! trazenje (kx,ky,kz) u ireducibilnoj zoni
    integer,       intent(in) :: ik
    integer,       intent(in) :: NkI, Nsymm
    real(kind=dp), intent(in) :: eps
    real(kind=dp), intent(in) :: kx, ky, kz
    real(kind=dp), intent(in) :: kI(:,:)
    real(kind=dp), intent(in) :: RI(:,:,:)
    integer      , intent(out):: R1, K1

    integer       :: i, j
    integer       :: it
    real(kind=dp) :: K11, K22, K33

    it = 1
    if (ik <= NkI) then
      R1 = 1
      K1 = ik
      it = 2
    else
      symmetry_loop: do  i = 2, Nsymm
        K11 = RI(i,1,1)*kx + RI(i,1,2)*ky + RI(i,1,3)*kz
        K22 = RI(i,2,1)*kx + RI(i,2,2)*ky + RI(i,2,3)*kz
        K33 = RI(i,3,1)*kx + RI(i,3,2)*ky + RI(i,3,3)*kz
        do  j = 1,NkI
          if (      abs(K11-kI(1,j)) <= eps &
              .and. abs(K22-kI(2,j)) <= eps &
              .and. abs(K33-kI(3,j)) <= eps ) then
            it = 2
            R1 = i
            K1 = j
            EXIT symmetry_loop
          end if
        end do
      end do symmetry_loop
    end if
    if (it == 1) then
      print*,'Can not find wave vector K=',ik, 'in I.B.Z.'
      STOP
    end if
    ! print *, 'R1:',R1,'K1:', K1
end subroutine findKinIBZ


subroutine findKQinBZ(KQx, KQy, KQz, eps, Nsymm, NkI, Ntot, NG, ktot, kI, RI, G, iG0, R2, K2)
  ! trazenje (KQx,KQy) prvo u FBZ a onda u IBZ
  integer,       intent(in)  :: Nsymm, NkI, Ntot, NG
  real(kind=dp), intent(in)  :: eps
  real(kind=dp), intent(in)  :: KQx, KQy, KQz
  real(kind=dp), intent(in)  :: ktot(:,:)
  real(kind=dp), intent(in)  :: kI(:,:)
  real(kind=dp), intent(in)  :: G(:,:)
  real(kind=dp), intent(in)  :: RI(:,:,:)
  integer,       intent(out) :: iG0, R2, K2

  integer :: it
  integer :: iG, jk, i, j

  it = 1
  iG_loop: do  iG = 1,NG
    k_loop_FBZ : do  jk = 1, Ntot
      if ( abs(KQx-G(1,iG)-ktot(1,jk)) <= eps .and. &
           abs(KQy-G(2,iG)-ktot(2,jk)) <= eps .and. &
           abs(KQz-G(3,iG)-ktot(3,jk)) <= eps ) then
        it = 2
        iG0 = iG
        symm_loop: do  i = 1, Nsymm
          K11 = sum(RI(i,1,1:3) * ktot(1:3,jk) )
          K22 = sum(RI(i,2,1:3) * ktot(1:3,jk) )
          K33 = sum(RI(i,3,1:3) * ktot(1:3,jk) )
          k_loop_IBZ: do  j = 1, NkI
            if ( abs(K11-kI(1,j)) <= eps .and. &
                 abs(K22-kI(2,j)) <= eps .and. &
                 abs(K33-kI(3,j)) <= eps ) then
              it = 3
              R2 = i
              K2 = j
              EXIT iG_loop
            end if
          end do k_loop_IBZ
        end do symm_loop
      end if
    end do k_loop_FBZ
  end do iG_loop  

  if (it == 1) then
    print*,'Can not find wave vector K+Q=',ik,'+',iq, 'in FBZ.'
    STOP
  else if (it == 2) then
    print*,'Can not find wave vector K+Q=',ik,'+',iq, 'in IBZ.'
    STOP
  end if

end subroutine findKQinBZ

  subroutine genFBZ(Nk,NkI,Nsymm,eps,kI,R,Ntot,ktot)
    ! Pomocu operacija tockaste grupe i vektora iz I.B.Z. 
    ! generira sve (medjusobno razlicite!!!) valne vektore u 1.B.Z.
    integer,       intent(in)  :: Nk, NkI, Nsymm
    real(kind=dp), intent(in)  :: eps 
    real(kind=dp), intent(in)  :: kI(:,:)
    real(kind=dp), intent(in)  :: R(:,:,:)

    integer,       intent(out) :: Ntot      ! ukupno jedinstvenih k-tocaka u FBZ
    real(kind=dp), intent(out) :: ktot(:,:) ! jedinstvene k-tocake u FBZ


    integer       :: it, jk
    integer       :: i
    integer       :: ik, lk
    integer       :: n, m
    real(kind=dp) :: k(3,Nk)


    jk = 0
    Ntot = 0 ! Total number of different points ''ktot'' inside 1.B.Z

    symm_loop: do  i = 1, Nsymm    ! loop over all symmetries
      k_loop_IBZ: do  ik = 1, NkI  ! loop over k points in IBZ
        it = 1
        jk = jk + 1
        do  n = 1, 3   ! loop over x,y,z 
          k(n,jk) = 0.0
          do  m = 1, 3 ! loop over x,y,z
            k(n,jk) = k(n,jk) + R(i,n,m)*kI(m,ik) ! kreira nove k tocke u BZ pomocu simetrije
          end do
        end do 

        if (jk > 1) then
          do  lk = 1, jk-1
            ! provjera jeli tocka razlicita od neke vec prije kreirane
            if ( abs(k(1,jk)-k(1,lk)) <= eps .and. &
                 abs(k(2,jk)-k(2,lk)) <= eps .and. &
                 abs(k(3,jk)-k(3,lk)) <= eps ) then 
              it = 2
            end if
          end do
        end if

        if (it == 1) then ! ne postoji dodaj ju
          ! !$omp atomic
          Ntot = Ntot+1
          ktot(1:3,Ntot) = k(1:3,jk)
        end if

      end do k_loop_IBZ
    end do symm_loop

    ! deallocate(k)

    ! output da vidimo kako izgleda FBZ
    open(887,FILE='fbz_check.dat',status='new')
    do  i = 1,Ntot
      write(887,*) ktot(1,i), ktot(2,i)  
    end do
    close(887)

  end subroutine genFBZ

  subroutine checkFBZintegration(Nband,NkI,Nsymm,Ntot,eps,kI,RI,Efermi,E,NelQE,Nel)
    ! Provjeri je li broj el. u FBZ (Nel) odgovara stvarnom broju el. u jed. cel. (NelQE)
    integer,       intent(in)  :: NelQE
    integer,       intent(in)  :: NkI, Nsymm, Nband, Ntot
    real(kind=dp), intent(in)  :: eps, Efermi
    real(kind=dp), intent(in)  :: kI(:,:)
    real(kind=dp), intent(in)  :: RI(:,:,:)
    real(kind=dp), intent(in)  :: E(:,:)
    real(kind=dp), intent(out) :: Nel

    integer       :: it, ik, n, i, j, K1
    real(kind=dp) :: kx,ky,kz

    Nel = 0 
    k_loop_FBZ : do  ik = 1,Ntot
      kx = ktot(1,ik)
      ky = ktot(2,ik)
      kz = ktot(3,ik)
      band_loop: do  n = 1, Nband
        if (n == 1) then
            it = 1
          if (ik <= NkI) then
            K1 = ik
            it = 2
          else
            symm_loop: do  i = 2, Nsymm
              K11 = RI(i,1,1)*kx + RI(i,1,2)*ky + RI(i,1,3)*kz
              K22 = RI(i,2,1)*kx + RI(i,2,2)*ky + RI(i,2,3)*kz
              K33 = RI(i,3,1)*kx + RI(i,3,2)*ky + RI(i,3,3)*kz
              k_loop_IBZ: do  j = 1, NkI
                if ( abs(K11-kI(1,j)) <= eps .and. &
                     abs(K22-kI(2,j)) <= eps .and. &
                     abs(K33-kI(3,j)) <= eps ) then
                  it = 2
                  K1 = j
                  ! zbroji broj el. u prvoj vrpci
                  if (E(K1,n) < Efermi) then 
                    Nel = Nel + 1.0
                    ! print *,'Nel',Nel,'band:',n
                  end if  
                  cycle band_loop
                end if
              end do k_loop_IBZ
            end do symm_loop
          end if
          if (it == 1) then
            print*,'Can not find wave vector K=',ik, 'in I.B.Z.'
            STOP
          end if
        end if
        
        ! zbroji broj el. u preostalim vrpcama
        if (E(K1,n) < Efermi) then 
          Nel = Nel + 1.0
          ! print *,'Nel',Nel,'band:',n
        end if  
    
      end do band_loop
    end do k_loop_FBZ
    Nel = 2.0*Nel / Ntot ! zbroji za en. manje od fermijeve
    
  end subroutine checkFBZintegration

subroutine genGlfandParity(lf,Ecut,NG,Gcar,G,Nlf,Nlfd,parG,Glf)
  ! Generate Reciprocal vectors for crystal local field 
  ! effects calculations in array Glf(3,Nlf)

  integer,          intent(in)  :: lf, NG, Nlfd
  real(kind=dp),    intent(in)  :: Ecut
  real(kind=dp),    intent(in)  :: Gcar
  real(kind=dp),    intent(in)  :: G(:,:)
  integer,          intent(out) :: Nlf
  integer,          intent(out) :: parG(:)
  real(kind=dp),    intent(out) :: Glf(:,:)

  integer       :: iG
  real(kind=dp) :: Eref

  Nlf = 0
  if (lf == 1) then
    do iG = 1, NG
      ! local field efekti samo u okomitom smjeru
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
          end if
        end if
      end if
    end do
  else
    do  iG = 1, NG
      ! local field efekti samo u svim smjerovima
      Eref = Gcar**2*sum(G(1:3,iG)**2) / 2.0
      if (Eref <= Ecut) then
        Nlf = Nlf+1
        Glf(:,Nlf) = G(1:3,iG)
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
    STOP
  end if

  ! neven debug
  ! print *,'Eref',Eref
  ! print *,'Glf(1:3,1:5)',Glf(1:3,1:5)
end subroutine genGlfandParity

subroutine genMnmK1K2(jump, eps, Nlf, iG0, NG1, NG2, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2)
  ! Konstrukcija stupca matricnih elementa nabojnih vrhova MnmK1K2(G) 
  integer,          intent(in)    :: iG0, Nlf, NG1, NG2
  integer,          intent(in)    :: R1,R2
  real(kind=dp),    intent(in)    :: eps
  real(kind=dp),    intent(in)    :: R(:,:,:)
  real(kind=dp),    intent(in)    :: RI(:,:,:)
  real(kind=dp),    intent(in)    :: Glf(:,:)
  real(kind=dp),    intent(in)    :: G(:,:)     ! polje valnih vektora G u recp. prost. za wfn.
  complex(kind=dp), intent(in)    :: C1(:)
  complex(kind=dp), intent(in)    :: C2(:)
  integer,          intent(inout) :: jump
  integer,          intent(inout) :: Gfast(:)
  complex(kind=dp), intent(out)   :: MnmK1K2(:)

  integer       :: iG, iG1, iG2
  integer       :: iGfast
  real(kind=dp) :: K11, K22, K33
  real(kind=dp) :: Gxx1,Gyy1,Gzz1
  real(kind=dp) :: Gxx2,Gyy2,Gzz2


  iGfast = 0
  MnmK1K2(1:Nlf) = cmplx(0.0) ! nabojni vrhovi
  do  iG = 1,Nlf ! suma po lokalnim fieldovima kojih ima Nlf
    do  iG1 = 1,NG1 ! vito zamjenjeno NGd sa NG1
      iGfast = iGfast + 1
      Gxx1 = G(1,iG1)
      Gyy1 = G(2,iG1)
      Gzz1 = G(3,iG1)
      K11 = R(R1,1,1)*Gxx1 + R(R1,1,2)*Gyy1 + R(R1,1,3)*Gzz1
      K22 = R(R1,2,1)*Gxx1 + R(R1,2,2)*Gyy1 + R(R1,2,3)*Gzz1
      K33 = R(R1,3,1)*Gxx1 + R(R1,3,2)*Gyy1 + R(R1,3,3)*Gzz1
      K11 = K11 + Glf(1,iG)
      K22 = K22 + Glf(2,iG)
      K33 = K33 + Glf(3,iG)
      K11 = K11 + G(1,iG0)
      K22 = K22 + G(2,iG0)
      K33 = K33 + G(3,iG0)
      Gxx1 = RI(R2,1,1)*K11 + RI(R2,1,2)*K22 + RI(R2,1,3)*K33
      Gyy1 = RI(R2,2,1)*K11 + RI(R2,2,2)*K22 + RI(R2,2,3)*K33
      Gzz1 = RI(R2,3,1)*K11 + RI(R2,3,2)*K22 + RI(R2,3,3)*K33
      ! !$omp critical(jump_operation)
      if (jump == 1) then
        ! !$omp single
        iG2_loop: do  iG2 = 1,NG2
          Gfast(iGfast) = NG2+1
          Gxx2 = G(1,iG2)
          Gyy2 = G(2,iG2)
          Gzz2 = G(3,iG2)
          if ( abs(Gxx2-Gxx1) < eps .and. &
               abs(Gyy2-Gyy1) < eps .and. &
               abs(Gzz2-Gzz1) < eps) then
              Gfast(iGfast) = iG2
              ! goto 1111
              EXIT iG2_loop
          end if
        end do iG2_loop
        ! !$omp end single
      end if
      ! !$omp end critical(jump_operation)
      ! 1111              continue
      iG2 = Gfast(iGfast)
      if (iG2 <= NG2) then
        MnmK1K2(iG) = MnmK1K2(iG) + conjg(C1(iG1))*C2(iG2)
        ! print *,'condition satisfied',MnmK1K2(iG)
      end if
    end do
  end do
  jump = 2
end subroutine genMnmK1K2

  subroutine genChi0(io,no,Nlf,domega,S0,Chi0)
    implicit none
    integer,          intent(in)    :: io
    integer,          intent(in)    :: no 
    integer,          intent(in)    :: Nlf
    real(kind=dp),    intent(in)    :: domega
    real(kind=dp),    intent(in)    :: S0(:,:,:)
    complex(kind=dp), intent(out)   :: Chi0(:,:)
    
    integer       :: jo, iG, jG
    real(kind=dp) :: oi, oj, fact
    real(kind=dp) :: ReChi0, ImChi0
   
    oi = (io-1)*domega
    Chi0(1:Nlf,1:Nlf) = cmplx(0.0,0.0)
    do  iG = 1,Nlf
      do  jG = 1,Nlf
        ReChi0 = 0.0
        ! static limit
        if (io == 1) then
          do  jo = 2,no
            oj = (jo-1)*domega
            fact = domega/oj
            ! analticki trikovi za integriranje 
            if (jo == 2) then 
              fact = 3.0/2.0
            else if (jo == no) then
              fact = 0.5*domega/oj
            end if
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
          end do
          ReChi0 = -2.0 * ReChi0
        else if (io == 2) then
          do  jo = 1,no
            oj = (jo-1)*domega
            if (jo /= io) then
              fact = domega/(oi-oj)
            else if (jo == 1) then
              fact = 1.0
            else if (jo == 2) then
              fact = 0.0
            else if (jo == 3) then
              fact = -3.0/2.0
            else if (jo == no) then
              fact = 0.5*domega/(oi-oj)
            end if
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            if (jo == 1 .or. jo == no) then
              fact=0.5*domega/(oi+oj)
            end if
            ReChi0 = ReChi0 - fact*S0(jo,iG,jG)
          end do
        else if (io == (no-1)) then
          do  jo = 1,no
            oj = (jo-1)*domega
            if (jo /= io) then
              fact = domega/(oi-oj)
            else if (jo == 1) then
              fact = 0.5*domega/(oi-oj)
            else if (jo == (no-2)) then
              fact = 3.0/2.0
            else if (jo == (no-1)) then
              fact = 0.0
            else if (jo == no) then
              fact = -1.0
            end if
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            if (jo == 1 .or. jo == no) then
              fact = 0.5*domega/(oi+oj)
            end if
            ReChi0 = ReChi0 - fact*S0(jo,iG,jG)
          end do
        else
          do  jo = 1,no
            oj = (jo-1)*domega
            if (jo /= io) then
              fact = domega/(oi-oj)
            else if (jo == 1) then
              fact=0.5*domega/(oi-oj)
            else if (jo == (io-1)) then
              fact=3.0/2.0
            else if (jo == io) then
              fact=0.0
            else if (jo == (io+1)) then
              fact=-3.0/2.0
            else if (jo == no) then
              fact = 0.5*domega/(oi-oj)
            end if
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            if (jo == 1 .or. jo == no) then
              fact = 0.5*domega/(oi+oj)
            end if
            ReChi0 = ReChi0 - fact*S0(jo,iG,jG)
          end do
        end if

        ImChi0 = -pi*S0(io,iG,jG)
        Chi0(iG,jG) = cmplx(ReChi0,ImChi0)

        ! kraj po iG,jG
      end do
    end do
  end subroutine genChi0


  subroutine genW1W2(io,no,domega,S0,W1,W2)
    implicit none
    integer,          intent(in)    :: io
    integer,          intent(in)    :: no 
    real(kind=dp),    intent(in)    :: domega
    real(kind=dp),    intent(in)    :: S0(:,:,:)
    real(kind=dp),    intent(out)   :: W1, W2   

    integer        :: jo  
    real(kind=dp)  :: fact 

    W1 = 0.0
    W2 = 0.0
    ! static limit
    if (io == 1) then
      do  jo = 2,no
        oj = (jo-1)*domega
        fact = domega/oj
        
        if (jo == 2) then 
          fact = 3.0/2.0
        else if (jo == no) then
          fact = 0.5*domega/oj
        end if

        W1 = W1 - fact*S0(jo,iG,jG)
      end do

      W2 = -W1

    else if (io == 2) then
      do  jo = 1,no
        oj= (jo-1)*domega
        
        if (jo /= io) then
          fact = domega/(oi-oj)
        else if (jo == 1) then
          fact = 1.0
        else if (jo == 2) then
          fact = 0.0
        else if (jo == 3) then
          fact = -3.0/2.0
        else if (jo == no) then
          fact = 0.5*domega/(oi-oj)
        end if             
        
        W1 = W1 + fact*S0(jo,iG,jG)
        fact = domega/(oi+oj)
        
        if (jo == 1 .or. jo == no) then
          fact = 0.5*domega/(oi+oj)
        end if

        W2 = W2 + fact*S0(jo,iG,jG)
      end do

    else if (io == (no-1)) then
      do  jo = 1,no
        oj = (jo-1)*domega

        if (jo /= io) then 
          fact = domega/(oi-oj)
        else if (jo == 1) then
          fact = 0.5*domega / (oi-oj)
        else if (jo == (no-2)) then
          fact = 3.0/2.0
        else if (jo == (no-1)) then
          fact = 0.0
        else if (jo == no) then
          fact = -1.0
        end if

        W1 = W1 + fact*S0(jo,iG,jG)
        fact = domega / (oi+oj)
        
        if (jo == 1 .or. jo == no) then
          fact = 0.5*domega / (oi+oj)
        end if

        W2 = W2 + fact*S0(jo,iG,jG)
      end do

    else
      do  jo = 1,no
        oj = (jo-1)*domega
        !vito - promjena if if if u else if
        if (jo /= io) then 
          fact = domega/(oi-oj)
        else if (jo == 1) then
          fact = 0.5*domega/(oi-oj)
        else if (jo == (io-1)) then
          fact = 3.0/2.0
        else if (jo == io) then
          fact = 0.0
        else if (jo == (io+1)) then
          fact = -3.0/2.0
        else if (jo == no) then
          fact = 0.5*domega/(oi-oj)
        end if
        W1 = W1 + fact*S0(jo,iG,jG)
        fact = domega/(oi+oj)
        
        if (jo == 1 .or. jo == no) then
          fact = 0.5*domega/(oi+oj)
        end if

        W2 = W2 + fact*S0(jo,iG,jG)
      end do
    end if
  end subroutine genW1W2

  subroutine genV(eps,qx,qy,Nlf,parG,Glf,GlfV,V)
    implicit none
    integer,        intent(in)    :: Nlf
    integer,        intent(in)    :: parG(:)
    real(kind=dp),  intent(in)    :: eps
    real(kind=dp),  intent(in)    :: qx,qy
    real(kind=dp),  intent(in)    :: Glf(:,:)
    real(kind=dp),  intent(in)    :: GlfV(:,:)
    real(kind=dp),  intent(inout) :: V(:,:)


    integer       :: iG, jG
    real(kind=dp) :: Gabs

    do  iG = 1,Nlf
      Gabs = sqrt( (qx+GlfV(1,iG))**2 + (qy+GlfV(2,iG))**2 )
      if (Gabs == 0.0) then 
        Gabs = eps
      end if
      do  jG = 1,Nlf
        V(iG,jG) = 0.0
        if (Glf(1,jG) == Glf(1,iG) ) then
          if (Glf(2,jG) == Glf(2,iG) ) then
            V(iG,jG) = 4.0*pi*(1.0-exp(-Gabs*c0)) / (Gabs*c0)
            V(iG,jG) = V(iG,jG)*( Gabs**2 - GlfV(3,iG)*GlfV(3,jG) )
            V(iG,jG) = V(iG,jG)/( Gabs**2 + GlfV(3,iG)**2 )
            V(iG,jG) = V(iG,jG)/( Gabs**2 + GlfV(3,jG)**2 )
            V(iG,jG) = -real(parG(iG)) * real(parG(jG)) * V(iG,jG) ! dble converted to real
            if (Glf(3,jG) == Glf(3,iG)) then
              V(iG,jG) = 4.0*pi / ( Gabs**2 + GlfV(3,iG)**2 ) + V(iG,jG)
            end if
          end if
        end if
      end do
    end do
  end subroutine genV

  subroutine genDielectricEpsilon(Nlf,Chi0,V,diel_epsilon)
    ! Epsilon (GG')  = I - V(GG')Chi0
    implicit none
    integer,          intent(in)    :: Nlf
    real(kind=dp),    intent(in)    :: V(:,:)    
    complex(kind=dp), intent(in)    :: Chi0(:,:)
    complex(kind=dp), intent(inout) :: diel_epsilon(:,:)

    integer :: iG, jG, kG
    complex(kind=dp) :: Imat(Nlf,Nlf)

    Imat(1:Nlf,1:Nlf) = cmplx(0.0,0.0)
    do  iG = 1,Nlf
      Imat(iG,iG) = cmplx(1.0,0.0)
    end do

    do  iG = 1,Nlf
      do  jG = 1,Nlf
        diel_epsilon(iG,jG) = Imat(iG,jG)
        do  kG = 1,Nlf
          diel_epsilon(iG,jG) = diel_epsilon(iG,jG) - Chi0(iG,kG)*V(kG,jG)
        end do
      end do
    end do
  end subroutine genDielectricEpsilon

  subroutine genChi(Nlf,diel_epsilon,Chi0,Chi)
    implicit none
    integer,          intent(in)     :: Nlf
    complex(kind=dp), intent(in)     :: diel_epsilon(:,:)
    complex(kind=dp), intent(in)     :: Chi0(:,:)
    complex(kind=dp), intent(inout)  :: Chi(:,:)

    integer :: iG, jG, kG

    Chi(1:Nlf,1:Nlf) = cmplx(0.0,0.0)
    do  iG = 1,Nlf
      do  jG = 1,Nlf
        do  kG = 1,Nlf
          Chi(iG,jG) = Chi(iG,jG) + diel_epsilon(iG,kG)*Chi0(kG,jG)
        end do
      end do
    end do
  end subroutine genChi

  subroutine genWT(io,Nlf,V,Chi,WT)
    !  SCREENED COULOMB INTERACTION W^T_GG'(Q,\omega)
    implicit none
    integer,          intent(in)    :: io
    integer,          intent(in)    :: Nlf
    real(kind=dp),    intent(in)    :: V(:,:)
    complex(kind=dp), intent(in)    :: Chi(:,:)
    complex(kind=dp), intent(inout) :: WT(:,:,:)

    integer :: iG, jG, kG1, kG2

    WT(io,1:Nlf,1:Nlf) = cmplx(0.0,0.0)
    do  iG = 1,Nlf
      do  jG = 1,Nlf
        do  kG1 = 1,Nlf
          do  kG2 = 1,Nlf
            WT(io,iG,jG) = WT(io,iG,jG) + V(iG,kG1)*Chi(kG1,kG2)*V(kG2,jG)
          end do
        end do
        WT(io,iG,jG) = V(iG,jG) + WT(io,iG,jG)
      end do
    end do
  end subroutine genWT

  subroutine loadKC(path,KC)
    character(len=100), intent(in)  :: path
    real(kind=dp),      intent(out) :: KC(3,3)
    
    integer :: j
    integer :: ios, lno=0
    character (len=35) :: tag, buffer

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

  subroutine loadG(KC,parG,G)
    ! Reading the reciprocal vectors in crystal coordinates and transformation
    ! in Cartesian cordinates.
    implicit none
    real(kind=dp),  intent(in)    :: KC(3,3)
    integer,        intent(inout) :: parG(:) ! paritet svakog valnog vektora G
    real(kind=dp),  intent(inout) :: G(:,:)  ! polje valnih vektora G u recp. prost. za wfn.

    integer :: ios1, ios2, lno
    integer :: n, m
    integer :: Gi(3)

    open(200,FILE='gvectors.xml',status='old',err=200,iostat=ios1)
    print *, 'File gvectors.xml oppened successfully.'
    lno = 0

    do  i=1,8
      read(200,*,err=201,iostat=ios2,end=202) ! dummy
      lno = lno + 1
    end do

    G(1:3,1:NG) = 0.0
    do  iG = 1,NG
      read(200,'(i10,i11,i11) ',err=201,iostat=ios2,end=202) Gi(1),Gi(2),Gi(3)
      lno = lno +1
      if (iG == 1) then
        if (Gi(1) /= 0 .or. Gi(2) /= 0 .or. Gi(3) /= 0) then
          print*,'*********************************'
          print*,'WARRNING!, G vectors input is wrong!!'
          print*,'G(1) is not (0,0,0)!!'
          STOP
        end if
      end if
    ! transformation in cart.coord (also!, after this all G components are in 2pi/a0 units)
      do n = 1,3
        do m = 1,3
          G(n,iG) = G(n,iG)+KC(n,m)*real(Gi(m)) ! DBLE converted to real
        end do
      end do
      parG(iG)=Gi(3)
    end do
  close(200)

  goto 5000
  200 write(*,*) 'error cant read file id 20, ist=',ist10
  201   write(*,*) '201 buffer1 read. Error reading line ',lno10+1,', iostat = ',ist11
  202   write(*,*) '202 buffer1 read. Number of lines read = ',lno10
  5000 continue 

  end subroutine loadG

  subroutine loadCs(thread_id,ik,savedir,K1,K2,n,m,NGd,jump, Nlf, iG0, eps, R, RI, G, Gfast, MnmK1K2)
    ! WARNING THIS SUBROUTINE IS BUGGY and therefore not used, MnmK1K2 gets assigned nonensensical values.

    ! Load Fourier coefficients (C1 spin up and C2 spin down) in the wfn. expansion 
    ! for a given (K1, K2, n, m) then generate MnmK1K2 for the loaded coefficients
    implicit none
    integer,            intent(in)    :: thread_id, ik
    integer,            intent(in)    :: K1, K2, n, m
    integer,            intent(in)    :: NGd
    character(len=100), intent(in)    :: savedir

    ! MnmK1K2 vars
    integer,            intent(in)    :: iG0, Nlf
    real(kind=dp),      intent(in)    :: eps
    real(kind=dp),      intent(in)    :: R(:,:,:)
    real(kind=dp),      intent(in)    :: RI(:,:,:)
    real(kind=dp),      intent(in)    :: G(:,:)    
    integer,            intent(inout) :: jump
    integer,            intent(inout) :: Gfast(:)
    complex(kind=dp),   intent(inout) :: MnmK1K2(:)
    

    integer                        :: NG1, NG2
    character(len=iotk_attlenx)    :: attr
    character(len=100)             :: pathk1, pathk2, bandn, bandm
    complex(kind=dp), allocatable  :: C1(:) ! pointer changed to allocatable (less memory efficent)
    complex(kind=dp), allocatable  :: C2(:)

    !$omp critical(pathk_read)
    ! otvara save/K.000x/evc.dat u atributu <evc band> ispod CnK(G) koef.
    call paths(savedir,K1,K2,n,m,pathk1,pathk2,bandn,bandm) 
    ! u ovom dijelu programa se iscitava iz binarnih fileova evc.dat za
    ! fiksni K1,K2, i vrpce n i m
    !Otvaranje atribute za INFO
    ! print *,'pathk1',pathk1
    call iotk_open_read(10+ik,pathk1)
    call iotk_scan_empty(10+ik,"INFO",attr=attr)
    call iotk_scan_attr(attr,"igwx",NG1)
    ! Alociranje polja C1
    allocate (C1(NG1))
    ! Ucitavanje podataka iza evc.n
    call iotk_scan_dat(10+ik,bandn,C1)
    call iotk_close_read(10+ik)
    !  Otvaranje atribute za INFO
    ! print *,'pathk2',pathk2
    call iotk_open_read(10+ik,pathk2)
    call iotk_scan_empty(10+ik,"INFO",attr=attr)
    call iotk_scan_attr(attr,"igwx",NG2)
    ! Alociranje polja C2
    allocate (C2(NG2))
    ! Ucitavanje podataka iza evc.m
    call iotk_scan_dat(10+ik,bandm,C2)
    call iotk_close_read(10+ik)
    !$omp end critical(pathk_read)

    if (NGd > NG1) then
      write(*,*) 'NGd is bigger than NG1=',NG1
      STOP
    else if (NGd > NG2) then
      write(*,*) 'NGd is bigger than NG2=',NG2
      STOP
    end if
    
    ! Konstrukcija stupca matricnih elementa nabojnih vrhova MnmK1K2(G)      
    call genMnmK1K2(jump, eps, Nlf, iG0, NG1, NG2, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2)

    ! neven debug
    ! print *,MnmK1K2(1:5)
    ! print *, jump, Nlf, iG0, NG1, NG2, eps, R(5,3,3), RI(5,3,3), G(3,3), Gfast(1:5),C1(4),C2(4)
    deallocate(C1)
    deallocate(C2)
  end subroutine loadCS

subroutine loadCsQE6( ik, ibnd, savedir, evc, igwx )
    ! read_a_wfc(ibnd, filename, evc, ik, xk, nbnd, ispin, npol, gamma_only, ngw, igwx )
    ! read QE 6.0 and greater, wfn coefficeints
    ! use iso_fortran_env, ONLY: DP=> REAL64
    implicit none 
    character (len=*), intent(in)               :: savedir
    integer,           intent(in)               :: ik, ibnd
    integer,           intent(out)              :: igwx 
    complex(DP),       intent(out), allocatable :: evc(:)

    character (len=300) :: path 

    integer  :: nbnd, ispin, npol, ngw, i, ik2
    real(dp) :: xk(3)
    ! integer  :: dummy_int   
    logical  :: gamma_only 
    integer  :: ios,iuni
    real(dp) :: scalef
    real(dp) :: b1(3), b2(3), b3(3) !, dummy_real 

    character(len=3)   :: str1 = 'wfc'
    character(len=4)   :: str3 ='.dat'
    character(len=100) :: ik_str
    write(ik_str,'(I10)') ik

    path = trim(savedir)//str1//trim(adjustl(ik_str))//str3

    iuni = 10 + ik
    
    open(unit = iuni, file = trim(adjustl(path)), form = 'unformatted', status = 'old', iostat=ios) 
    ! read(iuni) ik2, xk, ispin, gamma_only, scalef
    read(iuni) ! skip line
    read(iuni) ngw, igwx, npol, nbnd
    read(iuni) b1, b2, b3 

    ! avoid reading miller indices of G vectors below E_cut for this kpoint 
    ! if needed allocate  an integer array of dims (1:3,1:igwx) 

    read(iuni) ! dummy_int
    allocate (evc(npol*igwx))
    if ( ibnd > nbnd) then 
       print '("looking for band nr. ",I7," but there are only ",I7," bands in the file")',ibnd, nbnd
       stop
    end if 
    do i = 1, nbnd 
       if ( i == ibnd ) then 
          read(iuni) evc(1:npol*igwx) 
          exit 
       else 
          read(iuni) ! dummy_real
       end if 
    end do 
    close(iuni) 
    end subroutine loadCsQE6


  subroutine loadkIandE(path, NkI, Nband, Nocc, kI, E)
    implicit none
    integer,            intent(in)    :: NkI
    integer,            intent(in)    :: Nband, Nocc
    character(len=100), intent(in)    :: path
    real(kind=dp),      intent(inout) :: kI(:,:)
    real(kind=dp),      intent(inout) :: E(:,:)

    integer :: ios
    real(kind=dp),    parameter :: Hartree = 2.0D0*13.6056923D0

    open(400,FILE=path,status='old',err=500,iostat=ios) 
    do  ik = 1,NkI
      if (ik == 1) then
        read(400,*) 
      end if
      read(400,'(10X,3F10.3) ') kI(1,ik),kI(2,ik),kI(3,ik)
      read(400,'(10F8.4) ') (E(ik,i),i=1,Nband)
    end do
    close(400)
      
    goto 400
    500 write(*,*) 'Cannot open BAND file. iostat = ',ios
    stop
    400 continue
    
    ! konverzija en. u Hartree
    E(1:NkI,1:Nband) = E(1:NkI,1:Nband)/Hartree
    ! scissor operator, ispravlja/shifta DFT gap (na 1eV u ovom slucaju)
    E(1:NkI,Nocc+1:Nband) = E(1:NkI,Nocc+1:Nband) + 1.0/Hartree

  end subroutine loadkIandE

  subroutine writeInfo(qx, qy, qz, Gcar,Nsymm,Nlf, Ntot, NkI, Nband, eta, T, Nel, NelQE )
    implicit none
    integer       , intent(in) :: NelQE, Nsymm, Nlf, Ntot, NkI, Nband
    real(kind=dp) , intent(in) :: qx,qy,qz
    real(kind=dp) , intent(in) :: eta, T, Gcar
    real(kind=dp) , intent(in) :: Nel

    integer       :: ios
    real(kind=dp) :: absq, error
    real(kind=dp), parameter :: Hartree = 2.0D0*13.6056923D0

    absq = sqrt(qx**2+qy**2+qz**2)

    open(55,FILE='Info', err=700, iostat=ios)
    write(55,*) '***************General***********************'
    write(55,*) ''
    write(55,*) 'Number of point symmetry operation is',Nsymm
    write(55,'(a25,3f10.4,a5) ') 'Wave vector (qx,qy,qz)=(',qx*Gcar,qy*Gcar, qz*Gcar,') a.u.'
    write(55,'(a25,f8.4,a5) ') '|(qx,qy,qz)|=',absq*Gcar,'a.u.'
    if (lf == 1)write(55,*) 'Local field effcts in z-dir'
    if (lf == 3)write(55,*) 'Local field in all xyz-dir'
    write(55,*) 'Number of local field vectors is',Nlf
    write(55,*) 'Number of different K vectors in 1.B.Z. is',Ntot
    write(55,*) 'Number of K vectors in I.B.Z. is',  NkI
    write(55,*) 'Number of bands is               ', Nband
    write(55,'(a25,f8.4,a5) ') 'Eta damping is ',    eta*Hartree*1000.0,'meV'
    write(55,'(a25,f8.4,a5) ') 'Temperature is  ',   T*Hartree*1000.0,'meV'
    write(55,*) ''
    write(55,*) '************* Checking 1BZ integration*******'
    write(55,*) ''
    write(55,'(a40,f7.4) ') 'Number of electrons(1BZ integration)=',Nel
    write(55,*) 'Number of electrons(unit cell)=',NelQE
    error = abs((NelQE-Nel)/NelQE)
    write(55,'(a25,f8.4,a5) ') 'Relative error=',error*100.0,'%'
    if (error > 0.05) then
      write(55,*) 'WARRNING!!-1BZ INTEGRATION IS BAD!.'
      STOP
    end if
    close(55)
    
    goto 600
    700 write(*,*) 'Cannot open Info file. iostat = ',ios
    stop
    600 continue

  end subroutine writeInfo

  subroutine writeWT_Qi(iq,Nlf,domega,WT)
    ! ispis time ordered zasjenjene 
    ! kulonske interakcije W_GG'^T(Q,\omega)
    implicit none
    integer,          intent(in) :: iq, Nlf
    real(kind=dp),    intent(in) :: domega
    complex(kind=dp), intent(in) :: WT(:,:,:)

    integer :: nord, iG, jG, io
    real(kind=dp) :: o
    character(len=100) :: dato 

    dato = 'W_Qi'
    nord = index(dato,'i', back =.false.)
    if (iq < 10) then
      write(dato(nord:nord),'(i1) ')iq
    else if (iq >= 10 .and. iq < 100) then
      write(dato(nord:nord+1),'(i2) ')iq
    else
      write(dato(nord:nord+2),'(i3) ')iq
    end if
    
    open(74,FILE=dato)
    do  io = 1,1
      o = (io-1)*domega
      write(74,*) 'omega=',o,'Hartree'
      write(74,'(10F15.5) ')((WT(io,iG,jG),jG=1,Nlf),iG=1,Nlf)
    end do 
    close(74)
  end subroutine writeWT_Qi

subroutine writeKramKron_Qi(iq,qx, qy, qz, Gcar, KKS, SKK, G0, WT, V)
  ! ispisuje info file za provjeru Kramers–Kronig relacija
  integer,          intent(in) :: iq
  real(kind=dp),    intent(in) :: qx, qy, qz
  real(kind=dp),    intent(in) :: KKS, SKK
  real(kind=dp),    intent(in) :: Gcar
  real(kind=dp),    intent(in) :: V(:,:)
  complex(kind=dp), intent(in) :: G0
  complex(kind=dp), intent(in) :: WT(:,:,:) 

  integer            :: nord, ios
  real(kind=dp)      :: absq
  character(len=100) :: dato

  absq = sqrt(qx**2+qy**2+qz**2)

  dato = 'Kramers-Kron_Qi'
  nord = INDEX(dato,'i', back =.false.)
  if (iq < 10) then
    write(dato(nord:nord),'(i1) ') iq
  else if (iq >= 10 .and. iq < 100) then
    write(dato(nord:nord+1),'(i2) ') iq
  else
    write(dato(nord:nord+2),'(i3) ') iq
  end if
  

  open(33,FILE=dato, err=900, iostat=ios)
  write(33,'(a25,3f10.4,a5) ') 'Wave vector (qx,qy,qz)=(',qx*Gcar,qy*Gcar, qz*Gcar,') a.u.'
  write(33,'(a25,f8.4,a5) ') '|(qx,qy,qz)|=',absq*Gcar,'a.u.'
  write(33,*) 'int(WindKK-Wind)^2 =  ',KKS
  write(33,*) 'int(Wind)^2 =  ',SKK
  write(33,*) '****************************************'
  write(33,*) 'Kramers–Kronig relation relative error'
  write(33,'(5X,f7.2,a2) ') 100.0*abs(KKS/SKK), '%'
  write(33,*) 'Usporedba Gamma i WT'
  write(33,*) 'real[Gamma(o=0,1,1)]=',real(G0)
  write(33,*) 'real[WT(o=0,1,1)]/2=', real(WT(1,1,1)-V(1,1))/2.0
  close(33)

  goto 800
  900 write(*,*) 'Cannot open Kramers-Kron_Qi file. iostat = ',ios
  stop
  800 continue
end subroutine writeKramKron_Qi

end program surface_loss
