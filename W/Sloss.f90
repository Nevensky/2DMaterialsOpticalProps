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
!integer, parameter :: sp = real32
!integer, parameter :: dp = real64
!integer, parameter :: qp = real128

! integer, parameter :: sp = SELECTED_REAL_KIND(p=6)
! integer, parameter :: dp = SELECTED_REAL_KIND(p=15)
! integer, parameter :: qp = SELECTED_REAL_KIND(p=33)

! NOT used, read directly from QE
! integer :: nMPx,nMPy,nMPz
! parameter(nMPx=201,nMPy=201,nMPz=1)


integer :: debugCount= 0.0


integer :: ik,i,j,jk,it,lk,Ntot,iG0,Nsymm,iq, &
           io,n,m,iG,R1,K1,R2,K2, Nlf,NG1,   &
           NG2,iG1,iG2,jG,kG,jo,jump,loss,   &
           iGfast,ikmin,lf,kG1,kG2,nord

integer, parameter :: NkI = 6835 ! number of wave vectors in IBZ
integer, parameter :: Nband = 60 ! number of bands
integer, parameter :: NelQE = 18 ! Number of electrons(unit cell)
integer, parameter :: Nk = 48*NkI ! number of wave vectors in FBZ with no symmetry 
integer, parameter :: NGd = 4000 ! number of coefficients CG shulod be less than minimum number of coefficients all over all evc.n files ... moglo bi se dinamicki alocirati 
integer, parameter :: NG = 44121! zasto 8000 ?? total number of G vectors  
integer, parameter :: no = 2001 ! broj frekvencija
integer, parameter :: nq = 2 ! broj valnih vektora tu je 2 jer je rucno paralelizirano!
integer, parameter :: Nlfd = 50 ! dimenzija polja za local field zasto prozivoljno 50, ne moze se znati unaprijed

! file i/o debug
integer :: ist,ist2,ist9,ist10,ist11,ist12
integer :: lno,lno2,lno9,lno10,lno11,lno12

! scalars
real(kind=sp) :: kx,ky,kz
real(kind=sp) :: KQx,KQy,KQz
real(kind=sp) :: qgx,qgy,qgz
real(kind=sp) :: omin,omax
real(kind=sp) :: qx,qy,qz
real(kind=sp) :: kmin
real(kind=sp) :: domega
real(kind=sp) :: o
real(kind=sp) :: K11,K22,K33
real(kind=sp) :: Lor ! lorentzian
real(kind=sp) :: De 
real(kind=sp) :: Gabs ! 
real(kind=sp) :: kref ! trazi najmanju k-tocku sampliranu u MP meshu u kojem se moze izracunati ILS
real(kind=sp) :: Eref 
real(kind=sp) :: Gxx1,Gyy1,Gzz1
real(kind=sp) :: Gxx2,Gyy2,Gzz2
real(kind=sp) :: fact
real(kind=sp) :: oi,oj
real(kind=sp) :: ImChi0, ReChi0
! real(kind=sp) :: q ! ne koristi se
! real(kind=sp) :: qmax ! ne koristi se
real(kind=sp) :: Nel ! Number of electrons(1BZ integration)
real(kind=sp) :: absq
real(kind=sp) :: error
real(kind=sp) :: W1,W2
real(kind=sp) :: ImW
real(kind=sp) :: Wind
real(kind=sp) :: W2KK
real(kind=sp) :: KKS = 0.0
real(kind=sp) :: SKK = 0.0
real(kind=sp) :: WindKK
real(kind=sp) :: krefM

real(kind=sp), parameter :: pi = 4.D0*ATAN(1.D0)
real(kind=sp), parameter :: eV = 1.602176487D-19
real(kind=sp), parameter :: Hartree = 2.0D0*13.6056923D0
real(kind=sp), parameter :: Planck = 6.626196D-34
real(kind=sp), parameter :: three = 3.0d0 

real(kind=sp), parameter :: Efermi = 0.5554/Hartree ! Fermijeva en. u eV
real(kind=sp), parameter :: a0 = 5.9715 ! unit cell parameter in parallel direction in a.u.  
real(kind=sp), parameter :: c0 = 29.8575
real(kind=sp), parameter :: Gcar = 2.0*pi/a0 ! unit cell norm.
real(kind=sp), parameter :: eps = 1.0D-3 ! 1.0D-4 threshold
real(kind=sp), parameter :: T = 0.01/Hartree ! temperature in eV
real(kind=sp), parameter :: eta = 0.05/Hartree ! damping i\eta
real(kind=sp), parameter :: Ecut = 0.0 ! cutoff energy in Hartree for crystal local field calculations , for Ecut=0 S matrix is a scalar
real(kind=sp), parameter :: Vcell = 922.0586 ! unit-cell volume in a.u.^3 
real(kind=sp), parameter :: aBohr = 0.5291772D0 ! unit cell parameter in perpendicular direction in a.u. (z-separation between supercells)   



complex(kind=sp), parameter :: rone = cmplx(1.0,0.0)
complex(kind=sp), parameter :: czero = cmplx(0.0,0.0)
complex(kind=sp), parameter :: ione = cmplx(0.0,1.0)


complex(kind=sp) :: G0

! scalar arrays
integer, dimension(Nlfd*NGd) :: Gfast ! pomocna funkcija
integer, dimension(3) :: Gi ! pomocna funkcija
integer, dimension(NG) :: parG ! paritet svakog valnog vektora

! multidim arrays
real(kind=sp), dimension(3,NkI) :: kI
real(kind=sp), dimension(NkI,Nband) :: E ! vl. vr. danog k-i i band-i
real(kind=sp), dimension(48,3,3) :: R ! matr. simetrijskih op.
real(kind=sp), dimension(48,3,3) :: RI ! inverz od R
real(kind=sp), dimension(3,Nk) :: k
real(kind=sp), dimension(3,NK) :: ktot ! ukupno jedinstvenih k-tocaka u FBZ
real(kind=sp), dimension(3,NG) :: G ! polje valnih vektora G u recp. prost. za wfn.
real(kind=sp), dimension(Nlfd,Nlfd) :: V ! matr. gole coulomb. int.
real(kind=sp), dimension(no,Nlfd,Nlfd) :: S0 ! korelacijska matrica
real(kind=sp), dimension(3,3) :: KC ! pomocna funkcija
real(kind=sp), dimension(3,Nlfd) :: Glf ! local field effect polje valnih vekt. u rec. prost.
real(kind=sp), dimension(3,Nlfd) :: GlfV ! local field effect za nescreenanu (golu) int. V


! Nldf dimenziju ne znamo a priori zapravo, trebalo bi staviti sve te matrice allocatable i 
! naknadno ih alocirati

complex(kind=sp), dimension(Nlfd,Nlfd) :: Imat ! jedinicna matr.
complex(kind=sp), dimension(Nlfd,Nlfd) :: diel_epsilon ! Epsilon (GG')  = I - V(GG')Chi0
complex(kind=sp), dimension(Nlfd,Nlfd) :: Chi ! (eq. 2.88 nakon invertiranja) ;oprez bio je double precision

complex(kind=sp), dimension(Nlfd) :: MnmK1K2 ! nabojni vrhovi
complex(kind=sp), dimension(Nlfd,Nlfd) :: Chi0 ! (eq. 2.89)
complex(kind=sp), dimension(no,Nlfd,Nlfd) :: WT ! time ordered RPA screened coulomb int. (eq. 2.93)
complex(kind=sp), dimension(Nlfd,Nlfd) :: Gammap ! omega>0 ,eq....(skripta 5) \sum_{q,m} \int \dd omega' S(\omega')/{(\omega-\omega'-e_{k+q,m} +i\eta}) za GW se koristi se za ovaj dio 
complex(kind=sp), dimension(Nlfd,Nlfd) :: Gammam ! omega<0

character (len=100) :: bandn,bandm,nis,pathk1,pathk2,dato,root,outdir,path,fajl
character (len=35) :: tag,buffer

complex(kind=sp), pointer, dimension(:) :: C1,C2 ! Fourierovi koef. u razvoju wfn. iz QE 


! OpenMP vars
integer,parameter :: Nthreads = 4
integer :: thread_id


! MKL matrix inversion vars
integer :: info_trf, info_tri
integer,allocatable :: ipiv(:)
integer :: lwork
integer,allocatable :: work(:)


! real :: TES
! integer :: zora
! integer :: niba=100
! !$omp parallel shared(TES) private(zora)
! !$omp do reduction(+:TES)
! do zora=1,Niba
!    TES = TES+zora
! end do
! !$omp end do 
! !$omp master
! print *, 'TES:',TES
! !$omp end master
! !$omp end parallel 


! BRAVAIS LATTICE PARAMETERS

!     bravais-lattice index     =            4
!     lattice parameter (alat)  =       5.9715  a.u.
!     unit-cell volume          =     922.0586 (a.u.)^3
!     number of atoms/cell      =            3
!     number of atomic types    =            2
!     number of electrons       =        18.00
!     number of Kohn-Sham states=           13
!     kinetic-energy cutoff     =      50.0000  Ry
!     charge density cutoff     =     200.0000  Ry
!     convergence threshold     =      1.0E-06
!     mixing beta               =       0.7000
!     number of iterations used =            8  plain     mixing
!     Exchange-correlation      = SLA-PZ-NOGX-NOGC ( 1  1  0  0 0 0)

!     celldm(1)=   5.971535  celldm(2)=   0.000000  celldm(3)=   5.000000
!     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

!     crystal axes: (cart. coord. in units of alat)
!               a(1) = (   1.000000   0.000000   0.000000 )
!               a(2) = (  -0.500000   0.866025   0.000000 )
!               a(3) = (   0.000000   0.000000   5.000000 )

!     reciprocal axes: (cart. coord. in units 2 pi/alat)
!               b(1) = (  1.000000  0.577350 -0.000000 )
!               b(2) = (  0.000000  1.154701  0.000000 )
!               b(3) = (  0.000000 -0.000000  0.200000 )


!             QUANTUM ESSPRESSO IMPUTS:
root='../MoS2_201X201'
outdir='../MoS2_201X201/tmp'

lf = 1 ! crystal local field effect included in z for lf=1 or in x,y,z direction lf=3
loss = 1
jump = 1 ! za 1 preskace trazenje wfn. u IBZ za sve bands m i n
omin = 1.0D-5 ! raspon frekvencija u Hartreeima
omax = 2.0D0
domega = (omax-omin)/(no-1)



!           call FOR POINT GROUP TRANSFORMATIONS
!           Point group transformations are in Cartesian coordinate

call PointR(root,Nsymm,R,RI)
print *,"PointR done."


!           Upis valnih vektora iz irreducibilne Brillouinove
!           zone i pripadnih energijskih nivoa iz filea '****.band'.
!           wave vectors are in Cartesian coordinate



fajl='/MoS2.band'
path=trim(root)//trim(fajl)
open(40,FILE=path,status='old',err=400,iostat=ist9)

do  ik = 1,NkI
  if (ik == 1) then
    read(40,*) 
  end if
  read(40,'(10X,3F10.3) ') kI(1,ik),kI(2,ik),kI(3,ik)
  read(40,'(10F8.4) ') (E(ik,i),i=1,Nband)
end do
close(40)

! 10          FORMAT(10F8.4)
! 20          FORMAT(10X,f10.3,f10.3,f10.3)
GO TO 500
400 write(*,*) '100 cannot open file. iostat = ',ist9
500 continue


do  ik=1,NkI
  do  i=1,Nband
    E(ik,i) = E(ik,i)/Hartree
    if (i >= 10) then ! scissor operator, ispravlja/shifta DFT gap (na 1eV u ovom slucaju)
      E(ik,i) = E(ik,i) + 1.0/Hartree
    end if
  end do
end do





!            generator 1.B.Z.
!            Dio programa koji pomocu operacija tockaste grupe i vektora iz
!            I.B.Z. generira sve (MEDJUSOBNO RAZLICITE!!!) V. vektore u 1.B.Z.
!            Ntot-Tot number of different points ''ktot'' inside 1.B.Z


! !$omp parallel
! !$omp do
! do i=1,Nsymm
!   print *,i
! end do
! !$omp end do
! !$omp end parallel

! print *, 'DEBUG: entering parallel region'
! !$omp parallel shared(k,ktot,Ntot)

jk = 0
Ntot = 0
! !$omp do  reduction(+:Ntot)
do  i = 1,Nsymm ! loop over No. symmetries
  do  ik = 1,NkI  ! loop over k points in IBZ
    it = 1
    jk = jk+1
    do  n = 1,3 ! loop over kx,ky,kz 
      k(n,jk) = 0.0
      do  m = 1,3 ! loop over x,y,z
        k(n,jk) = k(n,jk) + R(i,n,m)*kI(m,ik) ! kreira nove k tocke u BZ pomocu simetrije
      end do
    end do
    if (jk > 1) then
      do  lk=1,jk-1
        ! vito - maknut nest if if if
        ! if (abs(k(1,jk)-k(1,lk)) <= eps) then
        !   if (abs(k(2,jk)-k(2,lk)) <= eps) then
        !     if (abs(k(3,jk)-k(3,lk)) <= eps) then
        !       it=2
        !     end if
        !   end if
        ! end if
        if ( abs(k(1,jk)-k(1,lk)) <= eps .and. &
             abs(k(2,jk)-k(2,lk)) <= eps .and. &
             abs(k(3,jk)-k(3,lk)) <= eps ) then ! je li razlicita tocka od neke prije vec kreirane
          it=2
        end if
      end do
    end if
    if (it == 1) then ! ne postoji dodaj ju
      Ntot=Ntot+1
      ktot(1:3,Ntot)=k(1:3,jk)
      ! ktot(1,Ntot)=k(1,jk)
      ! ktot(2,Ntot)=k(2,jk)
      ! ktot(3,Ntot)=k(3,jk)
    end if
  end do
end do
! !$omp end do
! !$omp end parallel
! print *, 'DEBUG: exiting parallel region'

! Checking 1BZ integration
! provjeri je li broj el. u FBZ (Nel) odgovara stvarnom broju el. u jed. cel. (NelQE)
Nel = 0 
k_loop_FBZ : do  ik = 1,Ntot
  kx = ktot(1,ik)
  ky = ktot(2,ik)
  kz = ktot(3,ik)
  band_loop: do  n = 1,Nband
    if (n == 1) then
        it = 1
      if (ik <= NkI) then
        K1 = ik
        it = 2
      else
        symm_loop: do  i = 2,Nsymm
          K11 = RI(i,1,1)*kx + RI(i,1,2)*ky + RI(i,1,3)*kz
          K22 = RI(i,2,1)*kx + RI(i,2,2)*ky + RI(i,2,3)*kz
          K33 = RI(i,3,1)*kx + RI(i,3,2)*ky + RI(i,3,3)*kz
          k_loop_IBZ: do  j = 1,NkI
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


open(887,FILE='fbz_check.dat',status='new')
do  i = 1,Ntot
  write(887,*) ktot(1,i),ktot(2,i)  ! output da vidimo kako izgleda FBZ
end do
close(887)



! KC transformation matrix from rec.cryst. axes to cart.koord.
! If G' is vector in rec.cryst. axes then a=KC*a' is vector in cart. axes


fajl='/MoS2.sc.out'
path=TRIM(root)//TRIM(fajl)
tag='     reciprocal axes: (cart. coord.'
open(30,FILE=path,status='old')
do  i = 1,100000
  read(30,'(a) ') buffer
  lno=lno+1
  if (buffer == tag) then
    do  j=1,3
      read(30,'(23X,3F10.3) ',err=10001,iostat=ist,end=20001) KC(1,j), KC(2,j), KC(3,j)
    end do
    EXIT
  end if
end do
! 70          FORMAT(23X,3F10.3)
close(30)


goto 8000
10001   write(*,*) 'Error reading line ',lno+1,', iostat = ',ist
20001   write(*,*) 'Number of lines read = ',lno
8000 continue


! Reading the reciprocal vectors in crystal coordinates and transformation
! in Cartesian cordinates.
open(20,FILE='gvectors.xml',status='old',err=200,iostat=ist10)
lno10 =0
print *, 'File gvectors.dat oppened successfully.'
do  i=1,8
  read(20,*,err=201,iostat=ist11,end=202) nis
  ! print *,nis
  lno10 = lno10 +1
end do
G = 0.0 ! vito - premjesteno iz n,m loopa
do  iG = 1,NG
  read(20,'(i10,i11,i11) ',err=201,iostat=ist11,end=202) Gi(1),Gi(2),Gi(3)
  lno10 = lno10 +1
  ! print *, lno10
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
    ! G(n,iG) = 0.0
    do m = 1,3
      G(n,iG) = G(n,iG)+KC(n,m)*real(Gi(m)) ! DBLE converted to real
    end do
  end do
  parG(iG)=Gi(3)
end do
! 100         FORMAT(i10,i11,i11)
close(20)


GO TO 5000
200 write(*,*) 'error cant read file id 20, ist=',ist10
201   write(*,*) '201 buffer1 read. Error reading line ',lno10+1,', iostat = ',ist11
202   write(*,*) '202 buffer1 read. Number of lines read = ',lno10
5000 continue


! Reciprocal vectors for crystal local field effects calculations in array ''Glf(3,Nlf) ''

Nlf = 0
if (lf == 1) then
  do iG = 1,NG
    ! local field efekti samo u okomitom smjeru
    if (G(1,iG) == 0.0 .and. G(2,iG) == 0.0) then
      Eref = Gcar**2 * G(3,iG)**2 / 2.0
      if (Eref <= Ecut) then
        Nlf = Nlf + 1
        Glf(1:2,Nlf) = 0.0
        ! Glf(1,Nlf) = 0.0
        ! Glf(2,Nlf) = 0.0
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
  do  iG = 1,NG
    ! local field efekti samo u svim smjerovima
    Eref = Gcar**2 *( G(1,iG)**2 + G(2,iG)**2 + G(3,iG)**2 ) / 2.0
    if (Eref <= Ecut) then
      Nlf = Nlf+1
      Glf(:,Nlf) = G(:,iG)
      ! Glf(1,Nlf) = G(1,iG)
      ! Glf(2,Nlf) = G(2,iG)
      ! Glf(3,Nlf) = G(3,iG)
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


! MKL matrix inversion vars
allocate(ipiv(MAX(1,MIN(Nlf, Nlf))))
lwork = Nlf
allocate(work(Nlf))

! neven debug
print *,'Eref',Eref
print *,'Glf(1:3,1:5)',Glf(1:3,1:5)


!             IBZ q LOOP STARTS HERE!!!
! iq=0 ne moze biti nula, opticki racun
! iq=2 do iq=...cutoff transfer q vektor!
! ikmin = min. valni vektor u BZ svi veci su visekratnici tog minimalnog
do  iq = 42,42 ! 42,61
  
!             searching min. q=(qx,qy,qz) in Gamma - M direction
  kmin = 1.0
  Ntot_loop: do  i = 1,Ntot ! loop over different k-points in FBZ
    ! kref = sqrt(sum(ktot(1:3,i)**2))
    kref=sqrt(ktot(1,i)*ktot(1,i)+ktot(2,i)*ktot(2,i)+ktot(3,i)*ktot(3,i)) 
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

  qx = (iq-1)*ktot(1,ikmin)
  qy = (iq-1)*ktot(2,ikmin)
  qz = (iq-1)*ktot(3,ikmin)
  absq = sqrt(qx**2 + qy**2 + qz**2)
  
!             Info file
  
  open(55,FILE='Info')
  write(55,*) '***************General***********************'
  write(55,*) ''
  write(55,*) 'Number of point symmetry operation is',Nsymm
  write(55,'(a25,3f10.4,a5) ') 'Wave vector (qx,qy,qz)=(',qx*Gcar,qy*Gcar, qz*Gcar,') a.u.'
  write(55,'(a25,f8.4,a5) ') '|(qx,qy,qz)|=',absq*Gcar,'a.u.'
  if (lf == 1)write(55,*) 'Local field effcts in z-dir'
  if (lf == 3)write(55,*) 'Local field in all xyz-dir'
  write(55,*) 'Number of local field vectors is',Nlf
  write(55,*) 'Number of different K vectors in 1.B.Z. is',Ntot
  write(55,*) 'Number of K vectors in I.B.Z. is',NkI
  write(55,*) 'Number of bands is               ',Nband
  write(55,'(a25,f8.4,a5) ') 'Eta damping is ',eta*Hartree*1000.0,'meV'
  write(55,'(a25,f8.4,a5) ') 'Temperature is  ',T*Hartree*1000.0,'meV'
  write(55,*) ''
  write(55,*) '************* Checking 1BZ integration*******'
  write(55,*) ''
  write(55,'(a40,f7.4) ') 'Number of electrons(1BZ integration)=',Nel
  write(55,*) 'Number of electrons(unit cell)=',NelQE
  error = abs((NelQE-Nel)/NelQE)
  write(55,'(a25,f8.4,a5) ') 'Relative error=',error*100.0,'%'
  if (error > 0.05) then
    write(55,*) 'WARRNING!!-1BZ INTEGRATION IS BAD!.'
  end if
  close(55)
  

  ! vito stari format inicializacije matrice
  ! do  io = 1,no
    ! do  iG = 1,Nlf
      ! do  jG=1,Nlf
        ! S0(io,iG,jG) = czero
      ! end do
    ! end do
  ! end do
  ! S0(1:no,1:Nlf,1:Nlf) = czero
  S0 = czero
  
  
!              1.B.Z  LOOP STARTS HERE !!!!

print *, 'DEBUG: entering parallel region'
!$omp parallel shared(S0,Gfast) ! default(private) num_threads(Nthreads)
! neven debug
! thread_id =  omp_get_thread_num()
! print *, 'thread id:',thread_id
!$omp do ! reduction(-:S0)
  do ik=1,Ntot   
! k_loop_FBZ_2nd:  

    ! neven maknuto zbog openmp
    ! open(122,FILE='status')
    ! write(122,*) 'iq=',iq
    ! write(122,*) 'ik=',ik
    ! close(122)
    
    
    kx = ktot(1,ik)
    ky = ktot(2,ik)
    kz = ktot(3,ik)
    
!  trazenje (kx,ky,kz) u ireducibilnoj zoni
    
    it = 1
    if (ik <= NkI) then
      R1=1
      K1=ik
      it=2
    else
      symmetry_loop: do  i = 2,Nsymm
        K11 = RI(i,1,1)*kx + RI(i,1,2)*ky+RI(i,1,3)*kz
        K22 = RI(i,2,1)*kx + RI(i,2,2)*ky+RI(i,2,3)*kz
        K33 = RI(i,3,1)*kx + RI(i,3,2)*ky+RI(i,3,3)*kz
        do  j = 1,NkI
          ! vito - smnanjenje if uvjeta
          ! if (abs(K11-kI(1,j)) <= eps) then
          !   if (abs(K22-kI(2,j)) <= eps) then
          !     if (abs(K33-kI(3,j)) <= eps) then
          !       it=2
          !       R1=i
          !       K1=j
          !       ! GO TO 5222
          !       EXIT symmetry_loop
          !     end if
          !   end if
          ! end if
          if (      abs(K11-kI(1,j)) <= eps &
              .and. abs(K22-kI(2,j)) <= eps &
              .and. abs(K33-kI(3,j)) <= eps ) then
            it=2
            R1=i
            K1=j
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
    ! 5222           continue
    
    
    it = 1
    KQx = kx + qx
    KQy = ky + qy
    KQz = kz + qz

    !$omp critical(printWaveVector)
    thread_id =  omp_get_thread_num()
    print *,'thread id:',thread_id,'ik: ',ik
    print *, 'KQx,KQy,KQz:',KQx,KQy,KQz
    !$omp end critical(printWaveVector)

!   trazenje (KQx,KQy) prvo u 1.B.Z a onda u I.B.Z.
    iG_loop: do  iG = 1,NG
      do  jk = 1,Ntot
        ! vito smanjenje If if if loop
        if ( abs(KQx-G(1,iG)-ktot(1,jk)) <= eps .and. &
             abs(KQy-G(2,iG)-ktot(2,jk)) <= eps .and. &
             abs(KQz-G(3,iG)-ktot(3,jk)) <= eps ) then
          it=2
          iG0 = iG
          do  i = 1,Nsymm
            ! K11 = RI(i,1,1)*ktot(1,jk) + RI(i,1,2)*ktot(2,jk) +  &
            !     RI(i,1,3)*ktot(3,jk)
            ! K22 = RI(i,2,1)*ktot(1,jk) + RI(i,2,2)*ktot(2,jk) +  &
            !     RI(i,2,3)*ktot(3,jk)
            ! K33 = RI(i,3,1)*ktot(1,jk) + RI(i,3,2)*ktot(2,jk) +  &
            !     RI(i,3,3)*ktot(3,jk)
            K11 = sum(RI(i,1,1:3) * ktot(1:3,jk) )
            K22 = sum(RI(i,2,1:3) * ktot(1:3,jk) )
            K33 = sum(RI(i,3,1:3) * ktot(1:3,jk) )
            do  j = 1,NkI
              if ( abs(K11-kI(1,j)) <= eps .and. &
                   abs(K22-kI(2,j)) <= eps .and. &
                   abs(K33-kI(3,j)) <= eps ) then
                it = 3
                R2 = i
                K2 = j
                EXIT iG_loop
                ! GO TO 2111
              end if
            end do
          end do
        end if
      end do
    end do iG_loop  

    ! vito prepravljeno maknuti ifovi
    ! ig_loop: do  iG = 1,NG
    !   do  jk = 1,Ntot
    !     ! vito smanjenje If if if loop
    !     if (abs(KQx-G(1,iG)-ktot(1,jk)) <= eps) then
    !       if (abs(KQy-G(2,iG)-ktot(2,jk)) <= eps) then
    !         if (abs(KQz-G(3,iG)-ktot(3,jk)) <= eps) then
    !           it=2
    !           iG0 = iG
    !           do  i = 1,Nsymm
    !             K11=RI(i,1,1)*ktot(1,jk)+RI(i,1,2)*ktot(2,jk)+  &
    !                 RI(i,1,3)*ktot(3,jk)
    !             K22=RI(i,2,1)*ktot(1,jk)+RI(i,2,2)*ktot(2,jk)+  &
    !                 RI(i,2,3)*ktot(3,jk)
    !             K33=RI(i,3,1)*ktot(1,jk)+RI(i,3,2)*ktot(2,jk)+  &
    !                 RI(i,3,3)*ktot(3,jk)
    !             do  j=1,NkI
    !               if (abs(K11-kI(1,j)) <= eps) then
    !                 if (abs(K22-kI(2,j)) <= eps) then
    !                   if (abs(K33-kI(3,j)) <= eps) then
    !                     it = 3
    !                     R2 = i
    !                     K2 = j
    !                     EXIT ig_loop
    !                     ! GO TO 2111
    !                   end if
    !                 end if
    !               end if
    !             end do
    !           end do
    !         end if
    !       end if
    !     end if
    !   end do
    ! end do ig_loop
    ! 2111                 continue
    
    
    if (it == 1) then
      print*,'Can not find wave vector K+Q=',ik,'+',iq, 'in 1.B.Z.'
      STOP
    else if (it == 2) then
      print*,'Can not find wave vector K+Q=',ik,'+',iq, 'in I.B.Z.'
      STOP
    end if
    
    
!              R1-integer, redni broj point operacije R1 u transformaciji ''K=R1*K1''.
!              K1-integer, redni broj valnog vektora K1 u transformaciji ''K=R1*K1''.
!              iG0 i R2-integeri, redni broj vektora reciprocne restke G0 i point operacije R2 u transformaciji ''K+Q=G0+R2*K2''.
!              K2-integer, redni broj valnog vektora K2 u transformaciji  ''K+Q=G0+R2*K2''.
    
    
!              petlje po vrpcama n i m
    
    bands_n_loop: do  n = 1,9 ! filled bands loop
      bands_m_loop: do  m = 10,Nband ! empty bands loop
        
        
        
        
        ! otvara save/K.000x/evc.dat u atributu <evc band> ispod CnK(G) koef.
        
        call paths(outdir,K1,K2,n,m,pathk1,pathk2,bandn,bandm) 
        
!         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za
!         fiksni K1,K2,n i m
        
        !$omp critical(pathk1_read)
!               Otvaranje atribute za INFO
        call iotk_open_read(10+ik,pathk1)
        call iotk_scan_empty(10+ik,"INFO",attr=attr)
        call iotk_scan_attr(attr,"igwx",NG1)
!               Alociranje polja C1
        allocate (C1(NG1))
!               Ucitavanje podataka iza evc.n
        call iotk_scan_dat(10+ik,bandn,C1)
        call iotk_close_read(10+ik)
        !$omp end critical(pathk1_read)

        !$omp critical(pathk2_read)
!               Otvaranje atribute za INFO
        call iotk_open_read(10+ik,pathk2)
        call iotk_scan_empty(10+ik,"INFO",attr=attr)
        call iotk_scan_attr(attr,"igwx",NG2)
!               Alociranje polja C2
        allocate (C2(NG2))
!               Ucitavanje podataka iza evc.m
        call iotk_scan_dat(10+ik,bandm,C2)
        call iotk_close_read(10+ik)
        !$omp end critical(pathk2_read)
        
!                Konstrukcija stupca matricnih elementa MnmK1K2(G)
        
        ! vito maknut check        
        ! if (NGd > NG1) then
        !   write(*,*) 'NGd is bigger than NG1=',NG1
        !   STOP
        ! else if (NGd > NG2) then
        !   write(*,*) 'NGd is bigger than NG2=',NG2
        !   STOP
        ! end if
        
        
!       matrix elements
        ! print *,'GOT HERE 1'
        iGfast = 0
        MnmK1K2 = czero ! nabojni vrhovi
        do  iG = 1,Nlf ! suma po lokalnim fieldovima kojih ima Nlf
          ! MnmK1K2(iG) = czero
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
            !$omp critical(jump_operation)
            if (jump == 1) then
              iG2_loop: do  iG2 = 1,NG2
                Gfast(iGfast) = NG2+1
                Gxx2 = G(1,iG2)
                Gyy2 = G(2,iG2)
                Gzz2 = G(3,iG2)
                ! vito if if if loop smanjen
                ! if (abs(Gxx2-Gxx1) < eps) then
                !   if (abs(Gyy2-Gyy1) < eps) then
                !     if (abs(Gzz2-Gzz1) < eps) then
                !       Gfast(iGfast) = iG2
                !       ! GO TO 1111
                !       EXIT ig2_loop
                !     end if
                !   end if
                ! end if
                if ( abs(Gxx2-Gxx1) < eps .and. &
                     abs(Gyy2-Gyy1) < eps .and. &
                     abs(Gzz2-Gzz1) < eps) then
                    Gfast(iGfast) = iG2
                    ! GO TO 1111
                    EXIT iG2_loop
                end if
              end do iG2_loop
            end if
            !$omp end critical(jump_operation)
            ! 1111              continue
            iG2 = Gfast(iGfast)
            ! ovo bi trebalo prouciti zasto ovako radi
            if (iG2 <= NG2) then
              MnmK1K2(iG) = MnmK1K2(iG) + conjg(C1(iG1))*C2(iG2)
            end if
          end do
        end do
        jump = 2
        
        ! print *,'GOT HERE 2'
!       omega loop
        do  io = 1,no
          o = (io-1)*domega
          De = o + E(K1,n) - E(K2,m) 
          Lor = -eta/(De**2 + eta**2) ! ovo bi analticki bila delta funkcija imag. dio od 1/De
          if (abs(Lor) >= 1.0D-3/eta) then
            do  iG = 1,Nlf
              do  jG = 1,Nlf
                S0(io,iG,jG) = S0(io,iG,jG) -  &
                    2.0*Lor*MnmK1K2(iG)*conjg(MnmK1K2(jG)) / (pi*Ntot*Vcell)
              end do
            end do
          end if
        ! neven debug  
        ! debugCount = debugCount +1
        ! if (debugCount>9000000 .and. debugCount<10000000) then
        ! write(118,*) io*Hartree,S0(io,1,1)
        ! end if 
        end do

        deallocate(C1)
        deallocate(C2)  
                
      end do bands_m_loop ! end of m do loop
      

    end do bands_n_loop ! end of n do loop
    jump=1

  end do !k_loop_FBZ_2nd !  end of 1.B.Z do loop
 !$omp end do
 !$omp end parallel
 print *, 'DEBUG: exiting parallel region'
  
! STOP

! Puting (qx,qy,qz) and Glf in cartesian coordinates
  
  qx = Gcar*qx ! convert qx *2p/a0
  qy = Gcar*qy
  qz = Gcar*qz
  
  do  iG = 1,Nlf
    GlfV(:,iG) = Gcar*Glf(:,iG)
  !   GlfV(1,iG) = Gcar*Glf(1,iG)
  !   GlfV(2,iG) = Gcar*Glf(2,iG)
  !   GlfV(3,iG) = Gcar*Glf(3,iG)
  end do
 

  
! Kramers-Kroning relacije
  
! new sum over omega
  omega_loop_A: do  io = 1,no-1
  ! print*,io
    oi = (io-1)*domega
    do  iG = 1,Nlf
      do  jG = 1,Nlf
        ReChi0 = 0.0
!       static limit
        if (io == 1) then
          do  jo = 2,no
            oj = (jo-1)*domega
            fact = domega/oj
            ! analticki trikovi za integriranje 
            if (jo == 2) then 
              fact = 3.0/2.0
            else if (jo == no) then
              fact = 0.5*domega/oj
            ! else
              ! print *,  "WARNING jo loop condition not satisfied."
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
            ! else
              ! print *,  "WARNING jo loop condition not satisfied."
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
            ! else 
              ! print *,  "WARNING jo loop condition not satisfied."
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
            ! else 
              ! print *,  "WARNING jo loop condition not satisfied."
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
        ! neven debug 
        ! print *,"Chi(",iG,jG,')=',Chi0(iG,jG)
        ! neven debug
        ! if (io==1) then
          ! write(18,*,action='write',position='append') ReChi0
          ! write(28,*,action='write',position='append') ImChi0
        ! end if

!                kraj po iG,jG
      end do
    end do
    ! neven debug 
    ! write(18,*) io*Hartree, ReChi0(1,1)
    
    
!                Calculation of the ''Chi''  by matrix invertion
    
    
!                MATRIX V(G,G')
    
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

    

    ! do  iG = 1,Nlf
    !   do jG = 1,Nlf
    !     Imat(iG,iG) = czero
    !   enddo
    !   Imat(iG,iG) = rone
    ! end do
    
    
    Imat = czero
    do  iG = 1,Nlf
      Imat(iG,iG) = rone
    end do

    do  iG = 1,Nlf
      do  jG = 1,Nlf
        diel_epsilon(iG,jG) = Imat(iG,jG)
        do  kG = 1,Nlf
          diel_epsilon(iG,jG) = diel_epsilon(iG,jG) - Chi0(iG,kG)*V(kG,jG)
        end do
      end do
    end do
    
    ! neven debug
    ! write(38,*) real(diel_epsilon(1,1)),aimag(diel_epsilon(1,1))
    
!  invertiranje matrice ''diel_epsilon = 1-Chi_0*V''
    
    ! call gjel(diel_epsilon,Nlf,Nlfd,Imat,Nlf,Nlfd)
    call sgetrf( Nlf,Nlfd, diel_epsilon, Nlf, ipiv, info_trf)
    call sgetri( Nlf, diel_epsilon, Nlf, ipiv, work, lwork, info_tri )
    Chi = czero
    do  iG = 1,Nlf
      do  jG = 1,Nlf
        do  kG = 1,Nlf
          Chi(iG,jG) = Chi(iG,jG) + diel_epsilon(iG,kG)*Chi0(kG,jG)
        end do
      end do
    end do
    
    
!  SCREENED COULOMB INTERACTION W^T_GG'(Q,\omega)
    WT(io,:,:) = czero
    do  iG = 1,Nlf
      do  jG = 1,Nlf
        ! WT(io,iG,jG) = czero
        do  kG1 = 1,Nlf
          do  kG2 = 1,Nlf
            WT(io,iG,jG) = WT(io,iG,jG) + V(iG,kG1)*Chi(kG1,kG2)*V(kG2,jG)
          end do
        end do
        WT(io,iG,jG) = V(iG,jG) + WT(io,iG,jG)
      end do
    end do

    write(20008,*) oi*Hartree,aimag(WT(io,1,1))
    write(10008,*) oi*Hartree,real(WT(io,1,1))
! kraj nove petlje po omega
  end do omega_loop_A
  
  ! neven debug
  check_WT : do io=1,no-1 
    write(88,*) real(WT(io,1,1))
    write(98,*) aimag(WT(io,1,1))
  end do check_WT

! neve debug openmp maknuto 
      ! ispis time ordered zasjenjene kulonske interakcije W_GG'^T(Q,\omega)
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
  ! 44              FORMAT(10F15.5)

  ! vito stari format  
  do  io = 1,no-1
    do  iG = 1,Nlf
      do  jG = 1,Nlf
        S0(io,iG,jG) = -(1.0/pi)*aimag(WT(io,iG,jG))
      end do
    end do
  end do
  ! S0(1:no-1,1:Nlf,1:Nlf) = -(1.0/pi)*aimag(WT(1:no-1,1:Nlf,1:Nlf))
  


  ! podatci za GW sve na dalje

  KKS = 0.0
  SKK = 0.0
  
  
! new sum over omega
  omega_loop_B: do  io = 1,no-1
    ! print*,io
    oi = (io-1)*domega
    do  iG = 1,Nlf
      do  jG = 1,Nlf
        W1 = 0.0
        W2 = 0.0
!      static limit
        if (io == 1) then

          do  jo = 2,no
            oj = (jo-1)*domega
            fact = domega/oj
            ! vito - promjena if if if u else if
            if (jo == 2) then 
              fact = 3.0/2.0
            else if (jo == no) then
              fact = 0.5*domega/oj
            ! else
              ! print *,  "WARNING jo loop condition not satisfied."
            end if
            W1 = W1 - fact*S0(jo,iG,jG)
          end do
          W2 = -W1

        else if (io == 2) then

          do  jo = 1,no
            oj= (jo-1)*domega
            !vito - promjena if if if u else if
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
            ! else
              ! print *,  "WARNING jo loop condition not satisfied."
            end if             
            W1 = W1 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            
            if (jo == 1 .or. jo == no) then
              fact = 0.5*domega/(oi+oj)
            end if
            W2 = W2 + fact*S0(jo,iG,jG)
          end do

        else if (io == (no-1)) then

          do  jo=1,no
            oj= (jo-1)*domega
            !vito - promjena if if if u else if
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
            ! else
              ! print *,  "WARNING jo loop condition not satisfied."
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
            oj= (jo-1)*domega
            !vito - promjena if if if u else if
            if (jo /= io) then 
              fact=domega/(oi-oj)
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
            ! else
              ! print *,  "WARNING jo loop condition not satisfied."
            end if
            W1 = W1 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            
            if (jo == 1 .or. jo == no) then
              fact = 0.5*domega/(oi+oj)
            end if
            W2 = W2 + fact*S0(jo,iG,jG)
          
          end do
        end if
        
        ImW = -pi*S0(io,iG,jG)
        ! stvari vezane uz GW...
        Gammap(iG,jG) = cmplx(W1,ImW)
        Gammam(iG,jG) = cmplx(-W2,0.0)
        if (iG == 1 .and. jG == 1) then
          W2KK = W2
          ! neven debug
          ! print *, 'W2KK<--- W2:',W2
        end if
        if (io == 1) then
          G0 = Gammap(1,1)
          ! neven debug
          ! print *,'ImW(io=1,iG,jG):',ImW
          ! print *, 'G0<--- Gammap(1,1):',Gammap(1,1)
        end if
        

!                kraj po iG,jG
      end do
    end do
    
!                Provjera KK relacija
    Wind = real(WT(io,1,1)-V(1,1))
    WindKK = real(Gammap(1,1)) - W2KK
    ! neven debug
    ! print *, 'WT: ',WT(io,1,1),' V(11): ',V(1,1) 
    ! print *, 'Wind: ',Wind
    ! print *, 'Gammap(1,1): ',Gammap(1,1),'W2KK: ',W2KK
    ! print *,'WindKK: ',WindKK
    fact = domega
    if (io == 1 .or. io == no-1) then
      fact = 0.5*domega
    end if
    KKS = KKS + fact*(WindKK-Wind)*(WindKK-Wind)
    SKK = SKK + fact*Wind*Wind
    
    ! neven debug
    ! print *,'KKS',KKS
    ! print *, 'SKK',SKK
    
!                kraj nove petlje po omega
  end do omega_loop_B
  ! vito greska?? ovaj fajl je vec zatvoren
  ! close(74)
  

  ! neven debug openmp maknutno
  
  dato = 'Kramers-Kron_Qi'
  nord = INDEX(dato,'i', back =.false.)
  if (iq < 10) then
    write(dato(nord:nord),'(i1) ') iq
  else if (iq >= 10 .and. iq < 100) then
    write(dato(nord:nord+1),'(i2) ') iq
  else
    write(dato(nord:nord+2),'(i3) ') iq
  end if
  
  !$omp master
  open(33+ik,FILE=dato)
  write(33+ik,'(a25,3f10.4,a5) ') 'Wave vector (qx,qy,qz)=(',qx*Gcar,qy*Gcar, qz*Gcar,') a.u.'
  write(33+ik,'(a25,f8.4,a5) ') '|(qx,qy,qz)|=',absq*Gcar,'a.u.'
  write(33+ik,*) 'int(WindKK-Wind)^2 =  ',KKS
  write(33+ik,*) 'int(Wind)^2 =  ',SKK
  write(33+ik,*) '****************************************'
  write(33+ik,*) 'Kramers–Kronig relation relative error'
  write(33+ik,'(5X,f7.2,a2) ') 100.0*abs(KKS/SKK), '%'
  ! 78               FORMAT(a23,f10.5)
  ! 79               FORMAT(a16,f10.5)
  ! 80               FORMAT(5X,f7.2,a2)
  write(33+ik,*) 'Usporedba Gamma i WT'
  write(33+ik,*) 'real[Gamma(o=0,1,1)]=',real(G0)
  write(33+ik,*) 'real[WT(o=0,1,1)]/2=', real(WT(1,1,1)-V(1,1))/2.0
  close(33+ik)
  !$omp end master
  
  
! kraj po q
end do

! MKL matrix inversion vars
deallocate(ipiv)
deallocate(work)

end program surface_loss
