PROGRAM surface_loss
 
!---------------------------------------------------------------------------
!        PROGRAM FOR ab initio-SURFACE LOSS CALCULATION FOR LAYERED SYSTEMS
!        USING SUPERCELL METHOD
!---------------------------------------------------------------------------
!        Quantum Esspresso:
!        verbosity           = 'high'
!       WARNING: VALID ONLY FOR NORM-CONSERVING PSEUDOPOTENTIALS
!---------------------------------------------------------------------------


!USE ISO_Fortran_env ! precision kind, fortran 2008
USE iotk_module
USE ModPointR
implicit none

CHARACTER (LEN=iotk_attlenx) :: attr
LOGICAL :: found

! single / double precision
!INTEGER, PARAMETER :: sp = real32
!INTEGER, PARAMETER :: dp = real64
!INTEGER, PARAMETER :: qp = real128

! INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(p=6)
! INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(p=15)
! INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(p=33)

! NOT used, read directly from QE
! integer :: nMPx,nMPy,nMPz
! parameter(nMPx=201,nMPy=201,nMPz=1)

INTEGER :: ik,i,j,jk,it,lk,Ntot,iG0,Nsymm,iq, &
           io,n,m,iG,R1,K1,R2,K2, Nlf,NG1,   &
           NG2,iG1,iG2,jG,kG,jo,jump,loss,   &
           iGfast,ikmin,lf,kG1,kG2,nord

INTEGER, PARAMETER :: NkI = 6835 ! number of wave vectors in IBZ
INTEGER, PARAMETER :: Nband = 60 ! number of bands
INTEGER, PARAMETER :: NelQE = 18 ! Number of electrons(unit cell)
INTEGER, PARAMETER :: Nk = 48*NkI ! number of wave vectors in FBZ with no symmetry 
INTEGER, PARAMETER :: NGd = 4000 ! number of coefficients CG shulod be less than minimum number of coefficients all over all evc.n files ... moglo bi se dinamicki alocirati 
INTEGER, PARAMETER :: NG = 44121! zasto 8000 ?? total number of G vectors  
INTEGER, PARAMETER :: no = 2001 ! broj frekvencija
INTEGER, PARAMETER :: nq = 2 ! broj valnih vektora tu je 2 jer je rucno paralelizirano!
INTEGER, PARAMETER :: Nlfd = 50 ! dimenzija polja za local field zasto prozivoljno 50, ne moze se znati unaprijed

! file i/o debug
INTEGER :: ist,ist2,ist9,ist10,ist11,ist12
INTEGER :: lno,lno2,lno9,lno10,lno11,lno12

! scalars
REAL(kind=sp) :: kx,ky,kz
REAL(kind=sp) :: KQx,KQy,KQz
REAL(kind=sp) :: qgx,qgy,qgz
REAL(kind=sp) :: omin,omax
REAL(kind=sp) :: qx,qy,qz
REAL(kind=sp) :: kmin
REAL(kind=sp) :: domega
REAL(kind=sp) :: o
REAL(kind=sp) :: K11,K22,K33
REAL(kind=sp) :: Lor ! lorentzian
REAL(kind=sp) :: De 
REAL(kind=sp) :: Gabs ! 
REAL(kind=sp) :: kref ! trazi najmanju k-tocku sampliranu u MP meshu u kojem se moze izracunati ILS
REAL(kind=sp) :: Eref 
REAL(kind=sp) :: Gxx1,Gyy1,Gzz1
REAL(kind=sp) :: Gxx2,Gyy2,Gzz2
REAL(kind=sp) :: fact
REAL(kind=sp) :: oi,oj
REAL(kind=sp) :: ImChi0, ReChi0
! REAL(kind=sp) :: q ! ne koristi se
! REAL(kind=sp) :: qmax ! ne koristi se
REAL(kind=sp) :: Nel ! Number of electrons(1BZ integration)
REAL(kind=sp) :: absq
REAL(kind=sp) :: error
REAL(kind=sp) :: W1,W2
REAL(kind=sp) :: ImW
REAL(kind=sp) :: Wind
REAL(kind=sp) :: W2KK
REAL(kind=sp) :: KKS = 0.0
REAL(kind=sp) :: SKK = 0.0
REAL(kind=sp) :: WindKK
REAL(kind=sp) :: krefM

REAL(kind=dp), PARAMETER :: pi = 4.D0*ATAN(1.D0)
REAL(kind=dp), PARAMETER :: eV = 1.602176487D-19
REAL(kind=dp), PARAMETER :: Hartree = 2.0D0*13.6056923D0
REAL(kind=dp), PARAMETER :: Planck = 6.626196D-34
REAL(kind=dp), PARAMETER :: three = 3.0d0 

REAL(kind=sp), PARAMETER :: Efermi = 0.5554/Hartree ! Fermijeva en. u eV
REAL(kind=sp), PARAMETER :: a0 = 5.9715 ! unit cell parameter in parallel direction in a.u.  
REAL(kind=sp), PARAMETER :: c0 = 29.8575
REAL(kind=sp), PARAMETER :: Gcar = 2.0*pi/a0 ! unit cell norm.
REAL(kind=sp), PARAMETER :: eps = 1.0D-4 ! threshold
REAL(kind=sp), PARAMETER :: T = 0.01/Hartree ! temperature in eV
REAL(kind=sp), PARAMETER :: eta = 0.05/Hartree ! damping i\eta
REAL(kind=sp), PARAMETER :: Ecut = 0.0 ! cutoff energy for crystal local field calculations  
REAL(kind=sp), PARAMETER :: Vcell = 922.0586 ! unit-cell volume in a.u.^3 
REAL(kind=sp), PARAMETER :: aBohr = 0.5291772D0 ! unit cell parameter in perpendicular direction in a.u. (z-separation between supercells)   



COMPLEX(kind=dp), PARAMETER :: ione = CMPLX(1.0,0.0)
COMPLEX(kind=dp), PARAMETER :: czero = CMPLX(0.0,0.0)
COMPLEX(kind=dp), PARAMETER :: rone = CMPLX(0.0,1.0)


COMPLEX(kind=sp) :: G0

! scalar arrays
INTEGER, DIMENSION(Nlfd*NGd) :: Gfast ! pomocna funkcija
INTEGER, DIMENSION(3) :: Gi ! pomocna funkcija
INTEGER, DIMENSION(NG) :: parG ! paritet svakog valnog vektora

! multidim arrays
REAL(kind=sp), DIMENSION(3,NkI) :: kI
REAL(kind=sp), DIMENSION(NkI,Nband) :: E ! vl. vr. danog k-i i band-i
REAL(kind=sp), DIMENSION(48,3,3) :: R ! matr. simetrijskih op.
REAL(kind=sp), DIMENSION(48,3,3) :: RI ! inverz od R
REAL(kind=sp), DIMENSION(3,Nk) :: k
REAL(kind=sp), DIMENSION(3,NK) :: ktot ! ukupno jedinstvenih k-tocaka u FBZ
REAL(kind=sp), DIMENSION(3,NG) :: G ! polje valnih vektora G u recp. prost. za wfn.
REAL(kind=sp), DIMENSION(Nlfd,Nlfd) :: V ! matr. gole coulomb. int.
REAL(kind=sp), DIMENSION(no,Nlfd,Nlfd) :: S0 ! korelacijska matrica
REAL(kind=sp), DIMENSION(3,3) :: KC ! pomocna funkcija
REAL(kind=sp), DIMENSION(3,Nlfd) :: Glf ! local field effect polje valnih vekt. u rec. prost.
REAL(kind=sp), DIMENSION(3,Nlfd) :: GlfV ! local field effect za nescreenanu (golu) int. V


! Nldf dimenziju ne znamo a priori zapravo, trebalo bi staviti sve te matrice allocatable i 
! naknadno ih alocirati

COMPLEX(kind=dp), DIMENSION(Nlfd,Nlfd) :: Imat ! jedinicna matr.
COMPLEX(kind=dp), DIMENSION(Nlfd,Nlfd) :: diel_epsilon ! Epsilon (GG')  = I - V(GG')Chi0
COMPLEX(kind=sp), DIMENSION(Nlfd,Nlfd) :: Chi ! (eq. 2.88 nakon invertiranja) ;oprez bio je double precision

COMPLEX(kind=sp), DIMENSION(Nlfd) :: MnmK1K2 ! nabojni vrhovi
COMPLEX(kind=sp), DIMENSION(Nlfd,Nlfd) :: Chi0 ! (eq. 2.89)
COMPLEX(kind=sp), DIMENSION(no,Nlfd,Nlfd) :: WT ! time ordered RPA screened coulomb int. (eq. 2.93)
COMPLEX(kind=sp), DIMENSION(Nlfd,Nlfd) :: Gammap ! omega>0 ,eq....(skripta 5) \sum_{q,m} \int \dd omega' S(\omega')/{(\omega-\omega'-e_{k+q,m} +i\eta}) za GW se koristi se za ovaj dio 
COMPLEX(kind=sp), DIMENSION(Nlfd,Nlfd) :: Gammam ! omega<0

CHARACTER (LEN=100) :: bandn,bandm,nis,pathk1,pathk2,dato,root,outdir,path,fajl
CHARACTER (LEN=35) :: tag,buffer

COMPLEX(kind=sp), POINTER, DIMENSION(:) :: C1,C2 ! Fourierovi koef. u razvoju wfn. iz QE 




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



!           CALL FOR POINT GROUP TRANSFORMATIONS
!           Point group transformations are in Cartesian coordinate

CALL PointR(root,Nsymm,R,RI)
PRINT *,"PointR done."


!           Upis valnih vektora iz irreducibilne Brillouinove
!           zone i pripadnih energijskih nivoa iz filea '****.band'.
!           wave vectors are in Cartesian coordinate



fajl='/MoS2.band'
path=trim(root)//trim(fajl)
OPEN(40,FILE=path,status='old',err=400,iostat=ist9)

DO  ik = 1,NkI
  IF(ik == 1) THEN
    READ(40,*) 
  END IF
  READ(40,'(10X,3F10.3)') kI(1,ik),kI(2,ik),kI(3,ik)
  READ(40,'(10F8.4)') (E(ik,i),i=1,Nband)
END DO
CLOSE(40)

! 10          FORMAT(10F8.4)
! 20          FORMAT(10X,f10.3,f10.3,f10.3)
GO TO 500
400 write(*,*) '100 cannot open file. iostat = ',ist9
500 CONTINUE

DO  ik=1,NkI
  DO  i=1,Nband
    E(ik,i) = E(ik,i)/Hartree
    IF(i >= 10) THEN ! scissor operator, ispravlja/shifta DFT gap (na 1eV u ovom slucaju)
      E(ik,i) = E(ik,i) + 1.0/Hartree
    END IF
  END DO
END DO





!            generator 1.B.Z.
!            Dio programa koji pomocu operacija tockaste grupe i vektora iz
!            I.B.Z. generira sve (MEDJUSOBNO RAZLICITE!!!) V. vektore u 1.B.Z.
!            Ntot-Tot number of different points ''ktot'' inside 1.B.Z



jk = 0
Ntot = 0
DO  i = 1,Nsymm ! loop over No. symmetries
  DO  ik = 1,NkI  ! loop over k points in IBZ
    it = 1
    jk = jk+1
    DO  n = 1,3 ! loop over kx,ky,kz 
      k(n,jk) = 0.0
      DO  m = 1,3 ! loop over x,y,z
        k(n,jk) = k(n,jk) + R(i,n,m)*kI(m,ik) ! kreira nove k tocke u BZ pomocu simetrije
      END DO
    END DO
    IF(jk > 1) THEN
      DO  lk=1,jk-1
        ! vito - maknut nest if if if
        ! IF(ABS(k(1,jk)-k(1,lk)) <= eps) THEN
        !   IF(ABS(k(2,jk)-k(2,lk)) <= eps) THEN
        !     IF(ABS(k(3,jk)-k(3,lk)) <= eps) THEN
        !       it=2
        !     END IF
        !   END IF
        ! END IF
        IF( ABS(k(1,jk)-k(1,lk)) <= eps .AND. &
            ABS(k(2,jk)-k(2,lk)) <= eps .AND. &
            ABS(k(3,jk)-k(3,lk)) <= eps ) THEN ! je li razlicita tocka od neke prije vec kreirane
          it=2
        END IF
      END DO
    END IF
    IF(it == 1) THEN ! ne postoji dodaj ju
      Ntot=Ntot+1
      ktot(1,Ntot)=k(1,jk)
      ktot(2,Ntot)=k(2,jk)
      ktot(3,Ntot)=k(3,jk)
    END IF
  END DO
END DO


!             Checking 1BZ integration
! NEVEN PROVJERI LOGIKU JOS JEDNOM
Nel = 0 ! provjeri je li broj el. u FBZ odgovara stvarnom broju el. u jed. cel. NelQE
k_loop_FBZ : DO  ik = 1,Ntot
  kx = ktot(1,ik)
  ky = ktot(2,ik)
  kz = ktot(3,ik)
  band_loop: DO  n = 1,Nband
    IF(n == 1) THEN
        it = 1
!      END IF
      IF(ik <= NkI) THEN
        K1 = ik
        it = 2
        ! dodano sa linije goto 5022
        ! IF(E(K1,n) < Efermi) THEN 
        !   Nel = Nel + 1.0
        ! END IF
      ELSE
        symm_loop: DO  i = 2,Nsymm
          K11 = RI(i,1,1)*kx + RI(i,1,2)*ky + RI(i,1,3)*kz
          K22 = RI(i,2,1)*kx + RI(i,2,2)*ky + RI(i,2,3)*kz
          K33 = RI(i,3,1)*kx + RI(i,3,2)*ky + RI(i,3,3)*kz
          k_loop_IBZ: DO  j = 1,NkI
            ! vito - if if if loop
            ! IF(ABS(K11-kI(1,j)) <= eps) THEN
            !   IF(ABS(K22-kI(2,j)) <= eps) THEN
            !     IF(ABS(K33-kI(3,j)) <= eps) THEN
            IF ( ABS(K11-kI(1,j)) <= eps .AND. &
                 ABS(K22-kI(2,j)) <= eps .AND. &
                 ABS(K33-kI(3,j)) <= eps ) THEN
              PRINT *,'FOUND IBZ k-vec:',j
              it = 2
              K1 = j
              ! GO TO 5022
              ! dodano sa linije goto 5022
              IF(E(K1,n) < Efermi) THEN 
                Nel = Nel + 1.0
                PRINT *,Nel
              END IF
              CYCLE band_loop
            END IF
          END DO k_loop_IBZ
        END DO symm_loop
      END IF
        IF(it == 1) THEN
          PRINT*,'Can not find wave vector K=',ik, 'in I.B.Z.'
          STOP
        END IF
        ! 5022          CONTINUE
    END IF
    ! vito - vjerojatno nepotrebno , te
    ! IF(E(K1,n) < Efermi) THEN 
    !   Nel = Nel + 1.0
    ! END IF
  END DO band_loop
END DO k_loop_FBZ
Nel = 2.0*Nel / Ntot ! zbroji za en. manje od fermijeve

STOP

OPEN(887,FILE='fbz_check.dat',status='new')
DO  i = 1,Ntot
  WRITE(887,*) ktot(1,i),ktot(2,i)  ! output da vidimo kako izgleda FBZ
END DO
CLOSE(887)



! KC transformation matrix from rec.cryst. axes to cart.koord.
! If G' is vector in rec.cryst. axes then a=KC*a' is vector in cart. axes


fajl='/MoS2.sc.out'
path=TRIM(root)//TRIM(fajl)
tag='     reciprocal axes: (cart. coord.'
OPEN(30,FILE=path,status='old')
DO  i = 1,100000
  READ(30,'(a)') buffer
  lno=lno+1
  IF(buffer == tag) THEN
    DO  j=1,3
      READ(30,'(23X,3F10.3)',err=10001,iostat=ist,end=20001) KC(1,j), KC(2,j), KC(3,j)
    END DO
    EXIT
  END IF
END DO
! 70          FORMAT(23X,3F10.3)
10001   write(*,*)'Error reading line ',lno+1,', iostat = ',ist
20001   write(*,*)'Number of lines read = ',lno
CLOSE(30)


! Reading the reciprocal vectors in crystal coordinates and transformation
! in Cartesian cordinates.
OPEN(20,FILE='gvectors.xml',status='old',err=200,iostat=ist10)
lno10 =0
PRINT *, 'File gvectors.dat oppened successfully.'
DO  i=1,8
  READ(20,*,err=201,iostat=ist11,end=202) nis
  PRINT *,nis
  lno10 = lno10 +1
END DO
G = 0.0 ! vito - premjesteno iz n,m loopa
DO  iG = 1,NG
  READ(20,'(i10,i11,i11)',err=201,iostat=ist11,end=202) Gi(1),Gi(2),Gi(3)
  lno10 = lno10 +1
  PRINT *, lno10
  IF(iG == 1) THEN
    IF(Gi(1) /= 0 .OR. Gi(2) /= 0 .OR. Gi(3) /= 0) THEN
      PRINT*,'*********************************'
      PRINT*,'WARRNING!, G vectors input is wrong!!'
      PRINT*,'G(1) is not (0,0,0)!!'
      STOP
    END IF
  END IF
! transformation in cart.coord (also!, after this all G components are in 2pi/a0 units)
  DO n = 1,3
    ! G(n,iG) = 0.0
    DO m = 1,3
      G(n,iG) = G(n,iG)+KC(n,m)*REAL(Gi(m)) ! DBLE converted to REAL
    END DO
  END DO
  parG(iG)=Gi(3)
END DO
! 100         FORMAT(i10,i11,i11)
CLOSE(20)

GO TO 5000
200 write(*,*) 'error cant read file id 20, ist=',ist10
201   write(*,*)'201 buffer1 read. Error reading line ',lno10+1,', iostat = ',ist11
202   write(*,*)'202 buffer1 read. Number of lines read = ',lno10
5000 CONTINUE


! Reciprocal vectors for crystal local field effects calculations in array ''Glf(3,Nlf)''

Nlf = 0
IF(lf == 1) THEN
  DO iG = 1,NG
    ! local field efekti samo u okomitom smjeru
    IF(G(1,iG) == 0.0 .AND. G(2,iG) == 0.0) THEN
      Eref = Gcar**2 * G(3,iG)**2 / 2.0
      IF(Eref <= Ecut) THEN
        Nlf = Nlf+1
        Glf(1,Nlf) = 0.0
        Glf(2,Nlf) = 0.0
        Glf(3,Nlf) = G(3,iG)
        IF( (parG(iG)/2)*2 == parG(iG) ) THEN
          parG(Nlf) = 1
        ELSE
          parG(Nlf) = -1
        END IF
      END IF
    END IF
  END DO
ELSE
  DO  iG = 1,NG
    ! local field efekti samo u svim smjerovima
    Eref = Gcar**2 *( G(1,iG)**2 + G(2,iG)**2 + G(3,iG)**2 ) / 2.0
    IF(Eref <= Ecut) THEN
      Nlf = Nlf+1
      Glf(:,Nlf) = G(:,iG)
      ! Glf(1,Nlf) = G(1,iG)
      ! Glf(2,Nlf) = G(2,iG)
      ! Glf(3,Nlf) = G(3,iG)
      IF( (parG(iG)/2)*2 == parG(iG) ) THEN
        parG(Nlf) = 1
      ELSE
        parG(Nlf) = -1
      END IF
    END IF
  END DO
END IF
IF(Nlf > Nlfd) THEN
  PRINT*,'Nlf is bigger than Nlfd'
  STOP
END IF




!             IBZ q LOOP STARTS HERE!!!
! iq=0 ne moze biti nula, opticki racun
! iq=2 do iq=...cutoff transfer q vektor!
! ikmin = min. valni vektor u BZ svi veci su visekratnici tog minimalnog
DO  iq = 42,61
  
!             searching min. q=(qx,qy,qz) in GM direction
  kmin = 1.0
  ntot_loop: DO  i = 1,Ntot ! loop over different k-points in FBZ
    kref = SQRT(ktot(1,i)*ktot(1,i) + ktot(2,i)*ktot(2,i) + ktot(3,i)*ktot(3,i))
    IF(kref == 0.0) THEN
      EXIT Ntot_loop
    ELSE IF(kref < kmin) THEN
      kmin = kref
      ikmin = i
      krefM = kmin
    END IF
  END DO ntot_loop
  
  
  qx = (iq-1)*ktot(1,ikmin)
  qy = (iq-1)*ktot(2,ikmin)
  qz = (iq-1)*ktot(3,ikmin)
  
  absq = SQRT(qx*qx + qy*qy + qz*qz)
  
!             Info file
  
  OPEN(55,FILE='Info')
  WRITE(55,*)'***************General***********************'
  WRITE(55,*)''
  WRITE(55,*)'Number of point symmetry operation is',Nsymm
  WRITE(55,'(a25,3f10.4,a5)')'Wave vector (qx,qy,qz)=(',qx*Gcar,qy*Gcar, qz*Gcar,') a.u.'
  WRITE(55,'(a25,f8.4,a5)')'|(qx,qy,qz)|=',absq*Gcar,'a.u.'
  IF(lf == 1)WRITE(55,*)'Local field effcts in z-dir'
  IF(lf == 3)WRITE(55,*)'Local field in all xyz-dir'
  WRITE(55,*)'Number of local field vectors is',Nlf
  WRITE(55,*)'Number of different K vectors in 1.B.Z. is',Ntot
  WRITE(55,*)'Number of K vectors in I.B.Z. is',NkI
  WRITE(55,*)'Number of bands is               ',Nband
  WRITE(55,'(a25,f8.4,a5)')'Eta damping is ',eta*Hartree*1000.0,'meV'
  WRITE(55,'(a25,f8.4,a5)')'Temperature is  ',T*Hartree*1000.0,'meV'
  WRITE(55,*)''
  WRITE(55,*)'************* Checking 1BZ integration*******'
  WRITE(55,*)''
  WRITE(55,'(a40,f7.4)')'Number of electrons(1BZ integration)=',Nel
  WRITE(55,*)'Number of electrons(unit cell)=',NelQE
  error=ABS((NelQE-Nel)/NelQE)
  WRITE(55,'(a25,f8.4,a5)')'Relative error=',error*100.0,'%'
  IF(error > 0.05) THEN
    WRITE(55,*)'WARRNING!!-1BZ INTEGRATION IS BAD!.'
  END IF
  CLOSE(55)
  ! 88         FORMAT(a25,3f10.4,a5)
  ! 99         FORMAT(a25,f8.4,a5)
  ! 12         FORMAT(a40,f7.4)
  
  
  DO  io = 1,no
    DO  iG = 1,Nlf
      DO  jG=1,Nlf
        S0(io,iG,jG) = czero
      END DO
    END DO
  END DO
  
  
!              1.B.Z  LOOP STARTS HERE !!!!
  
  DO  ik=1,Ntot
    
    OPEN(122,FILE='status')
    WRITE(122,*) 'iq=',iq
    WRITE(122,*) 'ik=',ik
    CLOSE(122)
    
    
    kx = ktot(1,ik)
    ky = ktot(2,ik)
    kz = ktot(3,ik)
    
!  trazenje (kx,ky,kz) u ireducibilnoj zoni
    
    it = 1
    IF(ik <= NkI) THEN
      R1=1
      K1=ik
      it=2
    ELSE
      symmetry_loop: DO  i = 2,Nsymm
        K11 = RI(i,1,1)*kx + RI(i,1,2)*ky+RI(i,1,3)*kz
        K22 = RI(i,2,1)*kx + RI(i,2,2)*ky+RI(i,2,3)*kz
        K33 = RI(i,3,1)*kx + RI(i,3,2)*ky+RI(i,3,3)*kz
        DO  j = 1,NkI
          ! vito - smnanjenje if uvjeta
          ! IF(ABS(K11-kI(1,j)) <= eps) THEN
          !   IF(ABS(K22-kI(2,j)) <= eps) THEN
          !     IF(ABS(K33-kI(3,j)) <= eps) THEN
          !       it=2
          !       R1=i
          !       K1=j
          !       ! GO TO 5222
          !       EXIT symmetry_loop
          !     END IF
          !   END IF
          ! END IF
          IF(       ABS(K11-kI(1,j)) <= eps &
              .and. ABS(K22-kI(2,j)) <= eps &
              .and. ABS(K33-kI(3,j)) <= eps ) THEN
            it=2
            R1=i
            K1=j
            EXIT symmetry_loop
          END IF
        END DO
      END DO symmetry_loop
    END IF
    IF(it == 1) THEN
      PRINT*,'Can not find wave vector K=',ik, 'in I.B.Z.'
      STOP
    END IF
    ! 5222           CONTINUE
    
    
    it = 1
    KQx = kx + qx
    KQy = ky + qy
    KQz = kz + qz
    
!   trazenje (KQx,KQy) prvo u 1.B.Z a onda u I.B.Z.
    iG_loop: DO  iG = 1,NG
      DO  jk = 1,Ntot
        ! vito smanjenje If if if loop
        IF(ABS(KQx-G(1,iG)-ktot(1,jk)) <= eps .and. &
           ABS(KQy-G(2,iG)-ktot(2,jk)) <= eps .and. &
           ABS(KQz-G(3,iG)-ktot(3,jk)) <= eps ) THEN
          it=2
          iG0 = iG
          DO  i = 1,Nsymm
            K11 = RI(i,1,1)*ktot(1,jk)+RI(i,1,2)*ktot(2,jk)+  &
                RI(i,1,3)*ktot(3,jk)
            K22 = RI(i,2,1)*ktot(1,jk)+RI(i,2,2)*ktot(2,jk)+  &
                RI(i,2,3)*ktot(3,jk)
            K33 = RI(i,3,1)*ktot(1,jk)+RI(i,3,2)*ktot(2,jk)+  &
                RI(i,3,3)*ktot(3,jk)
            DO  j = 1,NkI
              IF(ABS(K11-kI(1,j)) <= eps .and. &
                 ABS(K22-kI(2,j)) <= eps .and. &
                 ABS(K33-kI(3,j)) <= eps ) THEN
                it = 3
                R2 = i
                K2 = j
                EXIT iG_loop
                ! GO TO 2111
              END IF
            END DO
          END DO
        END IF
      END DO
    END DO iG_loop  

    ! vito prepravljeno maknuti ifovi
    ! ig_loop: DO  iG = 1,NG
    !   DO  jk = 1,Ntot
    !     ! vito smanjenje If if if loop
    !     IF(ABS(KQx-G(1,iG)-ktot(1,jk)) <= eps) THEN
    !       IF(ABS(KQy-G(2,iG)-ktot(2,jk)) <= eps) THEN
    !         IF(ABS(KQz-G(3,iG)-ktot(3,jk)) <= eps) THEN
    !           it=2
    !           iG0 = iG
    !           DO  i = 1,Nsymm
    !             K11=RI(i,1,1)*ktot(1,jk)+RI(i,1,2)*ktot(2,jk)+  &
    !                 RI(i,1,3)*ktot(3,jk)
    !             K22=RI(i,2,1)*ktot(1,jk)+RI(i,2,2)*ktot(2,jk)+  &
    !                 RI(i,2,3)*ktot(3,jk)
    !             K33=RI(i,3,1)*ktot(1,jk)+RI(i,3,2)*ktot(2,jk)+  &
    !                 RI(i,3,3)*ktot(3,jk)
    !             DO  j=1,NkI
    !               IF(ABS(K11-kI(1,j)) <= eps) THEN
    !                 IF(ABS(K22-kI(2,j)) <= eps) THEN
    !                   IF(ABS(K33-kI(3,j)) <= eps) THEN
    !                     it = 3
    !                     R2 = i
    !                     K2 = j
    !                     EXIT ig_loop
    !                     ! GO TO 2111
    !                   END IF
    !                 END IF
    !               END IF
    !             END DO
    !           END DO
    !         END IF
    !       END IF
    !     END IF
    !   END DO
    ! END DO ig_loop
    ! 2111                 CONTINUE
    
    
    IF(it == 1) THEN
      PRINT*,'Can not find wave vector K+Q=',ik,'+',iq, 'in 1.B.Z.'
      STOP
    ELSE IF(it == 2) THEN
      PRINT*,'Can not find wave vector K+Q=',ik,'+',iq, 'in I.B.Z.'
      STOP
    END IF
    
    
!              R1-integer, redni broj point operacije R1 u transformaciji ''K=R1*K1''.
!              K1-integer, redni broj valnog vektora K1 u transformaciji ''K=R1*K1''.
!              iG0 i R2-integeri, redni broj vektora reciprocne restke G0 i point operacije R2 u transformaciji ''K+Q=G0+R2*K2''.
!              K2-integer, redni broj valnog vektora K2 u transformaciji  ''K+Q=G0+R2*K2''.
    
    
!              petlje po vrpcama n i m
    
    DO  n = 1,9 ! filled bands loop
      DO  m = 10,Nband ! empty bands loop
        
        
        
        
        ! otvara save/K.000x/evc.dat u atributu <evc band> ispod CnK(G) koef.
        
        CALL paths(outdir,K1,K2,n,m,pathk1,pathk2,bandn,bandm) 
        
!         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za
!         fiksni K1,K2,n i m
        
!               Otvaranje atribute za INFO
        CALL iotk_open_read(10,pathk1)
        CALL iotk_scan_empty(10,"INFO",attr=attr)
        CALL iotk_scan_attr(attr,"igwx",NG1)
!               Alociranje polja C1
        allocate (C1(NG1))
!               Ucitavanje podataka iza evc.n
        CALL iotk_scan_dat(10,bandn,C1)
        CALL iotk_close_read(10)
!               Otvaranje atribute za INFO
        CALL iotk_open_read(10,pathk2)
        CALL iotk_scan_empty(10,"INFO",attr=attr)
        CALL iotk_scan_attr(attr,"igwx",NG2)
!               Alociranje polja C2
        allocate (C2(NG2))
!               Ucitavanje podataka iza evc.m
        CALL iotk_scan_dat(10,bandm,C2)
        CALL iotk_close_read(10)
        
!                Konstrukcija stupca matricnih elementa MnmK1K2(G)
        
        
        IF(NGd > NG1) THEN
          WRITE(*,*)'NGd is bigger than NG1=',NG1
          STOP
        ELSE IF(NGd > NG2) THEN
          WRITE(*,*)'NGd is bigger than NG2=',NG2
          STOP
        END IF
        
        
!       matrix elements
        iGfast = 0
        MnmK1K2 = czero ! nabojni vrhovi
        DO  iG = 1,Nlf ! suma po lokalnim fieldovima kojih ima Nlf
          ! MnmK1K2(iG) = czero
          DO  iG1 = 1,NGd
            iGfast = iGfast+1
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
            IF(jump == 1) THEN
              iG2_loop: DO  iG2 = 1,NG2
                Gfast(iGfast) = NG2+1
                Gxx2 = G(1,iG2)
                Gyy2 = G(2,iG2)
                Gzz2 = G(3,iG2)
                ! vito if if if loop smanjen
                ! IF(ABS(Gxx2-Gxx1) < eps) THEN
                !   IF(ABS(Gyy2-Gyy1) < eps) THEN
                !     IF(ABS(Gzz2-Gzz1) < eps) THEN
                !       Gfast(iGfast) = iG2
                !       ! GO TO 1111
                !       EXIT ig2_loop
                !     END IF
                !   END IF
                ! END IF
                IF( ABS(Gxx2-Gxx1) < eps .AND. &
                    ABS(Gyy2-Gyy1) < eps .AND. &
                    ABS(Gzz2-Gzz1) < eps) THEN
                    Gfast(iGfast) = iG2
                    ! GO TO 1111
                    EXIT iG2_loop
                END IF
              END DO iG2_loop
            END IF
            ! 1111              CONTINUE
            iG2 = Gfast(iGfast)
            ! ovo bi trebalo prouciti zasto ovako radi
            IF(iG2 <= NG2) THEN
              MnmK1K2(iG) = MnmK1K2(iG) + CONJG(C1(iG1))*C2(iG2)
            END IF
          END DO
        END DO
        jump = 2
        
        
!       omega loop
        DO  io = 1,no
          o = (io-1)*domega
          De = o + E(K1,n) - E(K2,m) 
          Lor = -eta/(De**2 + eta**2) ! ovo bi analticki bila delta funkcija imag. dio od 1/De
          IF(ABS(Lor) >= 1.0D-3/eta) THEN
            DO  iG = 1,Nlf
              DO  jG = 1,Nlf
                S0(io,iG,jG) = S0(io,iG,jG) -  &
                    2.0*Lor*MnmK1K2(iG)*CONJG(MnmK1K2(jG)) / (pi*Ntot*Vcell)
              END DO
            END DO
          END IF
          
          
        END DO
        
        deallocate(C1)
        deallocate(C2)  
                
      END DO ! end of m do loop
      

    END DO ! end of n do loop
    jump=1

  END DO !  end of 1.B.Z do loop
  
  
  
! Puting (qx,qy,qz) and Glf in cartesian coordinates
  
  qx = Gcar*qx ! convert qx *2p/a0
  qy = Gcar*qy
  qz = Gcar*qz
  
  GlfV(:,iG) = Gcar*Glf(:,iG)
  ! DO  iG = 1,Nlf
  !   GlfV(1,iG) = Gcar*Glf(1,iG)
  !   GlfV(2,iG) = Gcar*Glf(2,iG)
  !   GlfV(3,iG) = Gcar*Glf(3,iG)
  ! END DO
 

  
! Kramers-Kroning relacije
  
! ew sum over omega
  DO  io = 1,no-1
  ! print*,io
    oi = (io-1)*domega
    DO  iG = 1,Nlf
      DO  jG = 1,Nlf
        ReChi0 = 0.0
!       static limit
        IF(io == 1) THEN
          DO  jo = 2,no
            oj = (jo-1)*domega
            fact = domega/oj
            ! analticki trikovi za integriranje 
            IF(jo == 2) THEN 
              fact = 3.0/2.0
            ELSE IF(jo == no) THEN
              fact = 0.5*domega/oj
            ELSE
              PRINT *,  "WARNING jo loop condition not satisfied."
            END IF
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
          END DO
          ReChi0 = -2.0*ReChi0
        ELSE IF(io == 2) THEN
          DO  jo = 1,no
            oj = (jo-1)*domega
            IF(jo /= io) THEN
              fact = domega/(oi-oj)
            ELSE IF(jo == 1) THEN
              fact = 1.0
            ELSE IF(jo == 2) THEN
              fact = 0.0
            ELSE IF(jo == 3) THEN
              fact = -3.0/2.0
            ELSE IF(jo == no) THEN
              fact = 0.5*domega/(oi-oj)
            ELSE
              PRINT *,  "WARNING jo loop condition not satisfied."
            END IF
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            IF(jo == 1.OR.jo == no) THEN
              fact=0.5*domega/(oi+oj)
            END IF
            ReChi0 = ReChi0 - fact*S0(jo,iG,jG)
          END DO
        ELSE IF(io == (no-1)) THEN
          DO  jo = 1,no
            oj = (jo-1)*domega
            IF(jo /= io) THEN
              fact = domega/(oi-oj)
            ELSE IF(jo == 1) THEN
              fact = 0.5*domega/(oi-oj)
            ELSE IF(jo == (no-2)) THEN
              fact = 3.0/2.0
            ELSE IF(jo == (no-1)) THEN
              fact = 0.0
            ELSE IF(jo == no) THEN
              fact = -1.0
            ELSE 
              PRINT *,  "WARNING jo loop condition not satisfied."
            END IF
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            IF(jo == 1 .OR. jo == no) THEN
              fact = 0.5*domega/(oi+oj)
            END IF
            ReChi0 = ReChi0 - fact*S0(jo,iG,jG)
          END DO
        ELSE
          DO  jo = 1,no
            oj = (jo-1)*domega
            IF(jo /= io) THEN
              fact = domega/(oi-oj)
            ELSE IF(jo == 1) THEN
              fact=0.5*domega/(oi-oj)
            ELSE IF(jo == (io-1)) THEN
              fact=3.0/2.0
            ELSE IF(jo == io) THEN
              fact=0.0
            ELSE IF(jo == (io+1)) THEN
              fact=-3.0/2.0
            ELSE IF(jo == no) THEN
              fact = 0.5*domega/(oi-oj)
            ELSE 
              PRINT *,  "WARNING jo loop condition not satisfied."
            END IF
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            IF(jo == 1 .OR. jo == no) THEN
              fact = 0.5*domega/(oi+oj)
            END IF
            ReChi0 = ReChi0 - fact*S0(jo,iG,jG)
          END DO
        END IF
        
        ImChi0 = -pi*S0(io,iG,jG)
        Chi0(iG,jG) = CMPLX(ReChi0,ImChi0)
!                kraj po iG,jG
      END DO
    END DO
    
    
    
!                Calculation of the ''Chi''  by matrix invertion
    
    
!                MATRIX V(G,G')
    
    DO  iG = 1,Nlf
      Gabs = SQRT( (qx+GlfV(1,iG))*(qx+GlfV(1,iG)) +  &
                   (qy+GlfV(2,iG))*(qy+GlfV(2,iG)) )
      IF(Gabs == 0.0) THEN 
        Gabs = eps
      END IF
      DO  jG = 1,Nlf
        V(iG,jG) = 0.0
        IF(Glf(1,jG) == Glf(1,iG) ) THEN
          IF(Glf(2,jG) == Glf(2,iG) ) THEN
            V(iG,jG) = 4.0*pi*(1.0-EXP(-Gabs*c0)) / (Gabs*c0)
            V(iG,jG) = V(iG,jG)*( Gabs*Gabs - GlfV(3,iG)*GlfV(3,jG) )
            V(iG,jG) = V(iG,jG)/( Gabs*Gabs + GlfV(3,iG)*GlfV(3,iG) )
            V(iG,jG) = V(iG,jG)/( Gabs*Gabs + GlfV(3,jG)*GlfV(3,jG) )
            V(iG,jG) = -REAL(parG(iG))*REAL(parG(jG))*V(iG,jG) ! dble converted to real
            IF(Glf(3,jG) == Glf(3,iG)) THEN
              V(iG,jG) = 4.0*pi / ( Gabs*Gabs + GlfV(3,iG)*GlfV(3,iG) ) + V(iG,jG)
            END IF
          END IF
        END IF
      END DO
    END DO
    
    
    Imat = czero
    DO  iG = 1,Nlf
      Imat(iG,iG) = rone
    END DO
    
    DO  iG = 1,Nlf
      DO  jG = 1,Nlf
        diel_epsilon(iG,jG) = Imat(iG,jG)
        DO  kG = 1,Nlf
          diel_epsilon(iG,jG) = diel_epsilon(iG,jG) - Chi0(iG,kG)*V(kG,jG)
        END DO
      END DO
    END DO
    
    
!  invertiranje matrice ''diel_epsilon = 1-Chi_0*V''
    
    
    CALL gjel(diel_epsilon,Nlf,Nlfd,Imat,Nlf,Nlfd)
    
    Chi = czero
    DO  iG = 1,Nlf
      DO  jG = 1,Nlf
        DO  kG = 1,Nlf
          Chi(iG,jG) = Chi(iG,jG) + diel_epsilon(iG,kG)*Chi0(kG,jG)
        END DO
      END DO
    END DO
    
    
!  SCREENED COULOMB INTERACTION W^T_GG'(Q,\omega)
    
    DO  iG = 1,Nlf
      DO  jG = 1,Nlf
        WT(io,iG,jG) = czero
        DO  kG1 = 1,Nlf
          DO  kG2 = 1,Nlf
            WT(io,iG,jG) = WT(io,iG,jG) + V(iG,kG1)*Chi(kG1,kG2)*V(kG2,jG)
          END DO
        END DO
        WT(io,iG,jG) = V(iG,jG) + WT(io,iG,jG)
      END DO
    END DO
    
!                kraj nove petlje po omega
  END DO
  
  
! ispis time ordered zasjenjene kulonske interakcije W_GG'^T(Q,\omega)
  dato = 'W_Qi'
  nord = INDEX(dato,'i', back =.false.)
  IF(iq < 10) THEN
    WRITE(dato(nord:nord),'(i1)')iq
  ELSE IF(iq >= 10 .AND. iq < 100) THEN
    WRITE(dato(nord:nord+1),'(i2)')iq
  ELSE
    WRITE(dato(nord:nord+2),'(i3)')iq
  END IF
  
  OPEN(74,FILE=dato)
  DO  io = 1,1
    o = (io-1)*domega
    WRITE(74,*) 'omega=',o,'Hartree'
    WRITE(74,'(10F15.5)')((WT(io,iG,jG),jG=1,Nlf),iG=1,Nlf)
  END DO
  CLOSE(74)
  ! 44              FORMAT(10F15.5)
  
  DO  io = 1,no-1
    DO  iG = 1,Nlf
      DO  jG = 1,Nlf
        S0(io,iG,jG) = -(1.0/pi)*AIMAG(WT(io,iG,jG))
      END DO
    END DO
  END DO
  
  ! podatci za GW sve na dalje

  KKS = 0.0
  SKK = 0.0
  
  
! new sum over omega
  omega_loop: DO  io = 1,no-1
    ! print*,io
    oi = (io-1)*domega
    DO  iG = 1,Nlf
      DO  jG = 1,Nlf
        W1 = 0.0
        W2 = 0.0
!      static limit
        IF(io == 1) THEN

          DO  jo = 2,no
            oj = (jo-1)*domega
            fact = domega/oj
            ! vito - promjena if if if u else if
            IF(jo == 2) THEN 
              fact = 3.0/2.0
            ELSE IF(jo == no) THEN
              fact = 0.5*domega/oj
            ELSE
              PRINT *,  "WARNING jo loop condition not satisfied."
            END IF
            W1 = W1 - fact*S0(jo,iG,jG)
          END DO
          W2 = -W1

        ELSE IF(io == 2) THEN

          DO  jo = 1,no
            oj= (jo-1)*domega
            !vito - promjena if if if u else if
            IF(jo /= io) THEN
              fact = domega/(oi-oj)
            ELSE IF(jo == 1) THEN
              fact = 1.0
            ELSE IF(jo == 2) THEN
              fact = 0.0
            ELSE IF(jo == 3) THEN
              fact = -3.0/2.0
            ELSE IF(jo == no) THEN
              fact = 0.5*domega/(oi-oj)
            ELSE
              PRINT *,  "WARNING jo loop condition not satisfied."
            END IF             
            W1 = W1 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            IF(jo == 1 .OR. jo == no) THEN
              fact = 0.5*domega/(oi+oj)
            END IF
            W2 = W2 + fact*S0(jo,iG,jG)
          END DO

        ELSE IF(io == (no-1)) THEN

          DO  jo=1,no
            oj= (jo-1)*domega
            !vito - promjena if if if u else if
            IF(jo /= io) THEN 
              fact = domega/(oi-oj)
            ELSE IF(jo == 1) THEN
              fact = 0.5*domega / (oi-oj)
            ELSE IF(jo == (no-2)) THEN
              fact = 3.0/2.0
            ELSE IF(jo == (no-1)) THEN
              fact = 0.0
            ELSE IF(jo == no) THEN
              fact = -1.0
            ELSE
              PRINT *,  "WARNING jo loop condition not satisfied."
            END IF
            W1 = W1 + fact*S0(jo,iG,jG)
            fact = domega / (oi+oj)
            
            IF(jo == 1 .OR. jo == no) THEN
              fact = 0.5*domega / (oi+oj)
            END IF
            W2 = W2 + fact*S0(jo,iG,jG)
          
          END DO

        ELSE
          
          DO  jo = 1,no
            oj= (jo-1)*domega
            !vito - promjena if if if u else if
            IF(jo /= io) THEN 
              fact=domega/(oi-oj)
            ELSE IF(jo == 1) THEN
              fact = 0.5*domega/(oi-oj)
            ELSE IF(jo == (io-1)) THEN
              fact = 3.0/2.0
            ELSE IF(jo == io) THEN
              fact = 0.0
            ELSE IF(jo == (io+1)) THEN
              fact = -3.0/2.0
            ELSE IF(jo == no) THEN
              fact = 0.5*domega/(oi-oj)
            ELSE
              PRINT *,  "WARNING jo loop condition not satisfied."
            END IF
            W1 = W1 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            
            IF(jo == 1 .OR. jo == no) THEN
              fact = 0.5*domega/(oi+oj)
            END IF
            W2 = W2 + fact*S0(jo,iG,jG)
          
          END DO
        END IF
        
        ImW = -pi*S0(io,iG,jG)
        ! stvari vezane u GW...
        Gammap(iG,jG) = CMPLX(W1,ImW)
        Gammam(iG,jG) = CMPLX(-W2,0.0)
        IF(iG == 1 .AND. jG == 1) THEN
          W2KK = W2
          IF(io == 1) THEN
            G0 = Gammap(1,1)
          END IF
        END IF

!                kraj po iG,jG
      END DO
    END DO
    
!                Provjera KK relacija
    Wind = REAL(WT(io,1,1)-V(1,1))
    WindKK = REAL(Gammap(1,1)) - W2KK
    fact = domega
    IF(io == 1 .OR. io == no-1) THEN
      fact = 0.5*domega
    END IF
    KKS = KKS + fact*(WindKK-Wind)*(WindKK-Wind)
    SKK = SKK + fact*Wind*Wind
    
    
!                kraj nove petlje po omega
  END DO omega_loop
  CLOSE(74)
  
  
  dato = 'Kramers-Kron_Qi'
  nord = INDEX(dato,'i', back =.false.)
  IF(iq < 10) THEN
    WRITE(dato(nord:nord),'(i1)') iq
  ELSE IF(iq >= 10.AND.iq < 100) THEN
    WRITE(dato(nord:nord+1),'(i2)') iq
  ELSE
    WRITE(dato(nord:nord+2),'(i3)') iq
  END IF
  
  OPEN(33,FILE=dato)
  WRITE(33,'(a25,3f10.4,a5)') 'Wave vector (qx,qy,qz)=(',qx*Gcar,qy*Gcar, qz*Gcar,') a.u.'
  WRITE(33,'(a25,f8.4,a5)') '|(qx,qy,qz)|=',absq*Gcar,'a.u.'
  WRITE(33,*) 'int(WindKK-Wind)^2 =  ',KKS
  WRITE(33,*) 'int(Wind)^2 =  ',SKK
  WRITE(33,*) '****************************************'
  WRITE(33,*) 'Kramers–Kronig relation relative error'
  WRITE(33,'(5X,f7.2,a2)') 100.0*ABS(KKS/SKK), '%'
  ! 78               FORMAT(a23,f10.5)
  ! 79               FORMAT(a16,f10.5)
  ! 80               FORMAT(5X,f7.2,a2)
  WRITE(33,*) 'Usporedba Gamma i WT'
  WRITE(33,*) 'real[Gamma(o=0,1,1)]=',REAL(G0)
  WRITE(33,*) 'real[WT(o=0,1,1)]/2=', REAL(WT(1,1,1)-V(1,1))/2.0
  CLOSE(33)
  
  
! kraj po q
END DO

END PROGRAM surface_loss
