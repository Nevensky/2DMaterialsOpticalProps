PROGRAM surface_loss
 
! Code converted using TO_F90 by Alan Miller
! Date: 2020-06-17  Time: 12:20:53

!        PROGRAM FOR ab initio-SURFACE LOSS CALCULATION FOR LAYERED SYSTEMS
!        USING SUPERCELL METHOD


!        NkI -number of wave vectors in irreducible B. zone
!        Ntot-total number of the mutually different wave vector-program generates this number
!        Nband-number of the bands
!        NG-total number of G vectors
!        NGd-number of coefficients CG shulod me less than minimum number of coefficients all over all evc.n files
!        nMPx*nMPy*nMPz-Monkhorest-Pack sampling
!        Ef-Fermi energy
!        T-temperature in eV
!        nq-number of wave vectors in W calculation
!        Gama-Damping parameter in eV
!        Vcell-unit-cell volume in a.u.^3
!        a0-unit cell parameter in parallel direction in a.u.
!        c0-unit cell parameter in perpendicular direction in a.u. (z-separation between supercells)
!        Ecut-- crystal local field cut-off for W
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Quantum Esspresso:
!        verbosity           = 'high'
!        VALID JUST FOR NORM-CONSERVING PSEUDOPOTENTIALS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




use iotk_module
implicitnone
CHARACTER (LEN=1) :: iotk_attlenx) :: attr
LOGICAL :: :: found

INTEGER :: nki,nband,ik,i,nk,j,jk,it,lk,ntot,nsim,ik1,ik2,  &
    ng,io,no,iq,nq,nmpx,nmpy,nmpz,n,m,ig,r1,k1,r2,k2,ig01,ig02,  &
    nlf,ng1,ng2,ngd,ig1,ig2,jg,nlfd,kg,jo,jump,ni,mi,  &
    MOD,! mod = 1 -> racunaj 4-point polarizabilnost, mod = 2 izracunaj ladder struja-struja odzvinu funkciju i dodaj ju bubble  &
    igfast,ikmin,nelqe,lf,kg1,kg2,nord,io1,io2, nkbz, ! krug oko K i K' u BZ  &
    nkd,! dimenzija Fockov kernela, koliko ima k-ova u FBZ pol, ! polarizacija  &
    spin, ! spin nval, !broj valentnih vrpci bez LS vezanja  &
    nvalls, ! broj valentnih vrpci sa LS vezanjem
PARAMETER(nmpx=51,nmpy=51,nmpz=1,nki=460,nband=40,nelqe=26,  &
    nk=48*nki,ngd=4000,ng=8000,no=401,nq=121,nlfd=50,nkd=2700)


!        skalars
REAL*8 a0,c0,eps,kx,ky,kz,omin,omax,sgn,struja,  &
    qx,qy,qz,kmin,domega,o,ef,k11,k22,k33,t,lor,de,gabs,kref,  &
    eref,ecut,gxx1,gyy1,gzz1,gxx2,gyy2,gzz2,fact,vcell,oi,oj,  &
    q,gcar,abohr,zero,nel,absq,error,o1,o2,gama,  &
    qmax,!maximum transfer wave vector around K point..  ne smiju se preklapati krugovi oko k i k'  &
    kpointx,kpointy,kpointz,  &
    qmin, ! minimum transfer wave vector in WT calculation, q po kojem je sampliran wt(q,omega)  &
    vq, ! 2*pi/Q gw,alpha,epsq,  &
    valley ! faktor 2x jer imamo K i K'
doubleprecision pi,three,hartree,planck,ev,eta
PARAMETER(hartree=2.0D0*13.6056923D0,ef=1.5703/hartree,  &
    a0=5.9715,c0=29.8575,pi=3.141592654D0,gcar=2.0*pi/a0,  &
    eps=1.0D-4,t=0.01/hartree,eta=0.05/hartree,ev=1.602176487D-19,  &
    planck=6.626196D-34,ecut=0.0,vcell=922.0586, abohr=0.5291772D0,zero=0.0)
DOUBLE COMPLEX a,ione,czero,rone,em,
&    g0,w, ! pomocni za G i WT
&    l0i,l0j, ! vezano uz propagator slobodnih el.-supljina para (no scattering)
&    j1,j2 ! strujni vrhovi po polarizaciji i valnom vekt.


!        arrays
INTEGER :: gfast,gi,parg
REAL*8 ki,e,r,ri,k,ktot,g,glf,kc, kbz, ! krugovi u FBZ,  &
    delta ! Delta(ik1)=E(K1,nvalLS-1)-E(K1,nvalLS+1)
COMPLEX*8 mnmk1k2,wt, fockk, ! Fockov Kernel (tj. matr. el. od W)  &
    xi0, ! prepis jedn. lladd, !  4-point polarizabilnosti L^ladd  &
    current, ! strujni vrhovi  &
    chi_ladd, ! zatvaranje ladder dijagram a strujnim vrhovima  &
    chi_para ! Paramagnetska struja-struja polarizabilost  = obicni RPA bubble  (bez w)
DOUBLE COMPLEX fock, ! Fockov Kernel (tj. matr. el. od W) UNIT
DIMENSION ki(3,nki),e(nki,nband),r(48,3,3),ri(48,3,3),  &
    k(3,nk),ktot(3,nk),g(3,ng),glf(3,nlfd),mnmk1k2(2,nlfd),  &
    gfast(nlfd*ngd),kc(3,3),gi(3),wt(nq,nlfd,nlfd),parg(ng),  &
    kbz(3,nk),fockk(nkd,nkd),delta(nkd),fock(nkd,nkd),  &
    UNIT(nkd,nkd),xi0(nkd,nkd),lladd(nkd,nkd),current(2,nkd),  &
    chi_ladd(2),chi_para(2)


CHARACTER (LEN=100) :: bandn,bandm,nis,pathk1,pathk2,dato,  &
    root,path,fajl,root1,root2
CHARACTER (LEN=35) :: tag,buffer

INTEGER :: :: ngw, igwx, nbnd, nk1
COMPLEX(:: kind = 8),pointer, DIMENSION(:) :: c1,c2


rone=DCMPLX(1.0,0.0)
czero=DCMPLX(0.0,0.0)
ione=DCMPLX(0.0,1.0)


! BRAVAIS LATTICE PARAMETERS

!     bravais-lattice index     =            4
!     lattice parameter (alat)  =       5.9715  a.u.
!     unit-cell volume          =     922.0586 (a.u.)^3
!     number of atoms/cell      =            3
!     number of atomic types    =            2
!     *******************************OVO JE BROJ ELEKTRONA KAD IMAS FULL REL. PSEUDOPOTENTIAL
!     number of electrons       =        26.00
!     **************************************************************************************
!     number of Kohn-Sham states=           34
!     kinetic-energy cutoff     =      50.0000  Ry
!     charge density cutoff     =     200.0000  Ry
!     convergence threshold     =      1.0E-06
!     mixing beta               =       0.7000
!     number of iterations used =            8  plain     mixing
!     Exchange-correlation      = PBE ( 1  4  3  4 0 0)
!     Non magnetic calculation with spin-orbit


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
root='/home/vito/PROJECTS/MoS2-BSE/MoS2_51x51/'
!             Staticka zasjenjena kulonska interakcija WT je smjestena u folderu
root1='/home/vito/PROJECTS/MoS2-BSE/W/'



!             Crystal local field effects are included in z direction lf=1
!             Crystal local field effects are included in x,y,z direction lf=3
!             If mod=1 Fock kernel calculation
!             If mod=2 Ladder ireducible current-current polarisability calculation
!             SAMO ZA MOD=2
!             For spin up put spin=1
!             For spin down put spin=2

MOD=2
spin=1
lf=1
jump=1
three=3.0D0
omin=1.0D-5
omax=4.0/hartree ! 0eV do 4eV jer iza ladder term nije dobar
domega=(omax-omin)/(no-1)


!             maximum transfer wave vector around K point..
!             qmax=sqrt(3.0)/3.0
qmax=0.3
IF(qmax <= (1.0/3.0))valley=2.0
IF(qmax > (1.0/3.0))valley=1.0
!             minimum transfer wave vector in WT calculation
qmin=0.006
!             broj valentnih vrpci bez LS vezanja
nval=9
!             broj valentnih vrpci sa LS vezanjem (full-rel. PP)
nvalls=26
!             band gap correction (in eV)
gw=1.0 ! rucna korekcija, scissor operator



!           CALL FOR POINT GROUP TRANSFORMATIONS
!           Point group transformations are in Cartesian coordinate

CALL pointr(root,nsim,r,ri)



!           Upis valnih vektora iz irreducibilne Brillouinove
!           zone i pripadnih energijskih nivoa iz filea '****.band'.
!           wave vectors are in Cartesian coordinate


fajl='MoS2_LS.band'
path=trim(root)//trim(fajl)
OPEN(1,FILE=path)
DO  ik=1,nki
  IF(ik == 1)READ(1,*)nis
  READ (1,20)ki(1,ik),ki(2,ik),ki(3,ik)
  READ (1,10)(e(ik,i),i=1,nband)
END DO
CLOSE(1)
10          FORMAT(10F8.4)
20          FORMAT(10X,f10.3,f10.3,f10.3)


DO  ik=1,nki
  DO  i=1,nband
    e(ik,i)=e(ik,i)/hartree
    IF(i > nvalls)e(ik,i)=e(ik,i)+gw/hartree
  END DO
END DO

PRINT*,
PRINT*,
PRINT*,'Be careful MoS2 gap should be GW corrected'
PRINT*,'Band structure is full relativistic'
PRINT*,




!            generator 1.B.Z.
!            Dio programa koji pomocu operacija tockaste grupe i vektora iz
!            I.B.Z. generira sve (MEDJUSOBNO RAZLICITE!!!) v. vektore u 1.B.Z.
!            Ntot-Tot number of different points ''ktot'' inside 1.B.Z



jk=0
ntot=0
DO  i=1,nsim
  DO  ik=1,nki
    it=1
    jk=jk+1
    DO  n=1,3
      k(n,jk)=zero
      DO  m=1,3
        k(n,jk)=k(n,jk)+r(i,n,m)*ki(m,ik)
      END DO
    END DO
    IF(jk > 1)THEN
      DO  lk=1,jk-1
        IF(ABS(k(1,jk)-k(1,lk)) <= eps)THEN
          IF(ABS(k(2,jk)-k(2,lk)) <= eps)THEN
            IF(ABS(k(3,jk)-k(3,lk)) <= eps)THEN
              it=2
            END IF
          END IF
        END IF
      END DO
    END IF
    IF(it == 1)THEN
      ntot=ntot+1
      ktot(1,ntot)=k(1,jk)
      ktot(2,ntot)=k(2,jk)
      ktot(3,ntot)=k(3,jk)
    END IF
  END DO
END DO


!             Checking 1BZ integration
nel=0
DO  ik=1,ntot
  kx=ktot(1,ik)
  ky=ktot(2,ik)
  kz=ktot(3,ik)
  DO  n=1,nband
    IF(n == 1)THEN
      it=1
      IF(ik <= nki)THEN
        k1=ik
        it=2
      ELSE
        DO  i=2,nsim
          k11=ri(i,1,1)*kx+ri(i,1,2)*ky+ri(i,1,3)*kz
          k22=ri(i,2,1)*kx+ri(i,2,2)*ky+ri(i,2,3)*kz
          k33=ri(i,3,1)*kx+ri(i,3,2)*ky+ri(i,3,3)*kz
          DO  j=1,nki
            IF(DABS(k11-ki(1,j)) <= eps)THEN
              IF(DABS(k22-ki(2,j)) <= eps)THEN
                IF(DABS(k33-ki(3,j)) <= eps)THEN
                  it=2
                  k1=j
                  GO TO 5022
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF
      IF(it == 1)THEN
        PRINT*,'Can not find wave vector K=',ik, 'in I.B.Z.'
        STOP
      END IF
      5022          CONTINUE
    END IF
    IF(e(k1,n) < ef)nel=nel+1.0
  END DO
END DO
nel=nel/ntot

DO  i=1,ntot
  WRITE(887,*)ktot(1,i),ktot(2,i)
END DO


epsq=SQRT(4.0*pi*c0/(vcell*ntot))


!           KC transformation matrix from rec.cryst. axes to cart.koord.
!           If g' is vector in rec.cryst. axes then a=KC*a' is vector in cart. axes


fajl='/MoS2.sc.out'
path=trim(root)//trim(fajl)
tag='     reciprocal axes: (cart. coord.'
OPEN(1,FILE=path)
DO  i=1,100000
  READ(1,'(a)')buffer
  IF(buffer == tag)THEN
    DO  j=1,3
      READ(1,70)kc(1,j),kc(2,j),kc(3,j)
    END DO
    GO TO 998
  END IF
END DO
70          FORMAT(23X,3F10.3)
998         CONTINUE
CLOSE(1)


!           Reading the reciprocal vectors in crystal coordinates and transformation
!           in Cartezi cordinates.
OPEN (1,FILE='gvectors.dat')
DO  i=1,8
  READ(1,*)nis
END DO
DO  ig=1,ng
  READ(1,100)gi(1),gi(2),gi(3)
  IF(ig == 1)THEN
    IF(gi(1) /= 0.OR.gi(2) /= 0.OR.gi(3) /= 0)THEN
      PRINT*,'*********************************'
      PRINT*,'WARRNING!, G vectors input is wrong!!'
      PRINT*,'G(1) is not (0,0,0)!!'
      STOP
    END IF
  END IF
!           transformation in cart.coord (also!, after this all G components are in 2pi/a0 units)
  DO   n=1,3
    g(n,ig)=zero
    DO   m=1,3
      g(n,ig)=g(n,ig)+kc(n,m)*DBLE(gi(m))
    END DO
  END DO
  parg(ig)=gi(3)
END DO
100         FORMAT(i10,i11,i11)
CLOSE(1)



!            Reciprocal vectors for crystal local field effects calculations in array ''Glf(3,Nlf)''


nlf=0
IF(lf == 1)THEN
  DO  ig=1,ng
    IF(g(1,ig) == 0.0.AND.g(2,ig) == 0.0)THEN
      eref=gcar*gcar*g(3,ig)*g(3,ig)/2.0
      IF(eref <= ecut)THEN
        nlf=nlf+1
        glf(1,nlf)=0.0
        glf(2,nlf)=0.0
        glf(3,nlf)=g(3,ig)
        IF((parg(ig)/2)*2 == parg(ig))THEN
          parg(nlf)=1
        ELSE
          parg(nlf)=-1
        END IF
      END IF
    END IF
  END DO
ELSE
  DO  ig=1,ng
    eref=gcar*gcar*(g(1,ig)*g(1,ig)+g(2,ig)*g(2,ig)+ g(3,ig)*g(3,ig))/2.0
    IF(eref <= ecut)THEN
      nlf=nlf+1
      glf(1,nlf)=g(1,ig)
      glf(2,nlf)=g(2,ig)
      glf(3,nlf)=g(3,ig)
      IF((parg(ig)/2)*2 == parg(ig))THEN
        parg(nlf)=1
      ELSE
        parg(nlf)=-1
      END IF
    END IF
  END DO
END IF
IF(nlf > nlfd)THEN
  PRINT*,'Nlf is bigger than Nlfd'
  STOP
END IF



!             GENERIRANJE male BZ OKO 'K' TOCKE
!             K-point in Cartesi cordintes in  alatt 0.333333  0.577350  0.000000
kpointx=0.333333
kpointy=0.577350
kpointz=0.000000
nkbz=0


!            generiranje 'male' BZ

jk=0
nkbz=0
DO  i=1,nsim
  DO  ik=1,nki
    IF(ABS(ki(2,ik)) <= qmax)THEN
      it=1
      jk=jk+1
      DO  n=1,3
        k(n,jk)=zero
        DO  m=1,3
          k(n,jk)=k(n,jk)+r(i,n,m)*ki(m,ik)
        END DO
      END DO
      IF(jk > 1)THEN
        DO  lk=1,jk-1
          IF(ABS(k(1,jk)-k(1,lk)) <= eps)THEN
            IF(ABS(k(2,jk)-k(2,lk)) <= eps)THEN
              IF(ABS(k(3,jk)-k(3,lk)) <= eps)THEN
                it=2
              END IF
            END IF
          END IF
        END DO
      END IF
      IF(it == 1)THEN
        nkbz=nkbz+1
        kbz(1,nkbz)=k(1,jk)
        kbz(2,nkbz)=k(2,jk)
        kbz(3,nkbz)=k(3,jk)
      END IF
    END IF
  END DO
END DO

!            Pomak 'male BZ' iz 0 0 0 u K tocku


DO  i=1,nkbz
  kbz(1,i)=kbz(1,i)+kpointx
  kbz(2,i)=kbz(2,i)+kpointy
END DO

!            PROVJERA DA LI IMA ISTIH TOCAKA u OKLICI K tocke
DO  i=1,nkbz
  kx=kbz(1,i)
  ky=kbz(2,i)
  kz=kbz(3,i)
  DO  j=1,nkbz
    IF(j /= i)THEN
      IF(ABS(kx-kbz(1,j)) <= eps)THEN
        IF(ABS(ky-kbz(2,j)) <= eps)THEN
          IF(ABS(kz-kbz(3,j)) <= eps)THEN
            PRINT*,'Postoje jednake k tocke u okolini K tocke'
            STOP
          END IF
        END IF
      END IF
    END IF
  END DO
END DO




DO  i=1,nkbz
  WRITE(889,*)kbz(1,i),kbz(2,i)
END DO






!             Info file

OPEN(55,FILE='Info')
WRITE(55,*)'***************General***********************'
WRITE(55,*)''
WRITE(55,*)'Number of point symmetry operation is',nsim
IF(lf == 1)WRITE(55,*)'Local field effcts in z-dir'
IF(lf == 3)WRITE(55,*)'Local field in all xyz-dir'
WRITE(55,*)'Number of local field vectors is',nlf
WRITE(55,*)'Number of different K vectors in 1.B.Z. is',ntot
WRITE(55,*)'Number of K vectors in I.B.Z. is',nki
WRITE(55,*)'Number of K points around K-point is',nkbz
WRITE(55,*)'Number of bands is               ',nband
WRITE(55,*)''
WRITE(55,*)'************* Checking 1BZ integration*******'
WRITE(55,*)''
WRITE(55,12)'Number of electrons(1BZ integration)=',nel
WRITE(55,*)'Number of electrons(unit cell)=',nelqe
error=ABS((nelqe-nel)/nelqe)
WRITE(55,99)'Relative error=',error*100.0,'%'
IF(error > 0.05)THEN
  WRITE(55,*)'WARRNING!!-1BZ INTEGRATION IS BAD!.'
END IF
CLOSE(55)
88         FORMAT(a25,3F10.4,a5)
99         FORMAT(a25,f8.4,a5)
12         FORMAT(a40,f7.4)


IF(MOD == 2)GO TO 8889



!             Upis time ordered zasjenjene kulonske interakcije W_GG'^T(Q,\omega)
DO  iq=2,nq
  q=(iq-1)*qmin
  vq=(2.0*pi)/q
  dato='W_Qi'
  nord=INDEX(dato,'i', back =.false.)
  IF(iq < 10)THEN
    WRITE(dato(nord:nord),'(i1)')iq
  ELSE IF(iq >= 10.AND.iq < 100)THEN
    WRITE(dato(nord:nord+1),'(i2)')iq
  ELSE
    WRITE(dato(nord:nord+2),'(i3)')iq
  END IF
  path=trim(root1)//trim(dato)
  
  OPEN(74,FILE=path)
  READ(74,*)nis
  READ(74,44)((wt(iq,ig,jg),jg=1,nlf),ig=1,nlf)
  CLOSE(74)
  44            FORMAT(10F15.5)
  
!             trazenje alpha
  IF(iq == 2)alpha=(c0*vq/REAL(wt(2,1,1))-1.0)/qmin
  
  
END DO





!             valni vektor K
DO  ik1=1,nkbz
  
  PRINT*,ik1
  
  kx=kbz(1,ik1)
  ky=kbz(2,ik1)
  kz=kbz(3,ik1)
  
  
  
!              trazenje K prvo u 1.B.Z a onda u I.B.Z.
  it=1
  DO  ig=1,ng
    DO  jk=1,ntot
      IF(DABS(kx-g(1,ig)-ktot(1,jk)) <= eps)THEN
        IF(DABS(ky-g(2,ig)-ktot(2,jk)) <= eps)THEN
          IF(DABS(kz-g(3,ig)-ktot(3,jk)) <= eps)THEN
            it=2
            ig01=ig
            DO  i=1,nsim
              k11=ri(i,1,1)*ktot(1,jk)+ri(i,1,2)*ktot(2,jk)+  &
                  ri(i,1,3)*ktot(3,jk)
              k22=ri(i,2,1)*ktot(1,jk)+ri(i,2,2)*ktot(2,jk)+  &
                  ri(i,2,3)*ktot(3,jk)
              k33=ri(i,3,1)*ktot(1,jk)+ri(i,3,2)*ktot(2,jk)+  &
                  ri(i,3,3)*ktot(3,jk)
              DO  j=1,nki
                IF(DABS(k11-ki(1,j)) <= eps)THEN
                  IF(DABS(k22-ki(2,j)) <= eps)THEN
                    IF(DABS(k33-ki(3,j)) <= eps)THEN
                      it=3
                      r1=i
                      k1=j
                      EXIT
                    END IF
                  END IF
                END IF
              END DO
            END DO
          END IF
        END IF
      END IF
    END DO
  END DO
  2001                 CONTINUE
  
  
  IF(it == 1)THEN
    PRINT*,'Can not find wave vector K=',ik1, 'in 1.B.Z.'
    STOP
  ELSE IF(it == 2)THEN
    PRINT*,'Can not find wave vector K=',ik1, 'in I.B.Z.'
    STOP
  END IF
  
  
  
  
!             valni vektor K'
  DO  ik2=1,nkbz
    
    
    kx=kbz(1,ik2)
    ky=kbz(2,ik2)
    kz=kbz(3,ik2)
    
    
    it=1
!              trazenje K' prvo u 1.B.Z a onda u I.B.Z.
    
    DO  ig=1,ng
      DO  jk=1,ntot
        IF(DABS(kx-g(1,ig)-ktot(1,jk)) <= eps)THEN
          IF(DABS(ky-g(2,ig)-ktot(2,jk)) <= eps)THEN
            IF(DABS(kz-g(3,ig)-ktot(3,jk)) <= eps)THEN
              it=2
              ig02=ig
              DO  i=1,nsim
                k11=ri(i,1,1)*ktot(1,jk)+ri(i,1,2)*ktot(2,jk)+  &
                    ri(i,1,3)*ktot(3,jk)
                k22=ri(i,2,1)*ktot(1,jk)+ri(i,2,2)*ktot(2,jk)+  &
                    ri(i,2,3)*ktot(3,jk)
                k33=ri(i,3,1)*ktot(1,jk)+ri(i,3,2)*ktot(2,jk)+  &
                    ri(i,3,3)*ktot(3,jk)
                DO  j=1,nki
                  IF(DABS(k11-ki(1,j)) <= eps)THEN
                    IF(DABS(k22-ki(2,j)) <= eps)THEN
                      IF(DABS(k33-ki(3,j)) <= eps)THEN
                        it=3
                        r2=i
                        k2=j
                        EXIT
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END IF
          END IF
        END IF
      END DO
    END DO
    2111                 CONTINUE
    
    
    IF(it == 1)THEN
      PRINT*,'Can not find wave vector K=',ik2, 'in 1.B.Z.'
      STOP
    ELSE IF(it == 2)THEN
      PRINT*,'Can not find wave vector K=',ik2, 'in I.B.Z.'
      STOP
    END IF
    
    
! transfer vektor q koji ulazi u W(q)
    qx=kbz(1,ik2)-kbz(1,ik1)
    qy=kbz(2,ik2)-kbz(2,ik1)
    q=gcar*SQRT(qx*qx+qy*qy)
    iq=q/qmin+1
    
    
    
    IF(iq > nq)GO TO 9999
    
    
!              ''K=G01+R1*K1''.
!              ''K'=G02+R2*K2''.
!              K1 i K2 su integeri koji predstavljaju valne vektore u IBZ pridruzene
!              valnim vektorima K i K'
!              iG01 i iG02 su integeri rec vec. G koji translatira K i K' u 1BZ
    
!             vrpca n
    DO   n=nval,nval+1
      
      IF(n == nval) ni=1
      IF(n == nval+1) ni=2
!             vrpca m
      DO  m=n,n
        
        
        
        CALL paths(root,k1,k2,n,m,pathk1,pathk2,bandn,bandm)
        
        
        
        
!         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za
!         fiksni K1,K2,n i m
        
!               Otvaranje atribute za INFO
        CALL iotk_open_read(10,pathk1)
        CALL iotk_scan_empty(10,"INFO",attr=attr)
        CALL iotk_scan_attr(attr,"igwx",ng1)
!               Alociranje polja C1
        allocate (c1(ng1))
!               Ucitavanje podataka iza evc.n
        CALL iotk_scan_dat(10,bandn,c1)
        CALL iotk_close_read(10)
!               Otvaranje atribute za INFO
        CALL iotk_open_read(10,pathk2)
        CALL iotk_scan_empty(10,"INFO",attr=attr)
        CALL iotk_scan_attr(attr,"igwx",ng2)
!               Alociranje polja C2
        allocate (c2(ng2))
!               Ucitavanje podataka iza evc.m
        CALL iotk_scan_dat(10,bandm,c2)
        CALL iotk_close_read(10)
        
!                Konstrukcija stupca matricnih elementa MnmK1K2(ni,mi,G)
        
        
        IF(ngd > ng1)THEN
          WRITE(*,*)'NGd is bigger than NG1=',ng1
          STOP
        ELSE IF(ngd > ng2)THEN
          WRITE(*,*)'NGd is bigger than NG2=',ng2
          STOP
        END IF
        
        
!                 matrix elements
        igfast=0
        DO  ig=1,nlf
          mnmk1k2(ni,ig)=czero
          DO  ig1=1,ngd
            igfast=igfast+1
            gxx1=g(1,ig1)
            gyy1=g(2,ig1)
            gzz1=g(3,ig1)
            k11=r(r1,1,1)*gxx1+r(r1,1,2)*gyy1+r(r1,1,3)*gzz1
            k22=r(r1,2,1)*gxx1+r(r1,2,2)*gyy1+r(r1,2,3)*gzz1
            k33=r(r1,3,1)*gxx1+r(r1,3,2)*gyy1+r(r1,3,3)*gzz1
            k11=k11+glf(1,ig)
            k22=k22+glf(2,ig)
            k33=k33+glf(3,ig)
            k11=k11+g(1,ig02)-g(1,ig01)
            k22=k22+g(2,ig02)-g(2,ig01)
            k33=k33+g(3,ig02)-g(3,ig01)
            gxx1=ri(r2,1,1)*k11+ri(r2,1,2)*k22+ri(r2,1,3)*k33
            gyy1=ri(r2,2,1)*k11+ri(r2,2,2)*k22+ri(r2,2,3)*k33
            gzz1=ri(r2,3,1)*k11+ri(r2,3,2)*k22+ri(r2,3,3)*k33
            IF(jump == 1)THEN
              ig2_loop: DO  ig2=1,ng2
                gfast(igfast)=ng2+1
                gxx2=g(1,ig2)
                gyy2=g(2,ig2)
                gzz2=g(3,ig2)
                IF(DABS(gxx2-gxx1) < eps)THEN
                  IF(DABS(gyy2-gyy1) < eps)THEN
                    IF(DABS(gzz2-gzz1) < eps)THEN
                      gfast(igfast)=ig2
                      ! GO TO 1111
                      exit ig2_loop
                    END IF
                  END IF
                END IF
              END DO
            END IF
            ! 1111              CONTINUE
            ig2=gfast(igfast)
            IF(ig2 <= ng2)THEN
              mnmk1k2(ni,ig)=mnmk1k2(ni,ig)+ CONJG(c1(ig1))*c2(ig2)
            END IF
          END DO
        END DO
        jump=2
        
        
        deallocate(c1)
        deallocate(c2)
        
!               end of m loop
      END DO
      
      
!               end of n loop
    END DO
    jump=1
    
    
    fockk(ik1,ik2)=czero
    
    IF(ik1 /= ik2)THEN
      DO  ig=1,nlf
        DO  jg=1,nlf
          w=wt(iq,ig,jg)
          fockk(ik1,ik2)=fockk(ik1,ik2)- w*CONJG(mnmk1k2(1,ig))*mnmk1k2(2,jg)
        END DO
      END DO
      
!               komponenta Q=0 korekcija
    ELSE IF(ik1 == ik2)THEN
      fockk(ik1,ik2)=-vcell*ntot*LOG(1.0+alpha*epsq)/alpha
    END IF
    7768            CONTINUE
    7769            CONTINUE
    
    
    9999            CONTINUE
!               end of K'  loop
  END DO
!               end of K loop
END DO



!              upisivanje Fock kernela i DeltaE
OPEN(234,FILE='FockK')
WRITE(234,44)((fockk(ik1,ik2),ik2=1,nkbz),ik1=1,nkbz)
CLOSE(234)
GO TO 9981



8889           CONTINUE

!*********************************************
!              Ovdje pocinje mod=2
!*********************************************


PRINT*,'Reading Fock kernel'

!              citanje Fock kernela i DeltaE
OPEN(234,FILE='FockK')
READ(234,44)((fockk(ik1,ik2),ik2=1,nkbz),ik1=1,nkbz)
CLOSE(234)


!               GENERIRANJE STRUJNIH VRHOVA jmu

DO  pol=1,2
  
  IF(pol == 1)PRINT*,'Calculation of j_x'
  IF(pol == 2)PRINT*,'Calculation of j_z'
  
  
  
!               valni vektor K
  DO  ik1=1,nkbz
    
    kx=kbz(1,ik1)
    ky=kbz(2,ik1)
    kz=kbz(3,ik1)
    
!               trazenje K prvo u 1.B.Z a onda u I.B.Z.
    it=1
    DO  ig=1,ng
      DO  jk=1,ntot
        IF(DABS(kx-g(1,ig)-ktot(1,jk)) <= eps)THEN
          IF(DABS(ky-g(2,ig)-ktot(2,jk)) <= eps)THEN
            IF(DABS(kz-g(3,ig)-ktot(3,jk)) <= eps)THEN
              it=2
              ig01=ig
              DO  i=1,nsim
                k11=ri(i,1,1)*ktot(1,jk)+ri(i,1,2)*ktot(2,jk)+  &
                    ri(i,1,3)*ktot(3,jk)
                k22=ri(i,2,1)*ktot(1,jk)+ri(i,2,2)*ktot(2,jk)+  &
                    ri(i,2,3)*ktot(3,jk)
                k33=ri(i,3,1)*ktot(1,jk)+ri(i,3,2)*ktot(2,jk)+  &
                    ri(i,3,3)*ktot(3,jk)
                DO  j=1,nki
                  IF(DABS(k11-ki(1,j)) <= eps)THEN
                    IF(DABS(k22-ki(2,j)) <= eps)THEN
                      IF(DABS(k33-ki(3,j)) <= eps)THEN
                        it=3
                        r1=i
                        k1=j
                        EXIT
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END IF
          END IF
        END IF
      END DO
    END DO
    2101                 CONTINUE
    
    
    IF(it == 1)THEN
      PRINT*,'Can not find wave vector K=',ik1, 'in 1.B.Z.'
      STOP
    ELSE IF(it == 2)THEN
      PRINT*,'Can not find wave vector K=',ik1, 'in I.B.Z.'
      STOP
    END IF
    
!             valni vektor K'
    DO  ik2=ik1,ik1
      
      kx=kbz(1,ik2)
      ky=kbz(2,ik2)
      kz=kbz(3,ik2)
      
      
      it=1
!              trazenje K' prvo u 1.B.Z a onda u I.B.Z.
      
      DO  ig=1,ng
        DO  jk=1,ntot
          IF(DABS(kx-g(1,ig)-ktot(1,jk)) <= eps)THEN
            IF(DABS(ky-g(2,ig)-ktot(2,jk)) <= eps)THEN
              IF(DABS(kz-g(3,ig)-ktot(3,jk)) <= eps)THEN
                it=2
                ig02=ig
                DO  i=1,nsim
                  k11=ri(i,1,1)*ktot(1,jk)+ri(i,1,2)*ktot(2,jk)+  &
                      ri(i,1,3)*ktot(3,jk)
                  k22=ri(i,2,1)*ktot(1,jk)+ri(i,2,2)*ktot(2,jk)+  &
                      ri(i,2,3)*ktot(3,jk)
                  k33=ri(i,3,1)*ktot(1,jk)+ri(i,3,2)*ktot(2,jk)+  &
                      ri(i,3,3)*ktot(3,jk)
                  DO  j=1,nki
                    IF(DABS(k11-ki(1,j)) <= eps)THEN
                      IF(DABS(k22-ki(2,j)) <= eps)THEN
                        IF(DABS(k33-ki(3,j)) <= eps)THEN
                          it=3
                          r2=i
                          k2=j
                          EXIT
                        END IF
                      END IF
                    END IF
                  END DO
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
      2411           CONTINUE
      
      
      IF(it == 1)THEN
        PRINT*,'Can not find wave vector K=',ik2, 'in 1.B.Z.'
        STOP
      ELSE IF(it == 2)THEN
        PRINT*,'Can not find wave vector K=',ik2, 'in I.B.Z.'
        STOP
      END IF
      
      
      
!              ''K=G01+R1*K1''.
!              ''K'=G02+R2*K2''.
!              K1 i K2 su integeri koji predstavljaju valne vektore u IBZ pridruzene
!              valnim vektorima K i K'
!              iG01 i iG02 su integeri rec vec. G koji translatira K i K' u 1BZ
      
      
      
      
      DO  n=nval,nval
        DO  m=nval+1,nval+1
          
          
          
          CALL paths(root,k1,k2,n,m,pathk1,pathk2,bandn,bandm)
          
          
          
!         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za
!         fiksni K1,K2,n i m
          
!               Otvaranje atribute za INFO
          CALL iotk_open_read(10,pathk1)
          CALL iotk_scan_empty(10,"INFO",attr=attr)
          CALL iotk_scan_attr(attr,"igwx",ng1)
!               Alociranje polja C1
          allocate (c1(ng1))
!               Ucitavanje podataka iza evc.n
          CALL iotk_scan_dat(10,bandn,c1)
          CALL iotk_close_read(10)
!               Otvaranje atribute za INFO
          CALL iotk_open_read(10,pathk2)
          CALL iotk_scan_empty(10,"INFO",attr=attr)
          CALL iotk_scan_attr(attr,"igwx",ng2)
!               Alociranje polja C2
          allocate (c2(ng2))
!               Ucitavanje podataka iza evc.m
          CALL iotk_scan_dat(10,bandm,c2)
          CALL iotk_close_read(10)
          
!                Konstrukcija stupca matricnih elementa MnmK1K2(G)
          
          
          IF(ngd > ng1)THEN
            WRITE(*,*)'NGd is bigger than NG1=',ng1
            STOP
          ELSE IF(ngd > ng2)THEN
            WRITE(*,*)'NGd is bigger than NG2=',ng2
            STOP
          END IF
          
          
!              matrix elements
          igfast=0
          DO  ig=1,1
            current(pol,ik1)=czero
            DO  ig1=1,ngd
              igfast=igfast+1
              gxx1=g(1,ig1)
              gyy1=g(2,ig1)
              gzz1=g(3,ig1)
              k11=r(r1,1,1)*gxx1+r(r1,1,2)*gyy1+r(r1,1,3)*gzz1
              k22=r(r1,2,1)*gxx1+r(r1,2,2)*gyy1+r(r1,2,3)*gzz1
              k33=r(r1,3,1)*gxx1+r(r1,3,2)*gyy1+r(r1,3,3)*gzz1
              IF(pol == 1)THEN
                struja=(2.0*kx+glf(1,ig)+2.0*k11-2.0*g(1,ig01))*gcar
              ELSE IF(pol == 2)THEN
                struja=(2.0*kz+glf(3,ig)+2.0*k33-2.0*g(3,ig01))*gcar
              END IF
              k11=k11+glf(1,ig)
              k22=k22+glf(2,ig)
              k33=k33+glf(3,ig)
              k11=k11+g(1,ig02)-g(1,ig01)
              k22=k22+g(2,ig02)-g(2,ig01)
              k33=k33+g(3,ig02)-g(3,ig01)
              gxx1=ri(r2,1,1)*k11+ri(r2,1,2)*k22+ri(r2,1,3)*k33
              gyy1=ri(r2,2,1)*k11+ri(r2,2,2)*k22+ri(r2,2,3)*k33
              gzz1=ri(r2,3,1)*k11+ri(r2,3,2)*k22+ri(r2,3,3)*k33
              IF(jump == 1)THEN
                DO  ig2=1,ng2
                  gfast(igfast)=ng2+1
                  gxx2=g(1,ig2)
                  gyy2=g(2,ig2)
                  gzz2=g(3,ig2)
                  IF(DABS(gxx2-gxx1) < eps)THEN
                    IF(DABS(gyy2-gyy1) < eps)THEN
                      IF(DABS(gzz2-gzz1) < eps)THEN
                        gfast(igfast)=ig2
                        GO TO 1011
                      END IF
                    END IF
                  END IF
                END DO
              END IF
              1011           CONTINUE
              ig2=gfast(igfast)
              IF(ig2 <= ng2)THEN
                current(pol,ik1)=current(pol,ik1)+  &
                    0.5D0*CONJG(c1(ig1))*struja*c2(ig2)
              END IF
!              kraj po iG1
            END DO
!              kraj po C.L.F. iG
          END DO
          jump=2
          
          
          
          IF(spin == 1)delta(ik1)=e(k1,nvalls-1)-e(k1,nvalls+1)
          IF(spin == 2)delta(ik1)=e(k1,nvalls)-e(k1,nvalls+2)
          
          
!              end of m do loop
        END DO
        
        
!              end of n do loop
      END DO
      jump=1
      
!              ond of B.Z do loop
    END DO
  END DO
  
!              kraj po polarizaciji
END DO


!               Rijesavanje BSE-FOCK jednadzbe


IF(spin == 1)THEN
  OPEN(334,FILE='Pi_ladder_up_x')
  OPEN(335,FILE='Pi_ladder_up_z')
ELSE IF(spin == 2)THEN
  OPEN(334,FILE='Pi_ladder_down_x')
  OPEN(335,FILE='Pi_ladder_down_z')
END IF


DO    io=1,no
  o=omin+(io-1)*domega
  PRINT*,io
  
  DO  ik1=1,nkbz
    DO  ik2=1,nkbz
      fock(ik1,ik2)=czero
      UNIT(ik1,ik2)=czero
    END DO
    UNIT(ik1,ik1)=rone
  END DO
  
  
!               Konstrukcija matrice L^0*Xi^F
!               Konstrukcija matrice Xi0=L^0*Xi^FL^0
  
  DO  ik1=1,nkbz
    l0i=1.0/(o+delta(ik1)+ione*eta)
    DO  ik2=1,nkbz
      l0j=1.0/(o+delta(ik2)+ione*eta)
      fock(ik1,ik2)=UNIT(ik1,ik2)- l0i*fockk(ik1,ik2)/(ntot*vcell)
      xi0(ik1,ik2)=l0i*fockk(ik1,ik2)*l0j
    END DO
  END DO
  
  
! inverzija Fock kernela
  CALL gjel(fock,nkbz,nkd,UNIT,nkbz,nkd)
  
  
!               konstrukcija 4-point polarizabilnosti L^ladd
  
  
  
  DO  ik1=1,nkbz
    DO  ik2=1,nkbz
      lladd(ik1,ik2)=czero
      DO  ik=1,nkbz
        lladd(ik1,ik2)=lladd(ik1,ik2)+fock(ik1,ik)*xi0(ik,ik2)
      END DO
    END DO
  END DO
  
  
  
!               Konstrukcija Ladder ireduciblne current-current polarisabilnosti chi_ladd_\mu\mu
  
  
  
  DO  pol=1,2
    chi_ladd(pol)=czero
    chi_para(pol)=czero
    DO  ik1=1,nkbz
      l0i=(o/delta(ik1))*(1.0/(o+delta(ik1)+ione*eta))
      DO  ik2=1,nkbz
        j1=current(pol,ik1)
        j2=CONJG(current(pol,ik2))
        chi_ladd(pol)=chi_ladd(pol)-  &
            valley*j1*lladd(ik1,ik2)*j2/(ntot*vcell*ntot*vcell)
      END DO
!               Paramagnetska struja-struja polarizabilost
      chi_para(pol)= chi_para(pol)+ valley*j1*l0i*CONJG(j1)/(ntot*vcell)
    END DO
  END DO
  
  WRITE(334,*)o*hartree,chi_ladd(1)
  WRITE(335,*)o*hartree,chi_ladd(2)
  
  WRITE(104,*)o*hartree,imag(chi_para(1)+chi_ladd(1))
  WRITE(105,*)o*hartree,imag(chi_para(1))
  WRITE(106,*)o*hartree,REAL(chi_para(1)+chi_ladd(1))
  WRITE(107,*)o*hartree,REAL(chi_para(1))
  
  WRITE(204,*)o*hartree,imag(chi_para(2)+chi_ladd(2))
  WRITE(205,*)o*hartree,imag(chi_para(2))
  WRITE(206,*)o*hartree,REAL(chi_para(2)+chi_ladd(2))
  WRITE(207,*)o*hartree,REAL(chi_para(2))
  
  
  
  
!               omega do loop
END DO



CLOSE(334)
CLOSE(335)
CLOSE(336)



9981            CONTINUE





END PROGRAM surface_loss

