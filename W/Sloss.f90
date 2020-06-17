PROGRAM surface_loss
 
! Code converted using TO_F90 by Alan Miller
! Date: 2020-06-17  Time: 13:22:14

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
!        Gama-Damping parameter in eV
!        Ecut-cutoff energy for crystal local field calculations
!        Vcell-unit-cell volume in a.u.^3
!        a0-unit cell parameter in parallel direction in a.u.
!        c0-unit cell parameter in perpendicular direction in a.u. (z-separation between supercells)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Quantum Esspresso:
!        verbosity           = 'high'
!        VALID JUST FOR NORM-CONSERVING PSEUDOPOTENTIALS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




use iotk_module
implicitnone
CHARACTER (LEN=1) :: iotk_attlenx) :: attr
LOGICAL :: :: found

INTEGER :: nki,nband,ik,i,nk,j,jk,it,lk,ntot,ig0,nsim,iq,  &
    ng,io,no,nq,nmpx,nmpy,nmpz,n,m,ig,r1,k1,r2,k2,  &
    nlf,ng1,ng2,ngd,ig1,ig2,jg,nlfd,kg,jo,jump,loss,  &
    igfast,ikmin,nelqe,lf,kg1,kg2,nord
PARAMETER(nmpx=201,nmpy=201,nmpz=1,nki=6835,nband=60,nelqe=18,  &
    nk=48*nki,ngd=4000,ng=8000,no=2001,nq=2,nlfd=50)
! no je broj frekvencija,  nq je broj valnih vektora tu je 2 jer je rucno paralelizirano!
! Nlf


!        skalars
REAL*8 a0,c0,eps,kx,ky,kz,kqx,kqy,kqz,qgx,qgy,qgz,omin,omax,  &
    qx,qy,qz,kmin,domega,o,ef,k11,k22,k33,t,lor,de,gabs,kref,  &
    eref,ecut,gxx1,gyy1,gzz1,gxx2,gyy2,gzz2,fact,vcell,oi,oj,  &
    imchi0,rechi0,q,gcar,abohr,zero,nel,absq,error,qmax,  &
    gama,w1,w2,imw,wind,w2kk,kks,skk,windkk,krefm
doubleprecision pi,three,hartree,planck,ev
PARAMETER(hartree=2.0D0*13.6056923D0,ef=0.5554/hartree,  &
    a0=5.9715,c0=29.8575,pi=3.141592654D0,gcar=2.0*pi/a0,  &
    eps=1.0D-4,t=0.01/hartree,gama=0.05/hartree,  &
    ev=1.602176487D-19,planck=6.626196D-34,ecut=0.0,  &
    vcell=922.0586,abohr=0.5291772D0,zero=0.0)
DOUBLE COMPLEX a,ione,czero,rone,em,g0

!        arrays
INTEGER :: gfast,gi,parg,kq0
REAL*8 ki,e,r,ri,k,ktot,g,glf,v,s0,kc,glfv
DOUBLE COMPLEX UNIT,epsilon,chi
COMPLEX*8 mnmk1k2,chi0,wt,gammap,gammam
DIMENSION ki(3,nki),e(nki,nband),r(48,3,3),ri(48,3,3), k(3,nk),ktot(3,nk),  &
    v(nlfd,nlfd),! matr. gole coulomb. int.  &
    g(3,ng), ! polje valnih vektora G u recp. prost. za wfn.  &
    glfv(3,nlfd), glf(3,nlfd) ! generiran G vekt. (0,0,z) za V odnosn za chi.  &
    mnmk1k2(nlfd), ! nabojni vrhovi UNIT(nlfd,nlfd), ! jedinična matrica  &
    chi0(nlfd,nlfd), ! (eq. 2.89)  &
    epsilon(nlfd,nlfd), ! Epsilon (GG')  = I - V(GG')Chi0  &
    chi(nlfd,nlfd),   ! (eq. 2.88 nakon invertiranja)  &
    s0(no,nlfd,nlfd), ! korelacijska matrica gfast(nlfd*ngd),  &
    kc(3,3),gi(3), ! pomocne funkcije parg(ng), ! paritet svakog valnog vektora  &
    wt(no,nlfd,nlfd), ! time ordered RPa screened coulomb int. (eq. 2.93)  &
    gammap(nlfd,nlfd),gammam(nlfd,nlfd), ! za GW ne koristi se za ovaj dio  &
    kq0(100)

CHARACTER (LEN=100) :: bandn, bandm,  &
    nis, pathk1,pathk2,  &
    dato, root,path,fajl
CHARACTER (LEN=35) :: tag,buffer

INTEGER :: :: ngw, igwx, nbnd, nk1
COMPLEX(:: kind = 8),pointer, DIMENSION(:) :: c1,c2


rone=DCMPLX(1.0,0.0)  ! real 1
czero=DCMPLX(0.0,0.0) ! complex 0
ione=DCMPLX(0.0,1.0) ! imag 1


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
root='/home/vito/PROJECTS/MoS2-BSE/MoS2_201X201'

!             Crystal local field effects are included in z direction lf=1
!             Crystal local field effects are included in x,y,z direction lf=3


lf=1 ! crystal local field effect included in z for lf=1
jump=1 ! za 1 preskace trazenje wfn. u IBZ za sve bands m i n
three=3.0D0 ! broj 3 haha
omin=1.0D-5 ! raspon frekvencija u Ha
omax=2.0D0
domega=(omax-omin)/(no-1)



!           CALL FOR POINT GROUP TRANSFORMATIONS
!           Point group transformations are in Cartesian coordinate

CALL pointr(root,nsim,r,ri)



!           Upis valnih vektora iz irreducibilne Brillouinove
!           zone i pripadnih energijskih nivoa iz filea '****.band'.
!           wave vectors are in Cartesian coordinate



fajl='/MoS2.band'
path=trim(root)//trim(fajl)
OPEN(1,FILE=path)
DO  ik=1,nki ! k vektora u IBZ
  IF(ik == 1)READ(1,*)nis ! preskakanje
  READ (1,20)ki(1,ik),ki(2,ik),ki(3,ik)
  READ (1,10)(e(ik,i),i=1,nband)
END DO
CLOSE(1)
10          FORMAT(10F8.4)
20          FORMAT(10X,f10.3,f10.3,f10.3)


DO  ik=1,nki
  DO  i=1,nband
    e(ik,i)=e(ik,i)/hartree
    IF(i >= 10)e(ik,i)=e(ik,i)+1.0/hartree ! scissor op. ispravljanje DFT gapa na 1ev ( u ovom slucaju)
  END DO
END DO





!            generator 1.B.Z.
!            Dio programa koji pomocu operacija tockaste grupe i vektora iz
!            I.B.Z. generira sve (MEDJUSOBNO RAZLICITE!!!) v. vektore u 1.B.Z.
!            Ntot-Tot number of different points ''ktot'' inside 1.B.Z



jk=0 ! k tocka u ...?
ntot=0
DO  i=1,nsim ! loop over No. symmetries
  DO  ik=1,nki  ! loop over k points in IBZ
    it=1 ! ?
    jk=jk+1
    DO  n=1,3  ! loop over kx,ky,kz
      k(n,jk)=zero
      DO  m=1,3 ! loop over x,y,z
        k(n,jk)=k(n,jk)+r(i,n,m)*ki(m,ik) ! kreira nove k tocke u BZ pomocu simetrije
      END DO
    END DO
    IF(jk > 1)THEN
      DO  lk=1,jk-1
        IF(ABS(k(1,jk)-k(1,lk)) <= eps)THEN  ! je li razlicita tocka od neke prije vec kreirane
          IF(ABS(k(2,jk)-k(2,lk)) <= eps)THEN
            IF(ABS(k(3,jk)-k(3,lk)) <= eps)THEN
              it=2         ! preskakanje tocke
            END IF
          END IF
        END IF
      END DO
    END IF
    IF(it == 1)THEN ! ne postoji dodaj ju
      ntot=ntot+1
      ktot(1,ntot)=k(1,jk)
      ktot(2,ntot)=k(2,jk)
      ktot(3,ntot)=k(3,jk)
    END IF
  END DO
END DO


!             Checking 1BZ integration
nel=0 ! provjeri je li broj el. u FBZ odgovara stvarnom broju el. u jed. cel. Nelqe
DO  ik=1,ntot
  kx=ktot(1,ik)
  ky=ktot(2,ik)
  kz=ktot(3,ik)
  DO  n=1,nband ! loop over bands
    IF(n == 1)THEN
      it=1
      IF(ik <= nki)THEN
        k1=ik
        it=2
      ELSE
        DO  i=2,nsim ! loop over no. symmetries
          k11=ri(i,1,1)*kx+ri(i,1,2)*ky+ri(i,1,3)*kz
          k22=ri(i,2,1)*kx+ri(i,2,2)*ky+ri(i,2,3)*kz
          k33=ri(i,3,1)*kx+ri(i,3,2)*ky+ri(i,3,3)*kz
          DO  j=1,nki ! loop over k u IBZ
            IF(DABS(k11-ki(1,j)) <= eps)THEN ! eps proizvoljno mali broj
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
    IF(e(k1,n) < ef)nel=nel+1.0 ! zbroji za en. manje od fermijeve
  END DO
END DO
nel=2.0*nel/ntot

DO  i=1,ntot
  WRITE(887,*)ktot(1,i),ktot(2,i) ! output da vidimo kako izgleda FBZ
END DO



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
  DO   n=1,3 ! loop over dim
    g(n,ig)=zero
    DO   m=1,3
      g(n,ig)=g(n,ig)+kc(n,m)*DBLE(gi(m))
    END DO
  END DO
  parg(ig)=gi(3)  ! odreduje paritet za svaki Gi
END DO
100         FORMAT(i10,i11,i11)
CLOSE(1)



!            Reciprocal vectors for crystal local field effects calculations in array ''Glf(3,Nlf)''
! na temelju cutoffa eliminar G-vektore koji su izvan cut-offa
nlf=0 ! broj novih local field G vektora
IF(lf == 1)THEN
  DO  ig=1,ng ! loop over all G-vectors za 3D LFE
    IF(g(1,ig) == 0.0.AND.g(2,ig) == 0.0)THEN
      eref=gcar*gcar*g(3,ig)*g(3,ig)/2.0 ! energijski prag
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
  DO  ig=1,ng ! loop over all G-vectors za 1D LFE (samo u z-smjeru)
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




!             IBZ q LOOP STARTS HERE!!!


! iq=0 ne moze biti nula, opticki racun
! iq=2 do iq=...cutoff transfer q vektor!
! ikmin = min. valni vektor u BZ svi veci su visekratnici tog minimalnog
DO  iq=42,61
  
!             searching min. q=(qx,qy,qz) in GM direction
  kmin=1.0
  DO  i=1,ntot
    kref=SQRT(ktot(1,i)*ktot(1,i)+ ktot(2,i)*ktot(2,i)+ktot(3,i)*ktot(3,i))
    IF(kref == zero)GO TO 970
    IF(kref < kmin)THEN
      kmin=kref
      ikmin=i
      krefm=kmin
    END IF
    970           CONTINUE
  END DO
  
  
  qx=(iq-1)*ktot(1,ikmin)
  qy=(iq-1)*ktot(2,ikmin)
  qz=(iq-1)*ktot(3,ikmin)
  
  absq=SQRT(qx*qx+qy*qy+qz*qz)
  
!             Info file
  
  OPEN(55,FILE='Info')
  WRITE(55,*)'***************General***********************'
  WRITE(55,*)''
  WRITE(55,*)'Number of point symmetry operation is',nsim
  WRITE(55,88)'Wave vector (qx,qy,qz)=(',qx*gcar,qy*gcar, qz*gcar,') a.u.'
  WRITE(55,99)'|(qx,qy,qz)|=',absq*gcar,'a.u.'
  IF(lf == 1)WRITE(55,*)'Local field effcts in z-dir'
  IF(lf == 3)WRITE(55,*)'Local field in all xyz-dir'
  WRITE(55,*)'Number of local field vectors is',nlf
  WRITE(55,*)'Number of different K vectors in 1.B.Z. is',ntot
  WRITE(55,*)'Number of K vectors in I.B.Z. is',nki
  WRITE(55,*)'Number of bands is               ',nband
  WRITE(55,99)'gama dumping is ',gama*hartree*1000.0,'meV'
  WRITE(55,99)'Temperature is  ',t*hartree*1000.0,'meV'
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
  
  
  DO  io=1,no
    DO  ig=1,nlf
      DO  jg=1,nlf
        s0(io,ig,jg)=czero ! korelacijska funkcija (inicijalizacija)
      END DO
    END DO
  END DO
  
  
!              1.B.Z  LOOP STARTS HERE !!!!
  
  DO  ik=1,ntot ! loop over k-points in FBZ
    
    OPEN(122,FILE='status')
    WRITE(122,*)'iq=',iq
    WRITE(122,*)'ik=',ik
    CLOSE(122)
    
    
    kx=ktot(1,ik)
    ky=ktot(2,ik)
    kz=ktot(3,ik)
    
!              trazenje (kx,ky,kz) u ireducibilnoj zoni
    
    it=1
    IF(ik <= nki)THEN
      r1=1
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
                r1=i ! pridruzena point group
                k1=j  ! trazeni vektor u IBZ
                GO TO 5222
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
    5222           CONTINUE
    
    
    it=1
    kqx=kx+qx
    kqy=ky+qy
    kqz=kz+qz
    
!              trazenje (KQx,KQy) prvo u 1.B.Z a onda u I.B.Z. (jer novo genrirani k+q more bit negdje vani u FBZ)
    
    DO  ig=1,ng
      DO  jk=1,ntot
        IF(DABS(kqx-g(1,ig)-ktot(1,jk)) <= eps)THEN
          IF(DABS(kqy-g(2,ig)-ktot(2,jk)) <= eps)THEN
            IF(DABS(kqz-g(3,ig)-ktot(3,jk)) <= eps)THEN
              it=2
              ig0=ig
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
                        r2=i  ! pridruzena point gruop
                        k2=j  ! pridruzen taj vektor pronadjen u IBZ
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
      PRINT*,'Can not find wave vector K+Q=',ik,'+',iq, 'in 1.B.Z.'
      STOP
    ELSE IF(it == 2)THEN
      PRINT*,'Can not find wave vector K+Q=',ik,'+',iq, 'in I.B.Z.'
      STOP
    END IF
    
    
!              R1-integer, redni broj point operacije R1 u transformaciji ''K=R1*K1''.
!              K1-integer, redni broj valnog vektora K1 u transformaciji ''K=R1*K1''.
!              iG0 i R2-integeri, redni broj vektora reciprocne restke G0 i point operacije R2 u transformaciji ''K+Q=G0+R2*K2''.
!              K2-integer, redni broj valnog vektora K2 u transformaciji  ''K+Q=G0+R2*K2''.
    
    
!              petlje po vrpcama n i m
    
    DO  n=1,9 ! loop over full bands , oprez kod metala mora ici malo iznad popunjene (n+1)
      DO  m=10,nband !  loop over empty bands
        
        
        
        
        
        CALL paths(root,k1,k2,n,m,pathk1,pathk2,bandn,bandm) ! fajl za popunit wfn. za K1 i k2 i pripadajuce vrpce n i m
        
        
!         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za
!         fiksni K1,K2,n i m
        
!               Otvaranje atribute za INFO
        CALL iotk_open_read(10,pathk1)
        CALL iotk_scan_empty(10,"INFO",attr=attr) ! nalazi <info>
        CALL iotk_scan_attr(attr,"igwx",ng1) ! koliko imamo G-vektora /
!               Alociranje polja C1
        allocate (c1(ng1)) ! u ovo trpamo fourierove koeficijent u razvoju wfn.
!               Ucitavanje podataka iza evc.n
        CALL iotk_scan_dat(10,bandn,c1) ! cita iz binranog oblika
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
        
        
!                 matrix elements
        igfast=0
        DO  ig=1,nlf
          mnmk1k2(ig)=czero
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
            k11=k11+g(1,ig0)
            k22=k22+g(2,ig0)
            k33=k33+g(3,ig0)
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
                      GO TO 1111
                    END IF
                  END IF
                END IF
              END DO
            END IF
            1111              CONTINUE
            ig2=gfast(igfast)
            IF(ig2 <= ng2)THEN
              mnmk1k2(ig)=mnmk1k2(ig)+ CONJG(c1(ig1))*c2(ig2)
            END IF
          END DO
        END DO
        jump=2
        
        
!                 omega loop
        DO  io=1,no
          o=(io-1)*domega
          de=o+e(k1,n)-e(k2,m)
          lor=-gama/(de*de+gama*gama)
          IF(DABS(lor) >= 1.0D-3/gama)THEN
            DO  ig=1,nlf
              DO  jg=1,nlf
                s0(io,ig,jg)=s0(io,ig,jg)-  &
                    2.0*lor*mnmk1k2(ig)*CONJG(mnmk1k2(jg))/ (pi*ntot*vcell)
              END DO
            END DO
          END IF
          
          
        END DO
        
        deallocate(c1)
        deallocate(c2)
        
        
        222             CONTINUE
        
        
        
!                end of m do loop
      END DO
      
!                end of n do loop
    END DO
    jump=1
    
    834              CONTINUE
!                ond of 1.B.Z do loop
  END DO
  
  
  
!               Puting (qx,qy,qz) and Glf in cartezi coordinate
  
  qx=gcar*qx
  qy=gcar*qy
  qz=gcar*qz
  
  DO  ig=1,nlf
    glfv(1,ig)=gcar*glf(1,ig)
    glfv(2,ig)=gcar*glf(2,ig)
    glfv(3,ig)=gcar*glf(3,ig)
  END DO
  
  
  
!               new sum over omega
  DO  io=1,no-1
!                 print*,io
    oi=(io-1)*domega
    DO  ig=1,nlf
      DO  jg=1,nlf
        rechi0=0.0
!                  static limit
        IF(io == 1)THEN
          DO  jo=2,no
            oj=(jo-1)*domega
            fact=domega/oj
            IF(jo == 2)fact=3.0/2.0
            IF(jo == no)fact=0.5*domega/oj
            rechi0=rechi0+fact*s0(jo,ig,jg)
          END DO
          rechi0=-2.0*rechi0
        ELSE IF(io == 2)THEN
          DO  jo=1,no
            oj=(jo-1)*domega
            IF(jo /= io)fact=domega/(oi-oj)
            IF(jo == 1)fact=1.0
            IF(jo == 2)fact=0.0
            IF(jo == 3)fact=-3.0/2.0
            IF(jo == no)fact=0.5*domega/(oi-oj)
            rechi0=rechi0+fact*s0(jo,ig,jg)
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            rechi0=rechi0-fact*s0(jo,ig,jg)
          END DO
        ELSE IF(io == (no-1))THEN
          DO  jo=1,no
            oj=(jo-1)*domega
            IF(jo /= io)fact=domega/(oi-oj)
            IF(jo == 1)fact=0.5*domega/(oi-oj)
            IF(jo == (no-2))fact=3.0/2.0
            IF(jo == (no-1))fact=0.0
            IF(jo == no)fact=-1.0
            rechi0=rechi0+fact*s0(jo,ig,jg)
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            rechi0=rechi0-fact*s0(jo,ig,jg)
          END DO
        ELSE
          DO  jo=1,no
            oj=(jo-1)*domega
            IF(jo /= io)fact=domega/(oi-oj)
            IF(jo == 1)fact=0.5*domega/(oi-oj)
            IF(jo == (io-1))fact=3.0/2.0
            IF(jo == io)fact=0.0
            IF(jo == (io+1))fact=-3.0/2.0
            IF(jo == no)fact=0.5*domega/(oi-oj)
            rechi0=rechi0+fact*s0(jo,ig,jg)
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            rechi0=rechi0-fact*s0(jo,ig,jg)
          END DO
        END IF
        
        imchi0=-pi*s0(io,ig,jg)
        chi0(ig,jg)=CMPLX(rechi0,imchi0)
        
        
!                kraj po iG,jG
      END DO
    END DO
    
    
    
!                Calculation of the ''Chi''  by matrix invertion
    
    
!                MATRIX V(G,G')
    
    DO  ig=1,nlf
      gabs=SQRT((qx+glfv(1,ig))*(qx+glfv(1,ig))+  &
          (qy+glfv(2,ig))*(qy+glfv(2,ig)))
      IF(gabs == 0.0)gabs=eps
      DO  jg=1,nlf
        v(ig,jg)=0.0
        IF(glf(1,jg) == glf(1,ig))THEN
          IF(glf(2,jg) == glf(2,ig))THEN
            v(ig,jg)=4.0*pi*(1.0-EXP(-gabs*c0))/(gabs*c0)
            v(ig,jg)=v(ig,jg)*(gabs*gabs-glfv(3,ig)*glfv(3,jg))
            v(ig,jg)=v(ig,jg)/(gabs*gabs+glfv(3,ig)*glfv(3,ig))
            v(ig,jg)=v(ig,jg)/(gabs*gabs+glfv(3,jg)*glfv(3,jg))
            v(ig,jg)=-DBLE(parg(ig))*DBLE(parg(jg))*v(ig,jg)
            IF(glf(3,jg) == glf(3,ig))THEN
              v(ig,jg)=4.0*pi/(gabs*gabs+glfv(3,ig)*glfv(3,ig))+ v(ig,jg)
            END IF
          END IF
        END IF
      END DO
    END DO
    
    
    
    DO  ig=1,nlf
      DO  jg=1,nlf
        UNIT(ig,jg)=czero
      END DO
      UNIT(ig,ig)=rone
    END DO
    
    DO  ig=1,nlf
      DO  jg=1,nlf
        epsilon(ig,jg)=UNIT(ig,jg)
        DO  kg=1,nlf
          epsilon(ig,jg)=epsilon(ig,jg)-chi0(ig,kg)*v(kg,jg)
        END DO
      END DO
    END DO
    
    
!                invertiranje matrice ''epsilon = 1-Chi_0*V''
    
    
    CALL gjel(epsilon,nlf,nlfd,UNIT,nlf,nlfd)
    
    
    DO  ig=1,nlf
      DO  jg=1,nlf
        chi(ig,jg)=czero
        DO  kg=1,nlf
          chi(ig,jg)=chi(ig,jg)+epsilon(ig,kg)*chi0(kg,jg)
        END DO
      END DO
    END DO
    
    
!                SCREENED COULOMB INTERACTION W^T_GG'(Q,\omega)
    
    DO  ig=1,nlf
      DO  jg=1,nlf
        wt(io,ig,jg)=czero
        DO  kg1=1,nlf
          DO  kg2=1,nlf
            wt(io,ig,jg)=wt(io,ig,jg)+ v(ig,kg1)*chi(kg1,kg2)*v(kg2,jg)
          END DO
        END DO
        wt(io,ig,jg)=v(ig,jg)+wt(io,ig,jg)
      END DO
    END DO
    
!                kraj nove petlje po omega
  END DO
  
  
!               ispis time ordered zasjenjene kulonske interakcije W_GG'^T(Q,\omega)
  dato='W_Qi'
  nord=INDEX(dato,'i', back =.false.)
  IF(iq < 10)THEN
    WRITE(dato(nord:nord),'(i1)')iq
  ELSE IF(iq >= 10.AND.iq < 100)THEN
    WRITE(dato(nord:nord+1),'(i2)')iq
  ELSE
    WRITE(dato(nord:nord+2),'(i3)')iq
  END IF
  
  OPEN(74,FILE=dato)
  DO  io=1,1
    o=(io-1)*domega
    WRITE(74,*)'omega=',o,'Hartree'
    WRITE(74,44)((wt(io,ig,jg),jg=1,nlf),ig=1,nlf)
  END DO
  CLOSE(74)
  44              FORMAT(10F15.5)
  
  DO  io=1,no-1
    DO  ig=1,nlf
      DO  jg=1,nlf
        s0(io,ig,jg)=-(1.0/pi)*imag(wt(io,ig,jg))
      END DO
    END DO
  END DO
  
  
  kks=zero
  skk=zero
  
  
!                new sum over omega
  DO  io=1,no-1
!                 print*,io
    oi=(io-1)*domega
    DO  ig=1,nlf
      DO  jg=1,nlf
        w1=0.0
        w2=0.0
!                static limit
        IF(io == 1)THEN
          DO  jo=2,no
            oj=(jo-1)*domega
            fact=domega/oj
            IF(jo == 2)fact=3.0/2.0
            IF(jo == no)fact=0.5*domega/oj
            w1=w1-fact*s0(jo,ig,jg)
          END DO
          w2=-w1
        ELSE IF(io == 2)THEN
          DO  jo=1,no
            oj=(jo-1)*domega
            IF(jo /= io)fact=domega/(oi-oj)
            IF(jo == 1)fact=1.0
            IF(jo == 2)fact=0.0
            IF(jo == 3)fact=-3.0/2.0
            IF(jo == no)fact=0.5*domega/(oi-oj)
            w1=w1+fact*s0(jo,ig,jg)
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            w2=w2+fact*s0(jo,ig,jg)
          END DO
        ELSE IF(io == (no-1))THEN
          DO  jo=1,no
            oj=(jo-1)*domega
            IF(jo /= io)fact=domega/(oi-oj)
            IF(jo == 1)fact=0.5*domega/(oi-oj)
            IF(jo == (no-2))fact=3.0/2.0
            IF(jo == (no-1))fact=0.0
            IF(jo == no)fact=-1.0
            w1=w1+fact*s0(jo,ig,jg)
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            w2=w2+fact*s0(jo,ig,jg)
          END DO
        ELSE
          DO  jo=1,no
            oj=(jo-1)*domega
            IF(jo /= io)fact=domega/(oi-oj)
            IF(jo == 1)fact=0.5*domega/(oi-oj)
            IF(jo == (io-1))fact=3.0/2.0
            IF(jo == io)fact=0.0
            IF(jo == (io+1))fact=-3.0/2.0
            IF(jo == no)fact=0.5*domega/(oi-oj)
            w1=w1+fact*s0(jo,ig,jg)
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            w2=w2+fact*s0(jo,ig,jg)
          END DO
        END IF
        
        imw=-pi*s0(io,ig,jg)
        gammap(ig,jg)=CMPLX(w1,imw)
        gammam(ig,jg)=CMPLX(-w2,0.0)
        IF(ig == 1.AND.jg == 1)w2kk=w2
        IF(ig == 1.AND.jg == 1.AND.io == 1)g0=gammap(1,1)
        
!                kraj po iG,jG
      END DO
    END DO
    
!                Provjera KK relacija
    wind=REAL(wt(io,1,1)-v(1,1))
    windkk=REAL(gammap(1,1))-w2kk
    fact=domega
    IF(io == 1.OR.io == no-1)fact=0.5*domega
    kks=kks+fact*(windkk-wind)*(windkk-wind)
    skk=skk+fact*wind*wind
    
    
!                kraj nove petlje po omega
  END DO
  CLOSE(74)
  
  
  dato='Kramers-Kron_Qi'
  nord=INDEX(dato,'i', back =.false.)
  IF(iq < 10)THEN
    WRITE(dato(nord:nord),'(i1)')iq
  ELSE IF(iq >= 10.AND.iq < 100)THEN
    WRITE(dato(nord:nord+1),'(i2)')iq
  ELSE
    WRITE(dato(nord:nord+2),'(i3)')iq
  END IF
  
  OPEN(33,FILE=dato)
  WRITE(33,88)'Wave vector (qx,qy,qz)=(',qx*gcar,qy*gcar, qz*gcar,') a.u.'
  WRITE(33,99)'|(qx,qy,qz)|=',absq*gcar,'a.u.'
  WRITE(33,*)'int(WindKK-Wind)^2 =  ',kks
  WRITE(33,*)'int(Wind)^2 =  ',skk
  WRITE(33,*)'****************************************'
  WRITE(33,*)'Kramers–Kronig relation relative error'
  WRITE(33,80)100.0*ABS(kks/skk), '%'
  78               FORMAT(a23,f10.5)
  79               FORMAT(a16,f10.5)
  80               FORMAT(5X,f7.2,a2)
  WRITE(33,*)'Usporedba Gamma i WT'
  WRITE(33,*)'real[Gamma(o=0,1,1)]=',REAL(g0)
  WRITE(33,*)'real[WT(o=0,1,1)]/2=', REAL(:: wt(1,1,1)-v(1,1))/2.0
  CLOSE(33)
  
  
!                kraj po q
END DO

END PROGRAM surface_loss

