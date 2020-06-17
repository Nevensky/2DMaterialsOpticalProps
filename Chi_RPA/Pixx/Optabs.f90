PROGRAM surface_current
 
! Code converted using TO_F90 by Alan Miller
! Date: 2020-06-17  Time: 13:00:14

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
!        Vcell-unit-cell volume in a.u.^-3
!        a0-unit cell parameter in  a.u.^-1


use iotk_module
implicitnone
CHARACTER (LEN=1) :: iotk_attlenx) :: attr
LOGICAL :: :: found

INTEGER :: nki,nband,ik,i,nk,j,jk,it,lk,ntot,ig0,nsim,  &
    ng,io,no,iq,nq,nmpx,nmpy,nmpz,n,m,ig,r1,k1,r2,k2,  &
    nlf,ng1,ng2,ngd,ig1,ig2,jg,nlfd,kg,jo,MOD,jump,  &
    igfast,ikmin,nelqe,pol,frac,lf,nocc
PARAMETER(nmpx=51,nmpy=51,nmpz=1,nki=6835,nband=60,nocc=9,  &
    nelqe=18,nk=48*nki,ngd=4000,ng=8000,no=5001,nq=1,nlfd=20)


!        skalars
REAL*8 a0,eps,kx,ky,kz,kqx,kqy,kqz,qgx,qgy,qgz,omin,omax,qx,  &
    qy,qz,kmin,domega,o,ef,k11,k22,k33,t,f1,f2,lor,de,gabs,  &
    kref,eref,ecut,gxx1,gyy1,gzz1,gxx2,gyy2,gzz2,fact,vcell,  &
    oi,oj,q,gcar,abohr,zero,nel,error,lossf,struja,pabs,  &
    gama_intra,gama_inter,expo1,expo2,gama,cond,c0,strujay,  &
    strujaz,imchi0,rechi0,dgw
doubleprecision pi,three,hartree,planck,ev
PARAMETER(hartree=2.0D0*13.6056923D0,ef=0.5554/hartree,  &
    a0=5.9715,c0=29.8575,pi=3.141592654D0,gcar=2.0*pi/a0,  &
    eps=1.0D-4,t=0.025/hartree,gama_intra=0.025/hartree,  &
    gama_inter=0.025/hartree,ev=1.602176487D-19,gama=1.0/137.0,  &
    planck=6.626196D-34,ecut=0.0,vcell=922.0586, abohr=0.5291772D0,zero=0.0)
DOUBLE COMPLEX ione,czero,rone,em

!        arrays
INTEGER :: gfast,gi
REAL*8 ki,e,r,ri,k,ktot,g,glf,v,kc,f,kk
DOUBLE COMPLEX UNIT,epsilon
COMPLEX*8 mnmk1k21,pi_dia,pi_tot,mnmk1k22,qeff,s0
DIMENSION ki(3,nki),e(nki,nband),r(48,3,3),ri(48,3,3),  &
    k(3,nk),ktot(3,nk),v(nlfd,nlfd),g(3,ng),  &
    glf(3,nlfd),mnmk1k21(nlfd),UNIT(nlfd,nlfd),  &
    epsilon(nlfd,nlfd),s0(-no:no,nlfd,nlfd),gfast(nlfd*ngd),  &
    kc(3,3),gi(3),f(48,3),mnmk1k22(nlfd),pi_dia(nlfd,nlfd),  &
    qeff(nlfd,nlfd),pi_tot(nlfd,nlfd),kk(3,6)

CHARACTER (LEN=100) :: bandn,bandm,nis,pathk1,pathk2,dato1,  &
    root,path,fajl,dato2,dato3,root1,root2
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



!        QUANTUM ESSPRESSO IMPUTS:
root1='/home/vito/PROJECTS/MoS2-BSE'
root2='/MoS2_201X201/'
root=trim(root1)//trim(root2)




!             For free spectral function calculation put mod=1
!             For current-current response function calculation put mod=2
!             For x polarization put pol=1
!             For y polarization put pol=2
!             For z polarization put pol=3
!             For mixed component yz put pol=4
!             Crystal local field effects are included in z direction lf=1
!             Crystal local field effects are included in x,y,z direction lf=3

MOD=2
lf=1
pol=1

!             CORRELATION FUNCTIONS, CURRENT-CURRENT RESPONSE FUNCTIONS and
!             EFFECTIVE CHARGE CARRIERS MATRIX OUTPUTS

dato1='Corrfun'
dato2='Qeff'
IF(pol == 1)dato3='Pi_RPA_xx'
IF(pol == 2)dato3='Pi_RPA_yy'
IF(pol == 3)dato3='Pi_RPA_zz'

!             scissors shift
dgw=1.0/hartree



jump=1
three=3.0D0
omin=1.0D-5
omax=(50.0/hartree+omin)
domega=(omax-omin)/(no-1)





!      CALL FOR POINT GROUP TRANSFORMATIONS


CALL pointr(root,nsim,r,ri,f)



!           Upis valnih vektora iz irreducibilne Brillouinove
!           zone i pripadnih energijskih nivoa iz filea '*****.band'.
!           wave vectors are in cart.koord.

fajl='/MoS2.band'
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
    IF(i >= nocc+1)e(ik,i)=e(ik,i)+dgw
  END DO
END DO


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
nel=2.0*nel/ntot

!             Checking the existence of fraction translation operations
frac=0
DO   i=1,nsim
  IF(f(i,1) /= zero.OR.f(i,2) /= zero.OR.f(i,3) /= zero)THEN
    frac=1
  END IF
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
  DO   n=1,3
    g(n,ig)=zero
    DO   m=1,3
      g(n,ig)=g(n,ig)+kc(n,m)*DBLE(gi(m))
    END DO
  END DO
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
    END IF
  END DO
END IF
IF(nlf > nlfd)THEN
  PRINT*,'Nlf is bigger than Nlfd'
  STOP
END IF



!             q LOOP STARTS HERE!!!

DO  iq=nq,nq
  
!             searching for the smalest 'optical' q
  kmin=1.0
  DO  i=1,ntot
    kref=SQRT(ktot(1,i)*ktot(1,i)+ ktot(2,i)*ktot(2,i)+ktot(3,i)*ktot(3,i))
    IF(kref == zero)GO TO 970
    IF(kref < kmin)THEN
      kmin=kref
      ikmin=i
    END IF
    970           CONTINUE
  END DO
  
  qx=(iq-1)*ktot(1,ikmin)
  qy=(iq-1)*ktot(2,ikmin)
  qz=(iq-1)*ktot(3,ikmin)
  
!             Info file
  
  OPEN(55,FILE='Info')
  WRITE(55,*)'***************General***********************'
  WRITE(55,*)' Currently we calculate         ---->',dato1
  WRITE(55,*)' Currently we calculate         ---->',dato2
  WRITE(55,*)''
  WRITE(55,*)'Number of point symmetry operation is',nsim
  IF(frac == 0)WRITE(55,*)'Fraction translation is not detected'
  IF(frac == 1)WRITE(55,*)'Fraction translation is detected'
  WRITE(55,88)'Wave vector (qx,qy,qz)=(',qx*gcar,qy*gcar, qz*gcar,') a.u.'
  WRITE(55,99)'|(qx,qy,qz)|=',SQRT(qx*qx+qy*qy+qz*qz)*gcar,'a.u.'
  IF(lf == 1)WRITE(55,*)'Local field effcts in z-dir'
  IF(lf == 3)WRITE(55,*)'Local field in all xyz-dir'
  WRITE(55,*)'Number of local field vectors is',nlf
  WRITE(55,*)'Number of different K vectors in 1.B.Z. is',ntot
  WRITE(55,*)'Number of K vectors in I.B.Z. is',nki
  WRITE(55,*)'Number of bands is               ',nband
  WRITE(55,99)'Gama_intra is  ',gama_intra*hartree*1000.0,'meV'
  WRITE(55,99)'Gama_inter is  ',gama_inter*hartree*1000.0,'meV'
  WRITE(55,99)'Temperature is      ',t*hartree*1000.0,'meV'
  WRITE(55,*)''
  WRITE(55,*)'-Im(Chi(io,G1,G2))/pi is in file---->',dato1
  WRITE(55,*)' Qeff complex matrix is in file ---->',dato2
  WRITE(55,*)' Pi_munu is in file            ---->',dato3
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
  88     FORMAT(a25,3F10.4,a5)
  99     FORMAT(a25,f7.3,a5)
  12     FORMAT(a40,f8.4)
  
  
  IF(MOD == 2)GO TO 888
  
  
  
  s0=czero
  qeff=czero
  
  OPEN(74,FILE=dato1)
  OPEN(75,FILE=dato2)
  
  
!              1.B.Z  LOOP STARTS HERE !!!!
  
  DO  ik=1,ntot
    
    OPEN(33,FILE='status')
    WRITE(33,*)ik
    CLOSE(33)
!                print*,ik
    
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
                r1=i
                k1=j
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
    
!              trazenje (KQx,KQy) prvo u 1.B.Z a onda u I.B.Z.
    
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
    
    
    
    DO  n=1,nband
      DO  m=1,nband
        
        
        
        expo1=EXP((e(k1,n)-ef)/t)
        expo2=EXP((e(k2,m)-ef)/t)
        f1=1.0/(expo1+1.0)
        f2=1.0/(expo2+1.0)
        f1=f1-f2
        
        IF((DABS(f1) >= 1.0D-3).OR.(n == m))THEN
          
          
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
          
          
!                 matrix elements
          igfast=0
          DO  ig=1,nlf
            mnmk1k21(ig)=czero
            mnmk1k22(ig)=czero
            DO  ig1=1,ngd
              igfast=igfast+1
              gxx1=g(1,ig1)
              gyy1=g(2,ig1)
              gzz1=g(3,ig1)
              k11=r(r1,1,1)*gxx1+r(r1,1,2)*gyy1+r(r1,1,3)*gzz1
              k22=r(r1,2,1)*gxx1+r(r1,2,2)*gyy1+r(r1,2,3)*gzz1
              k33=r(r1,3,1)*gxx1+r(r1,3,2)*gyy1+r(r1,3,3)*gzz1
              IF(pol == 1)struja=(qx+2.0*kx+glf(1,ig)+2.0*k11)*gcar
              IF(pol == 2)struja=(qy+2.0*ky+glf(2,ig)+2.0*k22)*gcar
              IF(pol == 3)struja=(qz+2.0*kz+glf(3,ig)+2.0*k33)*gcar
              IF(pol == 4)THEN
                strujay=(qy+2.0*ky+glf(2,ig)+2.0*k22)*gcar
                strujaz=(qz+2.0*kz+glf(3,ig)+2.0*k33)*gcar
              END IF
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
                IF(pol /= 4)THEN
                  mnmk1k21(ig)=mnmk1k21(ig)+  &
                      0.5D0*CONJG(c1(ig1))*struja*c2(ig2)
                  mnmk1k22(ig)=mnmk1k21(ig)
                ELSE
                  mnmk1k21(ig)=mnmk1k21(ig)+  &
                      0.5D0*CONJG(c1(ig1))*strujay*c2(ig2)
                  mnmk1k22(ig)=mnmk1k22(ig)+  &
                      0.5D0*CONJG(c1(ig1))*strujaz*c2(ig2)
                END IF
              END IF
!                 kraj po iG1
            END DO
!                 kraj po C.L.F. iG
          END DO
          jump=2
          
          
          IF(n /= m)THEN
!                  omega loop
            DO  io=-no,no
              o=io*domega
              de=o+e(k1,n)-e(k2,m)
              lor=gama_inter/(de*de+gama_inter*gama_inter)
              IF(DABS(lor) >= 1.0D-5/gama_inter)THEN
                DO  ig=1,nlf
                  DO  jg=1,nlf
!                    -1/pi*ImChi_munu-for Kramers Kronig
                    s0(io,ig,jg)=s0(io,ig,jg)-  &
                        2.0*f1*lor*mnmk1k21(ig)*CONJG(mnmk1k22(jg))/  &
                        (pi*ntot*vcell)
                  END DO
                END DO
              END IF
            END DO
          ELSE IF(n == m.AND.DABS(e(k1,n)-ef) <= 10.0*t)THEN
!                  Effective number of charge carriers (tensor)
            fact=expo1/((expo1+1.0D0)*(expo1+1.0D0))
            fact=-fact/t
            DO  ig=1,nlf
              DO  jg=1,nlf
                qeff(ig,jg)=qeff(ig,jg)  &
                    +2.0*fact*mnmk1k21(ig)*CONJG(mnmk1k22(jg)) /(ntot*vcell)
              END DO
            END DO
!               end of intra/interband if loop
          END IF
          
          deallocate(c1)
          deallocate(c2)
          
        END IF
!               end of the occupation if loop
        
        
        
!              end of m do loop
      END DO
      
      
!              end of n do loop
    END DO
    jump=1
    
    834            CONTINUE
!              ond of 1.B.Z do loop
  END DO
  
  
!              WRITING CORFUN S0_\mu\nu
!              WRITTING Q_eff_\mu\nu
  DO  io=-no,no
    o=io*domega
    WRITE(74,*)'omega=',o,'Hartree'
    WRITE(74,44)((s0(io,ig,jg),jg=1,nlf),ig=1,nlf)
  END DO
  WRITE(75,44)((qeff(ig,jg),jg=1,nlf),ig=1,nlf)
  44             FORMAT(10F15.10)
  CLOSE(74)
  CLOSE(75)
  
  
  IF(MOD == 1)GO TO 999
  
!               SECOND PART OF THE PROGRAM mod=2
!               Calculation of the matrix '''Pi_\mu\nu'' by using matrix ''S0_\mu\nu(G,G')'' and Kramers-Krroning relations
  
  888             CONTINUE
  
  OPEN(74,FILE=dato1)
  DO  io=-no,no
    READ(74,*)nis
    READ(74,44)((s0(io,ig,jg),jg=1,nlf),ig=1,nlf)
  END DO
  CLOSE(74)
  
  OPEN(75,FILE=dato2)
  READ(75,44)((qeff(ig,jg),jg=1,nlf),ig=1,nlf)
  503             CONTINUE
  CLOSE(75)
  
  
  
  
!               Puting (qx,qy,qz) and Glf in cartezi coordinate
  
  qx=gcar*qx
  qy=gcar*qy
  qz=gcar*qz
  
  
  DO  ig=1,nlf
    glf(1,ig)=gcar*glf(1,ig)
    glf(2,ig)=gcar*glf(2,ig)
    glf(3,ig)=gcar*glf(3,ig)
  END DO
  
  
  
  
  OPEN(77,FILE=dato3)
!               new sum over omega
  DO  io=1,no-1
    PRINT*,io
    oi=(io-1)*domega
    DO  ig=1,nlf
      DO  jg=1,nlf
!                Real part of the response function Re(Chi)
        rechi0=0.0
!                  static limit
        IF(io == 1)THEN
          DO  jo=2,no
            oj=(jo-1)*domega
            fact=domega/oj
            IF(jo == 2)fact=3.0/2.0
            IF(jo == no)fact=0.5*domega/oj
            rechi0=rechi0+ fact*(REAL(s0(-jo+1,ig,jg))-REAL(s0(jo-1,ig,jg)))
          END DO
        ELSE IF(io == 2)THEN
          DO  jo=1,no
            oj=(jo-1)*domega
            IF(jo /= io)fact=domega/(oi-oj)
            IF(jo == 1)fact=1.0
            IF(jo == 2)fact=0.0
            IF(jo == 3)fact=-3.0/2.0
            IF(jo == no)fact=0.5*domega/(oi-oj)
            rechi0=rechi0+fact*REAL(s0(jo-1,ig,jg))
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            rechi0=rechi0+fact*REAL(s0(-jo+1,ig,jg))
          END DO
        ELSE IF(io == (no-1))THEN
          DO  jo=1,no
            oj=(jo-1)*domega
            IF(jo /= io)fact=domega/(oi-oj)
            IF(jo == 1)fact=0.5*domega/(oi-oj)
            IF(jo == (no-2))fact=3.0/2.0
            IF(jo == (no-1))fact=0.0
            IF(jo == no)fact=-1.0
            rechi0=rechi0+fact*REAL(s0(jo-1,ig,jg))
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            rechi0=rechi0+fact*REAL(s0(-jo+1,ig,jg))
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
            rechi0=rechi0+fact*REAL(s0(jo-1,ig,jg))
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            rechi0=rechi0+fact*REAL(s0(-jo+1,ig,jg))
          END DO
        END IF
        rechi0=rechi0+pi*imag(s0(io-1,ig,jg))
        
!                Imaginary part of the response function Im(Chi)
        imchi0=0.0
!                  static limit
        IF(io == 1)THEN
          DO  jo=2,no
            oj=(jo-1)*domega
            fact=domega/oj
            IF(jo == 2)fact=3.0/2.0
            IF(jo == no)fact=0.5*domega/oj
            imchi0=imchi0+ fact*(imag(s0(-jo+1,ig,jg))-imag(s0(jo-1,ig,jg)))
          END DO
        ELSE IF(io == 2)THEN
          DO  jo=1,no
            oj=(jo-1)*domega
            IF(jo /= io)fact=domega/(oi-oj)
            IF(jo == 1)fact=1.0
            IF(jo == 2)fact=0.0
            IF(jo == 3)fact=-3.0/2.0
            IF(jo == no)fact=0.5*domega/(oi-oj)
            imchi0=imchi0+fact*imag(s0(jo-1,ig,jg))
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            imchi0=imchi0+fact*imag(s0(-jo+1,ig,jg))
          END DO
        ELSE IF(io == (no-1))THEN
          DO  jo=1,no
            oj=(jo-1)*domega
            IF(jo /= io)fact=domega/(oi-oj)
            IF(jo == 1)fact=0.5*domega/(oi-oj)
            IF(jo == (no-2))fact=3.0/2.0
            IF(jo == (no-1))fact=0.0
            IF(jo == no)fact=-1.0
            imchi0=imchi0+fact*imag(s0(jo-1,ig,jg))
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            imchi0=imchi0+fact*imag(s0(-jo+1,ig,jg))
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
            imchi0=imchi0+fact*imag(s0(jo-1,ig,jg))
            fact=domega/(oi+oj)
            IF(jo == 1.OR.jo == no)fact=0.5*domega/(oi+oj)
            imchi0=imchi0+fact*imag(s0(-jo+1,ig,jg))
          END DO
        END IF
        imchi0=imchi0-pi*REAL(s0(io-1,ig,jg))
        
        IF(io == 1)pi_dia(ig,jg)=-CMPLX(rechi0,zero)
        pi_tot(ig,jg)=CMPLX(rechi0,imchi0)
        pi_tot(ig,jg)=pi_tot(ig,jg)+pi_dia(ig,jg)
        
!                   dodavanje intraband clana
        
        pi_tot(ig,jg)=pi_tot(ig,jg) +qeff(ig,jg)*oi/(oi+ione*gama_intra)
        
        
!                kraj po iG,jG
      END DO
    END DO
    
!                WRITTING TOTAL RESPONSE FUNCTION Pi_xx
    
    WRITE(77,*)oi*hartree,pi_tot(1,1)
    
    
    
!                end new sum over omega
  END DO
  CLOSE(77)
  
  999              CONTINUE
!                kraj po q
END DO
END PROGRAM surface_current

