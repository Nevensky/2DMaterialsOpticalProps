PROGRAM surface_current
 
! Code converted using TO_F90 by Alan Miller
! Date: 2020-06-17  Time: 13:00:14

!        NkI -number of wave vectors in irreducible B. zone
!        Ntot-total number of the mutually different wave vector-program generates this number
!        Nband-number of the bands
!        NG-total number of G vectors
!        NGd-number of coefficients CG shulod me less than minimum number of coefficients all over all evc.n files
!        nMPx*nMPy*nMPz-Monkhorest-Pack sampling
!        Efermi-Fermi energy
!        T-temperature in eV
!        alpha-Damping parameter in eV
!        Ecut-cutoff energy for crystal local field calculations
!        Vcell-unit-cell volume in a.u.^-3
!        a0-unit cell parameter in  a.u.^-1


use OMP_lib
use iotk_module
use ModPointR

implicit none

character(len=iotk_attlenx) :: attr1, attr2

! mod = 1 tenzor korelacijske funkcije, mod = 2 struja-struja tenzor (Kramers-Krroning)

integer :: ik,i,j,jk,it,lk,Ntot,iG0,Nsymm,iq, &
           io,jo, n,m,iG,R1,K1,R2,K2, Nlf,NG1,   &
           NG2,iG1,iG2,jG,kG, jump,   &
           iGfast,ikmin,kG1,kG2

integer :: frac, mod

character(len=3) :: pol, lf
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
real(kind=dp),    parameter :: aBohr = 0.5291772d0
real(kind=dp),    parameter :: alpha = 1.0/137.0
complex(kind=dp), parameter :: rone  = cmplx(1.0,0.0)
complex(kind=dp), parameter :: czero = cmplx(0.0,0.0)
complex(kind=dp), parameter :: ione  = cmplx(0.0,1.0)


!        skalars
REAL*8 :: a0,eps,kx,ky,kz,kqx,kqy,kqz,qgx,qgy,qgz,omin,omax,qx,  &
    qy,qz,kmin,domega,o,Efermi,k11,k22,k33,t,f1,f2,lor,de,gabs,  &
    kref,eref,ecut,gxx1,gyy1,gzz1,gxx2,gyy2,gzz2,fact,vcell,  &
    oi,oj,q,gcar,abohr,Nel,error,lossf,struja,pabs,  &
    gama_intra,gama_inter,expo1,expo2,alpha,cond,c0,strujay,  &
    strujaz,ImChi0,ReChi0,dGW
double precision pi,Hartree,planck,ev
PARAMETER(Efermi = 0.5554/Hartree,  &
    a0 = 5.9715,c0 = 29.8575,gcar = 2.0*pi/a0,  &
    eps = 1.0D-4,t = 0.025/Hartree,gama_intra = 0.025/Hartree,  &
    gama_inter = 0.025/Hartree,ev = 1.602176487D-19,ecut = 0.0,vcell = 922.0586)
DOUBLE COMPLEX :: ione,czero,rone,em

!        arrays
INTEGER :: Gfast,gi
REAL*8 :: ki,E,R,RI,k,ktot,g,glf,v,KC,f,kk

COMPLEX*8 :: MnmK1K2,Pi_dia,Pi_tot,MnmK1K22, &
             Qeff, &! efektivni naboj, matrica kad imam LFE
             S0

DIMENSION ki(3,NkI),E(NkI,Nband),R(48,3,3),RI(48,3,3),  &
    k(3,nk),ktot(3,nk),v(Nlfd,Nlfd),g(3,NG),  &
    glf(3,Nlfd),MnmK1K2(Nlfd), &
    ,S0(-no:no,Nlfd,Nlfd),Gfast(Nlfd*NGd),  &
    KC(3,3),gi(3),f(48,3),MnmK1K22(Nlfd),Pi_dia(Nlfd,Nlfd),  &
    Qeff(Nlfd,Nlfd),Pi_tot(Nlfd,Nlfd),kk(3,6)

CHARACTER (LEN = 100) :: bandn,bandm,nis,pathk1,pathk2,dato1,  &
    root,path,fajl,dato2,dato3,root1,root2, outdir
CHARACTER (LEN = 35) :: tag,buffer

COMPLEX,pointer, DIMENSION(:) :: c1,c2


rone = cmplx(1.0,0.0)
czero = cmplx(0.0,0.0)
ione = cmplx(0.0,1.0)

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
root1='/Users/nevensky/'
root2='MoS2_201X201'
outdir='/Users/Nevensky/tmp'
root = trim(root1)//trim(root2)




!             For free spectral function calculation put mod = 1
!             For current-current response function calculation put mod = 2
!             For x polarization put pol = 1
!             For y polarization put pol = 2
!             For z polarization put pol = 3
!             For mixed component yz put pol = 4
!             Crystal local field effects are included in z direction lf = 1
!             Crystal local field effects are included in x,y,z direction lf = 3

MOD = 1 ! racun struja-struja tenzora preko Kramers Kronings
lf = 1
pol = 'xx' ! polrizacija u x-smjeru

!             CORRELATION FUNCTIONS, CURRENT-CURRENT RESPONSE FUNCTIONS and
!             EFFECTIVE CHARGE CARRIERS MATRIX OUTPUTS

dato1 ='Corrfun'
dato2 ='Qeff'
dato3 = 'Pi_RPA_'//adjustl(trim(pol))
! if(pol == 'xx')dato3='Pi_RPA_xx'
! if(pol == 'yy')dato3='Pi_RPA_yy'
! if(pol == 'zz')dato3='Pi_RPA_zz'

!             scissors shift
dGW = 1.0/Hartree



jump = 1

omin = 1.0D-5
omax=(50.0/Hartree + omin)
domega=(omax-omin)/(no-1)





!      call FOR POINT GROUP TRANSFORMATIONS



call PointR(root,Nsymm,R,RI)




!           Upis valnih vektora iz irreducibilne Brillouinove
!           zone i pripadnih energijskih nivoa iz filea '*****.band'.
!           wave vectors are in cart.koord.

fajl='/MoS2.band'
path = trim(root)//trim(fajl)
OPEN(1,FILE = path)
do  ik = 1,NkI
  if(ik == 1)READ(1,*) nis
  READ (1,20)ki(1,ik),ki(2,ik),ki(3,ik)
  READ (1,10)(E(ik,i),i = 1,Nband)
end do
close(1)
10          FORMAT(10F8.4)
20          FORMAT(10X,f10.3,f10.3,f10.3)




do  ik = 1,NkI
  do  i = 1,Nband
    E(ik,i)=E(ik,i)/Hartree
    if(i >= nocc + 1)E(ik,i)=E(ik,i)+dGW
  end do
end do


!            generator 1.B.Z.
!            Dio programa koji pomocu operacija tockaste grupe i vektora iz
!            I.B.Z. generira sve (MEDJUSOBNO RAZLICITE!!!) v. vektore u 1.B.Z.
!            Ntot-Tot number of different points ''ktot'' inside 1.B.Z


jk = 0
ntot = 0
do  i = 1,Nsymm
  do  ik = 1,NkI
    it = 1
    jk = jk + 1
    do  n = 1,3
      k(n,jk)=0.0
      do  m = 1,3
        k(n,jk)=k(n,jk)+R(i,n,m)*ki(m,ik)
      end do
    end do
    if(jk > 1)then
      do  lk = 1,jk-1
        if(ABS(k(1,jk)-k(1,lk)) <= eps)then
          if(ABS(k(2,jk)-k(2,lk)) <= eps)then
            if(ABS(k(3,jk)-k(3,lk)) <= eps)then
              it = 2
            end if
          end if
        end if
      end do
    end if
    if(it == 1)then
      ntot = ntot + 1
      ktot(1,ntot)=k(1,jk)
      ktot(2,ntot)=k(2,jk)
      ktot(3,ntot)=k(3,jk)
    end if
  end do
end do

!             Checking 1BZ integration
Nel = 0
do  ik = 1,ntot
  kx = ktot(1,ik)
  ky = ktot(2,ik)
  kz = ktot(3,ik)
  do  n = 1,Nband
    if(n == 1)then
      it = 1
      if(ik <= NkI)then
        K1 = ik
        it = 2
      else
        do  i = 2,Nsymm
          k11 = RI(i,1,1)*kx + RI(i,1,2)*ky + RI(i,1,3)*kz
          k22 = RI(i,2,1)*kx + RI(i,2,2)*ky + RI(i,2,3)*kz
          k33 = RI(i,3,1)*kx + RI(i,3,2)*ky + RI(i,3,3)*kz
          do  j = 1,NkI
            if(abs(k11-ki(1,j)) <= eps)then
              if(abs(k22-ki(2,j)) <= eps)then
                if(abs(k33-ki(3,j)) <= eps)then
                  it = 2
                  K1 = j
                  GO TO 5022
                end if
              end if
            end if
          end do
        end do
      end if
      if(it == 1)then
        PRINT*,'Can not find wave vector K=',ik, 'in I.B.Z.'
        stop
      end if
      5022          CONTINUE
    end if
    if(E(K1,n) < Efermi)Nel = Nel + 1.0
  end do
end do
Nel = 2.0*Nel/ntot


!           KC transformation matrix from rec.cryst. axes to cart.koord.
!           if g' is vector in rec.cryst. axes then a = KC*a' is vector in cart. axes


fajl='/MoS2.sc.out'
path = trim(root)//trim(fajl)
tag='     reciprocal axes: (cart. coord.'
OPEN(1,FILE = path)
do  i = 1,100000
  READ(1,'(a)')buffer
  if(buffer == tag)then
    do  j = 1,3
      READ(1,70)KC(1,j),KC(2,j),KC(3,j)
    end do
    GO TO 998
  end if
end do
70          FORMAT(23X,3F10.3)
998         CONTINUE
close(1)


!           Reading the reciprocal vectors in crystal coordinates and transformation
!           in Cartezi cordinates.
call loadG(KC,G)


!            Reciprocal vectors for crystal local field effects calculations in array ''Glf(3,Nlf)''
call genGlf(lf,Ecut,NG,Gcar,G,Nlf,Nlfd,Glf)



!          IBZ   q LOOP STARTS HERE!!!

do  iq = qmin,qmax ! nq = 1 u optickom smo limesu, dakle ne treba nam do loop po q
  
!             searching for the smalest 'optical' q
call findMinQ(Ntot, ktot, qx, qy, qz)

!             Info file
  
  OPEN(55,FILE='Info')
  write(55,*)'***************General***********************'
  write(55,*)' Currently we calculate         ---->',dato1
  write(55,*)' Currently we calculate         ---->',dato2
  write(55,*)''
  write(55,*)'Number of point symmetry operation is',Nsymm
  ! if(frac == 0)write(55,*)'Fraction translation is not detected'
  ! if(frac == 1)write(55,*)'Fraction translation is detected'
  write(55,'(A25,3F10.4,A5)') 'Wave vector (qx,qy,qz)=(',qx*gcar,qy*gcar, qz*gcar,') a.u.'
  write(55,'(A25,F7.3,A5)') '|(qx,qy,qz)|=',SQRT(qx*qx + qy*qy + qz*qz)*gcar,'a.u.'
  if(lf == 'x') write(55,*) 'Local field effcts in z-dir'
  if(lf == 'xyz')write(55,*) 'Local field in all xyz-dir'
  write(55,*)'Number of local field vectors is',Nlf
  write(55,*)'Number of different K vectors in 1.B.Z. is',ntot
  write(55,*)'Number of K vectors in I.B.Z. is',NkI
  write(55,*)'Number of bands is               ',Nband
  write(55,'(A25,F7.3,A5)') 'Gama_intra is  ',gama_intra*Hartree*1000.0,'meV'
  write(55,'(A25,F7.3,A5)') 'Gama_inter is  ',gama_inter*Hartree*1000.0,'meV'
  write(55,'(A25,F7.3,A5)') 'Temperature is      ',t*Hartree*1000.0,'meV'
  write(55,*)''
  write(55,*)'-Im(Chi(io,G1,G2))/pi is in file---->',dato1
  write(55,*)' Qeff complex matrix is in file ---->',dato2
  write(55,*)' Pi_munu is in file            ---->',dato3
  write(55,*)''
  write(55,*)'************* Checking 1BZ integration*******'
  write(55,*)''
  write(55,'(A40,F8.4)')'Number of electrons(1BZ integration)=',Nel
  write(55,*)'Number of electrons(unit cell)=',NelQE
  error = abs((NelQE-Nel)/NelQE)
  write(55,'(A25,F7.3,A5)') 'Relative error=',error*100.0,'%'
  if(error > 0.05)then
    write(55,*)'WARRNING!!-1BZ INTEGRATION IS BAD!.'
  end if
  close(55)
  
  
  if(MOD == 2)GO TO 888
  
  S0 = cmplx(0.0,0.0) 
  Qeff = cmplx(0.0,0.0)
  
  OPEN(74,FILE = dato1)
  OPEN(75,FILE = dato2)
  
  
!              1.B.Z  LOOP STARTS HERE !!!!
  
  do  ik = 1,ntot
    
    OPEN(33,FILE='status')
    write(33,*)ik
    close(33)
!                print*,ik
    
    kx = ktot(1,ik)
    ky = ktot(2,ik)
    kz = ktot(3,ik)
    
!              trazenje (kx,ky,kz) u ireducibilnoj zoni
    
    it = 1
    if(ik <= NkI)then
      R1 = 1
      K1 = ik
      it = 2
    else
      do  i = 2,Nsymm
        k11 = RI(i,1,1)*kx + RI(i,1,2)*ky + RI(i,1,3)*kz
        k22 = RI(i,2,1)*kx + RI(i,2,2)*ky + RI(i,2,3)*kz
        k33 = RI(i,3,1)*kx + RI(i,3,2)*ky + RI(i,3,3)*kz
        do  j = 1,NkI
          if(abs(k11-ki(1,j)) <= eps)then
            if(abs(k22-ki(2,j)) <= eps)then
              if(abs(k33-ki(3,j)) <= eps)then
                it = 2
                R1 = i
                K1 = j
                GO TO 5222
              end if
            end if
          end if
        end do
      end do
    end if
    if(it == 1)then
      PRINT*,'Can not find wave vector K=',ik, 'in I.B.Z.'
      stop
    end if
    5222           CONTINUE
    
    
    it = 1
    kqx = kx + qx
    kqy = ky + qy
    kqz = kz + qz
    
!              trazenje (KQx,KQy) prvo u 1.B.Z a onda u I.B.Z.
    
    do  iG = 1,NG
      do  jk = 1,ntot
        if(abs(kqx-g(1,iG)-ktot(1,jk)) <= eps)then
          if(abs(kqy-g(2,iG)-ktot(2,jk)) <= eps)then
            if(abs(kqz-g(3,iG)-ktot(3,jk)) <= eps)then
              it = 2
              iG0 = iG
              do  i = 1,Nsymm
                k11 = RI(i,1,1)*ktot(1,jk)+RI(i,1,2)*ktot(2,jk)+  &
                    RI(i,1,3)*ktot(3,jk)
                k22 = RI(i,2,1)*ktot(1,jk)+RI(i,2,2)*ktot(2,jk)+  &
                    RI(i,2,3)*ktot(3,jk)
                k33 = RI(i,3,1)*ktot(1,jk)+RI(i,3,2)*ktot(2,jk)+  &
                    RI(i,3,3)*ktot(3,jk)
                do  j = 1,NkI
                  if(abs(k11-ki(1,j)) <= eps)then
                    if(abs(k22-ki(2,j)) <= eps)then
                      if(abs(k33-ki(3,j)) <= eps)then
                        it = 3
                        R2 = i
                        K2 = j
                        EXIT
                      end if
                    end if
                  end if
                end do
              end do
            end if
          end if
        end if
      end do
    end do
    2111                 CONTINUE
    
    
    if(it == 1)then
      PRINT*,'Can not find wave vector K + Q=',ik,'+',iq, 'in 1.B.Z.'
      stop
    else if(it == 2)then
      PRINT*,'Can not find wave vector K + Q=',ik,'+',iq, 'in I.B.Z.'
      stop
    end if
    
    
!              R1-integer, redni broj point operacije R1 u transformaciji ''K = R1*K1''.
!              K1-integer, redni broj valnog vektora K1 u transformaciji ''K = R1*K1''.
!              iG0 i R2-integeri, redni broj vektora reciprocne restke G0 i point operacije R2 u transformaciji ''K + Q = G0 + R2*K2''.
!              K2-integer, redni broj valnog vektora K2 u transformaciji  ''K + Q = G0 + R2*K2''.
    
    
!              petlje po vrpcama n i m
    
    
    
    do  n = 1,Nband
      do  m = 1,Nband
        
        
        
        expo1 = exp((E(K1,n)-Efermi)/t)
        expo2 = exp((E(K2,m)-Efermi)/t)
        f1 = 1.0/(expo1 + 1.0)
        f2 = 1.0/(expo2 + 1.0)
        f1 = f1-f2
        
        if((abs(f1) >= 1.0D-3).OR.(n == m))then
          
          
          call paths(outdir,K1,K2,n,m,pathk1,pathk2,bandn,bandm)
          
          
          
!         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za
!         fiksni K1,K2,n i m
          
!               Otvaranje atribute za INFO
          call iotk_open_read(10,pathk1)
          call iotk_scan_empty(10,"INFO",attr = attr)
          call iotk_scan_attr(attr,"igwx",NG1)
!               Alociranje polja C1
          allocate (c1(NG1))
!               Ucitavanje podataka iza evc.n
          call iotk_scan_dat(10,bandn,c1)
          call iotk_close_read(10)
!               Otvaranje atribute za INFO
          call iotk_open_read(10,pathk2)
          call iotk_scan_empty(10,"INFO",attr = attr)
          call iotk_scan_attr(attr,"igwx",NG2)
!               Alociranje polja C2
          allocate (c2(NG2))
!               Ucitavanje podataka iza evc.m
          call iotk_scan_dat(10,bandm,c2)
          call iotk_close_read(10)
          
!                Konstrukcija stupca matricnih elementa MnmK1K2(G)
          
          
          if(NGd > NG1)then
            write(*,*)'NGd is bigger than NG1=',NG1
            stop
          else if(NGd > NG2)then
            write(*,*)'NGd is bigger than NG2=',NG2
            stop
          end if
          
          
!                 matrix elements
          iGfast = 0
          do  iG = 1,Nlf
            MnmK1K2(iG)=cmplx(0.0,0.0)
            MnmK1K22(iG)=cmplx(0.0,0.0)
            do  iG1 = 1,NGd
              iGfast = iGfast + 1
              gxx1 = g(1,iG1)
              gyy1 = g(2,iG1)
              gzz1 = g(3,iG1)
              k11 = R(R1,1,1)*gxx1 + R(R1,1,2)*gyy1 + R(R1,1,3)*gzz1
              k22 = R(R1,2,1)*gxx1 + R(R1,2,2)*gyy1 + R(R1,2,3)*gzz1
              k33 = R(R1,3,1)*gxx1 + R(R1,3,2)*gyy1 + R(R1,3,3)*gzz1
              if (pol == 'xx') then
                struja=(qx + 2.0*kx + glf(1,iG)+2.0*k11)*gcar
              ELSEIF (pol == 'yy') then
                struja=(qy + 2.0*ky + glf(2,iG)+2.0*k22)*gcar
              ELSEIF (pol == 'zz')
                struja=(qz + 2.0*kz + glf(3,iG)+2.0*k33)*gcar
              ELSEIF(pol == 'yz') then
                strujay=(qy + 2.0*ky + glf(2,iG)+2.0*k22)*gcar
                strujaz=(qz + 2.0*kz + glf(3,iG)+2.0*k33)*gcar
              end if
              k11 = k11 + glf(1,iG)
              k22 = k22 + glf(2,iG)
              k33 = k33 + glf(3,iG)
              k11 = k11 + g(1,iG0)
              k22 = k22 + g(2,iG0)
              k33 = k33 + g(3,iG0)
              gxx1 = RI(R2,1,1)*k11 + RI(R2,1,2)*k22 + RI(R2,1,3)*k33
              gyy1 = RI(R2,2,1)*k11 + RI(R2,2,2)*k22 + RI(R2,2,3)*k33
              gzz1 = RI(R2,3,1)*k11 + RI(R2,3,2)*k22 + RI(R2,3,3)*k33
              if(jump == 1)then
                do  iG2 = 1,NG2
                  Gfast(iGfast)=NG2 + 1
                  gxx2 = g(1,iG2)
                  gyy2 = g(2,iG2)
                  gzz2 = g(3,iG2)
                  if(abs(gxx2-gxx1) < eps)then
                    if(abs(gyy2-gyy1) < eps)then
                      if(abs(gzz2-gzz1) < eps)then
                        Gfast(iGfast)=iG2
                        GO TO 1111
                      end if
                    end if
                  end if
                end do
              end if
              1111              CONTINUE
              iG2 = Gfast(iGfast)
              if(iG2 <= NG2)then

                ! ako je polarazcija je tipa xx, yy, zz 
                if(pol == 'xx' .or. pol== 'yy' .or. pol == 'zz')then
                  MnmK1K2(iG) = MnmK1K2(iG) +  &
                      0.5D0*CONJG(c1(iG1))*struja*c2(iG2)
                  MnmK1K22(iG) = MnmK1K2(iG) ! strujni vrhovi su isti
                else ! ako  su miksani xy, xz,... onda izvrsi ovo

                  MnmK1K2(iG) = MnmK1K2(iG)  +  &
                      0.5D0*CONJG(c1(iG1))*strujay*c2(iG2)
                  MnmK1K22(iG)=MnmK1K22(iG) +  &
                      0.5D0*CONJG(c1(iG1))*strujaz*c2(iG2)
                end if
              end if
!                 kraj po iG1
            end do
!                 kraj po C.L.F. iG
          end do

          jump = 2 ! za svaki valni vektor q i dani k zapamti Gfast i za svaku vrpcu preskaci taj postupak

          
          
          if(n /= m)then
!                  omega loop
            do  io = -no,no
              o = io*domega
              de = o + E(K1,n) - E(K2,m)

              ! gama_inter je sirina interband prijelaza
              lor = gama_inter/(de*de + gama_inter*gama_inter)
              ! reze repove lorentziana lijevo i desno, pazljivo, minimum 1.0d-3, preporuceno 1.0d-5
              if(abs(lor) >= 1.0D-5/gama_inter)then
                do  iG = 1,Nlf
                  do  jg = 1,Nlf
!                    -1/pi*ImChi_munu-for Kramers Kronig
                    S0(io,iG,jg)=S0(io,iG,jg)-  &
                        2.0*f1*lor*MnmK1K2(iG)*CONJG(MnmK1K22(jg))/  &
                        (pi*ntot*vcell)
                  end do
                end do
              end if
            end do
          else if(n == m.AND.abs(E(K1,n)-Efermi) <= 10.0*t)then
!                  Effective number of charge carriers (tensor)
            fact = expo1/((expo1 + 1.0D0)*(expo1 + 1.0D0))
            fact=-fact/t
            do  iG = 1,Nlf
              do  jg = 1,Nlf

                ! izracun intraband korelacijske funkcije
                Qeff(iG,jg)=Qeff(iG,jg)  &
                    +2.0*fact*MnmK1K2(iG)*CONJG(MnmK1K22(jg)) /(ntot*vcell)
              end do
            end do
!               end of intra/interband if loop
          end if
          
          deallocate(c1)
          deallocate(c2)
          
        end if
!               end of the occupation if loop
        
        
        
!              end of m do loop
      end do
      
      
!              end of n do loop
    end do
    jump = 1
    
    834            CONTINUE
!              ond of 1.B.Z do loop
  end do
  
  
!              WRITING CORFUN S0_\mu\nu
!              WRITTING Q_eff_\mu\nu
print *,'got here, line 722'
  do  io = -no,no ! opskurni razlog za prosirenje raspona frekvencija na negativne da se korektno izracuna spektar kristala koji nemaju centar inverzije

    o = io*domega
    write(74,*)'omega=',o,'Hartree'
    write(74,'(10F15.10)')((S0(io,iG,jG),jG = 1,Nlf),iG = 1,Nlf)
  end do
  write(75,'(10F15.10)')((Qeff(iG,jg),jG = 1,Nlf),iG = 1,Nlf)
  close(74)
  close(75)
  
  
  if(MOD == 1)GO TO 999
  
!               SECOND PART OF THE PROGRAM mod = 2
!               Calculation of the matrix '''Pi_\mu\nu'' by using matrix ''S0_\mu\nu(G,G')'' and Kramers-Krroning relations
  
  888             CONTINUE 

  
  OPEN(74,FILE = dato1)
  do  io=-no,no
    READ(74,*)nis
    READ(74,'(10F15.10)')((S0(io,iG,jg),jg = 1,Nlf),iG = 1,Nlf)
  end do
  close(74)
  
  OPEN(75,FILE = dato2)
  READ(75,'(10F15.10)')((Qeff(iG,jg),jg = 1,Nlf),iG = 1,Nlf)
  503             CONTINUE
  close(75)
  
  
  
  
!               Puting (qx,qy,qz) and Glf in cartezi coordinate
  
  qx = gcar*qx
  qy = gcar*qy
  qz = gcar*qz
  
  
  do  iG = 1,Nlf
    glf(1,iG)=gcar*glf(1,iG)
    glf(2,iG)=gcar*glf(2,iG)
    glf(3,iG)=gcar*glf(3,iG)
  end do
  
  
  
  
  OPEN(77,FILE = dato3)
!               new sum over omega
  do  io = 1,no-1
    PRINT*,io
    oi=(io-1)*domega
    do  iG = 1,Nlf
      do  jg = 1,Nlf
!                Real part of the response function Re(Chi)
        ReChi0 = 0.0
!                  static limit
        if(io == 1)then
          do  jo = 2,no
            oj=(jo-1)*domega
            fact = domega/oj
            if(jo == 2)fact = 3.0/2.0
            if(jo == no)fact = 0.5*domega/oj
            ReChi0 = ReChi0+ fact*(REAL(S0(-jo + 1,iG,jg))-REAL(S0(jo-1,iG,jg)))
          end do
        else if(io == 2)then
          do  jo = 1,no
            oj=(jo-1)*domega
            if(jo /= io)fact = domega/(oi-oj)
            if(jo == 1)fact = 1.0
            if(jo == 2)fact = 0.0
            if(jo == 3)fact=-3.0/2.0
            if(jo == no)fact = 0.5*domega/(oi-oj)
            ReChi0 = ReChi0 + fact*REAL(S0(jo-1,iG,jg))
            fact = domega/(oi + oj)
            if(jo == 1.OR.jo == no)fact = 0.5*domega/(oi + oj)
            ReChi0 = ReChi0 + fact*REAL(S0(-jo + 1,iG,jg))
          end do
        else if(io == (no-1))then
          do  jo = 1,no
            oj=(jo-1)*domega
            if(jo /= io)fact = domega/(oi-oj)
            if(jo == 1)fact = 0.5*domega/(oi-oj)
            if(jo == (no-2))fact = 3.0/2.0
            if(jo == (no-1))fact = 0.0
            if(jo == no)fact=-1.0
            ReChi0 = ReChi0 + fact*REAL(S0(jo-1,iG,jg))
            fact = domega/(oi + oj)
            if(jo == 1.OR.jo == no)fact = 0.5*domega/(oi + oj)
            ReChi0 = ReChi0 + fact*REAL(S0(-jo + 1,iG,jg))
          end do
        else
          do  jo = 1,no
            oj=(jo-1)*domega
            if(jo /= io)fact = domega/(oi-oj)
            if(jo == 1)fact = 0.5*domega/(oi-oj)
            if(jo == (io-1))fact = 3.0/2.0
            if(jo == io)fact = 0.0
            if(jo == (io + 1))fact=-3.0/2.0
            if(jo == no)fact = 0.5*domega/(oi-oj)
            ReChi0 = ReChi0 + fact*REAL(S0(jo-1,iG,jg))
            fact = domega/(oi + oj)
            if(jo == 1.OR.jo == no)fact = 0.5*domega/(oi + oj)
            ReChi0 = ReChi0 + fact*REAL(S0(-jo + 1,iG,jg))
          end do
        end if
        ReChi0 = ReChi0 + pi*aimag(S0(io-1,iG,jg))
        
!                Imaginary part of the response function Im(Chi)
        ImChi0 = 0.0
!                  static limit
        if(io == 1)then
          do  jo = 2,no
            oj=(jo-1)*domega
            fact = domega/oj
            if(jo == 2)fact = 3.0/2.0
            if(jo == no)fact = 0.5*domega/oj
            ImChi0 = ImChi0 + fact*(aimag(S0(-jo + 1,iG,jg))-aimag(S0(jo-1,iG,jg)))
          end do
        else if(io == 2)then
          do  jo = 1,no
            oj=(jo-1)*domega
            if(jo /= io)fact = domega/(oi-oj)
            if(jo == 1)fact = 1.0
            if(jo == 2)fact = 0.0
            if(jo == 3)fact=-3.0/2.0
            if(jo == no)fact = 0.5*domega/(oi-oj)
            ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jg))
            fact = domega/(oi + oj)
            if(jo == 1.OR.jo == no)fact = 0.5*domega/(oi + oj)
            ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jg))
          end do
        else if(io == (no-1))then
          do  jo = 1,no
            oj=(jo-1)*domega
            if(jo /= io)fact = domega/(oi-oj)
            if(jo == 1)fact = 0.5*domega/(oi-oj)
            if(jo == (no-2))fact = 3.0/2.0
            if(jo == (no-1))fact = 0.0
            if(jo == no)fact=-1.0
            ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jg))
            fact = domega/(oi + oj)
            if(jo == 1.OR.jo == no)fact = 0.5*domega/(oi + oj)
            ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jg))
          end do
        else
          do  jo = 1,no
            oj=(jo-1)*domega
            if(jo /= io)fact = domega/(oi-oj)
            if(jo == 1)fact = 0.5*domega/(oi-oj)
            if(jo == (io-1))fact = 3.0/2.0
            if(jo == io)fact = 0.0
            if(jo == (io + 1))fact=-3.0/2.0
            if(jo == no)fact = 0.5*domega/(oi-oj)
            ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jg))
            fact = domega/(oi + oj)
            if(jo == 1.OR.jo == no)fact = 0.5*domega/(oi + oj)
            ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jg))
          end do
        end if

        ! ovaj dio je razlicit od Sloss, S0 je kompleksno polje
        ImChi0 = ImChi0-pi*REAL(S0(io-1,iG,jg))
        
        if(io == 1)Pi_dia(iG,jg)=-cmplx(ReChi0,0.0) ! diamagnetski doprinos ??
        Pi_tot(iG,jg)=cmplx(ReChi0,ImChi0) ! Pi_RPA = PiDIJAMAGNETSKI + PiPARAMAGNETSKI
        Pi_tot(iG,jg)=Pi_tot(iG,jg)+Pi_dia(iG,jg)
        
!                   dodavanje intraband clana
        
        Pi_tot(iG,jg)=Pi_tot(iG,jg) +Qeff(iG,jg)*oi/(oi + cmplx(0.0,1.0)*gama_intra)
        
        
!                kraj po iG,jG
      end do
    end do
    
!                WRITTING TOTAL RESPONSE FUNCTION Pi_xx
    
    write(77,*)oi*Hartree,Pi_tot(1,1)
    
    
    
!                end new sum over omega
  end do
  close(77)
  
  999              CONTINUE
!                kraj po q
end do


contains
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

  subroutine genGlf(lf,Ecut,NG,Gcar,G,Nlf,Nlfd,Glf)
    ! Generate Reciprocal vectors for crystal local field 
    ! effects calculations in array Glf(3,Nlf)
  
    character(len=*), intent(in)    :: lf
    integer,          intent(in)    :: NG, Nlfd
    real(kind=dp),    intent(in)    :: Ecut
    real(kind=dp),    intent(in)    :: Gcar
    real(kind=dp),    intent(in)    :: G(:,:)
    integer,          intent(out)   :: Nlf
    real(kind=dp),    intent(inout) :: Glf(:,:)
  
    integer       :: iG
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
            end if
          end if
        end if
      end do
    elseif (lf == 'xyz') then
      ! local field efekti samo u svim smjerovima (xyz)
      do  iG = 1, NG
        Eref = Gcar**2*sum(G(1:3,iG)**2) / 2.0
        if (Eref <= Ecut) then
          Nlf = Nlf+1
          Glf(1:3,Nlf) = G(1:3,iG)
          end if
        end if
      end do
    end if
    if (Nlf > Nlfd) then
      print*,'Nlf is bigger than Nlfd'
      stop
    end if
  
  end subroutine genGlf

  subroutine loadG(KC,G)
    ! Reading the reciprocal vectors in crystal coordinates and transformation
    ! in Cartesian cordinates.
    implicit none
    real(kind=dp),  intent(in)    :: KC(3,3)
    ! integer,        intent(inout) :: parG(:) ! paritet svakog valnog vektora G
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
      ! parG(iG)=Gi(3)
    end do
  close(200)

  goto 5000
  200 write(*,*) 'error cant read file id 20, ist=',ist10
  201   write(*,*) '201 buffer1 read. Error reading line ',lno10+1,', iostat = ',ist11
  202   write(*,*) '202 buffer1 read. Number of lines read = ',lno10
  5000 continue 


end PROGRAM surface_current

