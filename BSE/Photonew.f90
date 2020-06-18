PROGRAM foton
 
! Code converted using TO_F90 by Alan Miller
! Date: 2020-06-17  Time: 13:15:04

!       program for the calculation of the current-current reponse  function in MoS2
!       XY-plane unit cell parameter  a0
!       Z-dir unit cell parameter c0
!       gama=e^{2}/(\hbar c) --  fine-structure constant


implicitnone

INTEGER :: n,no,nq,
& noladd, ! broj frekvencijskih koraka u ladder odzivnoj funkciji (polarizabilnosti)
& nl, ! broj slojeva (layera)
& pol

!       Nlfd je minimalno 2 zbog 2X2 blok matrice za p mod
PARAMETER(no=5001,noladd=401,nq=401,nl=30)

DOUBLE PRECISION :: zero,one,two,three,four,six,pi,eta,ev,  &
    hartree,abohr,planck,kb
REAL*8 a0,c0,gama,domega,h

PARAMETER(zero=0.0D0,one=1.0D0,two=2.0D0,three=3.0D0,  &
    four=4.0D0,six=6.0D0,pi=3.141592654D0,a0=5.9715,c0=29.8575D0,  &
    kb=1.3806503D-23,ev=1.602176487D-19,hartree=2.0D0*13.6056923D0,  &
    abohr=0.5291772D0,planck=6.626196D-34, eta=0.00001/hartree,gama=1.0D0/137.0D0)


!     ...scalars...
INTEGER :: i,j,k,m,iq,io
REAL*8 q,o,omax,dq,omin
COMPLEX*8 pi_ladder_down,pi_ladder_up,pi_rpa,oi,d0,  &
    beta,dxx,dyy,pis,pip,czero,ione,ieta,cone,dyz,dzz
!     ...arrays...
COMPLEX*8 pixx,pizz,pi0
DOUBLE COMPLEX eps,UNIT
DIMENSION pixx(no),pizz(no),pi0(nl,nl),eps(nl,nl), UNIT(nl,nl)
CHARACTER (LEN=100) :: nis,dato,root,path,fajl,rootrpa,  &
    rootladddown,rootladdup


czero=CMPLX(zero,zero)
cone=CMPLX(one,zero)
ione=CMPLX(zero,one)
ieta=CMPLX(zero,eta)

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
root='/home/vito/PROJECTS/MoS2-BSE/MoS2_51x51'

!             RPA IRREDUCIBLE POLARIZABILICI Chi_RPA FOLDER:
rootrpa='/home/vito/PROJECTS/MoS2-BSE/Chi_RPA'

!             LADDER IRREDUCIBLE POLARIZABILICI Chi_ladd_down FOLDER:
rootladddown='/home/vito/PROJECTS/MoS2-BSE/Chi_ladd_down'

!             LADDER IRREDUCIBLE POLARIZABILICI Chi_ladd_up FOLDER:
rootladdup='/home/vito/PROJECTS/MoS2-BSE/Chi_ladd_up'


!            udaljenost medju slojevima u Ang
h=3.0/abohr


!            for s-mode pol=1
!            for p-mode pol=2

pol=1


dq=0.00001D0
omin=1.0D-5
omax=(50.0/hartree+omin)
domega=(omax-omin)/(no-1)

fajl='/Pixx/Pi_RPA_xx'
path=trim(rootrpa)//trim(fajl)
OPEN(31,FILE=path)
fajl='/Pi_ladder_down_x'
path=trim(rootladddown)//trim(fajl)
OPEN(32,FILE=path)
fajl='/Pi_ladder_up_x'
path=trim(rootladdup)//trim(fajl)
OPEN(33,FILE=path)
DO  io=1,noladd
  READ(31,*)o,pi_rpa
  IF(io <= noladd)THEN
    READ(32,*)o,pi_ladder_down
    READ(33,*)o,pi_ladder_up
  END IF
  pixx(io)=pi_rpa+pi_ladder_down+pi_ladder_up
END DO
CLOSE(31)
CLOSE(32)
CLOSE(33)


fajl='/Pizz/Pi_RPA_zz'
path=trim(rootrpa)//trim(fajl)
OPEN(31,FILE=path)
fajl='/Pi_ladder_down_z'
path=trim(rootladddown)//trim(fajl)
OPEN(32,FILE=path)
fajl='/Pi_ladder_up_z'
path=trim(rootladdup)//trim(fajl)
OPEN(33,FILE=path)
DO  io=1,noladd
  READ(31,*)o,pi_rpa
  IF(io <= noladd)THEN
    READ(32,*)o,pi_ladder_down
    READ(33,*)o,pi_ladder_up
  END IF
  pizz(io)=pi_rpa+pi_ladder_down+pi_ladder_up
END DO
CLOSE(31)
CLOSE(32)
CLOSE(33)




!*******************************************************************
!              POCETAK PETLJI PO Q I OMEGA
!******************************************************************
OPEN(3,FILE='plot')

!              Q loop strts here
q_loop: DO iq=1,nq
  q=(iq-1)*dq
  
  PRINT*,iq
  
!              omega loop starts here
  omega_loop: DO io=3,noladd
    o=omin+(io-1)*domega
    oi=CMPLX(o,eta)
    beta=CMPLX(gama*gama*oi*oi-q*q)
    beta=SQRT(beta)
    
    
    
    dxx=2.0D0*pi*ione*c0*gama*gama/beta
    dyy=2.0D0*pi*ione*c0*beta/(oi*oi)
    dyz=-2.0*pi*ione*q*c0/(oi*oi)
    dzz=2.0*pi*ione*q*q*c0/(beta*oi*oi)
    
    
    
!             S-MOD
    IF(pol == 1)d0=dxx
!             P-MOD
    IF(pol == 2)d0=dyy
    
    
    eps=czero
    DO  i=1,nl
      DO  j=1,nl
        ! za multilayere nalazi efektivni epsilon = relativna permitivnost
        eps(i,j)= -pixx(io)*d0*EXP(ione*h*beta*ABS(REAL(i)-REAL(j)))
        UNIT(i,j)=czero
        pi0(i,j)=czero
      END DO
      eps(i,i)=cone+eps(i,i)
      UNIT(i,i)=cone
      pi0(i,i)=pixx(io)
    END DO
    
    CALL gjel(eps,nl,nl,UNIT,nl,nl) ! invertiranje (dio Dysonove jedn.)
    
    pip=czero
    DO  i=1,nl
      ! Dysonova jedn.
      pip=pip+eps(1,i)*pi0(i,1)  
    END DO
    
    
!           optical conductivity in units  pi*e^2/2h s-mode
    IF(iq == 1)WRITE(100,*)o*hartree,REAL(-ione*4.0*c0*pip/oi)
!            if(iq.eq.1)write(200,*)o*Hartree,real(-ione*4.0*c0*Pip/oi)
    
    
    
    
    WRITE(3,*)10.0*q/abohr, o*hartree,REAL(-ione*4.0*c0*pip/oi)
    
    
    
!            kraj po omega
  END DO omega_loop
  
  
  
!            end of q loop
END DO q_loop

CLOSE(3)


END PROGRAM foton





