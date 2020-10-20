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

! calc = 1 tenzor korelacijske funkcije, calc = 2 struja-struja tenzor (Kramers-Krroning)

integer :: ik,i,j,jk,it,lk,Ntot,iG0,Nsymm,iq, &
           io,jo, n,m,iG,R1,K1,R2,K2, Nlf,NG1,   &
           NG2,iG1,iG2,jG,kG, jump,   &
           iGfast,ikmin,kG1,kG2

integer :: frac, calc

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
real(kind=dp),    parameter :: pi      = 4.D0*atan(1.D0)
real(kind=dp),    parameter :: eV      = 1.602176487D-19
real(kind=dp),    parameter :: Hartree = 2.0D0*13.6056923D0
real(kind=dp),    parameter :: Planck  = 6.626196D-34
real(kind=dp),    parameter :: aBohr   = 0.5291772d0
real(kind=dp),    parameter :: alpha   = 1.0/137.0
complex(kind=dp), parameter :: rone    = cmplx(1.0,0.0)
complex(kind=dp), parameter :: czero   = cmplx(0.0,0.0)
complex(kind=dp), parameter :: ione    = cmplx(0.0,1.0)


!        skalars
real*8 :: a0,eps,
kx,ky,kz,KQx,KQy,KQz,qgx,qgy,qgz,
omin,omax,
qx, qy,qz, kmin,domega,o,
Efermi,
k11,k22,k33,
T,
f1,f2,df,Lor,dE,Gabs,  &
    kref,eref,ecut,Gxx1,Gyy1,Gzz1,Gxx2,Gyy2,Gzz2,fact,Vcell,  &
    oi,oj,q,Gcar,abohr,Nel,error,lossf,struja,pabs,  &
    Gamma_intra,Gamma_inter,expo1,expo2,alpha,cond,c0,struja_y,  &
    struja_z,ImChi0,ReChi0,dGW
double precision pi,Hartree,planck,ev
PARAMETER(Efermi = 0.5554/Hartree,  &
    a0 = 5.9715,c0 = 29.8575,Gcar = 2.0*pi/a0,  &
    eps = 1.0D-4,T = 0.025/Hartree,
    Gamma_intra = 0.025/Hartree,  &
    Gamma_inter = 0.025/Hartree,
    ev = 1.602176487D-19,
    ecut = 0.0,
    Vcell = 922.0586)

Gamma_intra = Gamma_intra/Hartree
Gamma_inter = Gamma_inter/Hartree

DOUBLE COMPLEX :: ione,czero,rone

!        arrays
INTEGER :: Gfast,Gi
real*8 :: ki,E,R,RI,k,ktot,G,Glf,v,KC,kk

COMPLEX*8 :: MnmK1K2,Pi_dia,Pi_tot,MnmK1K22, &
             Qeff, &! efektivni naboj, matrica kad imam LFE
             S0

DIMENSION ki(3,NkI),E(NkI,Nband),R(48,3,3),RI(48,3,3),  &
    k(3,nk),ktot(3,nk),G(3,NG),  &
    Glf(3,Nlfd),MnmK1K2(Nlfd), &
    ,S0(-no:no,Nlfd,Nlfd),Gfast(Nlfd*NGd),  &
    KC(3,3),Gi(3),f(48,3),MnmK1K22(Nlfd),Pi_dia(Nlfd,Nlfd),  &
    Qeff(Nlfd,Nlfd),Pi_tot(Nlfd,Nlfd),kk(3,6)

CHARACTER (LEN = 100) :: bandn,bandm,nis,pathk1,pathk2,dato1,  &
    root,path,fajl,dato2,dato3,root1,root2, outdir
CHARACTER (LEN = 35) :: tag,buffer

COMPLEX,pointer, DIMENSION(:) :: C1,C2


rone = cmplx(1.0,0.0)
czero = cmplx(0.0,0.0)
ione = cmplx(0.0,1.0)

!        QUANTUM ESSPRESSO IMPUTS:
root1='/Users/nevensky/'
root2='MoS2_201X201'
outdir='/Users/Nevensky/tmp'
root = trim(root1)//trim(root2)




!  For free spectral function calculation put calc = 1
!  For current-current response function calculation put calc = 2
!  For mixed component yz put pol = 4
!  Crystal local field effects are included in z direction lf = 1
!  Crystal local field effects are included in x,y,z direction lf = 3

calc = 1 
lf = 1
pol = 'xx' ! polrizacija u x-smjeru

! CORRELATION FUNCTIONS, CURRENT-CURRENT RESPONSE FUNCTIONS and
! EFFECTIVE CHARGE CARRIERS MATRIX OUTPUTS

dato1 ='Corrfun'
dato2 ='Qeff'
dato3 = 'Pi_RPA_'//adjustl(trim(pol))


! scissors shift
dGW = dGW/Hartree


! nevne debug - prebaceno u config.in
! jump = 1
! omin = 1.0D-5
! omax = (50.0/Hartree + omin)


domega = (omax-omin)/(no-1)



! POINT GROUP TRANSFORMATIONS
path = trim(rundir)//trim(scf_file)
call PointR(root,Nsymm,R,RI)
print *,"status: PointR done."



! Upis valnih vektora iz irreducibilne Brillouinove
! zone i pripadnih enerGijskih nivoa iz filea '*****.band'.
! wave vectors are in cart.koord.

path=trim(rundir)//trim(band_file)
call loadkIandE(path, NkI, Nband, Nocc, kI, E)
print *,"status: kI and E loaded."


! generator 1.B.Z.
! Dio programa koji pomocu operacija tockaste grupe i vektora iz
! I.B.Z. generira sve (MEDJUSOBNO RAZLICITE!!!) v. vektore u 1.B.Z.
! Ntot-Tot number of different points ''ktot'' inside 1.B.Z
call genFBZ(Nk,NkI,Nsymm,eps,kI,R,Ntot,ktot)

! Checking 1BZ integration
call checkFBZintegration(Nband,NkI,Nsymm,Ntot,eps,kI,RI,Efermi,E,NelQE,Nel)


! KC transformation matrix from rec.cryst. axes to cart.koord.
! if G' is vector in rec.cryst. axes then a = KC*a' is vector in cart. axes
call loadKC(path,KC)


! Reading the reciprocal vectors in crystal coordinates and transformation
! in Cartesian cordinates.
call loadG(KC,G)


! Reciprocal vectors for crystal local field effects calculations in array ''Glf(3,Nlf)''
call genGlf(lf,Ecut,NG,Gcar,G,Nlf,Nlfd,Glf)



! IBZ   q LOOP STARTS HERE!!!

q_loop: do  iq = qmin,qmax ! nq = 1 u optickom smo limesu, dakle ne treba nam do loop po q
  
  ! searching for the smalest 'optical' q
  call findMinQ(Ntot, ktot, qx, qy, qz)

  ! Info file
  call writeInfo(lf, pol, qx, qy, qz, Gcar, Nsymm, Nlf, Ntot, NkI, Nband, T, Nel, NelQE,Gamma_intra, Gamma_inter)
  
  
  if (calc == 2) GO TO 888
  
  S0(-no:no,1:Nlf,1:Nlf) = cmplx(0.0,0.0)
  Qeff(1:Nlf,1:Nlf) = cmplx(0.0,0.0)
  
  
  ! 1.B.Z  LOOP STARTS HERE !!!!
  !$omp do
  do  ik = 1,Ntot
    ! vito debug
    ! open(33,FILE='status')
    ! write(33,*)ik
    ! close(33)
    
    kx = ktot(1,ik)
    ky = ktot(2,ik)
    kz = ktot(3,ik)
    
  ! trazenje (kx,ky,kz) u ireducibilnoj zoni
  call findKinIBZ(ik, NkI, Nsymm, eps, kx, ky, kz, RI, kI, R1, K1)
    
    KQx = kx + qx
    KQy = ky + qy
    KQz = kz + qz
    
  !  trazenje (KQx,KQy) prvo u 1.B.Z a onda u I.B.Z.
  call findKQinBZ(KQx, KQy, KQz, eps, Nsymm, NkI, Ntot, NG, ktot, kI, RI, G, iG0, R2, K2)
      
  ! R1-integer, redni broj point operacije R1 u transformaciji ''K = R1*K1''.
  ! K1-integer, redni broj valnog vektora K1 u transformaciji ''K = R1*K1''.
  ! iG0 i R2-integeri, redni broj vektora reciprocne restke G0 i point operacije R2 u transformaciji ''K + Q = G0 + R2*K2''.
  ! K2-integer, redni broj valnog vektora K2 u transformaciji  ''K + Q = G0 + R2*K2''.
        
    
    do  n = 1,Nband
      expo1 = exp((E(K1,n)-Efermi)/T)
      f1 = 1.0/(expo1 + 1.0)

      call paths(savedir,K1,n,pathk1,bandn)

      !$omp critical(loadCs_1)
      iuni1 = 20 + 2*thread_id + ik*100000 + n*100
      call paths(savedir,K1,n,pathk1,bandn)
      ! print *,'pathk1',pathk1
      call iotk_open_read(iuni1,pathk1)
      call iotk_scan_empty(iuni1,"INFO",attr=attr1) ! Otvaranje atribute za INFO
      call iotk_scan_attr(attr1,"igwx",NG1)
      
      allocate(C1(NG1))                            ! Alociranje polja C1
      call iotk_scan_dat(iuni1,bandn,C1)           ! Ucitavanje podataka iza evc.n
      call iotk_close_read(iuni1)
      !$omp end critical(loadCs_1)
      
      if (NGd > NG1) then
        write(*,*) 'NGd is bigger than NG1=',NG1
        stop
      end if

      do  m = 1,Nband

        expo2 = exp((E(K2,m)-Efermi)/T)        
        f2 = 1.0/(expo2 + 1.0)
        df = f1-f2
        
        occupation_if: if ( (abs(df) >= df_cut) .or. (n == m) ) then  
          ! call paths(outdir,K1,K2,n,m,pathk1,pathk2,bandn,bandm)
          call paths(savedir,K2,m,pathk2,bandm)
          
          !$omp critical(loadCs_2)
          iuni2 = 21 + (2*thread_id+1) + (2*ik+1)*100000 + (2*m+1)*100
          call pathsB(savedir,K2,m,pathk2,bandm)
          ! print *,'pathk2',pathk2
          call iotk_open_read(iuni2,pathk2)
          call iotk_scan_empty(iuni2,"INFO",attr=attr2)
          call iotk_scan_attr(attr2,"igwx",NG2)
          
          allocate(C2(NG2))                            ! Alociranje polja C1
          call iotk_scan_dat(iuni2,bandm,C2)           ! Ucitavanje podataka iza evc.n
          call iotk_close_read(iuni2)
          !$omp end critical(loadCs_2)
          
          
          if (NGd > NG2) then
            write(*,*) 'NGd is bigger than NG2=',NG2
            stop
          end if

          ! Konstrukcija stupca matricnih elementa MnmK1K2(G)
          call genStrujniVhrovi(jump, eps, Nlf, iG0, NG1, NG2, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2, MnmK1K22)
          
          deallocate(C2)

          intraORinterband_if: if (n /= m) then
          omega_loop: do io = -no,no
            o = io*domega
            dE = o + E(K1,n) - E(K2,m)
            Lor = Gamma_inter/(dE**2 + Gamma_inter**2) ! Gamma_inter je sirina interband prijelaza
            ! Reze repove Lorentziana lijevo i desno, pazljivo, minimum 1.0d-3, preporuceno 1.0d-5
            if (abs(Lor) >= Lor_cut/Gamma_inter) then
              do  iG = 1,Nlf
                do  jG = 1,Nlf
                  ! -1/pi*ImChi_munu -> for Kramers Kronig
                  S0(io,iG,jG) = S0(io,iG,jG) - 2.0*df*Lor*MnmK1K2(iG)*conjg(MnmK1K22(jG)) / (pi*Ntot*Vcell)
                end do
              end do
            end if
          end do omega_loop
          else if (n == m .and. abs(E(K1,n)-Efermi) <= 10.0*T) then
            ! Effective number of charge carriers (tensor)
            fact = expo1/( (expo1 + 1.0D0)*(expo1 + 1.0D0) )
            fact = -fact/T
            do  iG = 1,Nlf
              do  jG = 1,Nlf
                ! izracun intraband korelacijske funkcije
                Qeff(iG,jG) = Qeff(iG,jG) + 2.0*fact*MnmK1K2(iG)*conjg(MnmK1K22(jG)) / (Ntot*Vcell)
              end do
            end do
          end if intraORinterband_if
        end if occupation_if
      end do bands_m_loop
      deallocate(C1)
    end do bands_n_loop
    jump = 1
    
    834            CONTINUE

  end do ! end of FBZ do loop
  !$omp end do
  !$omp end parallel
  
  
  ! WRITING CORFUN S0_\mu\nu
  ! WRITTING Q_eff_\mu\nu
  print *,'got here, line 722'
  open(74,FILE = dato1)
  open(75,FILE = dato2)
  do  io = -no,no ! opskurni razlog za prosirenje raspona frekvencija na negativne da se korektno izracuna spektar kristala koji nemaju centar inverzije
    o = io*domega
    write(74,*)'omega=',o,'Hartree'
    write(74,'(10F15.10)')((S0(io,iG,jG),jG = 1,Nlf),iG = 1,Nlf)
  end do
  write(75,'(10F15.10)')((Qeff(iG,jG),jG = 1,Nlf),iG = 1,Nlf)
  close(74)
  close(75)
  
  
  if (calc == 1) GO TO 999
  
!               SECOND PART OF THE PROGRAM calc = 2
!               Calculation of the matrix '''Pi_\mu\nu'' by using matrix ''S0_\mu\nu(G,G')'' and Kramers-Krroning relations
  
  888             CONTINUE 

  
  open(74,FILE = dato1)
  do  io=-no,no
    READ(74,*) ! nis
    READ(74,'(10F15.10)')((S0(io,iG,jG),jG = 1,Nlf),iG = 1,Nlf)
  end do
  close(74)
  
  open(75,FILE = dato2)
  READ(75,'(10F15.10)')((Qeff(iG,jG),jG = 1,Nlf),iG = 1,Nlf)
  close(75)
  
  
  
  
  ! Convert (qx,qy,qz) and Glf to Cartesian coordinates
  
  qx = Gcar*qx
  qy = Gcar*qy
  qz = Gcar*qz
  
  do  iG = 1,Nlf
    Glf(1:3,iG)=Gcar*Glf(1:3,iG)
  end do
  
  
  
  
  open(77,FILE = dato3)
  ! new sum over omega
  omega_loop: do io = 1,no-1
    ! print *, io
    oi = (io-1)*domega

    iG_loop: do iG = 1,Nlf
      jG_loop: do jG = 1,Nlf

        call genReChi0(io,no,iG,jG,oi,domega,S0,Rechi0)
        ! ReChi0 = 0.0 ! real part of the response function Re(Chi)
        ! ! static limit
        ! if (io == 1) then
        !   do  jo = 2,no
        !     oj=(jo-1)*domega
        !     fact = domega/oj
        !     if (jo == 2)fact = 3.0/2.0
        !     if (jo == no)fact = 0.5*domega/oj
        !     ReChi0 = ReChi0 + fact*(real(S0(-jo + 1,iG,jG)) -real(S0(jo-1,iG,jG)))
        !   end do
        ! else if (io == 2) then
        !   do  jo = 1,no
        !     oj=(jo-1)*domega
        !     if (jo /= io) then
        !       fact = domega/(oi-oj)
        !     endif
        !     if (jo == 1) then
        !       fact = 1.0
        !     else if (jo == 2) then
        !       fact = 0.0
        !     else if (jo == 3) then
        !       fact=-3.0/2.0
        !     else if (jo == no) then
        !       fact = 0.5*domega/(oi-oj)
        !     end if
        !     ReChi0 = ReChi0 + fact*real(S0(jo-1,iG,jG))
        !     fact = domega/(oi + oj)
        !     if (jo == 1 .or. jo == no) then
        !       fact = 0.5*domega/(oi + oj)
        !     end if
        !     ReChi0 = ReChi0 + fact*real(S0(-jo + 1,iG,jG))
        !   end do
        ! else if (io == (no-1)) then
        !   do  jo = 1,no
        !     oj=(jo-1)*domega
        !     if (jo /= io) then
        !       fact = domega/(oi-oj)
        !     end if

        !     if (jo == 1) then
        !       fact = 0.5*domega/(oi-oj)
        !     else if (jo == (no-2)) then
        !       fact = 3.0/2.0
        !     else if (jo == (no-1)) then
        !       fact = 0.0
        !     else if (jo == no) then
        !       fact=-1.0
        !     end if

        !     ReChi0 = ReChi0 + fact*real(S0(jo-1,iG,jG))
        !     fact = domega/(oi + oj)

        !     if (jo == 1 .or. jo == no) then
        !       fact = 0.5*domega/(oi + oj)
        !     end if

        !     ReChi0 = ReChi0 + fact*real(S0(-jo + 1,iG,jG))
        !   end do
        ! else
        !   do  jo = 1,no
        !     oj = (jo-1)*domega
        !     if (jo /= io) then
        !       fact = domega/(oi-oj)
        !     end if

        !     if (jo == 1) then
        !       fact = 0.5*domega/(oi-oj)
        !     else if (jo == (io-1)) then
        !       fact = 3.0/2.0
        !     else if (jo == io) then
        !       fact = 0.0
        !     else if (jo == (io + 1)) then
        !       fact=-3.0/2.0
        !     else if (jo == no) then
        !       fact = 0.5*domega/(oi-oj)
        !     end if

        !     ReChi0 = ReChi0 + fact*real(S0(jo-1,iG,jG))
        !     fact = domega/(oi + oj)

        !     if (jo == 1 .or. jo == no) then
        !       fact = 0.5*domega/(oi + oj)
        !     end if

        !     ReChi0 = ReChi0 + fact*real(S0(-jo + 1,iG,jG))
        !   end do
        ! end if
        ! ReChi0 = ReChi0 + pi*aimag(S0(io-1,iG,jG))
        
        call genImChi0(io,no,iG,jG,oi,domega,S0,ImChi0)
        ! ImChi0 = 0.0 ! Imaginary part of the response function Im(Chi)
        ! ! static limit
        ! if (io == 1) then
        !   do  jo = 2,no
        !     oj = (jo-1)*domega
        !     fact = domega/oj
        !     if (jo == 2) then
        !       fact = 3.0/2.0
        !     else if (jo == no) then
        !       fact = 0.5*domega/oj
        !     end if
        !     ImChi0 = ImChi0 + fact*(aimag(S0(-jo + 1,iG,jG))-aimag(S0(jo-1,iG,jG)))
        !   end do
        ! else if (io == 2) then
        !   do  jo = 1,no
        !     oj = (jo-1)*domega
        !     if (jo /= io) then
        !       fact = domega/(oi-oj)
        !     end if
        !     if (jo == 1) then 
        !       fact = 1.0
        !     else if (jo == 2) then 
        !       fact = 0.0
        !     else if (jo == 3) then 
        !       fact=-3.0/2.0
        !     else if (jo == no) then
        !       fact = 0.5*domega/(oi-oj)
        !     end if

        !     ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jG))
        !     fact = domega/(oi + oj)

        !     if (jo == 1 .or. jo == no) then
        !       fact = 0.5*domega/(oi + oj)
        !     end if
        !     ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jG))
        !   end do
        ! else if (io == (no-1)) then
        !   do  jo = 1,no
        !     oj = (jo-1)*domega
        !     if (jo /= io) then
        !       fact = domega/(oi-oj) 
        !     end if
        !     if (jo == 1) then
        !       fact = 0.5*domega/(oi-oj)
        !     else if (jo == (no-2)) then
        !       fact = 3.0/2.0
        !     else if (jo == (no-1)) then
        !       fact = 0.0
        !     else if (jo == no) then
        !       fact=-1.0
        !     end if

        !     ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jG))
        !     fact = domega/(oi + oj)

        !     if (jo == 1 .or. jo == no) then
        !       fact = 0.5*domega/(oi + oj)
        !     end if

        !     ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jG))
        !   end do
        ! else
        !   do  jo = 1,no
        !     oj = (jo-1)*domega
        !     if (jo /= io) then
        !       fact = domega/(oi-oj)
        !     end if

        !     if (jo == 1) then
        !       fact = 0.5*domega/(oi-oj)
        !     else if (jo == (io-1)) then
        !       fact = 3.0/2.0
        !     else if (jo == io) then
        !       fact = 0.0
        !     else if (jo == (io + 1)) then
        !       fact=-3.0/2.0
        !     else if (jo == no) then
        !       fact = 0.5*domega/(oi-oj)
        !     end if

        !     ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jG))
        !     fact = domega/(oi + oj)

        !     if (jo == 1 .or. jo == no) then
        !       fact = 0.5*domega/(oi + oj)
        !     end if

        !     ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jG))
        !   end do
        ! end if
        
        ! ImChi0 = ImChi0 - pi*real(S0(io-1,iG,jG)) ! ovaj dio je razlicit od Sloss, S0 je kompleksno polje
        
        if (io == 1) then 
          Pi_dia(iG,jG) = -cmplx(ReChi0,0.0) ! neven debug: diamagnetski doprinos ??
        end if
        Pi_tot(iG,jG) = cmplx(ReChi0,ImChi0) ! Pi_RPA = PiDIJAMAGNETSKI + PiPARAMAGNETSKI
        Pi_tot(iG,jG) = Pi_tot(iG,jG) + Pi_dia(iG,jG)
        
        Pi_tot(iG,jG) = Pi_tot(iG,jG) + Qeff(iG,jG)*oi/(oi + cmplx(0.0,1.0)*Gamma_intra) ! dodavanje intraband clana
        
      end do jG_loop
    end do iG_loop
    
    ! WRITTING TOTAL RESPONSE FUNCTION Pi for a given polarization 'pol' to file
    write(77,*) oi*Hartree,Pi_tot(1,1)
    
    
    
!                end new sum over omega
  end do omega_loop
  close(77)
  
  999              CONTINUE
!                kraj po q
end do q_loop


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

  subroutine writeInfo(lf, pol, qx, qy, qz, Gcar, Nsymm, Nlf, Ntot, NkI, Nband, T, Nel, NelQE,Gamma_intra, Gamma_inter)
    implicit none
    integer      ,    intent(in) :: NelQE, Nsymm, Nlf, Ntot, NkI, Nband
    real(kind=dp),    intent(in) :: qx,qy,qz
    real(kind=dp),    intent(in) :: eta, T, Gcar
    real(kind=dp),    intent(in) :: Nel
    real(kind=dp),    intent(in) :: Gamma_inter, Gamma_intra
    character(len=*), intent(in) :: lf, pol

    integer       :: ios
    real(kind=dp) :: absq, error
    real(kind=dp), parameter :: Hartree = 2.0D0*13.6056923D0

    absq = sqrt(qx**2+qy**2+qz**2)

    open(55,FILE='Info', err=700, iostat=ios)
    write(55,*)'***************General***********************'
    write(55,*)' Currently we calculate         ---->',dato1
    write(55,*)' Currently we calculate         ---->',dato2
    write(55,*)''
    write(55,*)'Number of point symmetry operation is',Nsymm
    ! if (frac == 0)write(55,*)'Fraction translation is not detected'
    ! if (frac == 1)write(55,*)'Fraction translation is detected'
    write(55,'(A25,3F10.4,A5)') 'Wave vector (qx,qy,qz)=(',qx*Gcar,qy*Gcar, qz*Gcar,') a.u.'
    write(55,'(A25,F7.3,A5)') '|(qx,qy,qz)|=',absq*Gcar,'a.u.'
    write(55,*) 'Local field effcts in '//lf//'-dir'
    write(55,*) 'Polarization in '//pol//'-dir'
    write(55,*)'Number of local field vectors is',Nlf
    write(55,*)'Number of different K vectors in 1.B.Z. is',Ntot
    write(55,*)'Number of K vectors in I.B.Z. is',NkI
    write(55,*)'Number of bands is               ',Nband
    write(55,'(A25,F7.3,A5)') 'Gamma_intra is  ',Gamma_intra*Hartree*1000.0,'meV'
    write(55,'(A25,F7.3,A5)') 'Gamma_inter is  ',Gamma_inter*Hartree*1000.0,'meV'
    write(55,'(A25,F7.3,A5)') 'Temperature is      ',T*Hartree*1000.0,'meV'
    write(55,*)''
    write(55,*)'-Im(Chi(io,G1,G2))/pi is in file---->',dato1
    write(55,*)' Qeff complex matrix is in file ---->',dato2
    write(55,*)' Pi_munu is in file             ---->',dato3
    write(55,*)''
    write(55,*)'************* Checking 1BZ integration*******'
    write(55,*)''
    write(55,'(A40,F8.4)')'Number of electrons(1BZ integration)=',Nel
    write(55,*)'Number of electrons(unit cell)=',NelQE
    error = abs((NelQE-Nel)/NelQE)
    write(55,'(A25,F7.3,A5)') 'Relative error=',error*100.0,'%'
    if (error > 0.05) then
      write(55,*)'WARRNING!!-1BZ INTEGRATION IS BAD!.'
      stop
    end if
    close(55)    
    
    goto 600
    700 write(*,*) 'Cannot open Info file. iostat = ',ios
    stop
    600 continue

  end subroutine writeInfo

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
      stop
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
    stop
  else if (it == 2) then
    print*,'Can not find wave vector K+Q=',ik,'+',iq, 'in IBZ.'
    stop
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
            stop
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
    else if (lf == 'xyz') then
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
        if (Gi(1) /= 0  .or.  Gi(2) /= 0  .or.  Gi(3) /= 0) then
          print*,'*********************************'
          print*,'WARRNING!, G vectors input is wrong!!'
          print*,'G(1) is not (0,0,0)!!'
          stop
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


  subroutine genReChi0(io,no,iG,jG,oi,domega,S0,ReChi0)
    implicit none
    integer,          intent(in)  :: io, no
    integer,          intent(in)  :: iG, jG
    real(kind=dp),    intent(in)  :: oi, domega
    complex(kind=dp), intent(in)  :: S0(:,:,:)
    real(kind=dp),    intent(out) :: ReChi0

    integer :: jo
    real(kind=dp) :: oj, fact

    ReChi0 = 0.0 ! real part of the response function
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
        ReChi0 = ReChi0 + fact*( real(S0(-jo+1, iG, jG)) - real(S0(jo-1, iG, jG)) )
      end do
    else if (io == 2) then
      do  jo = 1,no
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj)
        endif
        if (jo == 1) then
          fact = 1.0
        else if (jo == 2) then
          fact = 0.0
        else if (jo == 3) then
          fact=-3.0/2.0
        else if (jo == no) then
          fact = 0.5*domega/(oi-oj)
        end if
        ReChi0 = ReChi0 + fact*real(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)
        if (jo == 1 .or. jo == no) then
          fact = 0.5*domega/(oi + oj)
        end if
        ReChi0 = ReChi0 + fact*real(S0(-jo + 1,iG,jG))
      end do
    else if (io == (no-1)) then
      do  jo = 1,no
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj)
        end if
        if (jo == 1) then
          fact = 0.5*domega/(oi-oj)
        else if (jo == (no-2)) then
          fact = 3.0/2.0
        else if (jo == (no-1)) then
          fact = 0.0
        else if (jo == no) then
          fact=-1.0
        end if
        ReChi0 = ReChi0 + fact*real(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)
        if (jo == 1 .or. jo == no) then
          fact = 0.5*domega/(oi + oj)
        end if
        ReChi0 = ReChi0 + fact*real(S0(-jo + 1,iG,jG))
      end do
    else
      do  jo = 1,no
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj)
        end if
        if (jo == 1) then
          fact = 0.5*domega/(oi-oj)
        else if (jo == (io-1)) then
          fact = 3.0/2.0
        else if (jo == io) then
          fact = 0.0
        else if (jo == (io + 1)) then
          fact=-3.0/2.0
        else if (jo == no) then
          fact = 0.5*domega/(oi-oj)
        end if
        ReChi0 = ReChi0 + fact*real(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)
        if (jo == 1 .or. jo == no) then
          fact = 0.5*domega/(oi + oj)
        end if
        ReChi0 = ReChi0 + fact*real(S0(-jo + 1,iG,jG))
      end do
    end if
    ReChi0 = ReChi0 + pi*aimag(S0(io-1,iG,jG))
    
  end subroutine genReChi0

  subroutine genImChi0(io,no,iG,jG,oi,domega,S0,ImChi0)
    implicit none
    integer,          intent(in)  :: io, no
    integer,          intent(in)  :: iG, jG
    real(kind=dp),    intent(in)  :: oi, domega
    complex(kind=dp), intent(in)  :: S0(:,:,:)
    real(kind=dp),    intent(out) :: ImChi0

    integer :: jo
    real(kind=dp) :: oj, fact

    ImChi0 = 0.0 ! Imaginary part of the response function Im(Chi)
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
        ImChi0 = ImChi0 + fact*(aimag(S0(-jo + 1,iG,jG)) - aimag(S0(jo-1,iG,jG)))
      end do
    else if (io == 2) then
      do  jo = 1,no
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj)
        end if
        if (jo == 1) then 
          fact = 1.0
        else if (jo == 2) then 
          fact = 0.0
        else if (jo == 3) then 
          fact=-3.0/2.0
        else if (jo == no) then
          fact = 0.5*domega/(oi-oj)
        end if

        ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)

        if (jo == 1 .or. jo == no) then
          fact = 0.5*domega/(oi + oj)
        end if
        ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jG))
      end do
    else if (io == (no-1)) then
      do  jo = 1,no
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj) 
        end if
        if (jo == 1) then
          fact = 0.5*domega/(oi-oj)
        else if (jo == (no-2)) then
          fact = 3.0/2.0
        else if (jo == (no-1)) then
          fact = 0.0
        else if (jo == no) then
          fact=-1.0
        end if

        ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)

        if (jo == 1 .or. jo == no) then
          fact = 0.5*domega/(oi + oj)
        end if

        ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jG))
      end do
    else
      do  jo = 1,no
        oj = (jo-1)*domega
        if (jo /= io) then
          fact = domega/(oi-oj)
        end if

        if (jo == 1) then
          fact = 0.5*domega/(oi-oj)
        else if (jo == (io-1)) then
          fact = 3.0/2.0
        else if (jo == io) then
          fact = 0.0
        else if (jo == (io + 1)) then
          fact=-3.0/2.0
        else if (jo == no) then
          fact = 0.5*domega/(oi-oj)
        end if

        ImChi0 = ImChi0 + fact*aimag(S0(jo-1,iG,jG))
        fact = domega/(oi + oj)

        if (jo == 1 .or. jo == no) then
          fact = 0.5*domega/(oi + oj)
        end if

        ImChi0 = ImChi0 + fact*aimag(S0(-jo + 1,iG,jG))
      end do
    end if
    
    ImChi0 = ImChi0 - pi*real(S0(io-1,iG,jG)) ! ovaj dio je razlicit od Sloss, S0 je kompleksno polje
    
  end subroutine genImChi0

  subroutine loadkIandE(path, NkI, Nband, Nocc, kI, dGW,E)
    implicit none
    integer,            intent(in)    :: NkI
    integer,            intent(in)    :: Nband, Nocc
    character(len=100), intent(in)    :: path
    real(kind=dp),      intent(in)    :: dGW(:,:)
    real(kind=dp),      intent(inout) :: kI(:,:)
    real(kind=dp),      intent(inout) :: E(:,:)

    integer :: ios, ik, i
    real(kind=dp),    parameter :: Hartree = 2.0D0*13.6056923D0

    open(400,FILE=path,status='old',err=500,iostat=ios) 
    do  ik = 1,NkI
      if (ik == 1) then
        read(400,*) 
      end if
      read(400,'(10X,3F10.3)') kI(1,ik),kI(2,ik),kI(3,ik)
      read(400,'(10F8.4)') (E(ik,i),i=1,Nband)
    end do
    close(400)
      
    goto 400
    500 write(*,*) 'Cannot open BAND file. iostat = ',ios
    stop
    400 continue
    
    ! konverzija en. u Hartree
    E(1:NkI,1:Nband) = E(1:NkI,1:Nband)/Hartree
    ! scissor operator, ispravlja/shifta DFT gap (na 1eV u ovom slucaju)
    E(1:NkI,Nocc+1:Nband) = E(1:NkI,Nocc+1:Nband) + dGW/Hartree

  end subroutine loadkIandE

subroutine genStrujniVhrovi(jump, eps, Nlf, iG0, NG1, NG2, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2, MnmK1K22)
  ! Konstrukcijamatricnih elementa strujnih vrhova MnmK1K2(iG) i MnmK1K2(iG) 
  implicit none
  integer,          intent(in)    :: iG0, Nlf, NG1,
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
  complex(kind=dp), intent(out)   :: MnmK1K22(:)

  integer          :: iG, iG1, iG2
  integer          :: iGfast
  real(kind=dp)    :: K11, K22, K33
  real(kind=dp)    :: Gxx1,Gyy1,Gzz1
  real(kind=dp)    :: Gxx2,Gyy2,Gzz2
  complex(kind=dp) :: struja, struja_y, struja_z

  iGfast = 0
  MnmK1K2(1:Nlf)  = cmplx(0.0,0.0)
  MnmK1K22(1:Nlf) = cmplx(0.0,0.0)
  iG_loop: do  iG = 1,Nlf
    iG1_loop: do iG1 = 1,NGd
      iGfast = iGfast + 1
      Gxx1 = G(1,iG1)
      Gyy1 = G(2,iG1)
      Gzz1 = G(3,iG1)
      ! neven debug -> provjeri tocnost
      ! k11 = R(R1,1,1)*Gxx1 + R(R1,1,2)*Gyy1 + R(R1,1,3)*Gzz1
      ! k22 = R(R1,2,1)*Gxx1 + R(R1,2,2)*Gyy1 + R(R1,2,3)*Gzz1
      ! k33 = R(R1,3,1)*Gxx1 + R(R1,3,2)*Gyy1 + R(R1,3,3)*Gzz1
      k11 = sum( R(R1,1,1:3)*G(1:3,iG1) )
      k22 = sum( R(R1,2,1:3)*G(1:3,iG1) )
      k33 = sum( R(R1,3,1:3)*G(1:3,iG1) )
      if (pol == 'xx') then
        struja = (qx + 2.0*kx + Glf(1,iG) + 2.0*k11)*Gcar
      else if (pol == 'yy') then
        struja = (qy + 2.0*ky + Glf(2,iG) + 2.0*k22)*Gcar
      else if (pol == 'zz') then
        struja = (qz + 2.0*kz + Glf(3,iG) + 2.0*k33)*Gcar
      else if (pol == 'yz') then
        struja_y = (qy + 2.0*ky + Glf(2,iG) + 2.0*k22)*Gcar
        struja_z = (qz + 2.0*kz + Glf(3,iG) + 2.0*k33)*Gcar
      end if
      k11 = k11 + Glf(1,iG) + G(1,iG0)
      k22 = k22 + Glf(2,iG) + G(2,iG0)
      k33 = k33 + Glf(3,iG) + G(3,iG0)
      
      Gxx1 = RI(R2,1,1)*k11 + RI(R2,1,2)*k22 + RI(R2,1,3)*k33
      Gyy1 = RI(R2,2,1)*k11 + RI(R2,2,2)*k22 + RI(R2,2,3)*k33
      Gzz1 = RI(R2,3,1)*k11 + RI(R2,3,2)*k22 + RI(R2,3,3)*k33
      if (jump == 1) then
        iG2_loop: do iG2 = 1,NG2
          Gfast(iGfast) = NG2 + 1
          Gxx2 = G(1,iG2)
          Gyy2 = G(2,iG2)
          Gzz2 = G(3,iG2)
          if (       abs(Gxx2-Gxx1) < eps &
              .and. abs(Gyy2-Gyy1) < eps &
              .and. abs(Gzz2-Gzz1) < eps ) then
            Gfast(iGfast) = iG2
            exit iG2_loop
          end if
        end do iG2_loop
      end if
      iG2 = Gfast(iGfast)
      if (iG2 <= NG2) then
        ! ako je polarazcija je tipa xx, yy, zz 
        if (pol == 'xx' .or. pol== 'yy' .or. pol == 'zz') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5D0*conjg(C1(iG1)) * struja * C2(iG2)
          MnmK1K22(iG) = MnmK1K2(iG) ! strujni vrhovi su isti
        else ! ako  su miksani xy, xz,...
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5D0*conjg(C1(iG1)) * struja_y * C2(iG2)
          MnmK1K22(iG) = MnmK1K22(iG) + 0.5D0*conjg(C1(iG1)) * struja_z * C2(iG2)
        end if
      end if
    end do iG1_loop
  end do iG_loop
  jump = 2 ! za svaki valni vektor q i dani k zapamti Gfast i za svaku vrpcu preskaci taj postupak   
    
end subroutine genStrujniVhrovi

end PROGRAM surface_current

