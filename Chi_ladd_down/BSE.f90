PROGRAM surface_loss


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


use OMP_lib
use ModPointR

implicit none

INTEGER :: nki,nband,ik,i,nk,j,jk,it,lk,ntot,Nsymm,ik1,ik2,  &
    ng,io,no,iq,nq,nmpx,nmpy,nmpz,n,m,ig,r1,k1,r2,k2,ig01,ig02,  &
    nlf,ng1,ng2,ngd,ig1,ig2,jg,nlfd,kg,jo,jump,ni,mi,nj,mj,MOD,  &
    igfast,ikmin,nelqe,lf,kg1,kg2,nord,io1,io2,nkbz,ref,nkd,  &
    pol,spin,Nocc,Nval_ls,nmax,mmax
PARAMETER(nmpx = 51,nmpy = 51,nmpz = 1,nki = 460,nband = 40,nelqe = 26,  &
    nk = 48*nki,ngd = 4000,ng = 8000,no = 401,nq = 121,nlfd = 50,nkd = 2700)


!        skalars
real*8 a0,c0,eps,kx,ky,kz,omin,omax,sgn,struja,  &
    qx,qy,qz,kmin,domega,o,ef,k11,k22,k33,t,lor,de,gabs,kref,  &
    eref,ecut,gxx1,gyy1,gzz1,gxx2,gyy2,gzz2,fact,vcell,oi,oj,  &
    q,gcar,abohr,zero,nel,absq,error,o1,o2,gama,vbare,qmax,  &
    estar,kpointx,kpointy,kpointz,qmin,vq,dGW,alpha,epsq, valley
doubleprecision pi,three,hartree,planck,ev,eta
PARAMETER(hartree = 2.0D0*13.6056923D0,ef = 1.5703/hartree,  &
    a0 = 5.9715,c0 = 29.8575,pi = 3.141592654D0,gcar = 2.0*pi/a0,  &
    eps = 1.0D-4,t = 0.01/hartree,eta = 0.05/hartree,ev = 1.602176487D-19,  &
    planck = 6.626196D-34,ecut = 0.0,vcell = 922.0586, abohr = 0.5291772D0,zero = 0.0)
DOUBLE COMPLEX a,ione,cmplx(0.0,0.0),rone,em,g0,w,l0i,l0j,j1,j2


!        arrays
INTEGER :: gfast,gi,parg
real*8 ki,e,r,ri,k,ktot,g,glf,kc,v,kbz,delta
COMPLEX*8 mnmk1k2,wt,fockk,xi0,lladd,current,chi_ladd, chi_para
DOUBLE COMPLEX fock,Imat
DIMENSION ki(3,nki),e(nki,nband),r(48,3,3),ri(48,3,3),  &
    k(3,nk),ktot(3,nk),g(3,ng),glf(3,nlfd),mnmk1k2(2,nlfd),  &
    gfast(nlfd*ngd),kc(3,3),gi(3),wt(nq,nlfd,nlfd),parg(ng),  &
    kbz(3,nk),fockk(nkd,nkd),delta(nkd),fock(nkd,nkd),  &
    Imat(nkd,nkd),xi0(nkd,nkd),lladd(nkd,nkd),current(2,nkd),  &
    chi_ladd(2),chi_para(2)


CHARACTER (LEN = 100) :: bandn,bandm,nis,pathk1,pathk2,dato,  &
    rundir,path,fajl,wdir,root2
CHARACTER (LEN = 35) :: tag,buffer

INTEGER :: :: ngw, igwx, nbnd, nk1
COMPLEX(:: kind = 8),pointer, DIMENSION(:) :: c1,c2


ione = DCMPLX(0.0,1.0)


!             QUANTUM ESSPRESSO IMPUTS:
rundir='/home/vito/PROJECTS/MoS2-BSE/MoS2_51x51/'
!             Staticka zasjenjena kulonska interakcija WT je smjestena u folderu
wdir='/home/vito/PROJECTS/MoS2-BSE/W/'



!             If mod = 1 Fock kernel calculation
!             If mod = 2 Ladder ireducible current-current polarisability calculation
!             SAMO ZA MOD = 2
!             For spin up put spin = 1
!             For spin down put spin = 2

MOD = 2
spin = 2
lf = 1
jump = 1
omin = 1.0D-5
omax = 4.0/hartree
domega=(omax-omin)/(no-1)


!             maximum transfer wave vector around K point..
!             qmax = sqrt(3.0)/3.0
qmax = 0.3
!             minimum transfer wave vector in WT calculation
qmin = 0.006
if ( qmax <= (1.0/3.0) ) then
  valley = 2.0
else
  valley = 1.0
endif


Nocc = 9 ! Number of valence bands (when there is no LS coupling)
Nval_ls = 26 ! Number of valence bands (with LS coupling, full-rel. PP)
dGW = 1.0  ! band gap correction (in eV)



! POINT GROUP TRANSFORMATIONS (in Cartesian coords.)
path = trim(rundir)//"/"//trim(scf_file)
call PointR(path,Nsymm,R,RI)
print *,"status: PointR done."



! OPREZ umjesto Nocc stavljeno Nval_ls
path=trim(rundir)//trim(band_file)
call loadkIandE(path, NkI, Nband, Nval_ls, kI, dGW,E)
print *,"status: kI and E loaded."


!            generator 1.B.Z.
!            Dio programa koji pomocu operacija tockaste grupe i vektora iz
!            I.B.Z. generira sve (MEDJUSOBNO RAZLICITE!!!) v. vektore u 1.B.Z.
!            Ntot-Tot number of different points ''ktot'' inside 1.B.Z



jk = 0
ntot = 0
do  i = 1,Nsymm
  do  ik = 1,nki
    it = 1
    jk = jk+1
    do  n = 1,3
      k(n,jk)=zero
      do  m = 1,3
        k(n,jk)=k(n,jk)+r(i,n,m)*ki(m,ik)
      end do
    end do
    if (jk > 1) then
      do  lk = 1,jk-1
        if (abs(k(1,jk)-k(1,lk)) <= eps) then
          if (abs(k(2,jk)-k(2,lk)) <= eps) then
            if (abs(k(3,jk)-k(3,lk)) <= eps) then
              it = 2
            end if
          end if
        end if
      end do
    end if
    if (it == 1) then
      ntot = ntot+1
      ktot(1,ntot)=k(1,jk)
      ktot(2,ntot)=k(2,jk)
      ktot(3,ntot)=k(3,jk)
    end if
  end do
end do


!             Checking 1BZ integration
nel = 0
do  ik = 1,ntot
  kx = ktot(1,ik)
  ky = ktot(2,ik)
  kz = ktot(3,ik)
  do  n = 1,nband
    if (n == 1) then
      it = 1
      if (ik <= nki) then
        k1 = ik
        it = 2
      else
        do  i = 2,Nsymm
          k11 = ri(i,1,1)*kx+ri(i,1,2)*ky+ri(i,1,3)*kz
          k22 = ri(i,2,1)*kx+ri(i,2,2)*ky+ri(i,2,3)*kz
          k33 = ri(i,3,1)*kx+ri(i,3,2)*ky+ri(i,3,3)*kz
          do  j = 1,nki
            if (DABS(k11-ki(1,j)) <= eps) then
              if (DABS(k22-ki(2,j)) <= eps) then
                if (DABS(k33-ki(3,j)) <= eps) then
                  it = 2
                  k1 = j
                  GO TO 5022
                end if
              end if
            end if
          end do
        end do
      end if
      if (it == 1) then
        print *,'Can not find wave vector K=',ik, 'in I.B.Z.'
        STOP
      end if
      5022          CONTINUE
    end if
    if (e(k1,n) < ef)nel = nel+1.0
  end do
end do
nel = nel/ntot

do  i = 1,ntot
  write(887,*)ktot(1,i),ktot(2,i)
end do


epsq = sqrt(4.0*pi*c0/(vcell*ntot))


!           KC transformation matrix from rec.cryst. axes to cart.koord.
!           If g' is vector in rec.cryst. axes then a = KC*a' is vector in cart. axes


fajl='/MoS2.sc.out'
path = trim(rundir)//trim(fajl)
tag='     reciprocal axes: (cart. coord.'
open(1,FILE = path)
do  i = 1,100000
  read(1,'(a)')buffer
  if (buffer == tag) then
    do  j = 1,3
      read(1,70)kc(1,j),kc(2,j),kc(3,j)
    end do
    GO TO 998
  end if
end do
70          FORMAT(23X,3F10.3)
998         CONTINUE
close(1)


!           Reading the reciprocal vectors in crystal coordinates and transformation
!           in Cartezi cordinates.
open (1,FILE='gvectors.dat')
do  i = 1,8
  read(1,*)nis
end do
do  ig = 1,ng
  read(1,100)gi(1),gi(2),gi(3)
  if (ig == 1) then
    if (gi(1) /= 0 .or. gi(2) /= 0 .or. gi(3) /= 0) then
      print *,'*********************************'
      print *,'WARRNING!, G vectors input is wrong!!'
      print *,'G(1) is not (0,0,0)!!'
      STOP
    end if
  end if
!           transformation in cart.coord (also!, after this all G components are in 2pi/a0 units)
  do   n = 1,3
    g(n,ig)=zero
    do   m = 1,3
      g(n,ig)=g(n,ig)+kc(n,m)*DBLE(gi(m))
    end do
  end do
  parg(ig)=gi(3)
end do
100         FORMAT(i10,i11,i11)
close(1)



!            Reciprocal vectors for crystal local field effects calculations in array ''Glf(3,Nlf)''


nlf = 0
if (lf == 1) then
  do  ig = 1,ng
    if (g(1,ig) == 0.0 .and. g(2,ig) == 0.0) then
      eref = gcar*gcar*g(3,ig)*g(3,ig)/2.0
      if (eref <= ecut) then
        nlf = nlf+1
        glf(1,nlf)=0.0
        glf(2,nlf)=0.0
        glf(3,nlf)=g(3,ig)
        IF((parg(ig)/2)*2 == parg(ig)) then
          parg(nlf)=1
        else
          parg(nlf)=-1
        end if
      end if
    end if
  end do
else
  do  ig = 1,ng
    eref = gcar*gcar*(g(1,ig)*g(1,ig)+g(2,ig)*g(2,ig)+ g(3,ig)*g(3,ig))/2.0
    if (eref <= ecut) then
      nlf = nlf+1
      glf(1,nlf)=g(1,ig)
      glf(2,nlf)=g(2,ig)
      glf(3,nlf)=g(3,ig)
      IF((parg(ig)/2)*2 == parg(ig)) then
        parg(nlf)=1
      else
        parg(nlf)=-1
      end if
    end if
  end do
end if
if (nlf > nlfd) then
  print *,'Nlf is bigger than Nlfd'
  STOP
end if



!             GENERIRANJE male BZ OKO 'K' TOCKE
!             K-point in Cartesi cordintes in  alatt 0.333333  0.577350  0.000000
kpointx = 0.333333
kpointy = 0.577350
kpointz = 0.000000
nkbz = 0


!            generiranje 'male' BZ

jk = 0
nkbz = 0
do  i = 1,Nsymm
  do  ik = 1,nki
    if (abs(ki(2,ik)) <= qmax) then
      it = 1
      jk = jk+1
      do  n = 1,3
        k(n,jk)=zero
        do  m = 1,3
          k(n,jk)=k(n,jk)+r(i,n,m)*ki(m,ik)
        end do
      end do
      if (jk > 1) then
        do  lk = 1,jk-1
          if (abs(k(1,jk)-k(1,lk)) <= eps) then
            if (abs(k(2,jk)-k(2,lk)) <= eps) then
              if (abs(k(3,jk)-k(3,lk)) <= eps) then
                it = 2
              end if
            end if
          end if
        end do
      end if
      if (it == 1) then
        nkbz = nkbz+1
        kbz(1,nkbz)=k(1,jk)
        kbz(2,nkbz)=k(2,jk)
        kbz(3,nkbz)=k(3,jk)
      end if
    end if
  end do
end do

!            Pomak 'male BZ' u K tocku


do  i = 1,nkbz
  kbz(1,i)=kbz(1,i)+kpointx
  kbz(2,i)=kbz(2,i)+kpointy
end do

!            PROVJERA DA LI IMA ISTIH TOCAKA u OKLICI K tocke
do  i = 1,nkbz
  kx = kbz(1,i)
  ky = kbz(2,i)
  kz = kbz(3,i)
  do  j = 1,nkbz
    if (j /= i) then
      if (abs(kx-kbz(1,j)) <= eps) then
        if (abs(ky-kbz(2,j)) <= eps) then
          if (abs(kz-kbz(3,j)) <= eps) then
            print *,'Postoje jednake k tocke u okolini K tocke'
            STOP
          end if
        end if
      end if
    end if
  end do
end do




do  i = 1,nkbz
  write(889,*)kbz(1,i),kbz(2,i)
end do






!             Info file

open(55,FILE='Info')
write(55,*)'***************General***********************'
write(55,*)''
write(55,*)'Number of point symmetry operation is',Nsymm
if (lf == 1)write(55,*)'Local field effcts in z-dir'
if (lf == 3)write(55,*)'Local field in all xyz-dir'
write(55,*)'Number of local field vectors is',nlf
write(55,*)'Number of different K vectors in 1.B.Z. is',ntot
write(55,*)'Number of K vectors in I.B.Z. is',nki
write(55,*)'Number of K points around K-point is',nkbz
write(55,*)'Number of bands is               ',nband
write(55,*)''
write(55,*)'************* Checking 1BZ integration*******'
write(55,*)''
write(55,12)'Number of electrons(1BZ integration)=',nel
write(55,*)'Number of electrons(unit cell)=',nelqe
error = abs((nelqe-nel)/nelqe)
write(55,99)'Relative error=',error*100.0,'%'
if (error > 0.05) then
  write(55,*)'WARRNING!!-1BZ INTEGRATION IS BAD!.'
end if
close(55)
88         FORMAT(a25,3F10.4,a5)
99         FORMAT(a25,f8.4,a5)
12         FORMAT(a40,f7.4)


if (MOD == 2)GO TO 8889


! Loading of time ordered screened Coulomb interaction $W_GG'^T(Q,\omega)$
do  iq = 2,nq
  q = (iq-1)*qmin
  vq = (2.0*pi)/q
  dato ='W_Qi'
  nord = INDEX(dato,'i', back =.false.)
  if (iq < 10) then
    write(dato(nord:nord),'(i1)') iq
  else if (iq >= 10 .and. iq < 100) then
    write(dato(nord:nord+1),'(i2)')iq
  else
    write(dato(nord:nord+2),'(i3)')iq
  end if
  path = trim(wdir)//trim(dato)
  
  open(74,FILE = path)
  read(74,*)nis
  read(74,'10F15.5')((wt(iq,ig,jg),jg = 1,nlf),ig = 1,nlf)
  close(74)
  
!             trazenje alpha
  if (iq == 2)alpha=(c0*vq/real(wt(2,1,1))-1.0)/qmin
  
  
end do





!             valni vektor K
do  ik1 = 1,nkbz
  
  print *,ik1
  
  kx = kbz(1,ik1)
  ky = kbz(2,ik1)
  kz = kbz(3,ik1)
  
  
  
!              trazenje K prvo u 1.B.Z a onda u I.B.Z.
  it = 1
  do  ig = 1,ng
    do  jk = 1,ntot
      if (DABS(kx-g(1,ig)-ktot(1,jk)) <= eps) then
        if (DABS(ky-g(2,ig)-ktot(2,jk)) <= eps) then
          if (DABS(kz-g(3,ig)-ktot(3,jk)) <= eps) then
            it = 2
            ig01 = ig
            do  i = 1,Nsymm
              k11 = ri(i,1,1)*ktot(1,jk)+ri(i,1,2)*ktot(2,jk)+  &
                  ri(i,1,3)*ktot(3,jk)
              k22 = ri(i,2,1)*ktot(1,jk)+ri(i,2,2)*ktot(2,jk)+  &
                  ri(i,2,3)*ktot(3,jk)
              k33 = ri(i,3,1)*ktot(1,jk)+ri(i,3,2)*ktot(2,jk)+  &
                  ri(i,3,3)*ktot(3,jk)
              do  j = 1,nki
                if (DABS(k11-ki(1,j)) <= eps) then
                  if (DABS(k22-ki(2,j)) <= eps) then
                    if (DABS(k33-ki(3,j)) <= eps) then
                      it = 3
                      r1 = i
                      k1 = j
                      exit
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
  2001                 CONTINUE
  
  
  if (it == 1) then
    print *,'Can not find wave vector K=',ik1, 'in 1.B.Z.'
    STOP
  else if (it == 2) then
    print *,'Can not find wave vector K=',ik1, 'in I.B.Z.'
    STOP
  end if
  
  
  
  
!             valni vektor K'
  do  ik2 = 1,nkbz
    
    
    kx = kbz(1,ik2)
    ky = kbz(2,ik2)
    kz = kbz(3,ik2)
    
    
    it = 1
!              trazenje K' prvo u 1.B.Z a onda u I.B.Z.
    
    do  ig = 1,ng
      do  jk = 1,ntot
        if (DABS(kx-g(1,ig)-ktot(1,jk)) <= eps) then
          if (DABS(ky-g(2,ig)-ktot(2,jk)) <= eps) then
            if (DABS(kz-g(3,ig)-ktot(3,jk)) <= eps) then
              it = 2
              ig02 = ig
              do  i = 1,Nsymm
                k11 = ri(i,1,1)*ktot(1,jk)+ri(i,1,2)*ktot(2,jk)+  &
                    ri(i,1,3)*ktot(3,jk)
                k22 = ri(i,2,1)*ktot(1,jk)+ri(i,2,2)*ktot(2,jk)+  &
                    ri(i,2,3)*ktot(3,jk)
                k33 = ri(i,3,1)*ktot(1,jk)+ri(i,3,2)*ktot(2,jk)+  &
                    ri(i,3,3)*ktot(3,jk)
                do  j = 1,nki
                  if (DABS(k11-ki(1,j)) <= eps) then
                    if (DABS(k22-ki(2,j)) <= eps) then
                      if (DABS(k33-ki(3,j)) <= eps) then
                        it = 3
                        r2 = i
                        k2 = j
                        exit
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
    
    
    if (it == 1) then
      print *,'Can not find wave vector K=',ik2, 'in 1.B.Z.'
      STOP
    else if (it == 2) then
      print *,'Can not find wave vector K=',ik2, 'in I.B.Z.'
      STOP
    end if
    
    
    
    qx = kbz(1,ik2)-kbz(1,ik1)
    qy = kbz(2,ik2)-kbz(2,ik1)
    q = gcar*sqrt(qx*qx+qy*qy)
    iq = q/qmin+1
    
    
    
    if (iq > nq)GO TO 9999
    
    
!              ''K = G01+R1*K1''.
!              ''K'=G02+R2*K2''.
!              K1 i K2 su integeri koji predstavljaju valne vektore u IBZ pridruzene
!              valnim vektorima K i K'
!              iG01 i iG02 su integeri rec vec. G koji translatira K i K' u 1BZ
    
!             vrpca n
    do   n = Nocc,Nocc+1
      
      if (n == Nocc) ni = 1
      if (n == Nocc+1) ni = 2
!             vrpca m
      do  m = n,n
        
        
        
        CALL paths(rundir,k1,k2,n,m,pathk1,pathk2,bandn,bandm)
        
        
        
        
!         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za
!         fiksni K1,K2,n i m
        
!               Otvaranje atribute za INFO
        CALL iotk_open_read(10,pathk1)
        CALL iotk_scan_empty(10,"INFO",attr = attr)
        CALL iotk_scan_attr(attr,"igwx",ng1)
!               Alociranje polja C1
        allocate (c1(ng1))
!               Ucitavanje podataka iza evc.n
        CALL iotk_scan_dat(10,bandn,c1)
        CALL iotk_close_read(10)
!               Otvaranje atribute za INFO
        CALL iotk_open_read(10,pathk2)
        CALL iotk_scan_empty(10,"INFO",attr = attr)
        CALL iotk_scan_attr(attr,"igwx",ng2)
!               Alociranje polja C2
        allocate (c2(ng2))
!               Ucitavanje podataka iza evc.m
        CALL iotk_scan_dat(10,bandm,c2)
        CALL iotk_close_read(10)
        
!                Konstrukcija stupca matricnih elementa MnmK1K2(ni,mi,G)
        
        
        if (ngd > ng1) then
          write(*,*)'NGd is bigger than NG1=',ng1
          STOP
        else if (ngd > ng2) then
          write(*,*)'NGd is bigger than NG2=',ng2
          STOP
        end if
        
        
!                 matrix elements
        igfast = 0
        do  ig = 1,nlf
          mnmk1k2(ni,ig)=cmplx(0.0,0.0)
          do  ig1 = 1,ngd
            igfast = igfast+1
            gxx1 = g(1,ig1)
            gyy1 = g(2,ig1)
            gzz1 = g(3,ig1)
            k11 = r(r1,1,1)*gxx1+r(r1,1,2)*gyy1+r(r1,1,3)*gzz1
            k22 = r(r1,2,1)*gxx1+r(r1,2,2)*gyy1+r(r1,2,3)*gzz1
            k33 = r(r1,3,1)*gxx1+r(r1,3,2)*gyy1+r(r1,3,3)*gzz1
            k11 = k11+glf(1,ig)
            k22 = k22+glf(2,ig)
            k33 = k33+glf(3,ig)
            k11 = k11+g(1,ig02)-g(1,ig01)
            k22 = k22+g(2,ig02)-g(2,ig01)
            k33 = k33+g(3,ig02)-g(3,ig01)
            gxx1 = ri(r2,1,1)*k11+ri(r2,1,2)*k22+ri(r2,1,3)*k33
            gyy1 = ri(r2,2,1)*k11+ri(r2,2,2)*k22+ri(r2,2,3)*k33
            gzz1 = ri(r2,3,1)*k11+ri(r2,3,2)*k22+ri(r2,3,3)*k33
            if (jump == 1) then
              do  ig2 = 1,ng2
                gfast(igfast)=ng2+1
                gxx2 = g(1,ig2)
                gyy2 = g(2,ig2)
                gzz2 = g(3,ig2)
                if (DABS(gxx2-gxx1) < eps) then
                  if (DABS(gyy2-gyy1) < eps) then
                    if (DABS(gzz2-gzz1) < eps) then
                      gfast(igfast)=ig2
                      GO TO 1111
                    end if
                  end if
                end if
              end do
            end if
            1111              CONTINUE
            ig2 = gfast(igfast)
            if (ig2 <= ng2) then
              mnmk1k2(ni,ig)=mnmk1k2(ni,ig)+ CONJG(c1(ig1))*c2(ig2)
            end if
          end do
        end do
        jump = 2
        
        
        deallocate(c1)
        deallocate(c2)
        
!               end of m loop
      end do
      
      
!               end of n loop
    end do
    jump = 1
    
    
    fockk(ik1,ik2)=cmplx(0.0,0.0)
    
    if (ik1 /= ik2) then
      do  ig = 1,nlf
        do  jg = 1,nlf
          w = wt(iq,ig,jg)
          fockk(ik1,ik2)=fockk(ik1,ik2)- w*CONJG(mnmk1k2(1,ig))*mnmk1k2(2,jg)
        end do
      end do
      
!               komponenta Q = 0
    else if (ik1 == ik2) then
      fockk(ik1,ik2)=-vcell*ntot*LOG(1.0+alpha*epsq)/alpha
    end if
    7768            CONTINUE
    7769            CONTINUE
    
    
    9999            CONTINUE
!               end of K'  loop
  end do
!               end of K loop
end do



!              upisivanje Fock kernela i DeltaE
open(234,FILE='FockK')
write(234,'10F15.5')((fockk(ik1,ik2),ik2 = 1,nkbz),ik1 = 1,nkbz)
close(234)
GO TO 9981



8889           CONTINUE

!*********************************************
!              Ovdje pocinje mod = 2
!*********************************************


print *,'Reading Fock kernel'

!              citanje Fock kernela i DeltaE
open(234,FILE='FockK')
read(234,'10F15.5')((fockk(ik1,ik2),ik2 = 1,nkbz),ik1 = 1,nkbz)
close(234)


!               GENERIRANJE STRUJNIH VRHOVA jmu

do  pol = 1,2
  
  if (pol == 1)print *,'Calculation of j_x'
  if (pol == 2)print *,'Calculation of j_z'
  
  
  
!               valni vektor K
  do  ik1 = 1,nkbz
    
    kx = kbz(1,ik1)
    ky = kbz(2,ik1)
    kz = kbz(3,ik1)
    
!               trazenje K prvo u 1.B.Z a onda u I.B.Z.
    it = 1
    do  ig = 1,ng
      do  jk = 1,ntot
        if (DABS(kx-g(1,ig)-ktot(1,jk)) <= eps) then
          if (DABS(ky-g(2,ig)-ktot(2,jk)) <= eps) then
            if (DABS(kz-g(3,ig)-ktot(3,jk)) <= eps) then
              it = 2
              ig01 = ig
              do  i = 1,Nsymm
                k11 = ri(i,1,1)*ktot(1,jk)+ri(i,1,2)*ktot(2,jk)+  &
                    ri(i,1,3)*ktot(3,jk)
                k22 = ri(i,2,1)*ktot(1,jk)+ri(i,2,2)*ktot(2,jk)+  &
                    ri(i,2,3)*ktot(3,jk)
                k33 = ri(i,3,1)*ktot(1,jk)+ri(i,3,2)*ktot(2,jk)+  &
                    ri(i,3,3)*ktot(3,jk)
                do  j = 1,nki
                  if (DABS(k11-ki(1,j)) <= eps) then
                    if (DABS(k22-ki(2,j)) <= eps) then
                      if (DABS(k33-ki(3,j)) <= eps) then
                        it = 3
                        r1 = i
                        k1 = j
                        exit
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
    2101                 CONTINUE
    
    
    if (it == 1) then
      print *,'Can not find wave vector K=',ik1, 'in 1.B.Z.'
      STOP
    else if (it == 2) then
      print *,'Can not find wave vector K=',ik1, 'in I.B.Z.'
      STOP
    end if
    
!             valni vektor K'
    do  ik2 = ik1,ik1
      
      kx = kbz(1,ik2)
      ky = kbz(2,ik2)
      kz = kbz(3,ik2)
      
      
      it = 1
!              trazenje K' prvo u 1.B.Z a onda u I.B.Z.
      
      do  ig = 1,ng
        do  jk = 1,ntot
          if (DABS(kx-g(1,ig)-ktot(1,jk)) <= eps) then
            if (DABS(ky-g(2,ig)-ktot(2,jk)) <= eps) then
              if (DABS(kz-g(3,ig)-ktot(3,jk)) <= eps) then
                it = 2
                ig02 = ig
                do  i = 1,Nsymm
                  k11 = ri(i,1,1)*ktot(1,jk)+ri(i,1,2)*ktot(2,jk)+  &
                      ri(i,1,3)*ktot(3,jk)
                  k22 = ri(i,2,1)*ktot(1,jk)+ri(i,2,2)*ktot(2,jk)+  &
                      ri(i,2,3)*ktot(3,jk)
                  k33 = ri(i,3,1)*ktot(1,jk)+ri(i,3,2)*ktot(2,jk)+  &
                      ri(i,3,3)*ktot(3,jk)
                  do  j = 1,nki
                    if (DABS(k11-ki(1,j)) <= eps) then
                      if (DABS(k22-ki(2,j)) <= eps) then
                        if (DABS(k33-ki(3,j)) <= eps) then
                          it = 3
                          r2 = i
                          k2 = j
                          exit
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
      2411           CONTINUE
      
      
      if (it == 1) then
        print *,'Can not find wave vector K=',ik2, 'in 1.B.Z.'
        STOP
      else if (it == 2) then
        print *,'Can not find wave vector K=',ik2, 'in I.B.Z.'
        STOP
      end if
      
      
      
!              ''K = G01+R1*K1''.
!              ''K'=G02+R2*K2''.
!              K1 i K2 su integeri koji predstavljaju valne vektore u IBZ pridruzene
!              valnim vektorima K i K'
!              iG01 i iG02 su integeri rec vec. G koji translatira K i K' u 1BZ
      
      
      
      
      do  n = Nocc,Nocc
        do  m = Nocc+1,Nocc+1
          
          
          
          CALL paths(rundir,k1,k2,n,m,pathk1,pathk2,bandn,bandm)
          
          
          
!         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za
!         fiksni K1,K2,n i m
          
!               Otvaranje atribute za INFO
          CALL iotk_open_read(10,pathk1)
          CALL iotk_scan_empty(10,"INFO",attr = attr)
          CALL iotk_scan_attr(attr,"igwx",ng1)
!               Alociranje polja C1
          allocate (c1(ng1))
!               Ucitavanje podataka iza evc.n
          CALL iotk_scan_dat(10,bandn,c1)
          CALL iotk_close_read(10)
!               Otvaranje atribute za INFO
          CALL iotk_open_read(10,pathk2)
          CALL iotk_scan_empty(10,"INFO",attr = attr)
          CALL iotk_scan_attr(attr,"igwx",ng2)
!               Alociranje polja C2
          allocate (c2(ng2))
!               Ucitavanje podataka iza evc.m
          CALL iotk_scan_dat(10,bandm,c2)
          CALL iotk_close_read(10)
          
!                Konstrukcija stupca matricnih elementa MnmK1K2(G)
          
          
          if (ngd > ng1) then
            write(*,*)'NGd is bigger than NG1=',ng1
            STOP
          else if (ngd > ng2) then
            write(*,*)'NGd is bigger than NG2=',ng2
            STOP
          end if
          
          
!              matrix elements
          igfast = 0
          do  ig = 1,1
            current(pol,ik1)=cmplx(0.0,0.0)
            do  ig1 = 1,ngd
              igfast = igfast+1
              gxx1 = g(1,ig1)
              gyy1 = g(2,ig1)
              gzz1 = g(3,ig1)
              k11 = r(r1,1,1)*gxx1+r(r1,1,2)*gyy1+r(r1,1,3)*gzz1
              k22 = r(r1,2,1)*gxx1+r(r1,2,2)*gyy1+r(r1,2,3)*gzz1
              k33 = r(r1,3,1)*gxx1+r(r1,3,2)*gyy1+r(r1,3,3)*gzz1
              if (pol == 1) then
                struja=(2.0*kx+glf(1,ig)+2.0*k11-2.0*g(1,ig01))*gcar
              else if (pol == 2) then
                struja=(2.0*kz+glf(3,ig)+2.0*k33-2.0*g(3,ig01))*gcar
              end if
              k11 = k11+glf(1,ig)
              k22 = k22+glf(2,ig)
              k33 = k33+glf(3,ig)
              k11 = k11+g(1,ig02)-g(1,ig01)
              k22 = k22+g(2,ig02)-g(2,ig01)
              k33 = k33+g(3,ig02)-g(3,ig01)
              gxx1 = ri(r2,1,1)*k11+ri(r2,1,2)*k22+ri(r2,1,3)*k33
              gyy1 = ri(r2,2,1)*k11+ri(r2,2,2)*k22+ri(r2,2,3)*k33
              gzz1 = ri(r2,3,1)*k11+ri(r2,3,2)*k22+ri(r2,3,3)*k33
              if (jump == 1) then
                do  ig2 = 1,ng2
                  gfast(igfast)=ng2+1
                  gxx2 = g(1,ig2)
                  gyy2 = g(2,ig2)
                  gzz2 = g(3,ig2)
                  if (DABS(gxx2-gxx1) < eps) then
                    if (DABS(gyy2-gyy1) < eps) then
                      if (DABS(gzz2-gzz1) < eps) then
                        gfast(igfast)=ig2
                        GO TO 1011
                      end if
                    end if
                  end if
                end do
              end if
              1011           CONTINUE
              ig2 = gfast(igfast)
              if (ig2 <= ng2) then
                current(pol,ik1)=current(pol,ik1)+  &
                    0.5D0*CONJG(c1(ig1))*struja*c2(ig2)
              end if
!              kraj po iG1
            end do
!              kraj po C.L.F. iG
          end do
          jump = 2
          
          
          
          if (spin == 1)delta(ik1)=e(k1,Nval_ls-1)-e(k1,Nval_ls+1)
          if (spin == 2)delta(ik1)=e(k1,Nval_ls)-e(k1,Nval_ls+2)
          
          
!              end of m do loop
        end do
        
        
!              end of n do loop
      end do
      jump = 1
      
!              ond of B.Z do loop
    end do
  end do
  
!              kraj po polarizaciji
end do


!               Rijesavanje BSE-FOCK jednadzbe


if (spin == 1) then
  open(334,FILE='Pi_ladder_up_x')
  open(335,FILE='Pi_ladder_up_z')
else if (spin == 2) then
  open(334,FILE='Pi_ladder_down_x')
  open(335,FILE='Pi_ladder_down_z')
end if


do    io = 1,no
  o = omin+(io-1)*domega
  print *,io
  
  do  ik1 = 1,nkbz
    do  ik2 = 1,nkbz
      fock(ik1,ik2)=cmplx(0.0,0.0)
      Imat(ik1,ik2)=cmplx(0.0,0.0)
    end do
    Imat(ik1,ik1)=cmplx(1.0,0.0)
  end do
  
  
!               Konstrukcija matrice L^0*Xi^F
!               Konstrukcija matrice Xi0 = L^0*Xi^FL^0
  
  do  ik1 = 1,nkbz
    l0i = 1.0/(o+delta(ik1)+ione*eta)
    do  ik2 = 1,nkbz
      l0j = 1.0/(o+delta(ik2)+ione*eta)
      fock(ik1,ik2)=Imat(ik1,ik2)- l0i*fockk(ik1,ik2)/(ntot*vcell)
      xi0(ik1,ik2)=l0i*fockk(ik1,ik2)*l0j
    end do
  end do
  
  
  
  CALL gjel(fock,nkbz,nkd,Imat,nkbz,nkd)
  
  
!               konstrukcija 4-point polarizabilnosti L^ladd
  
  
  
  do  ik1 = 1,nkbz
    do  ik2 = 1,nkbz
      lladd(ik1,ik2)=cmplx(0.0,0.0)
      do  ik = 1,nkbz
        lladd(ik1,ik2)=lladd(ik1,ik2)+fock(ik1,ik)*xi0(ik,ik2)
      end do
    end do
  end do
  
  
  
!               Konstrukcija Ladder ireduciblne current-current polarisabilnosti chi_ladd_\mu\mu
  
  
  
  do  pol = 1,2
    chi_ladd(pol)=cmplx(0.0,0.0)
    chi_para(pol)=cmplx(0.0,0.0)
    do  ik1 = 1,nkbz
      l0i=(o/delta(ik1))*(1.0/(o+delta(ik1)+ione*eta))
      do  ik2 = 1,nkbz
        j1 = current(pol,ik1)
        j2 = CONJG(current(pol,ik2))
        chi_ladd(pol)=chi_ladd(pol)-  &
            valley*j1*lladd(ik1,ik2)*j2/(ntot*vcell*ntot*vcell)
      end do
!               Paramagnetska struja-struja polarizabilost
      chi_para(pol)= chi_para(pol)+ valley*j1*l0i*CONJG(j1)/(ntot*vcell)
    end do
  end do
  
  write(334,*)o*hartree,chi_ladd(1)
  write(335,*)o*hartree,chi_ladd(2)
  
 write(104,*)o*hartree,aimag(chi_para(1)+chi_ladd(1))
 write(105,*)o*hartree,aimag(chi_para(1))
 write(106,*)o*hartree,real(chi_para(1)+chi_ladd(1))
 write(107,*)o*hartree,real(chi_para(1))
  
  write(204,*)o*hartree,aimag(chi_para(2)+chi_ladd(2))
  write(205,*)o*hartree,aimag(chi_para(2))
  write(206,*)o*hartree,real(chi_para(2)+chi_ladd(2))
  write(207,*)o*hartree,real(chi_para(2))
  
  
  
  
!               omega do loop
end do



close(334)
close(335)
close(336)



9981            CONTINUE


contains
  subroutine loadkIandE(path, NkI, Nband, Nocc, kI, dGW,E)
    implicit none
    integer,            intent(in)    :: NkI
    integer,            intent(in)    :: Nband, Nocc
    character(len=100), intent(in)    :: path
    real(kind=dp),      intent(in)    :: dGW
    real(kind=dp),      intent(inout) :: kI(:,:)
    real(kind=dp),      intent(inout) :: E(:,:)

    integer :: ios, ik, i
    real(kind=dp),    parameter :: Hartree = 2.0D0*13.6056923D0

    open(400,FILE=path,status='old',err=500,iostat=ios) 
    do  ik = 1,NkI
      if (ik == 1) then
        read(400,*) 
      end if
      read(400,'(10X,3F10.6)') kI(1,ik),kI(2,ik),kI(3,ik)
#ifdef __linux__
      read(400,'(10F9.4)') (E(ik,i),i=1,Nband)
#endif
#ifdef __APPLE__
      read(400,'(10F9.3)') (E(ik,i),i=1,Nband)
#endif
    end do
    close(400)
      
    goto 400
    500 write(*,*) 'Cannot open BAND file. iostat = ',ios
    stop
    400 continue
    
    ! konverzija en. u Hartree
    E(1:NkI,1:Nband) = E(1:NkI,1:Nband)/Hartree
    ! scissor operator, ispravlja/shifta DFT gap (na 1eV u ovom slucaju)
    E(1:NkI,Nocc+1:Nband) = E(1:NkI,Nocc+1:Nband) + dGW

  end subroutine loadkIandE



end program surface_loss

