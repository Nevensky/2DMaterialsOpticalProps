module brillouin_zone
  use iso_fortran_env, only: dp => real64
  implicit none

  public :: loadKiandE, genFBZ, checkFBZintegration, &
            & findKinIBZ, findKQinIBZ, findMinQ
  private

contains

  subroutine loadkIandE(path, NkI, Nband, Nval, kI, E, dGW)
    !! Loads all wavevectors (in Cartesiand coords.) in the ireducible BZ and
    !! corresponding eigen-energies form the Quantum Espresso .band files
    use constants, only: Hartree
    integer,          intent(in)           :: NkI
    integer,          intent(in)           :: Nband, Nval
    character(len=*), intent(in)           :: path
    real(kind=dp),    intent(inout)        :: kI(:,:)
    real(kind=dp),    intent(inout)        :: E(:,:)
    real(kind=dp),    intent(in), optional :: dGW

    integer :: iuni, ios, ik, i
    ! real(kind=dp),    parameter :: Hartree = 2.0D0*13.6056923D0
    real(kind=dp) :: dGW_ = 0.0_dp
    if (present(dGW)) dGW_ = dGW

    open(newunit=iuni,FILE=path,status='old',iostat=ios) 
    if (ios /=0) then
      print *, 'ERROR: Cannot open BAND file: ',path
      stop
    end if
    do  ik = 1,NkI
      if (ik == 1) then
        read(iuni,*) 
      end if
      read(iuni,'(10X,3F10.6)') kI(1,ik),kI(2,ik),kI(3,ik)
#ifdef __linux__
      read(iuni,'(10F9.4)') (E(ik,i),i=1,Nband)
#endif
#ifdef __APPLE__
      read(iuni,'(10F9.3)') (E(ik,i),i=1,Nband)
#endif
    end do
    close(iuni)

    ! convert energy to Hartree
    E(1:NkI,1:Nband) = E(1:NkI,1:Nband)/Hartree
    ! scissor operator, shifts the DFT bandgap
    E(1:NkI,Nval+1:Nband) = E(1:NkI,Nval+1:Nband) + dGW_

  end subroutine loadkIandE

  subroutine genFBZ(Nk,NkI,Nsym,eps,kI,R,Ntot,ktot,writeOutput)
    !! Generates all unique wavectors in the 1st BZ by applying 
    !! point group transformations on the reducible BZ
    integer,           intent(in)  :: Nk, NkI, Nsym ! No. of k-points, iredducible k-kpoints, symm. ops.
    real(kind=dp),     intent(in)  :: eps       ! threshold to distinguish whether k-points are the same
    real(kind=dp),     intent(in)  :: kI(:,:)   ! k-points in the irreducible BZ
    real(kind=dp),     intent(in)  :: R(:,:,:)  ! point group transformation matrices
    
    integer,           intent(out) :: Ntot      ! total No. of unique k-point in the 1st BZ
    real(kind=dp),     intent(out) :: ktot(:,:) ! unique k-points in the 1st BZ
    logical, optional, intent(in)  :: writeOutput ! write FBZ to file yes/no?

    logical       :: unique
    integer       :: iuni, ios
    integer       :: i, ik, jk, lk
    integer       :: n, m
    real(kind=dp) :: k(3,Nk) ! all non-unique k-points within the 1st BZ

    jk = 0
    Ntot = 0 

    symm_loop: do  i = 1, Nsym    ! loop over all symmetries
      k_loop_IBZ: do  ik = 1, NkI  ! loop over k points in IBZ
        unique = .true.
        jk = jk + 1
        do  n = 1, 3   ! loop over kx, ky, kz 
          k(n,jk) = 0.0_dp
          do  m = 1, 3 ! loop over kx, ky, kz
            k(n,jk) = k(n,jk) + R(i,n,m)*kI(m,ik) ! uses rotational symmetry to reconstrucrt FBZ from IBZ
          end do
        end do 

        if (jk > 1) then
          do  lk = 1, jk-1
            ! Check if the given k-point is unique (i.e. skips if it was already added)
            if ( all ( abs(k(1:3,jk)-k(1:3,lk)) <= eps ) ) then 
              unique = .false.
            end if
          end do
        end if

        if (unique) then ! if it is unique add it to ktot
          Ntot = Ntot+1
          ktot(1:3,Ntot) = k(1:3,jk)
        end if

      end do k_loop_IBZ
    end do symm_loop

    ! write all (kx,ky) 1st BZ k-points to file
    if (present(writeOutput)) then
      if (writeOutput) then
        open(newunit=iuni,iostat=ios,file='fbz_check.dat',action='write',status='new')
        if (ios/=0) then
          stop 'ERROR: Could not open fbz_check.dat file in genFBZ()'
        end if
        do  i = 1,Ntot
          write(iuni,*) ktot(1,i), ktot(2,i)  
        end do
        close(iuni)
      end if
    end if
  end subroutine genFBZ

    subroutine checkFBZintegration(Nband,NkI,Nsym,Ntot,eps,kI,ktot,RI,Efermi,E,NelQE,Nel,T)
    !! Checks if the No. of electrons in the 1st BZ (Nel) equals 
    !! the number of electrons in the unit cell as calculated by Quantum Espresso (NelQE)   
    use statistics, only: FermiDirac
    integer,       intent(in)  :: NelQE
    integer,       intent(in)  :: NkI, Nsym, Nband, Ntot
    real(kind=dp), intent(in)  :: eps, Efermi
    real(kind=dp), intent(in)  :: kI(:,:)
    real(kind=dp), intent(in)  :: ktot(:,:)  ! all unique k-points in the FBZ

    real(kind=dp), intent(in)  :: RI(:,:,:)
    real(kind=dp), intent(in)  :: E(:,:)
    real(kind=dp), intent(out) :: Nel
    real(kind=dp), optional, intent(in) :: T

    logical       :: found
    integer       :: ik, n, i, j, l, K1
    real(kind=dp) :: kx,ky,kz
    real(kind=dp) :: Ni ! Fermi-Dirac occupation
    ! real(kind=dp) :: K11, K22, K33
    real(kind=dp) :: K(3) ! => k_IBZ = R_inv x k_FBZ = R_inv x ktot

    Nel = 0 
    k_loop_FBZ : do  ik = 1,Ntot
      band_loop: do  n = 1, Nband
        if (n == 1) then
            found = .false.
            if (ik <= NkI) then
              K1 = ik
              found = .true.
            else
              symm_loop: do  i = 2, Nsym
                forall (l=1:3) K(l) = sum ( RI(i,l,1:3)*ktot(1:3,ik) )
                k_loop_IBZ: do  j = 1, NkI
                  if ( all ( abs(K(1:3)-kI(1:3,j)) <= eps ) ) then
                    found = .true.
                    K1 = j
                    cycle band_loop
                  end if
                end do k_loop_IBZ
              end do symm_loop
            end if

            if (.not. found) then
              print*,'Can not find wave vector K=',ik, 'in I.B.Z.'
              stop
            end if

        end if

        ! debug neven: the counting of electrons should be moved to a new subroutine

        ! sums electrons in the remeaining bands
        if (present(T)) then ! temperature is given use Fermi-Dirac statistics
          Ni = FermiDirac(E(K1,n), Efermi, T)
          Nel = Nel + Ni
        else ! temperature not given, assuming its is an isolator
          if (E(K1,n) < Efermi) then 
            Ni = 1.0_dp
          else
            Ni = 0.0_dp
          end if  
        end if
        Nel = Nel + Ni
        ! print *,'Nel:', Nel,' Ni:',Ni,' band: ',n
    
      end do band_loop
    end do k_loop_FBZ
    ! debug neven: remove 2xNel for a spin-orbit calculation?
    Nel = 2.0*Nel / Ntot ! divide by Ntot due to overcounting caused by the k_loop_FBZ

  end subroutine checkFBZintegration


  subroutine findKinIBZ(ik, NkI, Nsym, eps, kx, ky, kz, RI, kI, iR1, iK1)
    !! Finds k-point (kx,ky,kz) in the ireducible Brillouin zone (IBZ)
    integer,       intent(in)  :: ik
    integer,       intent(in)  :: NkI, Nsym
    real(kind=dp), intent(in)  :: eps
    real(kind=dp), intent(in)  :: kx, ky, kz
    real(kind=dp), intent(in)  :: kI(:,:)
    real(kind=dp), intent(in)  :: RI(:,:,:)
    integer      , intent(out) :: iR1, iK1 ! K = R1*K1

    integer       :: i, j, l
    logical       :: found
    real(kind=dp) :: K(3), k_fbz(3)

    k_fbz = kx
    k_fbz = ky
    k_fbz = kz

    found = .false.
    if (ik <= NkI) then
      iR1 = 1
      iK1 = ik
      found = .true.
    else
      symmetry_loop: do  i = 2, Nsym
        forall (l=1:3) K(l) = sum ( RI(i,l,1:3)*k_fbz(1:3) )
        ! K(1) = sum (RI(i,1,1:3)*k_fbz(1:3) ) ! Kx
        ! K(2) = sum (RI(i,2,1:3)*k_fbz(1:3) ) ! Ky
        ! K(3) = sum (RI(i,3,1:3)*k_fbz(1:3) ) ! Kz
        do  j = 1,NkI
          if ( all ( abs( K(1:3)-kI(1:3,j) ) <= eps ) ) then
            found = .true.
            iR1 = i
            iK1 = j
            EXIT symmetry_loop
          end if
        end do
      end do symmetry_loop
    end if
    if (.not. found) then
      print *,'Can not find wave vector K=',ik, 'in I.B.Z.'
      stop
    end if
    ! print *, 'iR1:',iR1,'iK1:', iK1
  end subroutine findKinIBZ


  subroutine findKQinIBZ(KQx, KQy, KQz, eps, Nsym, NkI, Ntot, NG, kI, ktot, RI, G, iG0, iR2, iK2)
    !! Finds the k-point (KQx,KQy,KQz) in the 1st. Brillouin zone (FBZ) and then the ireducible Brillouin zone (IBZ)
    integer,       intent(in)  :: Nsym, NkI, Ntot, NG
    real(kind=dp), intent(in)  :: eps
    real(kind=dp), intent(in)  :: KQx, KQy, KQz
    real(kind=dp), intent(in)  :: kI(:,:)
    real(kind=dp), intent(in)  :: ktot(:,:)
    real(kind=dp), intent(in)  :: G(:,:)
    real(kind=dp), intent(in)  :: RI(:,:,:)
    integer,       intent(out) :: iG0, iR2, iK2 ! G0=rec.lattice; K + Q = G0 + R2*K2
  
    integer       :: found
    integer       :: iG, jk, i, j, l
    real(kind=dp) :: K(3), KQ(3)

    KQ(1) = KQx
    KQ(2) = KQy
    KQ(3) = KQy
  
    found = 1 ! false
    iG_loop: do  iG = 1,NG
      k_loop_FBZ : do  jk = 1, Ntot
        if ( all (abs(KQ(1:3)-G(1:3,iG)-ktot(1:3,jk)) <= eps) ) then
          found = 2 ! true for FBZ
          iG0 = iG
          symm_loop: do  i = 1, Nsym
            forall (l=1:3) K(l) = sum( RI(i,l,1:3) * ktot(1:3,jk) )
            ! K(1) = sum(RI(i,1,1:3) * ktot(1:3,jk) )
            ! K(2) = sum(RI(i,2,1:3) * ktot(1:3,jk) )
            ! K(3) = sum(RI(i,3,1:3) * ktot(1:3,jk) )
            k_loop_IBZ: do  j = 1, NkI
              if ( all( abs(K(1:3)-kI(1:3,j)) <= eps) ) then
                found = 3 ! true for IBZ (and also for FBZ)
                iR2 = i
                iK2 = j
                EXIT iG_loop
              end if
            end do k_loop_IBZ
          end do symm_loop
        end if
      end do k_loop_FBZ
    end do iG_loop  
  
    if (found == 1) then
      ! print*,'Can not find wave vector K+Q=',ik,'+',iq, 'in FBZ.'
      print*,'Can not find wave vector K+Q in FBZ.'
      stop
    else if (found == 2) then
      ! print*,'Can not find wave vector K+Q=',ik,'+',iq, 'in IBZ.'
      print*,'Can not find wave vector K+Q=',iK2, 'in IBZ.'
      stop
    end if
  
  end subroutine findKQinIBZ

  subroutine findMinQ(iq, Ntot, ktot, qx, qy, qz)
    !! Searches for the mininimal wavevector \(\mathbf{q}=(qx,qy,qz)\) in the \(\Gamma \to M \) direction
    integer,       intent(in)  :: iq, Ntot
    real(kind=dp), intent(in)  :: ktot(:,:)
    real(kind=dp), intent(out) :: qx, qy, qz

    integer       :: i, ikmin
    real(kind=dp) :: kmin, kref, krefM ! , absq
    ! real(kind=dp) :: q(3) ! (qx,qy,qz)

    kmin = 1.0
    Ntot_loop: do  i = 1, Ntot ! loop over different k-points in FBZ
      kref = sqrt(sum(ktot(1:3,i)**2))
      ! neven debug
      ! print *,'i=',i,' kref: ',kref
      if (kref == 0.0) then
        CYCLE Ntot_loop
      else if (kref < kmin) then
        kmin = kref
        ikmin = i
        krefM = kmin
        ! neven debug
        print *,'i=',i,'kmin: ',kmin
      end if
    end do Ntot_loop
    ! neve debug
    ! print *,'ikmin=',ikmin,'kmin=',kmin,'ktot(1:3,ikmin)',ktot
    qx = (iq-1) * ktot(1,ikmin)
    qy = (iq-1) * ktot(2,ikmin)
    qz = (iq-1) * ktot(3,ikmin)
    ! absq = sqrt(qx**2 + qy**2 + qz**2)

    ! q(1:3) = (iq-1) * ktot(1:3,ikmin)
    ! absq = sqrt( sum( q(1:3)**2 ) )
  end subroutine findMinQ

end module brillouin_zone