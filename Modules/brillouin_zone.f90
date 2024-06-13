module brillouin_zone
  use iso_fortran_env, only: dp => real64
  use notifications, only: error
  use utility, only: int2str
  implicit none

  public :: loadKiandE, scissorE, genFBZ, checkFBZintegration, checkFBZintegration_new, &
            & findKinIBZ, findKQinIBZ, findMinQ, invertR, transformR, reshapeR, genFBZpath
  private

contains

  subroutine loadkIandE(path, NkI, Nband, Nval, kI, E, dGW)
    !! DEPRECATED: io_xml module stores E and kI differently
    !! Loads all wavevectors (in Cartesiand coords.) in the ireducible BZ and
    !! corresponding eigen-energies form the Quantum Espresso .band files
    use constants, only: Hartree
    integer,          intent(in)           :: NkI
    integer,          intent(in)           :: Nband, Nval
    character(len=*), intent(in)           :: path
    real(kind=dp),    intent(inout)        :: kI(:,:)
    real(kind=dp),    intent(inout)        :: E(:,:) !! DFT eigenenergies (NkI x Nband)
    real(kind=dp),    intent(in), optional :: dGW

    integer :: iuni, ios, ik, i
    ! real(kind=dp),    parameter :: Hartree = 2.0D0*13.6056923D0
    real(kind=dp) :: dGW_ = 0.0_dp
    if (present(dGW)) dGW_ = dGW

    open(newunit=iuni,FILE=path,status='old',iostat=ios) 
    if (ios /=0) then
      call error('Cannot open BAND file: '//path)
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

  subroutine scissorE(Nval, E, dGW)
    !! Shifts the band gap by \( \Delta E_\text{GW} \) for each k point
    use constants, only: Hartree
    integer,          intent(in)           :: Nval   !! No. valence bands
    real(kind=dp),    intent(inout)        :: E(:,:) !! eigenenergies (Nband x NkI)
    real(kind=dp),    intent(in), optional :: dGW    !! scissor band gap correction 

    integer :: Nband !! No. of bands
    integer :: NkI   !! No. of k-points in FBZ
    real(kind=dp) :: dGW_ = 0.0_dp

    if (present(dGW)) dGW_ = dGW
    Nband = size(E,1)
    NkI   = size(E,2)

    ! scissor operator, shifts the DFT bandgap
    E(Nval+1:Nband,1:NkI) = E(Nval+1:Nband,1:NkI) + dGW_
  end subroutine scissorE

  subroutine invertR(R, RI)
    !! Computes inverse of each rotational symmetry matrix
    use matrix_inverse, only: invert
    real(kind=dp), intent(inout) :: R(:,:,:)  !! point group rotational matrices (3 x 3 x Nrot)
    real(kind=dp), intent(out), allocatable, optional :: RI(:,:,:) !! rotational matrices inverted

    integer :: i, Nrot
    real(kind=dp), allocatable :: RI_(:,:,:)
    allocate(RI_(size(R,1),size(R,2),size(R,3)))
    RI_ = R

    Nrot = size(R,3)
    do i = 1,Nrot
      ! print *, 'DEBUG: irot=',i
      call invert(RI_(:,:,i))
    enddo

    ! do i=1,Nrot
    !   RI_(:,:,i) = transpose(RI_(:,:,i))
    ! end do

    if (present(RI)) then
      if (.not. allocated(RI)) allocate(RI(size(R,1),size(R,2),size(R,3)))
      RI = RI_ 
    else
      R = RI_
    end if

    deallocate(RI_)
  end subroutine invertR

  subroutine transformR(alat, a, b, R)
    real(kind=dp), intent(in) :: alat         !! lattice constant [Ry]
    real(kind=dp), intent(in) :: a(:,:)       !! direct lattice [Ry]
    real(kind=dp), intent(in) :: b(:,:)       !! reciprocal lattice [2pi/alat]
    real(kind=dp), intent(inout) :: R(:,:,:)  !! point group transformation (rotation) matrix (3 x 3 x Nrot)
    
    integer :: ir, Nrot
    Nrot = size(R,3)

    do ir = 1,Nrot
      ! neven debug: lattice converted to units [alat]
      R(:,:,ir) = matmul(R(:,:,ir),a/alat)
      R(:,:,ir) = transpose(R(:,:,ir))
      R(:,:,ir) = matmul(R(:,:,ir),b) 
      R(:,:,ir) = transpose(R(:,:,ir)) ! last step to make printout equal to QE scf output
      ! print *,'==R(:,:,',iq,')=='
      ! print *,R(:,:,iq)
    end do
  end subroutine transformR

  subroutine reshapeR(Nsym, R)
    !! Reshapes R(3x3xNrot) to R(3x3xNsym) because all rotations 
    !! above Nsym are not present in crystal
    integer,                    intent(in)    :: Nsym     !! No. of symmetries in crystal
    real(kind=dp), allocatable, intent(inout) :: R(:,:,:) !! Rotational point group transformation matrices
    
    integer :: Nrot
    real(kind=dp), allocatable :: Rtmp(:,:,:)
    Nrot = size(R,3)

    allocate(Rtmp(3,3,Nsym), source = R(:,:,1:Nsym))
    call move_alloc(Rtmp,R)

    ! print *,'Nrot =', Nrot,'Nsym=',Nsym
    ! print *, 'shape(R)=', shape(R)
  end subroutine reshapeR

  subroutine genFBZ(kI, R, ktot, Ntot, eps, writeOutput)
    !! Generates all unique wavectors in the 1st BZ by applying 
    !! point group transformations on the reducible BZ
    ! integer,                    intent(in)  :: Nsym
    real(kind=dp),              intent(in)  :: kI(:,:)      !! k-points in the IBZ (3 x NkI)
    real(kind=dp),              intent(in)  :: R(:,:,:)     !! point group transformation matrices (3 x 3 x Nrot)
    real(kind=dp), allocatable, intent(out) :: ktot(:,:)    !! unique k-points in the FBZ
    integer,       optional,    intent(out) :: Ntot         !! No. of unique k-point in the FBZ
    real(kind=dp), optional,    intent(in)  :: eps          !! threshold to distinguish whether two k-points are the same
    logical,       optional,    intent(in)  :: writeOutput  !! write FBZ to fbz_chek.dat file (yes/no)

    logical       :: unique    !! .true. if k-point is unique
    integer       :: iuni, ios !! file i/o vars
    integer       :: i, ik, jk, lk
    integer       :: n, m

    integer :: NkI   !! No. k-points in the IBZ
    integer :: Nsym  !! No. symmetry opperations
    integer :: Nk    !! No. of k-points in the FBZ
    integer :: Ntot_ !! No. of unique k-points in the FBZ
    real(kind=dp), allocatable :: k(:,:) !! all non-unique k-points within the 1st BZ (3xNk)
    real(kind=dp), allocatable :: tmp(:,:) !! holds ktot temporarily (3xNtot)
    real(kind=dp) :: eps_

    NkI = size(kI,2) 
    Nsym = size(R,3) ! careful, this is not true if R is not reshaped from 3x3xNrot to 3x3xNsym
    Nk = 48*NkI
    Ntot_ = 0

    eps_ = 1.0d-4 ! default thershold
    if (present(eps)) eps_ = eps

    allocate(k(3,Nk))
    if (.not. allocated(ktot)) allocate(ktot(3,Nk))

    jk = 0
    symm_loop: do  i = 1, Nsym     ! loop over all symmetries
      k_loop_IBZ: do  ik = 1, NkI  ! loop over k points in IBZ
        unique = .true.
        jk = jk + 1
        do  n = 1, 3   ! loop over kx, ky, kz 
          k(n,jk) = 0.0_dp
          do  m = 1, 3 ! loop over kx, ky, kz
            k(n,jk) = k(n,jk) + R(m,n,i)*kI(m,ik) ! uses rotational symmetry to reconstrucrt FBZ from IBZ
          end do
        end do 

        if (jk > 1) then
          do  lk = 1, jk-1
            ! Check if the given k-point is unique (i.e. skips if it was already added)
            if ( all ( abs(k(1:3,jk)-k(1:3,lk)) <= eps_ ) ) then 
              unique = .false.
            end if
          end do
        end if
        if (unique) then ! if it is unique add it to ktot
          Ntot_ = Ntot_+1
          ktot(1:3,Ntot_) = k(1:3,jk)
        end if

      end do k_loop_IBZ
    end do symm_loop

    if (present(Ntot)) Ntot = Ntot_

    ! print *, 'Ntot=',Ntot
    allocate(tmp(3,Ntot_),source=ktot(:,1:Ntot_)) ! allocates tmp with the 3xNtot shape
    call move_alloc(tmp,ktot) ! moves allocation & contents from tmp to ktot

    ! write all (kx,ky) surface FBZ k-points to file
    if (present(writeOutput)) then
      if (writeOutput) then
        open(newunit=iuni,iostat=ios,file='fbz_check.dat',action='write',status='new')
        if (ios/=0) then
          call error('Could not create new fbz_check.dat file in genFBZ()')
        end if
        do  i = 1,Ntot_
          write(iuni,*) ktot(1,i), ktot(2,i)
        end do
        close(iuni)
      end if
    end if

    if (allocated(k)) deallocate(k)
  end subroutine genFBZ

  ! subroutine genFBZ(Nk,NkI,Nsym,eps,kI,R,Ntot,ktot,writeOutput)
  !   !! Generates all unique wavectors in the 1st BZ by applying 
  !   !! point group transformations on the reducible BZ
  !   integer,           intent(in)  :: Nk, NkI, Nsym !! No. of k-points, iredducible k-kpoints, symm. ops.
  !   real(kind=dp),     intent(in)  :: eps           !! threshold to distinguish whether k-points are the same
  !   real(kind=dp),     intent(in)  :: kI(:,:)       !! k-points in the irreducible BZ
  !   real(kind=dp),     intent(in)  :: R(:,:,:)      !! point group transformation matrices
  !   integer,           intent(out) :: Ntot          !! total No. of unique k-point in the 1st BZ
  !   real(kind=dp),     intent(out) :: ktot(:,:)     !! unique k-points in the 1st BZ
  !   logical, optional, intent(in)  :: writeOutput   !! write FBZ to file yes/no?

  !   logical       :: unique
  !   integer       :: iuni, ios
  !   integer       :: i, ik, jk, lk
  !   integer       :: n, m
  !   real(kind=dp) :: k(3,Nk) ! all non-unique k-points within the 1st BZ

  !   jk = 0
  !   Ntot = 0 

  !   symm_loop: do  i = 1, Nsym    ! loop over all symmetries
  !     k_loop_IBZ: do  ik = 1, NkI  ! loop over k points in IBZ
  !       unique = .true.
  !       jk = jk + 1
  !       do  n = 1, 3   ! loop over kx, ky, kz 
  !         k(n,jk) = 0.0_dp
  !         do  m = 1, 3 ! loop over kx, ky, kz
  !           k(n,jk) = k(n,jk) + R(i,n,m)*kI(m,ik) ! uses rotational symmetry to reconstrucrt FBZ from IBZ
  !         end do
  !       end do 

  !       if (jk > 1) then
  !         do  lk = 1, jk-1
  !           ! Check if the given k-point is unique (i.e. skips if it was already added)
  !           if ( all ( abs(k(1:3,jk)-k(1:3,lk)) <= eps ) ) then 
  !             unique = .false.
  !           end if
  !         end do
  !       end if

  !       if (unique) then ! if it is unique add it to ktot
  !         Ntot = Ntot+1
  !         ktot(1:3,Ntot) = k(1:3,jk)
  !       end if

  !     end do k_loop_IBZ
  !   end do symm_loop

  !   ! write all (kx,ky) 1st BZ k-points to file
  !   if (present(writeOutput)) then
  !     if (writeOutput) then
  !       open(newunit=iuni,iostat=ios,file='fbz_check.dat',action='write',status='new')
  !       if (ios/=0) then
  !         call error('Could not open fbz_check.dat file in genFBZ()')
  !       end if
  !       do  i = 1,Ntot
  !         write(iuni,*) ktot(1,i), ktot(2,i)  
  !       end do
  !       close(iuni)
  !     end if
  !   end if
  ! end subroutine genFBZ

  subroutine checkFBZintegration(Nband,NkI,Nsym,Ntot,eps,kI,ktot,RI,Efermi,E,Nel,T)
    !! Checks if the No. of electrons in the 1st BZ (Nel) equals 
    !! the number of electrons in the unit cell as calculated by Quantum Espresso (NelQE)   
    use statistics, only: FermiDirac
    integer,       intent(in)  :: NkI, Nsym, Nband, Ntot
    real(kind=dp), intent(in)  :: Efermi     !! Fermi energy [Hartree]
    real(kind=dp), intent(in)  :: eps        !! threshold for two k-points to be considered equal
    real(kind=dp), intent(in)  :: kI(:,:)    !! k-points in the IBZ (3 x NkI)
    real(kind=dp), intent(in)  :: ktot(:,:)  !! all unique k-points in the FBZ

    real(kind=dp), intent(in)  :: RI(:,:,:)  !! inverse of point group rotantional matrices (3 x 3 x Nrot)
    real(kind=dp), intent(in)  :: E(:,:)     !! DFT eigenenergies [Hartree]
    real(kind=dp), intent(out) :: Nel        !! No. electrons computed in our code
    real(kind=dp), optional, intent(in) :: T !! electron temperature [Hartree]

    logical       :: found
    integer       :: ik, n, i, j, l, K1
    ! real(kind=dp) :: kx, ky, kz
    real(kind=dp) :: Ni !! Fermi-Dirac occupation
    ! real(kind=dp) :: K11, K22, K33
    real(kind=dp) :: K(3) !! => k_IBZ = R_inv x k_FBZ = R_inv x ktot

    Nel = 0.0_dp
    k_loop_FBZ : do  ik = 1,Ntot
      band_loop: do  n = 1, Nband
        if (n == 1) then
            found = .false.
            if (ik <= NkI) then
              K1 = ik
              found = .true.
            else
              symm_loop: do  i = 2, Nsym
                do l=1,3 
                  K(l) = sum ( RI(1:3,l,i)*ktot(1:3,ik) )
                end do
                k_loop_IBZ: do  j = 1, NkI
                  if ( all ( abs(K(1:3)-kI(1:3,j)) <= eps ) ) then
                    found = .true.
                    K1 = j
                    ! cycle band_loop ! WRONG, misses counting Ni for first band
                    goto 5022  ! why not exit sym_loop ?
                  end if
                end do k_loop_IBZ
              end do symm_loop
            end if
            5022 continue
            if (.not. found) then
              call error('Can not find wave vector K='//int2str(ik)//'in IBZ')
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

  subroutine checkFBZintegration_new(NelQE, Nel, kI, ktot, RI, E, Efermi, T, eps, spinorbit)
    !! Checks if the No. of electrons in the FBZ (Nel) equals 
    !! the number of electrons in the unit cell as calculated by Quantum Espresso (NelQE)   
    use statistics, only: FermiDirac
    real(kind=dp), intent(in)  :: NelQE        !! No. of electrons from QE calculation
    real(kind=dp), intent(in)  :: Efermi       !! Fermi energy
    real(kind=dp), intent(in)  :: kI(:,:)      !! k-points in the IBZ
    real(kind=dp), intent(in)  :: ktot(:,:)    !! all unique k-points in the FBZ
    real(kind=dp), intent(in)  :: RI(:,:,:)    !! inverted rotational symmetry matrices (3 x 3 x Nrot)
    real(kind=dp), intent(in)  :: E(:,:)       !! eigenenergies at each k-point (Nband x NkI)
    real(kind=dp), intent(out) :: Nel          !! No. of electrons
    real(kind=dp), optional, intent(in) :: T   !! electron temperature
    real(kind=dp), optional, intent(in) :: eps !! thershold for two k-points being the same
    logical,       optional, intent(in) :: spinorbit !! does the Quantum Espresso calculation include LS coupling?

    logical       :: found
    integer       :: ik, jk, n, i, l, K1
    integer       :: NkI, Nband, Ntot, Nsym
    ! real(kind=dp) :: kx,ky,kz
    real(kind=dp) :: Ni ! Fermi-Dirac occupation
    ! real(kind=dp) :: K11, K22, K33
    real(kind=dp) :: K(3) ! => k_IBZ = R_inv x k_FBZ = R_inv x ktot
    real(kind=dp) :: eps_
    real(kind=dp) :: eps_occupation = 0.05

    NkI = size(kI,2)
    Nsym = size(RI,3) ! careful, this only works if R has been reshaped from 3x3xNrot to 3x3xNsym
    Nband = size(E,1)
    Ntot = size(ktot,2)

    eps_ = 1.0d-4 ! default thershold
    if (present(eps)) eps_ = eps

    Nel = 0.0_dp
    k_loop_FBZ : do  ik = 1,Ntot
      found = .false.
      if (ik <= NkI) then
        K1 = ik
        found = .true.
      else
        symm_loop: do  i = 2, Nsym
          do l=1,3
            K(l) = sum ( RI(1:3,l,i)*ktot(1:3,ik) )
          end do
          k_loop_IBZ: do  jk = 1, NkI
            if ( all ( abs(K(1:3)-kI(1:3,jk)) <= eps_ ) ) then
              found = .true.
              K1 = jk ! IBZ index of a FBZ k-point
              exit symm_loop ! debug neven: is this ok?
            end if
          end do k_loop_IBZ
        end do symm_loop
      end if
      if (.not. found) then
        call error('Can not find wave vector iK='//int2str(ik)//'in IBZ')
      end if
      band_loop: do  n = 1, Nband
        ! calculate occupation of each band
        if (present(T)) then ! temperature is given use Fermi-Dirac statistics
          Ni = FermiDirac(E(n,K1), Efermi, T)
        else ! temperature not given, assuming the cyrstal is insulating
          if (E(n,K1) < Efermi) then 
            Ni = 1.0_dp
          else
            Ni = 0.0_dp
          end if  
        end if
        Nel = Nel + Ni ! occupation to sum of electrons
        ! print *,'DEBUG: Nel:', Nel,' Ni:',Ni,' band: ',n
      end do band_loop
    end do k_loop_FBZ
    
    Nel = 2*Nel / Ntot ! divide by Ntot due to overcounting caused by the k_loop_FBZ
    ! debug neven: remove 2xNel for a spin-orbit calculation?
    if (present(spinorbit)) then 
      if (spinorbit) Nel = Nel/2
    end if

    print *,'DEBUG: Nel: ',dble(Nel),'NelQE: ',NelQE
    write(*,'(A17,F6.2,A12)'), 'status: NelQE/Nel', 100.0*abs(NelQE-Nel)/NelQE,' % missmatch'
    if (abs(NelQE-dble(Nel)) > eps_occupation) then
      call error('Incorrect No. of electrons in FBZ.')
    end if

  end subroutine checkFBZintegration_new

  subroutine findKinIBZ(ik, kx, ky, kz, kI, RI, iR1, iK1, eps)
    !! Finds k-point (kx,ky,kz) in the ireducible Brillouin zone (IBZ)
    ! integer,       intent(in)  :: Nsym          !! No. of symmetry operations for the given crystal
    integer,       intent(in)  :: ik            !! index of k-point in IBZ
    real(kind=dp), intent(in)  :: kx, ky, kz    !! k-point in IBZ
    real(kind=dp), intent(in)  :: kI(:,:)       !! k-points in the IBZ (3 x NkI)
    real(kind=dp), intent(in)  :: RI(:,:,:)     !! inverted rotational symmetry matrices (3 x 3 x Nsym)
    integer      , intent(out) :: iR1           !! index of R1 in K = R1*K1
    integer      , intent(out) :: iK1           !! index of K1 in K = R1*K1
    real(kind=dp), intent(in), optional  :: eps !! thershold for two k-points being the same

    logical       :: found
    integer       :: i, j, l
    integer       :: NkI, Nsym
    real(kind=dp) :: K(3), k_fbz(3)
    real(kind=dp) :: eps_

    NkI = size(kI,2)
    Nsym = size(RI,3) ! careful, this only works if R has been reshaped from 3x3xNrot to 3x3xNsym

    k_fbz = kx
    k_fbz = ky
    k_fbz = kz

    eps_ = 1.0d-4 ! default thershold
    if (present(eps)) eps_ = eps

    found = .false.
    if (ik <= NkI) then
      iR1 = 1
      iK1 = ik
      found = .true. ! k was already in the IBZ
    else
      symmetry_loop: do  i = 2, Nsym
        do l=1,3
          K(l) = sum ( RI(1:3,l,i)*k_fbz(1:3) )
        end do
        ! K(1) = sum (RI(1:3,1,i)*k_fbz(1:3) ) ! Kx
        ! K(2) = sum (RI(1:3,2,i)*k_fbz(1:3) ) ! Ky
        ! K(3) = sum (RI(1:3,3,i)*k_fbz(1:3) ) ! Kz
        do  j = 1,NkI
          if ( all ( abs( K(1:3)-kI(1:3,j) ) <= eps_ ) ) then
            found = .true.
            iR1 = i
            iK1 = j
            exit symmetry_loop
          end if
        end do
      end do symmetry_loop
    end if
    if (.not. found) then
      call error('Can not find wave vector K='//int2str(ik)//'in IBZ')
    end if
    ! print *, 'iR1:',iR1,'iK1:', iK1
  end subroutine findKinIBZ


  subroutine findKQinIBZ(KQx, KQy, KQz, kI, ktot, RI, G, iG0, iR2, iK2, eps)
    !! Finds the k-point (KQx,KQy,KQz) in the 1st. Brillouin zone (FBZ) and then the ireducible Brillouin zone (IBZ)
    ! integer,       intent(in)  :: Nsym          !! No. of symmetry operations for the given crystal
    real(kind=dp), intent(in)  :: KQx, KQy, KQz !! K+Q wavevector possibly outside of FBZ
    real(kind=dp), intent(in)  :: kI(:,:)       !! k-points in the IBZ
    real(kind=dp), intent(in)  :: ktot(:,:)     !! k-points in the FBZ
    real(kind=dp), intent(in)  :: G(:,:)        !! G-vectors
    real(kind=dp), intent(in)  :: RI(:,:,:)     !! inverted rotational symmetry tensor (3 x 3 x Nrot)
    integer,       intent(out) :: iG0           !! index of the rec.lattice vector;
    integer,       intent(out) :: iR2           !! index of rotational sym. op. in the RI tensor for K + Q = G0 + R2*K2
    integer,       intent(out) :: iK2           !! index of k-point K2 in K + Q = G0 + R2*K2
    real(kind=dp), intent(in), optional  :: eps !! thershold for two k-points being the same
  
    logical       :: found_ibz, found_fbz
    integer       :: iG, jk, i, j, l
    integer       :: NkI, Ntot, NG, Nsym
    real(kind=dp) :: K(3), KQ(3)
    real(kind=dp) :: eps_

    Nsym = size(RI,3) ! careful this only works if R has been reshaped from 3x3xNrot to 3x3xNsym
    NkI  = size(kI,2)
    Ntot = size(ktot,2)
    NG   = size(G,2)

    KQ(1) = KQx
    KQ(2) = KQy
    KQ(3) = KQz

    eps_ = 1.0d-4 ! default thershold
    if (present(eps)) eps_ = eps
    found_fbz = .false.
    found_ibz = .false.
    iG_loop : do  iG = 1,NG
    write(104,*) G(1:3,iG)
      k_loop_FBZ : do  jk = 1, Ntot
        if ( all (abs(KQ(1:3)-G(1:3,iG)-ktot(1:3,jk)) <= eps_) ) then
          found_fbz = .true.
          iG0 = iG
          symm_loop: do  i = 1, Nsym
            do l=1,3
              K(l) = sum( RI(1:3,l,i) * ktot(1:3,jk) )
            end do
            k_loop_IBZ: do  j = 1, NkI
              if ( all( abs(K(1:3)-kI(1:3,j)) <= eps_) ) then
                found_ibz = .true. ! true for IBZ (and also for FBZ)
                iR2 = i
                iK2 = j
                exit iG_loop
              end if
            end do k_loop_IBZ
          end do symm_loop
        end if
      end do k_loop_FBZ
    end do iG_loop
  
    if (.not. found_fbz) then
      ! print*,'Can not find wave vector K+Q=',ik,'+',iq, 'in FBZ.'
      call error('Can not find wave vector K+Q in FBZ.')
    else if (.not. found_ibz) then
      ! print*,'Can not find wave vector K+Q=',ik,'+',iq, 'in IBZ.'
      call error('Can not find wave vector K+Q='//int2str(iK2)//'in IBZ.')
    end if
  
  end subroutine findKQinIBZ

  subroutine findMinQ(iq, ktot, q, absq)
    !! Searches for the mininimal wavevector \(\mathbf{q}=(qx,qy,qz)\) in the \(\Gamma \to M \) direction
    integer,       intent(in)            :: iq        !! q-vector index
    real(kind=dp), intent(in)            :: ktot(:,:) !! k-points in the FBZ (3 x Nk)
    real(kind=dp), intent(out)           :: q(3)      !! q-vector at index iq
    real(kind=dp), intent(out), optional :: absq      !! norm of q-vector at index iq

    integer       :: i, ikmin
    real(kind=dp) :: kmin, kref, krefM ! , absq
    integer       :: Ntot
    
    Ntot = size(ktot,2)

    kmin = 1.0_dp
    Ntot_loop: do  i = 1, Ntot ! loop over different k-points in FBZ
      kref = sqrt(sum(ktot(1:3,i)**2))
      ! neven debug
      ! print *,'i=',i,' kref: ',kref
      if (kref == 0.0_dp) then
        cycle Ntot_loop
      else if (kref < kmin) then
        kmin = kref
        ikmin = i
        krefM = kmin
        ! neven debug
        print *,'i=',i,'kmin: ',kmin,'ktot(ikmin):',ktot(1:3,ikmin)
      end if
    end do Ntot_loop
    ! neven debug
    ! print *,'ikmin=',ikmin,'kmin=',kmin,'ktot(1:3,ikmin)',ktot
    ! qx = (iq-1) * ktot(1,ikmin)
    ! qy = (iq-1) * ktot(2,ikmin)
    ! qz = (iq-1) * ktot(3,ikmin)
    ! absq = sqrt(qx**2 + qy**2 + qz**2)

    q(1:3) = (iq-1) * ktot(1:3,ikmin)
    if (present(absq)) absq = sqrt( sum( q(1:3)**2 ) )
  end subroutine findMinQ

  subroutine genFBZpath(sympts, ktot, eps, kgmkg, writeOutput)
    !! Generates high-symmetry path along G->M->K->G
    character(len=*),  intent(in)  :: sympts      !! labels of high-symmetry points
    real(kind=dp),     intent(in)  :: ktot(:,:)   !! k-points in the FBZ
    real(kind=dp),     intent(in)  :: eps         !! threshold for two k-points being the same
    integer,           intent(out) :: kgmkg(:)    !! k-point indices of high-symmetry path in the FBZ
    logical, optional, intent(in)  :: writeOutput !! write high-symmetry path to file yes/no

    logical :: writeOutput_ = .false.
    integer :: iuni
    integer :: i, ik, Ntot, Ngmkg

    Ntot = size(ktot,2)
    Ngmkg = 0

    if (present(writeOutput)) writeOutput_ = writeOutput

    if (sympts/='GMKG') then
      call error('Requested high-symmetry path not implemented.')
    end if

    ! G -> M
    print *,'G -> M'
    do  i = 1,Ntot
      if(ktot(1,i) == 0.0_dp .and. ktot(2,i) >= 0.0_dp) then
        Ngmkg = Ngmkg + 1
        kgmkg(Ngmkg)=i
        print *,'Ngmkg=',Ngmkg,'i_kI=',i
      end if
    end do

    ! M -> K
    print *,'M -> K'
    do i = 1,Ntot
      if (ktot(1,i) > 0.0_dp .and. ktot(2,i) > 0.0_dp) then
        if (abs(ktot(2,i)-1.0_dp/sqrt(3.0_dp)) < eps) then
          Ngmkg = Ngmkg + 1
          kgmkg(Ngmkg) = i
          print *,'Ngmkg=',Ngmkg,'i_kI=',i
        end if
      end if
    end do

    ! K -> G
    print *,'K->G'
    do i = Ntot,1,-1
      if (ktot(1,i) > 0.0_dp .AND. ktot(2,i) > 0.0_dp) then
        if (abs(ktot(2,i)/ktot(1,i)-sqrt(3.0_dp)) < eps) then
          if (abs(ktot(2,i)-1.0_dp/sqrt(3.0_dp)) > eps) then
            Ngmkg = Ngmkg + 1
            kgmkg(Ngmkg) = i
            print *,'Ngmkg=',Ngmkg,'i_kI=',i
          end if
        end if
      end if
    end do
    
    ! outputs 2D path in the FBZ
    if (writeOutput_) then
      open(newunit=iuni,file='highsymm_FBZ_path.dat',action='write',status='new')
      do i = 1,Ngmkg
        ik = kGMKG(i)
        write(iuni,*) i, ktot(1,ik), ktot(2,ik)
      enddo
      close(iuni)
    end if

  end subroutine genFBZpath

end module brillouin_zone