module brillouin_zone
  use iso_fortran_env, only: dp => real64
  implicit none

  public :: loadKiandE, genFBZ, checkFBZintegration, &
            & findKinIBZ, findKQinBZ, findMinQ
  private

contains

  subroutine loadkIandE(path, NkI, Nband, Nocc, kI, dGW,E)
    ! Loading of all wavevectors (in Cartesiand coords.) in the ireducible BZ and
    ! corresponding eigen-energies form the Quantum Espresso .band files
    integer,            intent(in)    :: NkI
    integer,            intent(in)    :: Nband, Nocc
    character(len=100), intent(in)    :: path
    real(kind=dp),      intent(in)    :: dGW
    real(kind=dp),      intent(inout) :: kI(:,:)
    real(kind=dp),      intent(inout) :: E(:,:)

    integer :: iuni, ios, ik, i
    real(kind=dp),    parameter :: Hartree = 2.0D0*13.6056923D0

    open(newunit=iuni,FILE=path,status='old',err=500,iostat=ios) 
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
      
    goto 400
    500 write(*,*) 'Cannot open BAND file. iostat = ',ios
    stop
    400 continue
    
    ! konverzija en. u Hartree
    E(1:NkI,1:Nband) = E(1:NkI,1:Nband)/Hartree
    ! scissor operator, ispravlja/shifta DFT gap (na 1eV u ovom slucaju)
    E(1:NkI,Nocc+1:Nband) = E(1:NkI,Nocc+1:Nband) + dGW

  end subroutine loadkIandE

  subroutine genFBZ(Nk,NkI,Nsymm,eps,kI,R,Ntot,ktot)
    ! Generates all unique wavectors in the 1st BZ by applying 
    ! point group transformations on the reducible BZ
    integer,       intent(in)  :: Nk, NkI, Nsymm ! No. of k-points, iredducible k-kpoints, symm. ops.
    real(kind=dp), intent(in)  :: eps       ! threshold to distinguish whether k-points are the same
    real(kind=dp), intent(in)  :: kI(:,:)   ! k-points in the irreducible BZ
    real(kind=dp), intent(in)  :: R(:,:,:)  ! point group transformation matrices

    integer,       intent(out) :: Ntot      ! total No. of unique k-point in the 1st BZ
    real(kind=dp), intent(out) :: ktot(:,:) ! unique k-points in the 1st BZ

    integer       :: iuni, ios
    integer       :: it, jk
    integer       :: i
    integer       :: ik, lk
    integer       :: n, m
    real(kind=dp) :: k(3,Nk) ! all non-unique k-points within the 1st BZ

    jk = 0
    Ntot = 0 

    symm_loop: do  i = 1, Nsymm    ! loop over all symmetries
      k_loop_IBZ: do  ik = 1, NkI  ! loop over k points in IBZ
        it = 1
        jk = jk + 1
        do  n = 1, 3   ! loop over kx, ky, kz 
          k(n,jk) = 0.0
          do  m = 1, 3 ! loop over kx, ky, kz
            k(n,jk) = k(n,jk) + R(i,n,m)*kI(m,ik) ! kreira nove k tocke u BZ pomocu simetrije
          end do
        end do 

        if (jk > 1) then
          do  lk = 1, jk-1
            ! Check if the given k-point is unique (skips if it was already added)
            if ( abs(k(1,jk)-k(1,lk)) <= eps .and. &
                 abs(k(2,jk)-k(2,lk)) <= eps .and. &
                 abs(k(3,jk)-k(3,lk)) <= eps ) then 
              it = 2
            end if
          end do
        end if

        if (it == 1) then ! its unique, add it to ktot
          Ntot = Ntot+1
          ktot(1:3,Ntot) = k(1:3,jk)
        end if

      end do k_loop_IBZ
    end do symm_loop

    ! write all 1st BZ k-points to file
    open(newunit=iuni,iostat=ios,file='fbz_check.dat')
    do  i = 1,Ntot
      write(ios,*) ktot(1,i), ktot(2,i)  
    end do
    close(ios)

  end subroutine genFBZ

    subroutine checkFBZintegration(Nband,NkI,Nsymm,Ntot,eps,kI,ktot, RI,Efermi,E,NelQE,Nel)
    ! Checks if the No. of electrons in the 1st BZ (Nel) equals 
    ! the number of electrons in the unit cell as calculated by Quantum Espresso (NelQE)
    integer,       intent(in)  :: NelQE
    integer,       intent(in)  :: NkI, Nsymm, Nband, Ntot
    real(kind=dp), intent(in)  :: eps, Efermi
    real(kind=dp), intent(in)  :: kI(:,:)
    real(kind=dp), intent(in)  :: ktot(:,:)  ! all unique k-points in the FBZ

    real(kind=dp), intent(in)  :: RI(:,:,:)
    real(kind=dp), intent(in)  :: E(:,:)
    real(kind=dp), intent(out) :: Nel

    integer       :: it, ik, n, i, j, K1
    real(kind=dp) :: kx,ky,kz
    real(kind=dp) :: K11, K22, K33

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
                  ! sums electrons in the first band
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
        
        ! sums electrons in the remeaining bands
        if (E(K1,n) < Efermi) then 
          Nel = Nel + 1.0
          ! print *,'Nel',Nel,'band:',n
        end if  
    
      end do band_loop
    end do k_loop_FBZ
    Nel = 2.0*Nel / Ntot ! sums electrons for en. smaller than Efermi
    
  end subroutine checkFBZintegration


  subroutine findKinIBZ(ik, NkI, Nsymm, eps, kx, ky, kz, RI, kI, iR1, iK1)
    ! Finds k-point (kx,ky,kz) in the ireducible BZ
    integer,       intent(in)  :: ik
    integer,       intent(in)  :: NkI, Nsymm
    real(kind=dp), intent(in)  :: eps
    real(kind=dp), intent(in)  :: kx, ky, kz
    real(kind=dp), intent(in)  :: kI(:,:)
    real(kind=dp), intent(in)  :: RI(:,:,:)
    integer      , intent(out) :: iR1, iK1 ! K = R1*K1

    integer       :: i, j
    integer       :: it
    real(kind=dp) :: K11, K22, K33

    it = 1
    if (ik <= NkI) then
      iR1 = 1
      iK1 = ik
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
            iR1 = i
            iK1 = j
            EXIT symmetry_loop
          end if
        end do
      end do symmetry_loop
    end if
    if (it == 1) then
      print *,'Can not find wave vector K=',ik, 'in I.B.Z.'
      stop
    end if
    ! print *, 'iR1:',iR1,'iK1:', iK1
  end subroutine findKinIBZ


  subroutine findKQinBZ(KQx, KQy, KQz, eps, Nsymm, NkI, Ntot, NG, kI, ktot, RI, G, iG0, iR2, iK2)
    ! Finds the k-point (KQx,KQy,KQz) in the 1st. BZ a then in the ireducible BZ
    integer,       intent(in)  :: Nsymm, NkI, Ntot, NG
    real(kind=dp), intent(in)  :: eps
    real(kind=dp), intent(in)  :: KQx, KQy, KQz
    real(kind=dp), intent(in)  :: kI(:,:)
    real(kind=dp), intent(in)  :: ktot(:,:)
    real(kind=dp), intent(in)  :: G(:,:)
    real(kind=dp), intent(in)  :: RI(:,:,:)
    integer,       intent(out) :: iG0, iR2, iK2 ! G0=rec.lattice; K + Q = G0 + R2*K2
  
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
                iR2 = i
                iK2 = j
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
  end subroutine findMinQ

end module brillouin_zone