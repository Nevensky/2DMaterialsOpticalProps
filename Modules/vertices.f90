module vertices
  use iso_fortran_env, only: dp => real64
  use notifications, only: error
  implicit none
  
  public :: genChargeVertices, genCurrentVertices, genChargeVertices_v2,genCurrentVertices_v2
  private

contains
subroutine genChargeVertices(n,m,K1,K2,Npw_k1, Npw_k2, jump, eps, iG0, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2)
  !! Generates charge vertex matrix elements \( \rho_{nm,\mathbf{k,k'}}(\mathbf{G}) \)
  integer,          intent(in)    :: Npw_k1, Npw_k2
  integer,          intent(in)    :: iG0, K1,K2, n, m
  integer,          intent(in)    :: R1, R2     !! indices of rotational symmetry needed to obtain \(\mathbf{k,k'}\)
  real(kind=dp),    intent(in)    :: eps        !! convergence cutoff between G-vectors
  real(kind=dp),    intent(in)    :: R(:,:,:)   !! rotation matrices for all crystal symmetries (3x3xNrot)
  real(kind=dp),    intent(in)    :: RI(:,:,:)  !! inverted rotation matrices for all crystal symmetries
  real(kind=dp),    intent(in)    :: Glf(:,:)   !! local field vectors
  real(kind=dp),    intent(in)    :: G(:,:)     !! G-vectors of the reciprocal. lattice
  complex(kind=dp), intent(in)    :: C1(:)      !! Fourier coefficients \( c_{n,\mathbf{k}}(\mathbf{G'}) \)
  complex(kind=dp), intent(in)    :: C2(:)      !! Fourier coefficients \( c_{m,\mathbf{k}'}(\mathbf{G''}) \)
  logical,          intent(inout) :: jump       !! if true skips recomputing \(G''\) indices
  integer,          intent(inout) :: Gfast(:)   !! arrary for storing indices of \(G'' \)
  complex(kind=dp), intent(out)   :: MnmK1K2(:) !! charge vertices \( \rho_{nm,\mathbf{k,k'}} \)

  integer :: i, iG, iG1, iG2
  integer :: iGfast
  integer :: Nlf, NG1, NG2
  ! real(kind=dp) :: K11, K22, K33
  ! real(kind=dp) :: Gxx1,Gyy1,Gzz1
  ! real(kind=dp) :: Gxx2,Gyy2,Gzz2
  real(kind=dp) :: Gprime(3), K(3)
  Nlf = size(MnmK1K2,1)
  ! NG1 = size(C1,1)  ! = npwx instead of npw(ik)
  ! NG2 = size(C2,1)  ! = npwx instead of npw(ik)
  ! NG2 = size(Gfast,1) ! larger than needed, requires extra if 
  NG1 = Npw_k1 ! No. of planewaves in C1
  NG2 = Npw_k2 ! No. of planewaves in C2

  iGfast = 0
  MnmK1K2(1:Nlf) = dcmplx(0.0_dp,0.0_dp) 
  do iG = 1, Nlf      ! sum over local field vectors
    do iG1 = 1, NG1   ! neven debug: NGd repalced by NG1
      iGfast = iGfast + 1

      K = 0
      do i=1,3
        ! K11, K22, K33
        K(i) = sum ( R(i,1:3,R1)*G(1:3,iG1) ) 
        K(i) = K(i) + Glf(i,iG) + G(i,iG0)
      end do

      Gprime = 0
      do i=1,3
        ! Gxx1, Gyy1, Gzz1
        Gprime(i) = sum ( RI(i,1:3,R2)*K(1:3) ) 
      end do

      if (jump == .false.) then
        iG2_loop: do iG2 = 1,NG2
          Gfast(iGfast) = NG2+1
          if ( all( abs(G(1:3,iG2) - Gprime(1:3)) < eps) ) then
              Gfast(iGfast) = iG2
              exit iG2_loop
          end if
        end do iG2_loop
      end if
      iG2 = Gfast(iGfast)

      if (iG2 <= NG2) then 
        ! print *,'shape C1 C2',shape(C1),shape(C2)
        ! print *,'iG=',iG,'iG1=',iG1,'iG2=',iG2,'npwk=',Npw_k
        MnmK1K2(iG) = MnmK1K2(iG) + conjg(C1(iG1))*C2(iG2)
      end if

    end do
  end do
  jump = .true.
  ! if (K1 == K2) then
  ! print *, 'n,m=',n,m
  ! print *, ' K1, K2=',K1, K2
  ! print *, ' MnmK1K2(1)=',MnmK1K2(1)
  ! ! stop
  ! end if
end subroutine genChargeVertices


subroutine genChargeVertices_v2(n,m,Npw_k1,Npw_k2,jump, eps, iG0, K1, K2, R1, R2, R, RI, Glf, G, igk_k, Gfast, C1, C2, MnmK1K2)
  !! Generates charge vertex matrix elements \( \rho_{nm,\mathbf{k,k'}}(\mathbf{G}) \)
  integer, intent(in) :: Npw_k1, Npw_k2, n, m
  integer,          intent(in)    :: iG0
  integer,          intent(in)    :: K1, K2     !! indices of K and K+Q k-points in the IBZ
  integer,          intent(in)    :: R1, R2     !! indices of rotational symmetry needed to obtain \(\mathbf{k,k'}\)
  integer,          intent(in)    :: igk_k(:,:) !! k+G to G correspondence (Npw x NkI)
  real(kind=dp),    intent(in)    :: eps        !! convergence cutoff between G-vectors
  real(kind=dp),    intent(in)    :: R(:,:,:)   !! rotation matrices for all crystal symmetries (3x3xNrot)
  real(kind=dp),    intent(in)    :: RI(:,:,:)  !! inverted rotation matrices for all crystal symmetries
  real(kind=dp),    intent(in)    :: Glf(:,:)   !! local field vectors
  real(kind=dp),    intent(in)    :: G(:,:)     !! G-vectors of the reciprocal. lattice
  complex(kind=dp), intent(in)    :: C1(:)      !! Fourier coefficients \( c_{n,\mathbf{k}}(\mathbf{G'}) \)
  complex(kind=dp), intent(in)    :: C2(:)      !! Fourier coefficients \( c_{m,\mathbf{k}'}(\mathbf{G''}) \)
  logical,          intent(inout) :: jump       !! if true skips recomputing \(G''\) indices
  integer,          intent(inout) :: Gfast(:)   !! arrary for storing indices of \(G'' \)
  complex(kind=dp), intent(out)   :: MnmK1K2(:) !! charge vertices \( \rho_{nm,\mathbf{k,k'}} \)

  integer :: i, iG, iG1, iG2
  integer :: iGfast
  integer :: Nlf, NG1, NG2
  real(kind=dp) :: Gprime(3), K(3)

  Nlf = size(MnmK1K2,1)
  ! NG1 = size(C1,1)
  ! NG2 = size(C2,1)
  ! NG2 = size(igk_k,1) ! = npwx (this prevents igk_k going out of bounds)
  NG1 = Npw_k1
  NG2 = Npw_k2

  ! print *,'NG1=',NG1,'Npw(ik)=',Npw_k

  iGfast = 0
  MnmK1K2(1:Nlf) = dcmplx(0.0_dp,0.0_dp) 
  iG_loop: do iG = 1, Nlf      ! suma po lokalnim fieldovima kojih ima Nlf
    iG1_loop: do iG1 = 1, NG1 ! neven debug: zamjenjen vitin NGd sa NG1 (je li ok?)
      iGfast = iGfast + 1

      ! print *, 'iG1=',iG1,'igk_k(iG1,K1)=',igk_k(iG1,K1)
      ! if (igk_k(iG1,K1)==0) then
      !   ! print *, 'iG1=',iG1,'K1=',K1,'NG1=',NG1,'igk_k=',igk_k(iG1,K1)
      !   ! print *,'DEBUG: skipeed iG1_loop'
      !   cycle iG_loop
      ! end if

      K = 0
      do i=1,3
        ! K11, K22, K33
        ! print *, G(1:3,igk_k(iG1,K1))
        ! K(i) = sum ( R(i,1:3,R1)*G(1:3,iG1) )
        ! G(1:3,iG1) --> G(1:3,igk_k(iG1,K1)
        K(i) = sum ( R(i,1:3,R1)*G(1:3,igk_k(iG1,K1)) )
        K(i) = K(i) + Glf(i,iG) + G(i,iG0)
      end do

      Gprime = 0
      do i=1,3
        ! Gxx1, Gyy1, Gzz1
        Gprime(i) = sum ( RI(i,1:3,R2)*K(1:3) ) 
      end do

      if (jump == .false.) then
        iG2_loop: do iG2 = 1,NG2
          ! if (igk_k(iG2,K2)==0) then
            ! print *, 'iG2=',iG2,'K2=',K2,'NG2=',NG2,'igk_k=',igk_k(iG2,K2)
            ! print *,'DEBUG: skipeed iG2_loop'
            ! stop
            ! cycle iG1_loop
          ! end if
          ! write(208,*) 'iG2=',iG2,'K2=',K2
          ! write(208,*) 'igk_k(iG2,K2)=',igk_k(iG2,K2)
          Gfast(iGfast) = NG2 + 1
          ! if ( all( abs(G(1:3,iG2) - Gprime(1:3)) < eps) ) then
          ! G(1:3,iG2) --> G(1:3,igk_k(iG2,K2)
          if ( all( abs(G(1:3,igk_k(iG2,K2)) - Gprime(1:3)) < eps) ) then
              Gfast(iGfast) = iG2
              exit iG2_loop
          end if
        end do iG2_loop
      end if
      iG2 = Gfast(iGfast)
      ! print *,'jump=',jump,'iG2=',iG2,'iGfast=',iGfast

      if (iG2==0) then
        print *, 'n,m,K1,K2,iG1,iG2=',n,m,K1,K2,iG1,iG2
        call error('iG2 is zero!')
      end if
      if (iG2 <= NG2) then ! .and. iG2/=0 --> if igk_k has zero elements
        ! print *,'shape c2,c2=',shape(C1),shape(C2)
        ! print *,'C1*=',conjg(C1(iG1))
        ! print *,'C2=',C2(iG2)
        ! print *,'mnmk1k2=',MnmK1K2(iG)
        MnmK1K2(iG) = MnmK1K2(iG) + conjg(C1(iG1))*C2(iG2)
      end if

    end do iG1_loop
  end do iG_loop
  ! print *, 'n,m=',n,m
  ! print *, ' K1, K2=',K1, K2
  ! print *, ' MnmK1K2(1)=',MnmK1K2(1)
  jump = .true.
end subroutine genChargeVertices_v2


subroutine genCurrentVertices(pol, jump, eps, Gcar, qx,qy,qz, kx,ky,kz, iG0, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2, MnmK1K22)
  !! Construct matrix elements for current vertices \( j^\mu_{n\mathbf{K},m\mathbf{K'}}(\mathbf{G}) \)
  !! stored in MnmK1K2(iG) i MnmK1K22(iG) matrices respectively
  character(len=3), intent(in)    :: pol
  integer,          intent(in)    :: iG0 !, Nlf, NG1, NG2, NGd
  integer,          intent(in)    :: R1,R2
  real(kind=dp),    intent(in)    :: eps
  real(kind=dp),    intent(in)    :: Gcar
  real(kind=dp),    intent(in)    :: qx,qy,qz !! transfer wavevector
  real(kind=dp),    intent(in)    :: kx,ky,kz !! initial wavevector 
  real(kind=dp),    intent(in)    :: R(:,:,:) !! point group rotational transformation matrices (3x3xNrot)
  real(kind=dp),    intent(in)    :: RI(:,:,:) !! inverse of point group rot. transform. mat. (3x3xNrot)
  real(kind=dp),    intent(in)    :: Glf(:,:) 
  real(kind=dp),    intent(in)    :: G(:,:)  !! reciprocal lattice vectors \( \mathbf{G}_i \) for the density
  complex(kind=dp), intent(in)    :: C1(:)   !! Fourier coefficients of Bloch states for a given band,k-point pair (NG)
  complex(kind=dp), intent(in)    :: C2(:)   !! Fourier coefficients of Bloch states for a given band,k-point pair (NG)
  logical,          intent(inout) :: jump
  integer,          intent(inout) :: Gfast(:)
  complex(kind=dp), intent(out)   :: MnmK1K2(:) !! charge vertices for 1st polarization component (NG)
  complex(kind=dp), intent(out), optional :: MnmK1K22(:)  !! charge vertices for 2nd polarization component (NG)

  integer          :: i, iG, iG1, iG2
  integer          :: iGfast
  integer          :: Nlf, NG1, NG2
  ! real(kind=dp)    :: K11, K22, K33
  ! real(kind=dp)    :: Gxx1,Gyy1,Gzz1
  ! real(kind=dp)    :: Gxx2,Gyy2,Gzz2
  real(kind=dp)    :: Gprime(3)
  ! real(kind=dp)    :: current, current_y, current_z
  real(kind=dp)    :: k(3), q(3), K_(3)
  real(kind=dp)    :: current_xyz(3)

  Nlf = size(MnmK1K2,1)
  NG1 = size(C1,1)
  NG2 = size(C2,1)

  k(1) = kx
  k(2) = ky
  k(3) = kz

  q(1) = qx
  q(2) = qy
  q(3) = qz
  
  iGfast = 0
  MnmK1K2(1:Nlf)  = dcmplx(0.0_dp,0.0_dp)
  if (present(MnmK1K22)) MnmK1K22(1:Nlf) = dcmplx(0.0_dp,0.0_dp)

  iG_loop: do  iG = 1, Nlf
    iG1_loop: do iG1 = 1,NG1 ! debug neven: zamjeneno NGd sa NG1
      iGfast = iGfast + 1

      do i=1,3
        K_(i) = sum( R(i,1:3,R1)*G(1:3,iG1) )
      end do

      current_xyz(1:3) = Gcar * ( q(1:3) + 2.0_dp*k(1:3) + Glf(1:3,iG) + 2.0*K_(1:3) )
      K_(1:3) = K_(1:3) + Glf(1:3,iG) + G(1:3,iG0)
      
      do i=1,3
        Gprime(i) = sum ( RI(i,1:3,R2)*K_(i) )
      end do
      
      if (jump == .false.) then
        iG2_loop: do iG2 = 1,NG2
          Gfast(iGfast) = NG2 + 1
          if ( all( abs(G(1:3,iG2)-Gprime(1:3)) < eps) ) then
            Gfast(iGfast) = iG2
            exit iG2_loop
          end if
        end do iG2_loop
      end if

      iG2 = Gfast(iGfast)

      if (iG2 <= NG2) then
        if (pol == 'xx') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5_dp*conjg(C1(iG1)) * current_xyz(1) * C2(iG2)
          if (present(MnmK1K22)) MnmK1K22(iG) = MnmK1K2(iG) ! current vertices are the same
        elseif (pol == 'yy') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5_dp*conjg(C1(iG1)) * current_xyz(2) * C2(iG2)
          if (present(MnmK1K22)) MnmK1K22(iG) = MnmK1K2(iG) ! current vertices are the same
        elseif (pol == 'zz') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5_dp*conjg(C1(iG1)) * current_xyz(3) * C2(iG2)
          if (present(MnmK1K22))  MnmK1K22(iG) = MnmK1K2(iG) ! current vertices are the same
        elseif (pol =='yz' .and. present(MnmK1K22)) then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5_dp*conjg(C1(iG1)) * current_xyz(2) * C2(iG2)
          MnmK1K22(iG) = MnmK1K22(iG) + 0.5_dp*conjg(C1(iG1)) * current_xyz(3) * C2(iG2)
        elseif (pol =='zy' .and. present(MnmK1K22)) then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5_dp*conjg(C1(iG1)) * current_xyz(3) * C2(iG2)
          MnmK1K22(iG) = MnmK1K22(iG) + 0.5_dp*conjg(C1(iG1)) * current_xyz(2) * C2(iG2)
        else
          call error('Specified mixed polarization component not supported.'//adjustl(trim(pol))//' not allowed.')
        end if
      end if
    end do iG1_loop
  end do iG_loop
  jump = .true. ! for each wave-vector q and its corresponding k, store Gfast and skip the procedure for each future band (outside this func.) 
    
end subroutine genCurrentVertices

subroutine genCurrentVertices_v2(pol, jump, eps, Gcar, qx, qy, qz, kx, ky, kz, iG0, K1, K2, R1, R2, R, RI, Glf, G, igk_k, Gfast, C1, C2, MnmK1K2, MnmK1K22)
  !! Construct matrix elements for current vertices \( j^\mu_{n\mathbf{K},m\mathbf{K'}}(\mathbf{G}) \)
  !! stored in MnmK1K2(iG) i MnmK1K22(iG) matrices respectively
  character(len=3), intent(in)    :: pol
  integer,          intent(in)    :: iG0 !, Nlf, NG1, NG2, NGd
  integer,          intent(in)    :: K1, K2     !! indices of K and K+Q k-points in the IBZ
  integer,          intent(in)    :: R1,R2
  real(kind=dp),    intent(in)    :: eps
  real(kind=dp),    intent(in)    :: Gcar
  real(kind=dp),    intent(in)    :: qx,qy,qz !! transfer wavevector
  real(kind=dp),    intent(in)    :: kx,ky,kz !! initial wavevector 
  integer,          intent(in)    :: igk_k(:,:) !! k+G to G correspondence (Npw x NkI)
  real(kind=dp),    intent(in)    :: R(:,:,:) !! point group rotational transformation matrices (3x3xNrot)
  real(kind=dp),    intent(in)    :: RI(:,:,:) !! inverse of point group rot. transform. mat. (3x3xNrot)
  real(kind=dp),    intent(in)    :: Glf(:,:) 
  real(kind=dp),    intent(in)    :: G(:,:)  !! reciprocal lattice vectors \( \mathbf{G}_i \) for the density
  complex(kind=dp), intent(in)    :: C1(:)   !! Fourier coefficients of Bloch states for a given band,k-point pair (NG)
  complex(kind=dp), intent(in)    :: C2(:)   !! Fourier coefficients of Bloch states for a given band,k-point pair (NG)
  logical,          intent(inout) :: jump
  integer,          intent(inout) :: Gfast(:)
  complex(kind=dp), intent(out)   :: MnmK1K2(:) !! charge vertices for 1st polarization component (NG)
  complex(kind=dp), intent(out), optional :: MnmK1K22(:)  !! charge vertices for 2nd polarization component (NG)

  integer          :: i, iG, iG1, iG2
  integer          :: iGfast
  integer          :: Nlf, NG1, NG2
  ! real(kind=dp)    :: K11, K22, K33
  ! real(kind=dp)    :: Gxx1,Gyy1,Gzz1
  ! real(kind=dp)    :: Gxx2,Gyy2,Gzz2
  real(kind=dp)    :: Gprime(3)
  ! real(kind=dp)    :: current, current_y, current_z
  real(kind=dp)    :: k(3), q(3), K_(3)
  real(kind=dp)    :: current_xyz(3)

  Nlf = size(MnmK1K2,1)
  NG1 = size(C1,1)
  NG2 = size(C2,1)

  k(1) = kx
  k(2) = ky
  k(3) = kz

  q(1) = qx
  q(2) = qy
  q(3) = qz
  
  iGfast = 0
  MnmK1K2(1:Nlf)  = dcmplx(0.0_dp,0.0_dp)
  if (present(MnmK1K22)) MnmK1K22(1:Nlf) = dcmplx(0.0_dp,0.0_dp)

  iG_loop: do  iG = 1, Nlf
    iG1_loop: do iG1 = 1,NG1 ! debug neven: zamjeneno NGd sa NG1
      iGfast = iGfast + 1

      do i=1,3
        K_(i) = sum( R(i,1:3,R1)*G(1:3,igk_k(iG1,K1)) )
      end do

      current_xyz(1:3) = Gcar * ( q(1:3) + 2.0*k(1:3) + Glf(1:3,iG) + 2.0*K_(1:3) )
      K_(1:3) = K_(1:3) + Glf(1:3,iG) + G(1:3,iG0)
      
      do i=1,3
        Gprime(i) = sum ( RI(i,1:3,R2)*K_(i) )
      end do
      
      if (jump == .false.) then
        iG2_loop: do iG2 = 1,NG2
          Gfast(iGfast) = NG2 + 1
          if ( all( abs(G(1:3,igk_k(iG2,K2))-Gprime(1:3)) < eps) ) then
            Gfast(iGfast) = iG2
            exit iG2_loop
          end if
        end do iG2_loop
      end if

      iG2 = Gfast(iGfast)

      if (iG2 <= NG2) then
        if (pol == 'xx') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5_dp*conjg(C1(iG1)) * current_xyz(1) * C2(iG2)
          if (present(MnmK1K22)) MnmK1K22(iG) = MnmK1K2(iG) ! current vertices are the same
        elseif (pol == 'yy') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5_dp*conjg(C1(iG1)) * current_xyz(2) * C2(iG2)
          if (present(MnmK1K22)) MnmK1K22(iG) = MnmK1K2(iG) ! current vertices are the same
        elseif (pol == 'zz') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5_dp*conjg(C1(iG1)) * current_xyz(3) * C2(iG2)
          if (present(MnmK1K22))  MnmK1K22(iG) = MnmK1K2(iG) ! current vertices are the same
        elseif (pol =='yz' .and. present(MnmK1K22)) then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5_dp*conjg(C1(iG1)) * current_xyz(2) * C2(iG2)
          MnmK1K22(iG) = MnmK1K22(iG) + 0.5_dp*conjg(C1(iG1)) * current_xyz(3) * C2(iG2)
        elseif (pol =='zy' .and. present(MnmK1K22)) then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5_dp*conjg(C1(iG1)) * current_xyz(3) * C2(iG2)
          MnmK1K22(iG) = MnmK1K22(iG) + 0.5_dp*conjg(C1(iG1)) * current_xyz(2) * C2(iG2)
        else
          call error('Specified mixed polarization component not supported.'//adjustl(trim(pol))//' not allowed.')
        end if
      end if
    end do iG1_loop
  end do iG_loop
  jump = .true. ! for each wave-vector q and its corresponding k, store Gfast and skip the procedure for each future band (outside this func.) 
    
end subroutine genCurrentVertices_v2

end module vertices
