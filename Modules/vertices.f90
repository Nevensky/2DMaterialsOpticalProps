module vertices
  use iso_fortran_env, only: dp => real64
  implicit none
  
  public :: genChargeVertices, genCurrentVertices
  private

contains
subroutine genChargeVertices(jump, eps, Nlf, iG0, NG1, NG2, NGd, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2)
  !! Generates charge vertex matrix elements \( \rho_{nm,\mathbf{k,k'}}(\mathbf{G}) \)
  integer,          intent(in)    :: iG0, Nlf, NG1, NG2, NGd
  integer,          intent(in)    :: R1, R2     !! indices of rotational symmetry needed to obtain \(\mathbf{k,k'}\)
  real(kind=dp),    intent(in)    :: eps        !! convergence cutoff between G-vectors
  real(kind=dp),    intent(in)    :: R(:,:,:)   !! rotation matrices for all crystal symmetries
  real(kind=dp),    intent(in)    :: RI(:,:,:)  !! inverted rotation matrices for all crystal symmetries
  real(kind=dp),    intent(in)    :: Glf(:,:)   !! local field vectors
  real(kind=dp),    intent(in)    :: G(:,:)     !! G-vectors of the reciprocal. lattice
  complex(kind=dp), intent(in)    :: C1(:)      !! Fourier coefficients \( c_{n,\mathbf{k}}(\mathbf{G'}) \)
  complex(kind=dp), intent(in)    :: C2(:)      !! Fourier coefficients \( c_{m,\mathbf{k}'}(\mathbf{G''}) \)
  logical,          intent(inout) :: jump       !! if true skips recomputing \(G''\) indices
  integer,          intent(inout) :: Gfast(:)   !! arrary for storing indices of \(G'' \)
  complex(kind=dp), intent(out)   :: MnmK1K2(:) !! charge vertices \( \rho_{nm,\mathbf{k,k'}} \)

  integer       :: i, iG, iG1, iG2
  integer       :: iGfast
  ! real(kind=dp) :: K11, K22, K33
  ! real(kind=dp) :: Gxx1,Gyy1,Gzz1
  ! real(kind=dp) :: Gxx2,Gyy2,Gzz2
  real(kind=dp) :: Gprime(3), K(3)

  iGfast = 0
  MnmK1K2(1:Nlf) = cmplx(0.0_dp,0.0_dp) 
  do  iG = 1,Nlf     ! suma po lokalnim fieldovima kojih ima Nlf
    do  iG1 = 1,NGd  ! vito zamjenjeno NGd sa NG1
      iGfast = iGfast + 1

      do i=1,3
        ! K11, K22, K33
        K(i) = sum ( R(R1,i,1:3)*G(1:3,iG1) ) 
        K(i) = K(i) + Glf(i,iG) + G(i,iG0)
      end do

      do i=1,3
        ! Gxx1, Gyy1, Gzz1
        Gprime(i) = sum ( RI(R2,i,1:3)*K(1:3) ) 
      end do

      if (jump == .false.) then
        iG2_loop: do  iG2 = 1,NG2
          Gfast(iGfast) = NG2+1
          if ( all( abs(G(1:3,iG2) - Gprime(1:3)) < eps) ) then
              Gfast(iGfast) = iG2
              exit iG2_loop
          end if
        end do iG2_loop
      end if

      iG2 = Gfast(iGfast)

      if (iG2 <= NG2) then
        MnmK1K2(iG) = MnmK1K2(iG) + conjg(C1(iG1))*C2(iG2)
      else
        print *, "iG = ",iG,"iG1 = ",iG1,"iG2 = ",iG2
        stop 'ERROR: iG2 <= NG2 in genChargeVertices()'
      end if

    end do
  end do
  jump = .true.
end subroutine genChargeVertices

subroutine genCurrentVertices(pol, jump, eps, Gcar, qx,qy,qz, kx,ky,kz, Nlf, iG0, NG1, NG2, NGd, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2, MnmK1K22)
  !! Construct matrix elements for current vertices \( j^\mu_{n\mathbf{K},m\mathbf{K'}}(\mathbf{G}) \)
  !! stored in MnmK1K2(iG) i MnmK1K22(iG) matrices respectively
  character(len=3), intent(in)    :: pol
  integer,          intent(in)    :: iG0, Nlf, NG1, NG2, NGd
  integer,          intent(in)    :: R1,R2
  real(kind=dp),    intent(in)    :: eps
  real(kind=dp),    intent(in)    :: Gcar
  real(kind=dp),    intent(in)    :: qx,qy,qz
  real(kind=dp),    intent(in)    :: kx,ky,kz
  real(kind=dp),    intent(in)    :: R(:,:,:)
  real(kind=dp),    intent(in)    :: RI(:,:,:)
  real(kind=dp),    intent(in)    :: Glf(:,:)
  real(kind=dp),    intent(in)    :: G(:,:)  ! polje valnih vektora G u recp. prost. za wfn.
  complex(kind=dp), intent(in)    :: C1(:)   ! dim NG
  complex(kind=dp), intent(in)    :: C2(:)   ! dim NG (bcs the input is dim iband x NG )
  logical,          intent(inout) :: jump
  integer,          intent(inout) :: Gfast(:)
  complex(kind=dp), intent(out)   :: MnmK1K2(:) 
  complex(kind=dp), intent(out)   :: MnmK1K22(:)

  integer          :: i, iG, iG1, iG2
  integer          :: iGfast
  ! real(kind=dp)    :: K11, K22, K33
  ! real(kind=dp)    :: Gxx1,Gyy1,Gzz1
  ! real(kind=dp)    :: Gxx2,Gyy2,Gzz2
  real(kind=dp)    :: Gprime(3)
  ! real(kind=dp)    :: current, current_y, current_z
  real(kind=dp)    :: k(3), q(3), K_(3)
  real(kind=dp)    :: current_xyz(3)

  k(1) = kx
  k(2) = ky
  k(3) = kz

  q(1) = qx
  q(2) = qy
  q(3) = qz
  
  iGfast = 0
  MnmK1K2(1:Nlf)  = cmplx(0.0_dp,0.0_dp)
  MnmK1K22(1:Nlf) = cmplx(0.0_dp,0.0_dp)

  iG_loop: do  iG = 1,Nlf
    iG1_loop: do iG1 = 1,NGd
      iGfast = iGfast + 1

      do i=1,3
        K_(i) = sum( R(R1,i,1:3)*G(1:3,iG1) )
      end do

      current_xyz(1:3) = Gcar * ( q(1:3) + 2.0*k(1:3) + Glf(1:3,iG) + 2.0*K_(1:3) )
      K_(1:3) = K_(1:3) + Glf(1:3,iG) + G(1:3,iG0)
      
      do i=1,3
        Gprime(i) = sum ( RI(R2,i,1:3)*K_(i) )
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
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5D0*conjg(C1(iG1)) * current_xyz(1) * C2(iG2)
          MnmK1K22(iG) = MnmK1K2(iG) ! current vertices are the same
        elseif (pol == 'yy') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5D0*conjg(C1(iG1)) * current_xyz(2) * C2(iG2)
          MnmK1K22(iG) = MnmK1K2(iG) ! current vertices are the same
        elseif (pol == 'zz') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5D0*conjg(C1(iG1)) * current_xyz(3) * C2(iG2)
          MnmK1K22(iG) = MnmK1K2(iG) ! current vertices are the same
        elseif (pol =='yz') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5D0*conjg(C1(iG1)) * current_xyz(2) * C2(iG2)
          MnmK1K22(iG) = MnmK1K22(iG) + 0.5D0*conjg(C1(iG1)) * current_xyz(3) * C2(iG2)
        elseif (pol =='zy') then
          MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5D0*conjg(C1(iG1)) * current_xyz(3) * C2(iG2)
          MnmK1K22(iG) = MnmK1K22(iG) + 0.5D0*conjg(C1(iG1)) * current_xyz(2) * C2(iG2)
        else
          print *,'WARNING Specified mixed polarization component not supported.'//adjustl(trim(pol))//' not allowed.'
          stop
        end if
      else
        print *, "iG = ",iG,"iG1 = ",iG1,"iG2 = ",iG2
        stop "ERROR: iG2 <= NG2 in genCurrentVertices()"
      end if
    end do iG1_loop
  end do iG_loop
  jump = .true. ! for each wave-vector q and its corresponding k, store Gfast and skip the procedure for each future band (outside this func.) 
    
end subroutine genCurrentVertices


end module vertices
