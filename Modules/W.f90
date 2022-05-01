module Wmod
  use iso_fortran_env, only: dp => real64
  implicit none
  
  public :: genChargeVertices
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

end module Wmod
