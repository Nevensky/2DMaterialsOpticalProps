module W
  use iso_fortran_env, only: dp => real64
  implicit none
  
  public :: genChargeVertices
  private

contains
subroutine genChargeVertices(jump, eps, Nlf, iG0, NG1, NG2, NGd, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2)
  ! Construction of charge vertex matrix elements MnmK1K2(G) 
  integer,          intent(in)    :: iG0, Nlf, NG1, NG2, NGd
  integer,          intent(in)    :: R1,R2
  real(kind=dp),    intent(in)    :: eps
  real(kind=dp),    intent(in)    :: R(:,:,:)
  real(kind=dp),    intent(in)    :: RI(:,:,:)
  real(kind=dp),    intent(in)    :: Glf(:,:)
  real(kind=dp),    intent(in)    :: G(:,:)     ! polje valnih vektora G u recp. prost. za wfn.
  complex(kind=dp), intent(in)    :: C1(:)
  complex(kind=dp), intent(in)    :: C2(:)
  logical,          intent(inout) :: jump
  integer,          intent(inout) :: Gfast(:)
  complex(kind=dp), intent(out)   :: MnmK1K2(:)

  integer       :: i, iG, iG1, iG2
  integer       :: iGfast
  ! real(kind=dp) :: K11, K22, K33
  ! real(kind=dp) :: Gxx1,Gyy1,Gzz1
  ! real(kind=dp) :: Gxx2,Gyy2,Gzz2
  real(kind=dp) :: Gprime(3), K(3)

  iGfast = 0
  MnmK1K2(1:Nlf) = cmplx(0.0_dp,0.0_dp) ! nabojni vrhovi
  do  iG = 1,Nlf     ! suma po lokalnim fieldovima kojih ima Nlf
    do  iG1 = 1,NGd  ! vito zamjenjeno NGd sa NG1
      iGfast = iGfast + 1

      do i=1,3
        ! K11, K22, K33
        K(i) = sum ( R(R1,i,1:3)*G(1:3,iG1) ) + Glf(i,iG) + G(i,iG0)
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

end module W