!! This modules reads the output of a completed Wannier90 calculation
module io_wannier90
  use iso_fortran_env, only : iostat_end, dp => real64
  implicit none
  public :: loadE, loadU
  private 

contains

  subroutine loadE(path)
    character(len=300), intent(in) :: path
    continue
  end subroutine loadE

  subroutine loadU(path, k, U)
    !! Reads Wannier90 .mat files containing the unitary transformation
    !! matrices from the Bloch to the Wannier basis at each k-point
    character(len=300), intent(in) :: path
    real(dp),    allocatable, intent(out) :: k(:,:)   ! k-points IBZ or FBZ??
    complex(dp), allocatable, intent(out) :: U(:,:,:) !! (Nk x Nwan x Nband)

    integer  :: iuni
    integer  :: Nk, Nwan, Nband
    integer  :: ik, iwan, n
    real(dp) :: ReU, ImU

    open(newunit=iuni, file=trim(path), action='read', status='old')
    read(iuni,*) ! skip line
    read(iuni,'(3I12)') Nk, Nwan, Nband
    allocate(k(3,Nk))
    allocate(U(Nk,Nwan,Nband))
    do ik=1,Nk
      read(iuni,*) ! skip line
      read(iuni,('(3F15.10)')) k(1:3,ik)
      print *,'ik:',ik
      print *, k(:,ik)
      print *, 'WARNING: loadU() possibly we have U(Nk x Nband x Nwan) column major order. DOUBLE CHECK!'
      do iwan=1,Nwan
        do n=1,Nband
          read(iuni,'(2F15.10)') ReU, ImU
          U(ik,iwan,n) = cmplx(ReU, ImU)
        end do
      end do
    end do
    close(iuni)

    ! print *, 'Nk: 'Nk,'Nwan: ', Nwan,'Nband:', Nband
    ! deallocate(k)
    ! deallocate(U)
  end subroutine loadU


end module io_wannier90