module matrix_inverse
!! Returns the inverse of a matrix calculated by finding the
!! LU decomposition. Depends on LAPACK, preferably Intel MKL.

   use lapack95, only: getrf, getri
   use iso_fortran_env, only: dp => real64, idp => int16

   implicit none

   public :: invert, invert_cmplx, invert_real
   private

   interface invert
   !! Inverts real or complex matrix (double precision)
      procedure invert_cmplx, invert_real
   end interface invert

contains
  subroutine invert_cmplx(A)
    !! Inverts complex matrix (double precision)
    complex(dp), dimension(:,:), intent(inout) :: A

    complex(dp), allocatable :: work(:)         ! work array for LAPACK
    integer                  :: ipiv(size(A,1)) ! pivot indices
    integer                  :: n, info, ios 
    integer(idp)             :: lwork           ! work array size
    ! integer(idp), parameter :: memmax = 60*sizeof(dp) 

    ! External procedures defined in LAPACK for double complex matrices
    external zgetrf
    external zgetri

    n = size(A,1) ! number of rows of A

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call zgetrf(n, n, A, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular! (zgetrf)'
    end if


    ! ZGETRI computes the inverse of a matrix using 
    ! the LU factorization computed by ZGETRF.

    ! Query optimal work array for best performance
    lwork = -1
    allocate(work(1))

    call zgetri(n, A, n, ipiv, work, lwork, info) ! query optimal lwork value
    lwork = work(1)
    deallocate(work)
    !if (sizeof(lwork)>lwork) lwork = memmax
    allocate(work(lwork),stat=ios)
    if (ios /= 0) then
      print *,'ERROR: Failed allocation of work array. (zgetri) ios=',ios
      stop
    end if
    ! print *, " sizeof(A) = ",sizeof(A)/1024.0**3,"GB"
    ! print *, " sizeof(work) = ",sizeof(work)/1024.0**3,"GB"

    ! Perform inverse with the optimal work array
    call zgetri(n, A, n, ipiv, work, lwork, info)

    if (info /= 0) then
       stop 'ERROR: Matrix inversion failed! (zgetri)'
    end if

    deallocate(work)
  end subroutine invert_cmplx

  subroutine invert_real(A)
    !! Inverts real matrix (double precision)
    real(dp), dimension(:,:), intent(inout) :: A

    real(dp), allocatable :: work(:)          ! work array for LAPACK
    integer               :: ipiv(size(A,1))  ! pivot indices
    integer               :: n, info, info2, info3, ios 
    integer(idp)          :: lwork            ! work array size
    ! integer(idp), parameter :: memmax = 30*sizeof(dp) 

    ! External procedures defined in LAPACK for double real matrices
    external dgetrf
    external dgetri

    n = size(A,1) ! number of rows of A

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(n, n, A, n, ipiv, info)

    if (info /= 0) then
       print *, 'Matrix is numerically singular! (dgetrf) info=', info
       stop
    end if

    ! DGETRI computes the inverse of a matrix using 
    ! the LU factorization computed by DGETRF.

    ! Query optimal work array for best performance
    allocate(work(1))
    call dgetri(n, A, n, ipiv, work, -1, info2) ! query optimal lwork value
    if (info2 /= 0) then
       print *, 'ERROR: Workspace query failed! (dgetri) info2=', info2
       stop
    end if

    lwork = work(1)
    deallocate(work)
    !if (sizeof(lwork)>lworkmax) lwork = max_available_memory
    allocate(work(lwork),stat=ios)
    if (ios /= 0) then
            print *,'ERROR: Failed allocation of work array. (dgetri) ios=',ios
            stop
    end if

    ! print *, " sizeof(A) = ",sizeof(A)/1024.0**3,"GB"
    ! print *, " sizeof(work) = ",sizeof(work)/1024.0**3,"GB"

    ! Perform inverse with the optimal work array
    call dgetri(n, A, n, ipiv, work, lwork, info3)

    if (info3 /= 0) then
       print *, 'ERROR: Matrix inversion failed! (dgetri) info=', info3
       stop
    end if

    deallocate(work)
  end subroutine invert_real

end module matrix_inverse
