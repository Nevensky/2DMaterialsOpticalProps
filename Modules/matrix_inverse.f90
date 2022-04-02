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
       stop 'Matrix is numerically singular!'
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
            stop 'ERROR: Failed allocation of work array. (zgetri)'
    end if
    ! print *, " sizeof(A) = ",sizeof(A)/1024.0**3,"GB"
    ! print *, " sizeof(work) = ",sizeof(work)/1024.0**3,"GB"

    ! Perform inverse with the optimal work array
    call zgetri(n, A, n, ipiv, work, lwork, info)

    if (info /= 0) then
       stop 'ERROR: Matrix inversion failed!'
    end if

    deallocate(work)
  end subroutine invert_cmplx

  subroutine invert_real(A)
    !! Inverts real matrix (double precision)
    real(dp), dimension(:,:), intent(inout) :: A

    real(dp), allocatable :: work(:)         ! work array for LAPACK
    integer               :: ipiv(size(A,1)) ! pivot indices
    integer               :: n, info, ios 
    integer(idp)          :: lwork           ! work array size
    ! integer(idp), parameter :: memmax = 30*sizeof(dp) 

    ! External procedures defined in LAPACK for double real matrices
    external dgetrf
    external dgetri

    n = size(A,1) ! number of rows of A

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(n, n, A, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if


    ! DGETRI computes the inverse of a matrix using 
    ! the LU factorization computed by DGETRF.

    ! Query optimal work array for best performance
    lwork = -1
    allocate(work(1))

    call dgetri(n, A, n, ipiv, work, lwork, info) ! query optimal lwork value
    lwork = work(1)
    deallocate(work)
    !if (sizeof(lwork)>lwork) lwork = memmax
    allocate(work(lwork),stat=ios)
    if (ios /= 0) then
            stop 'ERROR: Failed allocation of work array. (zgetri)'
    end if
    ! print *, " sizeof(A) = ",sizeof(A)/1024.0**3,"GB"
    ! print *, " sizeof(work) = ",sizeof(work)/1024.0**3,"GB"

    ! Perform inverse with the optimal work array
    call dgetri(n, A, n, ipiv, work, lwork, info)

    if (info /= 0) then
       stop 'ERROR: Matrix inversion failed!'
    end if

    deallocate(work)
  end subroutine invert_real

end module matrix_inverse
