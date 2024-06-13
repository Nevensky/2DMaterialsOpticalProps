module matrix_inverse
!! Returns the inverse of a matrix calculated by finding the
!! LU decomposition. Depends on LAPACK, preferably Intel MKL.

   use lapack95, only: getrf, getri
   use iso_fortran_env, only: dp => real64, idp => int16
   use notifications, only: error

   implicit none

   public :: invert, checkIdentity !, invert_cmplx, invert_real
   private

   interface invert
   !! Inverts real or complex matrix (double precision)
      procedure invert_cmplx, invert_real
   end interface invert

   interface checkIdentity
   !! Checks if real or complex matrix multiplied by its inverse 
   !! gives an identity matrix (double precision)
      procedure check_identity_cmplx, check_identity_real
   end interface checkIdentity

contains
  subroutine invert_cmplx(A)
    !! Inverts complex matrix (double precision)
    complex(dp), dimension(:,:), intent(inout) :: A

    complex(dp), allocatable :: work(:)         !! work array for LAPACK
    integer                  :: ipiv(size(A,1)) !! pivot indices
    integer                  :: n, info, ios 
    integer(idp)             :: lwork           !! work array size
    ! integer(idp), parameter :: memmax = 60*sizeof(dp) 

    ! External procedures defined in LAPACK for double complex matrices
    external zgetrf
    external zgetri

    n = size(A,1) ! number of rows of A

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call zgetrf(n, n, A, n, ipiv, info)

    if (info /= 0) then
       call error('Matrix is numerically singular! (zgetrf)')
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
      print *, 'ios=',ios
      call error('Failed allocation of work array. (zgetri)')
    end if
    ! print *, " sizeof(A) = ",sizeof(A)/1024.0**3,"GB"
    ! print *, " sizeof(work) = ",sizeof(work)/1024.0**3,"GB"

    ! Perform inverse with the optimal work array
    call zgetri(n, A, n, ipiv, work, lwork, info)

    if (info /= 0) then
       call error('Matrix inversion failed! (zgetri)')
    end if

    deallocate(work)
  end subroutine invert_cmplx

  subroutine invert_real(A)
    !! Inverts real matrix (double precision)
    real(dp), dimension(:,:), intent(inout) :: A

    real(dp), allocatable :: work(:)          !! work array for LAPACK
    integer               :: ipiv(size(A,1))  !! pivot indices
    integer               :: n, info, info2, info3, ios 
    integer(idp)          :: lwork            !! work array size
    ! integer(idp), parameter :: memmax = 30*sizeof(dp) 

    ! External procedures defined in LAPACK for double real matrices
    external dgetrf
    external dgetri

    n = size(A,1) ! number of rows of A

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(n, n, A, n, ipiv, info)

    if (info /= 0) then
       print *, 'info=',info
       call error('Matrix is numerically singular! (dgetrf)')
    end if

    ! DGETRI computes the inverse of a matrix using 
    ! the LU factorization computed by DGETRF.

    ! Query optimal work array for best performance
    allocate(work(1))
    call dgetri(n, A, n, ipiv, work, -1, info2) ! query optimal lwork value
    if (info2 /= 0) then
      print *, 'info2=', info2
      call error('Workspace query failed! (dgetri)')
    end if

    lwork = work(1)
    deallocate(work)
    !if (sizeof(lwork)>lworkmax) lwork = max_available_memory
    allocate(work(lwork),stat=ios)
    if (ios /= 0) then
            print *,'ios=',ios
            call error('Failed allocation of work array. (dgetri)')
    end if

    ! print *, " sizeof(A) = ",sizeof(A)/1024.0**3,"GB"
    ! print *, " sizeof(work) = ",sizeof(work)/1024.0**3,"GB"

    ! Perform inverse with the optimal work array
    call dgetri(n, A, n, ipiv, work, lwork, info3)

    if (info3 /= 0) then
       print *, 'info3=', info3
       call error('Matrix inversion failed! (dgetri)')
    end if

    deallocate(work)
  end subroutine invert_real

subroutine check_identity_cmplx(A,A_inverted,threshold_)
   !! Checks if a complex matrix was correctly inverted i.e. if A⁻¹ · A = 𝟙
   !! WARNING: should be double checked see:
   !! https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-0/gemm-001.html
   implicit none
   complex(kind=dp), intent(in)           :: A(:,:)
   complex(kind=dp), intent(in)           :: A_inverted(:,:)
   real(kind=dp),    intent(in), optional :: threshold_

   complex(kind=dp), allocatable :: Icheck(:,:)

   integer :: Nlf
   ! dgemm vars (inversion success check)
   integer          :: i, K
   real(kind=dp)    :: alpha = 1.0_dp
   real(kind=dp)    :: beta  = 0.0_dp
   complex(kind=dp) :: checkIdentity, checkIdentity2
   complex(kind=dp) :: traceIdentity
   real(kind=dp)    :: threshold = dble(10d-4) ! default value arbitrarily chosen

   if (present(threshold_)) threshold = threshold_

   Nlf = size(A,1)
   K = Nlf   
   allocate(Icheck(Nlf,Nlf))

   call zgemm('N','N', Nlf, Nlf, K, alpha, A_inverted, Nlf, A, K, beta, Icheck, Nlf)
   traceIdentity = sum( (/ ( abs(Icheck(i,i)), i=1, size(Icheck, 1)) /) )
   checkIdentity = sum(abs(Icheck)) - traceIdentity
   checkIdentity2 = dcmplx(Nlf,0.0_dp) - traceIdentity
   if ( dble(checkIdentity)>threshold .or. &
      & aimag(checkIdentity)>threshold .or. &
      & dble(checkIdentity2)>threshold .or. &
      & aimag(checkIdentity2)>threshold) then
     print *, 'Nlf: ', Nlf
     print *, 'size Ainv: ', size(A_inverted), 'size A:',size(A)
     print *, 'Re(check): ', dble(checkIdentity)
     print *, 'Im(check): ', aimag(checkIdentity)
     print *, 'Re(check2): ', dble(checkIdentity2)
     print *, 'Im(check2): ', aimag(checkIdentity2)
     call error('Matrix inversion failed. (A⁻¹ · A ≠ 𝟙)')
   endif

   deallocate(Icheck)
end subroutine check_identity_cmplx

subroutine check_identity_real(A,A_inverted,threshold_)
   !! Checks if a real matrix was correctly inverted i.e. if A⁻¹ · A = 𝟙
   !! WARNING: should be double checked see:
   !! https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-0/gemm-001.html
   implicit none
   real(kind=dp), intent(in)           :: A(:,:)
   real(kind=dp), intent(in)           :: A_inverted(:,:)
   real(kind=dp), intent(in), optional :: threshold_

   real(kind=dp), allocatable :: Icheck(:,:)

   integer :: Nlf
   ! dgemm vars (inversion success check)
   integer          :: i, K
   real(kind=dp)    :: alpha = 1.0_dp
   real(kind=dp)    :: beta  = 0.0_dp
   real(kind=dp)    :: checkIdentity, checkIdentity2
   real(kind=dp)    :: traceIdentity
   real(kind=dp)    :: threshold = dble(10d-4) ! default value arbitrarily chosen

   if (present(threshold_)) threshold = threshold_

   Nlf = size(A,1)
   K = Nlf   
   allocate(Icheck(Nlf,Nlf))

   call dgemm('N','N', Nlf, Nlf, K, alpha, A_inverted, Nlf, A, K, beta, Icheck, Nlf)
   traceIdentity = sum( (/ ( abs(Icheck(i,i)), i=1, size(Icheck, 1)) /) )
   checkIdentity = sum(abs(Icheck)) - traceIdentity
   checkIdentity2 = Nlf - traceIdentity
   if ( dble(checkIdentity)>threshold .or. &
      & dble(checkIdentity2)>threshold) then
     print *, 'Nlf: ', Nlf
     print *, 'size Ainv: ', size(A_inverted), 'size A:',size(A)
     print *, 'check: ', dble(checkIdentity)
     print *, 'check2: ', dble(checkIdentity2)
     call error('Matrix inversion failed. (A⁻¹ · A ≠ 𝟙)')
   endif

   deallocate(Icheck)
end subroutine check_identity_real


end module matrix_inverse
