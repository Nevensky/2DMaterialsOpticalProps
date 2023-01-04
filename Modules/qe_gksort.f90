module qe_gksort
use iso_fortran_env, only: dp=>real64
implicit none

public  :: init_igk
private :: hpsort_eps, gk_sort

contains

subroutine init_igk( nks, xk, npwx, ngm, g, gcutw, igk_k)
  !! Subroutine adapted from Quantum Espresso 6.5 distribution. (PW/src/pwcom.f90)
  !! Initialize indices \(\text{igk_k}\) and number of plane waves
  !! per k-point:
  !
  !! * \( (k_{ik} + G)_i = k_{ik} + G_\text{igk} \);
  !! * i = 1, \text{ngk}(\text{ik});
  !! * \text{igk} = \text{igk}_k(i,ik).

  integer,  intent(in)  :: nks, npwx, ngm
  real(dp), intent(in)  :: gcutw !! wave function cut-off = ecutwfc/tpiba**2 where tpiba=2pi/alat
  real(dp), intent(in)  :: g(3,ngm), xk(:,:)
  integer,  intent(out), allocatable :: igk_k(:,:) ! neven debug: ja dodao

  integer, allocatable  :: ngk(:)
  real(dp), allocatable :: gk(:)
  integer :: ik,ig

  if (.not. allocated(igk_k)) allocate( igk_k(npwx,nks) )
  allocate(gk(npwx))
  allocate(ngk(nks))
  igk_k(:,:) = 0
  ! ... The following loop must NOT be called more than once in a run
  ! ... or else there will be problems with variable-cell calculations
  do ik = 1, nks
      ! neven debug: npwx sam ja dodao
      ! neven debug: promijenio sam igk_k(1,ik) u igk_k(1:npwx,ik)
      call gk_sort( xk(1:3,ik), npwx, ngm, g, gcutw, ngk(ik), igk_k(1:npwx,ik), gk )
  end do
  !
  deallocate( gk )
end subroutine init_igk

subroutine gk_sort( k, npwx, ngm, g, ecut, ngk, igk, gk )
   !! Subroutine adapted from Quantum Espresso 6.5 distribution. (PW/src/gk_sort.f90)
   !! Sorts k+g in order of increasing magnitude, up to ecut.

   real(dp), intent(in)  :: k(3)      !! k point
   integer,  intent(in)  :: ngm       !! No. of G vectors
   real(dp), intent(in)  :: g(3,ngm)  !! Coordinates of G vectors
   real(dp), intent(in)  :: ecut      !! the cut-off energy
   integer,  intent(in ) :: npwx      !! Max. No. of PW for wavefunctions
   integer,  intent(out) :: ngk       !! No. of k+G vectors inside the "ecut sphere"
   integer,  intent(out) :: igk(npwx) !! the correspondence k+G <-> G
   real(dp), intent(out) :: gk(npwx)  !! the moduli of k+G
   
   ! local variables
   integer :: ng   ! counter on   G vectors
   integer :: nk   ! counter on k+G vectors
   real(dp) :: q   ! |k+G|^2
   real(dp) :: q2x ! upper bound for |G|
   real(dp), parameter :: eps8  = 1.0e-8_dp

   ! first we count the number of k+G vectors inside the cut-off sphere
   q2x = ( sqrt( sum(k(:)**2) ) + sqrt( ecut ) )**2

   ngk = 0
   igk(:) = 0
   gk (:) = 0.0_dp

   do ng = 1, ngm
      q = sum( ( k(:) + g(:,ng) )**2 )
      if ( q <= eps8 ) q = 0.0_DP
      ! ... here if |k+G|^2 <= Ecut
      if ( q <= ecut ) THEN
         ngk = ngk + 1
         if ( ngk > npwx ) stop 'ERROR: Array gk out-of-bounds (gk_sort).'

         gk(ngk) = q
         igk(ngk) = ng ! set the initial value of index array
      else
         ! if |G| > |k| + SQRT( Ecut )  stop search and order vectors
         if ( sum( g(:,ng)**2 ) > ( q2x + eps8 ) ) EXIT
      end if
   end do

   if ( ng > ngm ) stop 'ERROR: Unexpected exit from do-loop (gk_sort).'
   ! ... order vector gk keeping initial position in index
   call hpsort_eps( ngk, gk, igk, eps8 )
   ! ... now order true |k+G|
   do nk = 1, ngk
      gk(nk) = sum( (k(:) + g(:,igk(nk)) )**2 )
   end do
end subroutine gk_sort

subroutine hpsort_eps (n, ra, ind, eps)
  !! Subroutine adapted from Quantum Espresso 6.5 distribution. (Modules/sort.f90)
  !! sort an array ra(1:n) into ascending order using heapsort algorithm,
  !! and considering two elements being equal if their values differ
  !! for less than "eps".
  !! n is input, ra is replaced on output by its sorted rearrangement.
  !! create an index table (ind) by making an exchange in the index array
  !! whenever an exchange is made on the sorted data array (ra).
  !! in case of equal values in the data array (ra) the values in the
  !! index array (ind) are used to order the entries.
  !! if on input ind(1)  = 0 then indices are initialized in the routine,
  !! if on input ind(1) != 0 then indices are assumed to have been
  !!                initialized before entering the routine and these
  !!                indices are carried around during the sorting process
  !!
  !! no work space needed !
  !! free us from machine-dependent sorting-routines !
  !!
  !! adapted from Numerical Recipes pg. 329 (new edition)
  !-input/output variables
  integer, intent(in) :: n  
  integer, intent(inout) :: ind (*)  
  real(dp), intent(inout) :: ra (*)
  real(dp), intent(in) :: eps
  !-local variables
  integer :: i, ir, j, l, iind  
  real(dp) :: rra  
  ! initialize index array
  if (ind (1) == 0) then  
     do i = 1, n  
        ind (i) = i  
     enddo
  endif
  ! nothing to order
  if (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  
  ir = n  
  sorting: do 
    ! still in hiring phase
    if ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    else  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       if ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       endif
    endif
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    do while ( j .le. ir )  
       if ( j .lt. ir ) then  
          ! compare to better underling
          if ( abs(ra(j)-ra(j+1)).ge.eps ) then  
             if (ra(j).lt.ra(j+1)) j = j + 1
          else
             ! this means ra(j) == ra(j+1) within tolerance
             if (ind (j) .lt.ind (j + 1) ) j = j + 1
          endif
       endif
       ! demote rra
       if ( abs(rra - ra(j)).ge.eps ) then  
          if (rra.lt.ra(j)) then
             ra (i) = ra (j)  
             ind (i) = ind (j)  
             i = j  
             j = j + j  
          else
             ! set j to terminate do-while loop
             j = ir + 1  
          end if
       else
          !this means rra == ra(j) within tolerance
          ! demote rra
          if (iind.lt.ind (j) ) then
             ra (i) = ra (j)
             ind (i) = ind (j)
             i = j
             j = j + j
          else
             ! set j to terminate do-while loop
             j = ir + 1
          endif
       end if
    enddo
    ra (i) = rra  
    ind (i) = iind  

  end do sorting    
end subroutine hpsort_eps

end module qe_gksort
