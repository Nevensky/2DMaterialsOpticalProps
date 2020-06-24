! *==POINTR.spg  processed by SPAG 6.72Dc at 09:08 on  2 Jun 2020
! MODULE POINTR

use ISO_Fortran_env
INTEGER, PARAMETER :: sp = real32
INTEGER, PARAMETER :: dp = real64

IMPLICIT NONE
      SUBROUTINE POINTR(Root,Nsim,R,Ri)

        INTEGER :: i, nsim, is, n, m

        REAL(kind=sp), DIMENSION(48,3,3) :: R
        REAL(kind=sp), DIMENSION(48,3,3) :: RI

        REAL(kind=dp), DIMENSION(3,3) :: T
        REAL(kind=dp), DIMENSION(3,3) :: unit

        REAL(kind=sp), PARAMETER :: zero = 0.0
        REAL(kind=sp), PARAMETER :: one = 1.0

        REAL(kind=sp) :: x,y,z

        CHARACTER(len=11 ) :: buffer1,tag1
        CHARACTER(len=7  ) :: buffer2,tag2
        CHARACTER(len=100) :: root,fajl,path  
 
 
!       point group transformations R readed from 'MoS2.sc.out'
!       matrices R are in cart. coordinate system  because k
!       from IBZ are printed  in cart.coord.
 
 
      fajl = '/MoS2.sc.out'
      path = TRIM(Root)//TRIM(fajl)
      tag1 = '     atomic'
      tag2 = ' cryst.'
 
 
 
      OPEN (1,FILE=path)
 
 ! how many symmetries we have
      nsim_loop : DO i = 1 , 5000
         READ (1,'(A)') buffer1
         IF ( buffer1==tag1 ) THEN
            READ (1,'(X)')
            READ (1,'(X)')
            READ (1,*) Nsim
!            GOTO 100
            EXIT nsim_loop
         ENDIF
      END DO nsim_loop
 
      is = 0
      read_matrix_loop: DO i = 1 , 5000
         READ (1,'(a)') buffer2
         IF ( buffer2==tag2 ) THEN
            is = is + 1
            READ (1,'(X)')
            READ (1,'(X)')
            READ (1,'(X)')
            READ (1,'(19X,3F11.3)') x , y , z
            R(is,1,1) = x
            R(is,1,2) = y
            R(is,1,3) = z
            READ (1,'(19X,3F11.3)') x , y , z
            R(is,2,1) = x
            R(is,2,2) = y
            R(is,2,3) = z
            READ (1,'(19X,3F11.3)') x , y , z
            R(is,3,1) = x
            R(is,3,2) = y
            R(is,3,3) = z
            IF ( is==Nsim ) THEN
               EXIT read_matrix_loop
               !GOTO 200
            END IF
         END IF
      END DO read_matrix_loop
! 200  CLOSE (1)
      CLOSE (1)
 
 
!       INVERSION
      DO i = 1 , Nsim
         DO n = 1 , 3
            DO m = 1 , 3
               unit(n,m) = DCMPLX(zero,zero)
               t(n,m) = DCMPLX(R(i,n,m),zero)
            END DO
            unit(n,n) = DCMPLX(one,zero)
         END DO
         CALL GJEL(t,3,3,unit,3,3)
         DO n = 1 , 3
            DO m = 1 , 3
               Ri(i,n,m) = REAL(t(n,m))
            ENDDO
         ENDDO
      ENDDO
!99001 FORMAT (19X,3F11.3)
!99002 FORMAT (X)
 
 
 
      END
 
 ! END MODULE 
 
 
 
