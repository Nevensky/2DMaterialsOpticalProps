!*==POINTR.spg  processed by SPAG 6.72Dc at 09:08 on  2 Jun 2020
      SUBROUTINE POINTR(Root,Nsim,R,Ri)
 
      IMPLICIT NONE
!*--POINTR5
 
      INTEGER i , Nsim , is , n , m
      REAL*8 R , Ri , ZERO , ONE , x , y , z
      DOUBLE COMPLEX t , unit
      PARAMETER (ZERO=0.0,ONE=1.0)
      DIMENSION R(48,3,3) , t(3,3) , Ri(48,3,3) , unit(3,3)
      CHARACTER*11 , BUFFER1 , TAG1*0 
      CHARACTER*7 , BUFFER2 , TAG2*0 
      CHARACTER*100 Root , fajl , path
 
 
!       point group transformations R readed from 'MoS2.sc.out'
!       matrices R are in cart. coordinate system  because k
!       from IBZ are printed  in cart.coord.
 
 
      fajl = '/MoS2.sc.out'
      path = TRIM(Root)//TRIM(fajl)
      tag1 = '     atomic'
      tag2 = ' cryst.'
 
 
 
      OPEN (1,FILE=path)
 
      DO i = 1 , 5000
         READ (1,'(a)') buffer1
         IF ( buffer1==tag1 ) THEN
            READ (1,99002)
            READ (1,99002)
            READ (1,*) Nsim
            GOTO 100
         ENDIF
      ENDDO
 
 100  is = 0
      DO i = 1 , 5000
         READ (1,'(a)') buffer2
         IF ( buffer2==tag2 ) THEN
            is = is + 1
            READ (1,99002)
            READ (1,99002)
            READ (1,99002)
            READ (1,99001) x , y , z
            R(is,1,1) = x
            R(is,1,2) = y
            R(is,1,3) = z
            READ (1,99001) x , y , z
            R(is,2,1) = x
            R(is,2,2) = y
            R(is,2,3) = z
            READ (1,99001) x , y , z
            R(is,3,1) = x
            R(is,3,2) = y
            R(is,3,3) = z
            IF ( is==Nsim ) GOTO 200
         ENDIF
      ENDDO
 200  CLOSE (1)
 
 
!       INVERTION
      DO i = 1 , Nsim
         DO n = 1 , 3
            DO m = 1 , 3
               unit(n,m) = DCMPLX(ZERO,ZERO)
               t(n,m) = DCMPLX(R(i,n,m),ZERO)
            ENDDO
            unit(n,n) = DCMPLX(ONE,ZERO)
         ENDDO
         CALL GJEL(t,3,3,unit,3,3)
         DO n = 1 , 3
            DO m = 1 , 3
               Ri(i,n,m) = REAL(t(n,m))
            ENDDO
         ENDDO
      ENDDO
99001 FORMAT (19X,3F11.3)
99002 FORMAT (X)
 
 
 
      END
 
 
 
 
