SUBROUTINE pointr(root,nsim,r,ri)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2020-06-17  Time: 12:59:12

implicitnone


CHARACTER (LEN=100), INTENT(IN OUT)      :: root
INTEGER, INTENT(IN)                      :: nsim
REAL*8, INTENT(OUT)                      :: r(48,3,3)
REAL*8, INTENT(OUT)                      :: ri(48,3,3)
INTEGER :: i, is,n,m
REAL*8  x,y,z
DOUBLE COMPLEX t,UNIT
REAL*8, PARAMETER :: zero=0.0
REAL*8, PARAMETER :: one=1.0
DIMENSION  t(3,3), UNIT(3,3)
CHARACTER (LEN=11,buffer1,tag1) ::
CHARACTER (LEN=7,buffer2,tag2) ::
CHARACTER (LEN=100) :: fajl,path


!       point group transformations R readed from 'MoS2.sc.out'
!       matrices R are in cart. coordinate system  because k
!       from IBZ are printed  in cart.coord.


fajl='/MoS2.sc.out'
path=trim(root)//trim(fajl)
tag1='     atomic'
tag2=' cryst.'



OPEN(1,FILE=path)

DO  i=1,5000
  READ(1,'(a)')buffer1
  IF(buffer1 == tag1)THEN
    READ(1,80)
    READ(1,80)
    READ(1,*)nsim
    EXIT
  END IF
END DO
999     CONTINUE

is=0
DO  i=1,5000
  READ(1,'(a)')buffer2
  IF(buffer2 == tag2)THEN
    is=is+1
    READ(1,80)
    READ(1,80)
    READ(1,80)
    READ(1,70)x,y,z
    r(is,1,1)=x
    r(is,1,2)=y
    r(is,1,3)=z
    READ(1,70)x,y,z
    r(is,2,1)=x
    r(is,2,2)=y
    r(is,2,3)=z
    READ(1,70)x,y,z
    r(is,3,1)=x
    r(is,3,2)=y
    r(is,3,3)=z
    IF(is == nsim)EXIT
  END IF
END DO
998     CONTINUE
CLOSE(1)
70      FORMAT(19X,3F11.3)
80      FORMAT(x)


!       INVERTION
DO  i=1,nsim
  DO  n=1,3
    DO  m=1,3
      UNIT(n,m)=DCMPLX(zero,zero)
      t(n,m)=DCMPLX(r(i,n,m),zero)
    END DO
    UNIT(n,n)=DCMPLX(one,zero)
  END DO
  CALL gjel(t,3,3,UNIT,3,3)
  DO  n=1,3
    DO  m=1,3
      ri(i,n,m)=REAL(t(n,m))
    END DO
  END DO
END DO



RETURN
END SUBROUTINE pointr





