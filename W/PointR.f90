! *==POINTR.spg  processed by SPAG 6.72Dc at 09:08 on  2 Jun 2020
MODULE ModPointR
USE ISO_Fortran_env

IMPLICIT NONE

INTEGER, PARAMETER :: sp = real32
INTEGER, PARAMETER :: dp = real64

CONTAINS
   SUBROUTINE PointR(root,Nsim,R,Ri)

      INTEGER :: i, nsim, is, n, m

      INTEGER :: lskip

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
 

      INTEGER :: ist3,ist4,ist5,ist6,ist7
      INTEGER :: lno3=0,lno4=0,lno5=0
 
!    point group transformations R readed from 'MoS2.sc.out'
!    matrices R are in cart. coordinate system  because k
!    from IBZ are printed  in cart.coord.
 
 
   fajl = '/MoS2.sc.out'
   path = TRIM(root)//TRIM(fajl)
   tag1 = '     atomic'
   tag2 = ' cryst.'
 
 
 
   OPEN (1,FILE=path,status='old',err=100,iostat=ist6)
 
!  how many symmetries we have
   nsim_loop : DO i = 1 , 5000
      READ (1,'(A)',err=1000,iostat=ist5,end=2000) buffer1
      lno4 = lno4+1
      IF ( buffer1==tag1 ) THEN
         READ (1,'(X)',err=101,iostat=ist7,end=201)
         READ (1,'(X)',err=101,iostat=ist7,end=202)
         READ (1,'(X)',err=101,iostat=ist7,end=203)
         READ (1,*) Nsim
         PRINT *,Nsim
!         GOTO 100
         EXIT nsim_loop
      ENDIF
   END DO nsim_loop
 
   is = 0
   read_matrix_loop: DO i = 1 , 5000
   lno3 = lno3+1
      READ (1,'(a)',err=1001,iostat=ist3,end=2001) buffer2
      IF ( buffer2==tag2 ) THEN
        is = is + 1
        DO lskip = 1,3
        READ (1,*) 
        END DO
        ! READ (1,'(X)',err=1002,iostat=ist4,end=2002)
        ! READ (1,'(X)',err=1002,iostat=ist4,end=2002)
        ! READ (1,'(X)',err=1002,iostat=ist4,end=2002)
        READ (1,'(19X,3F11.3)',err=1003,iostat=ist4,end=2003) x , y , z
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
        END IF
      END IF
   END DO read_matrix_loop
   CLOSE (1)

101   write(*,*)'buffer1=tag1. Error reading line ',lno4+1,', iostat = ',ist7
201   write(*,*)'buffer1=tag1. Number of lines read = ',lno4
102   write(*,*)'buffer1=tag1. Error reading line ',lno4+1,', iostat = ',ist7
202   write(*,*)'buffer1=tag1. Number of lines read = ',lno4
103   write(*,*)'buffer1=tag1. Error reading line ',lno4+1,', iostat = ',ist7
203   write(*,*)'buffer1=tag1. Number of lines read = ',lno4
104   write(*,*)'buffer1=tag1. Error reading line ',lno4+1,', iostat = ',ist7
204   write(*,*)'buffer1=tag1. Number of lines read = ',lno4

100   write(*,*)'cannot open file. iostat = ',ist6
1000   write(*,*)'buffer1 read. Error reading line ',lno4+1,', iostat = ',ist3
2000   write(*,*)'buffer1 read. Number of lines read = ',lno4
1001   write(*,*)'Error reading line ',lno3+1,', iostat = ',ist3
2001   write(*,*)'Number of lines read = ',lno3
1002   write(*,*)'buffer2=tag2. Error reading line ',lno3+1,', iostat = ',ist4
2002   write(*,*)'buffer2=tag2. Number of lines read = ',lno3 
1003   write(*,*)'buffer2=tag2. reding xyz Error reading line ',lno3+1,', iostat = ',ist4
2003   write(*,*)'buffer2=tag2. reding xyz Number of lines read = ',lno3 
!    INVERSION
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
 
 
 
   END SUBROUTINE PointR
 
END MODULE ModPointR

 
 
