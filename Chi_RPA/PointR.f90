! *==POINTR.spg  processed by SPAG 6.72Dc at 09:08 on  2 Jun 2020
MODULE ModPointR
USE ISO_Fortran_env

IMPLICIT NONE

INTEGER, PARAMETER :: sp = real32
INTEGER, PARAMETER :: dp = real64

CONTAINS
   SUBROUTINE PointR(path,Nsymm,R,Ri)

      INTEGER :: i, Nsymm, is, n, m

      INTEGER :: lskip

      REAL(kind=dp), DIMENSION(48,3,3) :: R
      REAL(kind=dp), DIMENSION(48,3,3) :: RI

      COMPLEX(kind=dp), DIMENSION(3,3) :: T
      COMPLEX(kind=dp), DIMENSION(3,3) :: Imat

      REAL(kind=dp), PARAMETER :: zero = 0.0
      REAL(kind=dp), PARAMETER :: one = 1.0

      REAL(kind=dp) :: x,y,z

      CHARACTER(len=11 ) :: buffer1,tag1
      CHARACTER(len=7  ) :: buffer2,tag2
      CHARACTER(len=100) :: path  
 

      INTEGER :: ist3,ist4,ist5,ist6,ist7
      INTEGER :: lno3=0,lno4=0,lno5=0
 
!    point group transformations R readed from 'MoS2.sc.out'
!    matrices R are in cart. coordinate system  because k
!    from IBZ are printed  in cart.coord.
 
 
   tag1 = '     atomic'
   tag2 = ' cryst.'
 
 
 
   OPEN (1,FILE=path,status='old',err=100,iostat=ist6)
 
!  how many symmetries we have
   Nsymm_loop : DO i = 1 , 5000
      READ (1,'(A)',err=1000,iostat=ist5,end=2000) buffer1
      lno4 = lno4+1
      IF ( buffer1==tag1 ) THEN
         DO lskip = 1,3
            READ (1,*,err=101,iostat=ist7,end=201) 
         END DO
         READ (1,*) Nsymm
         PRINT *,"Nsymm = ",Nsymm
         EXIT Nsymm_loop
      ENDIF
   END DO Nsymm_loop
   PRINT *, 'Exited Nsymm_loop'

   is = 0
   load_R_loop: DO i = 1 , 5000
   ! lno3 = lno3+1
   ! PRINT *,'line =',lno3
      READ (1,'(a)',err=1001,iostat=ist3,end=2001) buffer2
      IF ( buffer2==tag2 ) THEN
        is = is + 1
        DO lskip = 1,3
         READ (1,*,err=1002,iostat=ist4,end=2002) 
        END DO
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
        ! PRINT *,'symm. op. number is=', is 
        ! WRITE(*,'(3F11.3,/,3F11.3,/,3F11.3)') R(is,:,:)
        IF ( is==Nsymm ) THEN
            EXIT load_R_loop
        END IF
      END IF
   END DO load_R_loop
   CLOSE (1)
   ! PRINT *, 'EXITED R matrix loading loop.'

GOTO 2004
101   write(*,*)'101 buffer1=tag1. Error reading line ',lno4+1,', iostat = ',ist7
201   write(*,*)'201 buffer1=tag1. Number of lines read = ',lno4
102   write(*,*)'102 buffer1=tag1. Error reading line ',lno4+1,', iostat = ',ist7
202   write(*,*)'202 buffer1=tag1. Number of lines read = ',lno4
103   write(*,*)'103 buffer1=tag1. Error reading line ',lno4+1,', iostat = ',ist7
203   write(*,*)'203 buffer1=tag1. Number of lines read = ',lno4
104   write(*,*)'104 buffer1=tag1. Error reading line ',lno4+1,', iostat = ',ist7
204   write(*,*)'204 buffer1=tag1. Number of lines read = ',lno4

100    write(*,*)'100 cannot open file. iostat = ',ist6
1000   write(*,*)'1000 buffer1 read. Error reading line ',lno4+1,', iostat = ',ist3
2000   write(*,*)'200 buffer1 read. Number of lines read = ',lno4
1001   write(*,*)'1001 Error reading line ',lno3+1,', iostat = ',ist3
2001   write(*,*)'2001 Number of lines read = ',lno3
1002   write(*,*)'1002 buffer2=tag2. Error reading line ',lno3+1,', iostat = ',ist4
2002   write(*,*)'2002 buffer2=tag2. Number of lines read = ',lno3 
1003   write(*,*)'1003 buffer2=tag2. reding xyz Error reading line ',lno3+1,', iostat = ',ist4
2003   write(*,*)'2003 buffer2=tag2. reding xyz Number of lines read = ',lno3 
2004 CONTINUE

!    INVERSION
   DO i = 1 , Nsymm
      DO n = 1 , 3
         DO m = 1 , 3
            Imat(n,m) = DCMPLX(zero,zero)
            T(n,m) = DCMPLX(R(i,n,m),zero)
         END DO
         Imat(n,n) = DCMPLX(one,zero)
      END DO
      ! PRINT *,'EXITED n loop'
      PRINT *,'symmetry op. matrix R'
      WRITE(*,'(A8,I3/3F11.4/3F11.4/3F11.4)'),'i_symm= ',i, R(i,:,:)
      ! PRINT *,'preparing to invert helper matrix T'
      ! WRITE(*,'(A8,I3/6F11.4/6F11.4/6F11.4)'),'i_symm= ',i, T
      CALL GJEL(T,3,3,Imat,3,3)
      PRINT *,'COMPLETED matrix inversion'
      DO n = 1 , 3
         DO m = 1 , 3
            Ri(i,n,m) = REAL(T(n,m))
         ENDDO
      ENDDO
      WRITE(*,'(A8,I3/3F11.4/3F11.4/3F11.4)'),'i_symm= ',i, Ri(i,:,:)
   ENDDO
 
 
 
   END SUBROUTINE PointR
 
END MODULE ModPointR

 
 
