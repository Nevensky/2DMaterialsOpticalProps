module ModPointR
use iso_fortran_env, only: dp => real64
use matrix_inverse,  only: invert_real

implicit none

public :: PointR
private

contains
   subroutine PointR(path, Nsym, R, Ri)
      !  Point group transformations R read from 
      !  Quantum Espresso SCF output file.
      !  Matrices R are in cartesian coordinates because k points
      !  from IBZ are printed in cart. coord.
      integer :: i, Nsym, is, n, m

      real(dp) :: R(48,3,3)
      real(dp) :: Ri(48,3,3)

      real(dp) :: x,y,z

      character(len=11 ) :: buffer1, tag1
      character(len=7  ) :: buffer2, tag2
      character(len=100) :: path  
 

      ! file reading debug vars
      integer :: iuni
      integer :: ist3, ist4, ist5, ist6, ist7
      integer :: lno3 = 0, lno4 = 0
      integer :: lskip ! for skipping lines 
 
      tag1 = '     atomic'
      tag2 = ' cryst.'
 
      open(newunit=iuni, file=path, status='old', err=100, iostat=ist6)

      ! count how many symmetries we have
      Nsym_loop : do i = 1 , 5000
         read (1,'(A)', err=1000, iostat=ist5, end=2000) buffer1
         lno4 = lno4+1
         if ( buffer1==tag1 ) THEN
            do lskip = 1,3
               read (1,*,err=101,iostat=ist7,end=201) 
            end do
            read (1,*) Nsym
            print *,"Nsym = ",Nsym
            EXIT Nsym_loop
         endif
      end do Nsym_loop
      print *, 'Exited Nsym_loop'

      is = 0
      load_R_loop: do i = 1 , 5000
      lno3 = lno3+1
         read (1,'(a)',err=1001,iostat=ist3,end=2001) buffer2
         if ( buffer2==tag2 ) THEN
           is = is + 1
           do lskip = 1,3
            read (1,*,err=1002,iostat=ist4,end=2002) 
           end do
           read (1,'(19X,3F11.3)',err=1003,iostat=ist4,end=2003) x , y , z
           R(is,1,1) = x
           R(is,1,2) = y
           R(is,1,3) = z
           read (1,'(19X,3F11.3)') x , y , z
           R(is,2,1) = x
           R(is,2,2) = y
           R(is,2,3) = z
           read (1,'(19X,3F11.3)') x , y , z
           R(is,3,1) = x
           R(is,3,2) = y
           R(is,3,3) = z
           ! PRINT *,'symm. op. number is=', is 
           ! WRITE(*,'(3F11.3,/,3F11.3,/,3F11.3)') R(is,:,:)
           if ( is==Nsym ) THEN
               EXIT load_R_loop
           end if
         end if
      end do load_R_loop
      close(iuni)
      ! PRINT *, 'EXITED R matrix reading loop.'

      ! Matrix inversion

      Ri = R ! intialize inverse matrix with same content as R
      do is = 1 , Nsym
         print *,'symmetry op. matrix R'
         write(*,'(A8,I3/3F11.4/3F11.4/3F11.4)'),'i_symm= ',is, R(is,:,:)
         call invert_real(Ri(is,:,:))
         print *,'inverted symmetry op. matrix Ri'
         write(*,'(A8,I3/3F11.4/3F11.4/3F11.4)'),'i_symm= ',is, Ri(is,:,:)
      enddo
 
 GOTO 2004
101   write(*,*)'101 Error reading line ',lno4+1,', iostat = ',ist7
201   write(*,*)'201 Number of lines read = ',lno4
102   write(*,*)'102 Error reading line ',lno4+1,', iostat = ',ist7
202   write(*,*)'202 Number of lines read = ',lno4
103   write(*,*)'103 Error reading line ',lno4+1,', iostat = ',ist7
203   write(*,*)'203 Number of lines read = ',lno4
104   write(*,*)'104 Error reading line ',lno4+1,', iostat = ',ist7
204   write(*,*)'204 Number of lines read = ',lno4

100    write(*,*)'100 cannot open file. iostat = ',ist6
1000   write(*,*)'1000 buffer1 read. Error reading line ',lno4+1,', iostat = ',ist5
2000   write(*,*)'200 buffer1 read. Number of lines read = ',lno4
1001   write(*,*)'1001 Error reading line ',lno3+1,', iostat = ',ist3
2001   write(*,*)'2001 Number of lines read = ',lno3
1002   write(*,*)'1002 Error reading line ',lno3+1,', iostat = ',ist4
2002   write(*,*)'2002 Number of lines read = ',lno3 
1003   write(*,*)'1003 reding xyz Error reading line ',lno3+1,', iostat = ',ist4
2003   write(*,*)'2003 reding xyz Number of lines read = ',lno3 
2004 CONTINUE
 
 
      end subroutine PointR
 
end module ModPointR

 
 
