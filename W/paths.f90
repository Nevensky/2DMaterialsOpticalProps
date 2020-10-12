      ! subroutine paths(savedir,k1,k2,n,m,pathk1,pathk2,bandn,bandm)
      ! implicit none
 
      ! integer k1 , k2 , n , m , nord
      ! character(len=100) savedir , folder , file , pathk1 , pathk2 ,  bandn , bandm
 
 
 
      ! file = '/evc.dat'
 
      ! folder = '/K0000i'
      ! nord = index(folder,'i',back=.false.)
      ! if ( k1<10 ) then
      !    write (folder(nord:nord),'(i1)') k1
      ! elseif ( k1>=10 .and. k1<100 ) then
      !    write (folder(nord-1:nord),'(i2)') k1
      ! elseif ( k1>=100 .and. k1<1000 ) then
      !    write (folder(nord-2:nord),'(i3)') k1
      ! elseif ( k1>=1000 .and. k1<10000 ) then
      !    write (folder(nord-3:nord),'(i4)') k1
      ! else
      !    write (folder(nord-4:nord),'(i5)') k1
      ! endif
      ! pathk1 = trim(savedir)//trim(folder)//trim(file)
 
      ! folder = '/K0000i'
      ! nord = index(folder,'i',back=.false.)
      ! if ( k2<10 ) then
      !    write (folder(nord:nord),'(i1)') k2
      ! elseif ( k2>=10 .and. k2<100 ) then
      !    write (folder(nord-1:nord),'(i2)') k2
      ! elseif ( k2>=100 .and. k2<1000 ) then
      !    write (folder(nord-2:nord),'(i3)') k2
      ! elseif ( k2>=1000 .and. k2<10000 ) then
      !    write (folder(nord-3:nord),'(i4)') k2
      ! else
      !    write (folder(nord-4:nord),'(i5)') k2
      ! endif
      ! pathk2 = trim(savedir)//trim(folder)//trim(file)
 
 
 
      ! bandn = 'evc.n'
      ! nord = index(bandn,'n',back=.false.)
      ! if ( n<10 ) then
      !    write (bandn(nord:nord),'(i1)') n
      ! elseif ( n>=10 .and. n<100 ) then
      !    write (bandn(nord:nord+1),'(i2)') n
      ! else
      !    write (bandn(nord:nord+2),'(i3)') n
      ! endif
      ! bandm = 'evc.n'
      ! nord = index(bandm,'n',back=.false.)
      ! if ( m<10 ) then
      !    write (bandm(nord:nord),'(i1)') m
      ! elseif ( m>=10 .and. m<100 ) then
      !    write (bandm(nord:nord+1),'(i2)') m
      ! else
      !    write (bandm(nord:nord+2),'(i3)') m
      ! endif
 
 
 
 
      ! end subroutine paths
 
 
subroutine paths(savedir,k,n,pathk,bandn)
   implicit none
   
   integer k , n, nord
   character(len=100) savedir , folder , file , pathk ,  bandn
   
   file = '/evc.dat'
   
   folder = '/K0000i'
   nord = index(folder,'i',back=.false.)
   if ( k<10 ) then
      write (folder(nord:nord),'(i1)') k
   elseif ( k>=10 .and. k<100 ) then
      write (folder(nord-1:nord),'(i2)') k
   elseif ( k>=100 .and. k<1000 ) then
      write (folder(nord-2:nord),'(i3)') k
   elseif ( k>=1000 .and. k<10000 ) then
      write (folder(nord-3:nord),'(i4)') k
   else
      write (folder(nord-4:nord),'(i5)') k
   endif
   pathk = trim(savedir)//trim(folder)//trim(file)
   
   bandn = 'evc.n'
   nord = index(bandn,'n',back=.false.)
   if ( n<10 ) then
      write (bandn(nord:nord),'(i1)') n
   elseif ( n>=10 .and. n<100 ) then
      write (bandn(nord:nord+1),'(i2)') n
   else
      write (bandn(nord:nord+2),'(i3)') n
   endif
end subroutine paths

subroutine pathsB(savedir,k,n,pathk,bandn)
   implicit none
   
   integer k , n, nord
   character(len=100) savedir , folder , file , pathk ,  bandn
   
   file = '/evc.dat'
   
   folder = '/K0000i'
   nord = index(folder,'i',back=.false.)
   if ( k<10 ) then
      write (folder(nord:nord),'(i1)') k
   elseif ( k>=10 .and. k<100 ) then
      write (folder(nord-1:nord),'(i2)') k
   elseif ( k>=100 .and. k<1000 ) then
      write (folder(nord-2:nord),'(i3)') k
   elseif ( k>=1000 .and. k<10000 ) then
      write (folder(nord-3:nord),'(i4)') k
   else
      write (folder(nord-4:nord),'(i5)') k
   endif
   pathk = trim(savedir)//trim(folder)//trim(file)
   
   bandn = 'evc.n'
   nord = index(bandn,'n',back=.false.)
   if ( n<10 ) then
      write (bandn(nord:nord),'(i1)') n
   elseif ( n>=10 .and. n<100 ) then
      write (bandn(nord:nord+1),'(i2)') n
   else
      write (bandn(nord:nord+2),'(i3)') n
   endif
end subroutine pathsB