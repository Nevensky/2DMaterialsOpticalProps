	subroutine gjel(a,n,np,b,m,mp)
c
c-----------------------------------------------------------------------
c
c	subroutine for solving a system of linear equation using
c	Gauss-Jordan elimination
c
c-----------------------------------------------------------------------
c
	implicitnone
c
	integer n,np,m,mp,nmax
	parameter(nmax=1000)
	double complex a,b
	dimension a(np,np),b(np,mp)
	integer i,icol,irow,j,k,l,ll,indxc,indxr,ipiv
	double precision big
	double complex dum,pivinv,czero
c       parameter(czero=cmplx(0.0,0.0))
	dimension indxc(nmax),indxr(nmax),ipiv(nmax)

                  


        czero=cmplx(0.0,0.0)
	do 11 j=1,n
  11    ipiv(j)=0
c
  	do 22 i=1,n
  	   big=0.0
  	   do 13 j=1,n
  	      if(ipiv(j).ne.1) then
  	         do 12 k=1,n
  	            if(ipiv(k).eq.0) then
  	               if(cdabs(a(j,k)).ge.big) then
  	                  big=cdabs(a(j,k))
  	                  irow=j
  	                  icol=k
  	               endif
  	            else if(ipiv(k).gt.1) then
  	               print*,'singular matrix in gjel'
                       stop
  	            endif
  12           continue
            endif
  13     continue
  	   ipiv(icol)=ipiv(icol)+1
  	   if(irow.ne.icol) then
  	      do 14 l=1,n
  	         dum=a(irow,l)
  	         a(irow,l)=a(icol,l)
  14	      a(icol,l)=dum
		do 15 l=1,m
		   dum=b(irow,l)
		   b(irow,l)=b(icol,l)
  15        b(icol,l)=dum
  	   endif
  	   indxr(i)=irow
  	   indxc(i)=icol
  	   if(a(icol,icol).eq.czero)then
           print*,'singular matrix in gjel'
           endif
      	   pivinv=dcmplx(1.0,0.0)/a(icol,icol)
  	   a(icol,icol)=dcmplx(1.0,0.0)
  	   do 16 l=1,n
  16     a(icol,l)=a(icol,l)*pivinv
  	   do 17 l=1,m
  17     b(icol,l)=b(icol,l)*pivinv
  	   do 21 ll=1,n
  	      if(ll.ne.icol) then
  	         dum=a(ll,icol)
  	         a(ll,icol)=czero
  	         do 18 l=1,n
  18		   a(ll,l)=a(ll,l)-a(icol,l)*dum
  	         do 19 l=1,m
  19		   b(ll,l)=b(ll,l)-b(icol,l)*dum
  		endif
  21     continue
  22  continue
c
	do 24 l=n,1,-1
	   if(indxr(l).ne.indxc(l)) then
	      do 23 k=1,n
	         dum=a(k,indxr(l))
	         a(k,indxr(l))=a(k,indxc(l))
  23        a(k,indxc(l))=dum
  	   endif
  24  continue
c
	return
	end
