!*==GJEL.spg  processed by SPAG 6.72Dc at 09:04 on  2 Jun 2020
      SUBROUTINE GJEL(A,N,Np,B,M,Mp)
!
!-----------------------------------------------------------------------
!
!	subroutine for solving a system of linear equation using
!	Gauss-Jordan elimination
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!*--GJEL12
!
      INTEGER N , Np , M , Mp , NMAX
      PARAMETER (NMAX=1000)
      DOUBLE COMPLEX A , B
      DIMENSION A(Np,Np) , B(Np,Mp)
      INTEGER i , icol , irow , j , k , l , ll , indxc , indxr , ipiv
      DOUBLE PRECISION big
      DOUBLE COMPLEX dum , pivinv , czero
!       parameter(czero=cmplx(0.0,0.0))
      DIMENSION indxc(NMAX) , indxr(NMAX) , ipiv(NMAX)
 
 
 
 
      czero = CMPLX(0.0,0.0)
      DO j = 1 , N
         ipiv(j) = 0
      ENDDO
!
      DO i = 1 , N
         big = 0.0
         DO j = 1 , N
            IF ( ipiv(j).NE.1 ) THEN
               DO k = 1 , N
                  IF ( ipiv(k).EQ.0 ) THEN
                     IF ( CDABS(A(j,k)).GE.big ) THEN
                        big = CDABS(A(j,k))
                        irow = j
                        icol = k
                     ENDIF
                  ELSEIF ( ipiv(k).GT.1 ) THEN
                     PRINT * , 'singular matrix in gjel'
                     STOP
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         ipiv(icol) = ipiv(icol) + 1
         IF ( irow.NE.icol ) THEN
            DO l = 1 , N
               dum = A(irow,l)
               A(irow,l) = A(icol,l)
               A(icol,l) = dum
            ENDDO
            DO l = 1 , M
               dum = B(irow,l)
               B(irow,l) = B(icol,l)
               B(icol,l) = dum
            ENDDO
         ENDIF
         indxr(i) = irow
         indxc(i) = icol
         IF ( A(icol,icol).EQ.czero ) PRINT * ,                         &
             &'singular matrix in gjel'
         pivinv = DCMPLX(1.0,0.0)/A(icol,icol)
         A(icol,icol) = DCMPLX(1.0,0.0)
         DO l = 1 , N
            A(icol,l) = A(icol,l)*pivinv
         ENDDO
         DO l = 1 , M
            B(icol,l) = B(icol,l)*pivinv
         ENDDO
         DO ll = 1 , N
            IF ( ll.NE.icol ) THEN
               dum = A(ll,icol)
               A(ll,icol) = czero
               DO l = 1 , N
                  A(ll,l) = A(ll,l) - A(icol,l)*dum
               ENDDO
               DO l = 1 , M
                  B(ll,l) = B(ll,l) - B(icol,l)*dum
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
      DO l = N , 1 , -1
         IF ( indxr(l).NE.indxc(l) ) THEN
            DO k = 1 , N
               dum = A(k,indxr(l))
               A(k,indxr(l)) = A(k,indxc(l))
               A(k,indxc(l)) = dum
            ENDDO
         ENDIF
      ENDDO
!
      END
