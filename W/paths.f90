!*==PATHS.spg  processed by SPAG 6.72Dc at 09:06 on  2 Jun 2020
      SUBROUTINE PATHS(Root1,K1,K2,N,M,Pathk1,Pathk2,Bandn,Bandm)
      IMPLICIT NONE
!*--PATHS4
 
      INTEGER K1 , K2 , N , M , nord
      CHARACTER(len=100) Root1 , root2 , folder , file , Pathk1 , Pathk2 ,  Bandn , Bandm
 
 
 
      root2 = '/MoS2.save'
      file = '/evc.dat'
 
      folder = '/K0000i'
      nord = INDEX(folder,'i',BACK=.FALSE.)
      IF ( K1<10 ) THEN
         WRITE (folder(nord:nord),'(i1)') K1
      ELSEIF ( K1>=10 .AND. K1<100 ) THEN
         WRITE (folder(nord-1:nord),'(i2)') K1
      ELSEIF ( K1>=100 .AND. K1<1000 ) THEN
         WRITE (folder(nord-2:nord),'(i3)') K1
      ELSEIF ( K1>=1000 .AND. K1<10000 ) THEN
         WRITE (folder(nord-3:nord),'(i4)') K1
      ELSE
         WRITE (folder(nord-4:nord),'(i5)') K1
      ENDIF
      Pathk1 = TRIM(Root1)//TRIM(root2)//TRIM(folder)//TRIM(file)
 
      folder = '/K0000i'
      nord = INDEX(folder,'i',BACK=.FALSE.)
      IF ( K2<10 ) THEN
         WRITE (folder(nord:nord),'(i1)') K2
      ELSEIF ( K2>=10 .AND. K2<100 ) THEN
         WRITE (folder(nord-1:nord),'(i2)') K2
      ELSEIF ( K2>=100 .AND. K2<1000 ) THEN
         WRITE (folder(nord-2:nord),'(i3)') K2
      ELSEIF ( K2>=1000 .AND. K2<10000 ) THEN
         WRITE (folder(nord-3:nord),'(i4)') K2
      ELSE
         WRITE (folder(nord-4:nord),'(i5)') K2
      ENDIF
      Pathk2 = TRIM(Root1)//TRIM(root2)//TRIM(folder)//TRIM(file)
 
 
 
      Bandn = 'evc.n'
      nord = INDEX(Bandn,'n',BACK=.FALSE.)
      IF ( N<10 ) THEN
         WRITE (Bandn(nord:nord),'(i1)') N
      ELSEIF ( N>=10 .AND. N<100 ) THEN
         WRITE (Bandn(nord:nord+1),'(i2)') N
      ELSE
         WRITE (Bandn(nord:nord+2),'(i3)') N
      ENDIF
      Bandm = 'evc.n'
      nord = INDEX(Bandm,'n',BACK=.FALSE.)
      IF ( M<10 ) THEN
         WRITE (Bandm(nord:nord),'(i1)') M
      ELSEIF ( M>=10 .AND. M<100 ) THEN
         WRITE (Bandm(nord:nord+1),'(i2)') M
      ELSE
         WRITE (Bandm(nord:nord+2),'(i3)') M
      ENDIF
 
 
 
 
      END
 
 
