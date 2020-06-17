SUBROUTINE paths(root1,k1,k2,n,m,pathk1,pathk2,bandn,bandm)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2020-06-17  Time: 13:03:20


CHARACTER (LEN=100), INTENT(IN OUT)      :: root1
INTEGER, INTENT(IN OUT)                  :: k1
INTEGER, INTENT(IN OUT)                  :: k2
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: m
CHARACTER (LEN=100), INTENT(OUT)         :: pathk1
CHARACTER (LEN=100), INTENT(OUT)         :: pathk2
CHARACTER (LEN=100), INTENT(OUT)         :: bandn
CHARACTER (LEN=100), INTENT(OUT)         :: bandm
INTEGER :: nord
CHARACTER (LEN=100) :: root2,folder,FILE,



root2='/MoS2.save'
FILE='/evc.dat'

folder='/K0000i'
nord=INDEX(folder,'i', back = .false.)
IF(k1 < 10)THEN
  WRITE (folder(nord:nord),'(i1)')k1
ELSE IF(k1 >= 10.AND.k1 < 100)THEN
  WRITE (folder(nord-1:nord),'(i2)')k1
ELSE IF(k1 >= 100.AND.k1 < 1000)THEN
  WRITE (folder(nord-2:nord),'(i3)')k1
ELSE IF(k1 >= 1000.AND.k1 < 10000)THEN
  WRITE (folder(nord-3:nord),'(i4)')k1
ELSE
  WRITE (folder(nord-4:nord),'(i5)')k1
END IF
pathk1=trim(root1)//trim(root2) //trim(folder)//trim(FILE)

folder='/K0000i'
nord=INDEX(folder,'i', back = .false.)
IF(k2 < 10)THEN
  WRITE (folder(nord:nord),'(i1)')k2
ELSE IF(k2 >= 10.AND.k2 < 100)THEN
  WRITE (folder(nord-1:nord),'(i2)')k2
ELSE IF(k2 >= 100.AND.k2 < 1000)THEN
  WRITE (folder(nord-2:nord),'(i3)')k2
ELSE IF(k2 >= 1000.AND.k2 < 10000)THEN
  WRITE (folder(nord-3:nord),'(i4)')k2
ELSE
  WRITE (folder(nord-4:nord),'(i5)')k2
END IF
pathk2=trim(root1)//trim(root2) //trim(folder)//trim(FILE)



bandn='evc.n'
nord=INDEX(bandn,'n', back = .false.)
IF(n < 10)THEN
  WRITE (bandn(nord:nord),'(i1)')n
ELSE IF(n >= 10.AND.n < 100)THEN
  WRITE (bandn(nord:nord+1),'(i2)')n
ELSE
  WRITE (bandn(nord:nord+2),'(i3)')n
END IF
bandm='evc.n'
nord=INDEX(bandm,'n', back = .false.)
IF(m < 10)THEN
  WRITE (bandm(nord:nord),'(i1)')m
ELSE IF(m >= 10.AND.m < 100)THEN
  WRITE (bandm(nord:nord+1),'(i2)')m
ELSE
  WRITE (bandm(nord:nord+2),'(i3)')m
END IF




RETURN
END SUBROUTINE paths


