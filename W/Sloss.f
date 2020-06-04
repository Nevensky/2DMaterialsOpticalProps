         program Surface_LOSS

c        PROGRAM FOR ab initio-SURFACE LOSS CALCULATION FOR LAYERED SYSTEMS
c        USING SUPERCELL METHOD  


c        NkI -number of wave vectors in irreducible B. zone
c        Ntot-total number of the mutually different wave vector-program generates this number  
c        Nband-number of the bands
c        NG-total number of G vectors  
c        NGd-number of coefficients CG shulod me less than minimum number of coefficients all over all evc.n files  
c        nMPx*nMPy*nMPz-Monkhorest-Pack sampling
c        Ef-Fermi energy 
c        T-temperature in eV 
c        Gama-Damping parameter in eV
c        Ecut-cutoff energy for crystal local field calculations  
c        Vcell-unit-cell volume in a.u.^3 
c        a0-unit cell parameter in parallel direction in a.u.  
c        c0-unit cell parameter in perpendicular direction in a.u. (z-separation between supercells)   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        Quantum Esspresso: 
c        verbosity           = 'high' 
c        VALID JUST FOR NORM-CONSERVING PSEUDOPOTENTIALS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 


        
         use iotk_module
         implicitnone      
         character(iotk_attlenx) :: attr
         logical :: found

         INTEGER NkI,Nband,ik,i,Nk,j,jk,it,lk,Ntot,iG0,nsim,iq,
     &   NG,io,no,nq,nMPx,nMPy,nMPz,n,m,iG,R1,K1,R2,K2,
     &   Nlf,NG1,NG2,NGd,iG1,iG2,jG,Nlfd,kG,jo,jump,loss,
     &   iGfast,ikmin,NelQE,lf,kG1,kG2,nord
         PARAMETER(nMPx=201,nMPy=201,nMPz=1,NkI=6835,Nband=60,NelQE=18,
     &   Nk=48*NkI,NGd=4000,NG=8000,no=2001,nq=2,Nlfd=50)
         ! no je broj frekvencija,  nq je broj valnih vektora tu je 2 jer je rucno paralelizirano!
         ! Nlf 
    

c        skalars
         REAL*8 a0,c0,eps,kx,ky,kz,KQx,KQy,KQz,qGx,qGy,qGz,omin,omax,
     &   qx,qy,qz,kmin,domega,o,Ef,K11,K22,K33,T,Lor,De,Gabs,kref,
     &   Eref,Ecut,Gxx1,Gyy1,Gzz1,Gxx2,Gyy2,Gzz2,fact,Vcell,oi,oj,
     &   ImChi0,ReChi0,Q,Gcar,aBohr,zero,Nel,absq,error,qmax,
     &   Gama,W1,W2,ImW,Wind,W2KK,KKS,SKK,WindKK,krefM 
         DOUBLEPRECISION pi,three,Hartree,planck,eV
         PARAMETER(Hartree=2.0d0*13.6056923d0,EF=0.5554/Hartree,
     &   a0=5.9715,c0=29.8575,pi=3.141592654d0,Gcar=2.0*pi/a0,
     &   eps=1.0d-4,T=0.01/Hartree,Gama=0.05/Hartree,
     &   eV=1.602176487d-19,planck=6.626196d-34,Ecut=0.0,
     &   Vcell=922.0586,aBohr=0.5291772d0,zero=0.0)
         DOUBLE COMPLEX a,ione,czero,rone,eM,G0
    
c        arrays
         INTEGER Gfast,Gi,parG,kQ0
         REAL*8 kI,E,R,RI,k,ktot,G,Glf,V,S0,KC,GlfV
         DOUBLE COMPLEX unit,epsilon,Chi
         COMPLEX*8 MnmK1K2,Chi0,WT,Gammap,Gammam 
         DIMENSION kI(3,NkI),E(NkI,Nband),R(48,3,3),RI(48,3,3),
     &   k(3,Nk),ktot(3,Nk),
     &   V(Nlfd,Nlfd),! matr. gole coulomb. int.
     &   G(3,NG), ! polje valnih vektora G u recp. prost. za wfn.
     &   GlfV(3,Nlfd), Glf(3,Nlfd) ! generiran G vekt. (0,0,z) za V odnosn za chi.
     &   MnmK1K2(Nlfd), ! nabojni vrhovi
     &   unit(Nlfd,Nlfd), ! jedinična matrica
     &   Chi0(Nlfd,Nlfd), ! (eq. 2.89)
     &   epsilon(Nlfd,Nlfd), ! Epsilon (GG')  = I - V(GG')Chi0
     &   Chi(Nlfd,Nlfd),   ! (eq. 2.88 nakon invertiranja)
     &   S0(no,Nlfd,Nlfd), ! korelacijska matrica
     &   Gfast(Nlfd*NGd),
     &   KC(3,3),Gi(3), ! pomocne funkcije
     &   parG(NG), ! paritet svakog valnog vektora
     &   WT(no,Nlfd,Nlfd), ! time ordered RPa screened coulomb int. (eq. 2.93)
     &   Gammap(Nlfd,Nlfd),Gammam(Nlfd,Nlfd), ! za GW ne koristi se za ovaj dio
     &   kQ0(100) 
 
         CHARACTER*100 bandn,
     &   bandm,
     &   nis,
     &   pathK1,pathK2, 
     &   dato,
     &   root,path,fajl
         CHARACTER*35 tag,buffer

         integer :: ngw, igwx, nbnd, nk1
         complex(kind = 8),pointer, dimension(:) :: C1,C2

         
          rone=dcmplx(1.0,0.0)  ! real 1
          czero=dcmplx(0.0,0.0) ! complex 0
          ione=dcmplx(0.0,1.0) ! imag 1


c BRAVAIS LATTICE PARAMETERS 

c     bravais-lattice index     =            4
c     lattice parameter (alat)  =       5.9715  a.u.
c     unit-cell volume          =     922.0586 (a.u.)^3
c     number of atoms/cell      =            3
c     number of atomic types    =            2
c     number of electrons       =        18.00
c     number of Kohn-Sham states=           13
c     kinetic-energy cutoff     =      50.0000  Ry
c     charge density cutoff     =     200.0000  Ry
c     convergence threshold     =      1.0E-06
c     mixing beta               =       0.7000
c     number of iterations used =            8  plain     mixing
c     Exchange-correlation      = SLA-PZ-NOGX-NOGC ( 1  1  0  0 0 0)

c     celldm(1)=   5.971535  celldm(2)=   0.000000  celldm(3)=   5.000000
c     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

c     crystal axes: (cart. coord. in units of alat)
c               a(1) = (   1.000000   0.000000   0.000000 )  
c               a(2) = (  -0.500000   0.866025   0.000000 )  
c               a(3) = (   0.000000   0.000000   5.000000 )  

c     reciprocal axes: (cart. coord. in units 2 pi/alat)
c               b(1) = (  1.000000  0.577350 -0.000000 )  
c               b(2) = (  0.000000  1.154701  0.000000 )  
c               b(3) = (  0.000000 -0.000000  0.200000 )  


c             QUANTUM ESSPRESSO IMPUTS:
              root='/home/vito/PROJECTS/MoS2-BSE/MoS2_201X201'
              
c             Crystal local field effects are included in z direction lf=1 
c             Crystal local field effects are included in x,y,z direction lf=3 


              lf=1 ! crystal local field effect included in z for lf=1
              jump=1 ! za 1 preskace trazenje wfn. u IBZ za sve bands m i n
              three=3.0d0 ! broj 3 haha
              omin=1.0d-5 ! raspon frekvencija u Ha
              omax=2.0d0
              domega=(omax-omin)/(no-1) 

         

c           CALL FOR POINT GROUP TRANSFORMATIONS
c           Point group transformations are in Cartesian coordinate 

            call PointR(root,nsim,R,RI)



c           Upis valnih vektora iz irreducibilne Brillouinove 
c           zone i pripadnih energijskih nivoa iz filea '****.band'. 
c           wave vectors are in Cartesian coordinate



            fajl='/MoS2.band'
            path=TRIM(root)//TRIM(fajl)
            open(1,FILE=path) 
            do 11 ik=1,NkI ! k vektora u IBZ
            if(ik.eq.1)READ(1,*)nis ! preskakanje
            READ (1,20)kI(1,ik),kI(2,ik),kI(3,ik)
            READ (1,10)(E(ik,i),i=1,Nband)
11          continue 
            CLOSE(1)
10          FORMAT(10F8.4)
20          FORMAT(10X,F10.3,F10.3,F10.3)  


            do 333 ik=1,NkI 
            do 444 i=1,Nband 
            E(ik,i)=E(ik,i)/Hartree
            if(i.ge.10)E(ik,i)=E(ik,i)+1.0/Hartree ! scissor op. ispravljanje DFT gapa na 1eV ( u ovom slucaju)
444         continue
333         continue

 



c            generator 1.B.Z. 
c            Dio programa koji pomocu operacija tockaste grupe i vektora iz 
c            I.B.Z. generira sve (MEDJUSOBNO RAZLICITE!!!) v. vektore u 1.B.Z.  
c            Ntot-Tot number of different points ''ktot'' inside 1.B.Z 

       

             jk=0 ! k tocka u ...?
             Ntot=0
             do 13 i=1,nsim ! loop over No. symmetries
             do 14 ik=1,NkI  ! loop over k points in IBZ
             it=1 ! ?
             jk=jk+1
             do 345 n=1,3  ! loop over kx,ky,kz 
             k(n,jk)=zero
             do 346 m=1,3 ! loop over x,y,z
             k(n,jk)=k(n,jk)+R(i,n,m)*kI(m,ik) ! kreira nove k tocke u BZ pomocu simetrije
346          continue                    
345          continue           
             if(jk.gt.1)then   
             do 17 lk=1,jk-1
             if(abs(k(1,jk)-k(1,lk)).le.eps)then  ! je li razlicita tocka od neke prije vec kreirane
             if(abs(k(2,jk)-k(2,lk)).le.eps)then      
             if(abs(k(3,jk)-k(3,lk)).le.eps)then    
             it=2         ! preskakanje tocke  
             endif
             endif
             endif
17           continue
             endif 
             if(it.eq.1)then ! ne postoji dodaj ju
             Ntot=Ntot+1  
             ktot(1,Ntot)=k(1,jk)
             ktot(2,Ntot)=k(2,jk)
             ktot(3,Ntot)=k(3,jk)
             endif
14           continue  
13           continue


c             Checking 1BZ integration
              Nel=0 ! provjeri je li broj el. u FBZ odgovara stvarnom broju el. u jed. cel. NelQE
              do 883 ik=1,Ntot
              kx=ktot(1,iK)
              ky=ktot(2,iK)
              kz=ktot(3,iK)
              do 442 n=1,Nband ! loop over bands
              if(n.eq.1)then
              it=1
              if(ik.le.nkI)then
              K1=ik
              it=2
              else
              do 631 i=2,nsim ! loop over no. symmetries
              K11=RI(i,1,1)*kx+RI(i,1,2)*ky+RI(i,1,3)*kz
              K22=RI(i,2,1)*kx+RI(i,2,2)*ky+RI(i,2,3)*kz
              K33=RI(i,3,1)*kx+RI(i,3,2)*ky+RI(i,3,3)*kz
              do 642 j=1,nkI ! loop over k u IBZ
              if(dabs(K11-KI(1,j)).le.eps)then ! eps proizvoljno mali broj
              if(dabs(K22-KI(2,j)).le.eps)then 
              if(dabs(K33-KI(3,j)).le.eps)then 
              it=2
              K1=j  
              goto 5022
              endif
              endif
              endif
642           continue
631           continue
              endif
              if(it.eq.1)then
              print*,'Can not find wave vector K=',iK,
     &        'in I.B.Z.'  
              stop 
              endif
5022          continue         
              endif
              if(E(k1,n).lt.Ef)Nel=Nel+1.0 ! zbroji za en. manje od fermijeve
442           continue
883           continue
              Nel=2.0*Nel/Ntot

              do 521 i=1,ntot
521           write(887,*)ktot(1,i),ktot(2,i) ! output da vidimo kako izgleda FBZ



c           KC transformation matrix from rec.cryst. axes to cart.koord.    
c           If g' is vector in rec.cryst. axes then a=KC*a' is vector in cart. axes 
           

            fajl='/MoS2.sc.out'
            path=TRIM(root)//TRIM(fajl)
            tag='     reciprocal axes: (cart. coord.' 
            open(1,FILE=path) 
            do 166 i=1,100000
            read(1,'(a)')buffer
            if(buffer.eq.tag)then 
            do 753 j=1,3       
753         read(1,70)KC(1,j),KC(2,j),KC(3,j)
            goto 998
            endif
166         continue
70          FORMAT(23X,3F10.3)  
998         continue
            close(1)

                   
c           Reading the reciprocal vectors in crystal coordinates and transformation 
c           in Cartezi cordinates.             
            OPEN (1,FILE='gvectors.dat')
            do 24 i=1,8
24          READ(1,*)nis
            do 25 iG=1,NG
            READ(1,100)Gi(1),Gi(2),Gi(3)
            if(iG.eq.1)then
            if(Gi(1).ne.0.or.Gi(2).ne.0.or.Gi(3).ne.0)then
            print*,'*********************************'
            print*,'WARRNING!, G vectors input is wrong!!'
            print*,'G(1) is not (0,0,0)!!'
            stop
            endif
            endif
c           transformation in cart.coord (also!, after this all G components are in 2pi/a0 units)
            do 776  n=1,3 ! loop over dim
            G(n,iG)=zero
            do 777  m=1,3
            G(n,iG)=G(n,iG)+KC(n,m)*dble(Gi(m))
777         continue
776         continue 
            parG(iG)=Gi(3)  ! odreduje paritet za svaki Gi
25          continue
100         FORMAT(I10,I11,I11)
            CLOSE(1) 


           
c            Reciprocal vectors for crystal local field effects calculations in array ''Glf(3,Nlf)'' 
! na temelju cutoffa eliminar G-vektore koji su izvan cut-offa
             Nlf=0 ! broj novih local field G vektora
             if(lf.eq.1)then
             do 813 iG=1,NG ! loop over all G-vectors za 3D LFE
             if(G(1,iG).eq.0.0.and.G(2,iG).eq.0.0)then
             Eref=Gcar*Gcar*G(3,iG)*G(3,iG)/2.0 ! energijski prag
             if(Eref.le.Ecut)then 
             Nlf=Nlf+1
             Glf(1,Nlf)=0.0
             Glf(2,Nlf)=0.0
             Glf(3,Nlf)=G(3,iG)
             if((parG(iG)/2)*2.eq.parG(iG))then
             parG(Nlf)=1    
             else           
             parG(Nlf)=-1
             endif 
             endif
             endif
813          continue
             else
             do 812 iG=1,NG ! loop over all G-vectors za 1D LFE (samo u z-smjeru)
             Eref=Gcar*Gcar*(G(1,iG)*G(1,iG)+G(2,iG)*G(2,iG)+
     &       G(3,iG)*G(3,iG))/2.0
             if(Eref.le.Ecut)then 
             Nlf=Nlf+1
             Glf(1,Nlf)=G(1,iG)
             Glf(2,Nlf)=G(2,iG)
             Glf(3,Nlf)=G(3,iG)
             if((parG(iG)/2)*2.eq.parG(iG))then
             parG(Nlf)=1    
             else           
             parG(Nlf)=-1
             endif 
             endif
812          continue
             endif              
             if(Nlf.gt.Nlfd)then
             print*,'Nlf is bigger than Nlfd'
             stop
             endif




c             IBZ q LOOP STARTS HERE!!! 


! iq=0 ne moze biti nula, opticki racun
! iq=2 do iq=...cutoff transfer q vektor!
! ikmin = min. valni vektor u BZ svi veci su visekratnici tog minimalnog
              do 801 iq=42,61 
  
c             searching min. q=(qx,qy,qz) in GM direction                
              kmin=1.0
              do 913 i=1,Ntot
              kref=sqrt(ktot(1,i)*ktot(1,i)+
     &        ktot(2,i)*ktot(2,i)+ktot(3,i)*ktot(3,i)) 
              if(kref.eq.zero)goto 970
              if(kref.lt.kmin)then
              kmin=kref      
              ikmin=i
              krefM=kmin
              endif
970           continue
913           continue
 

              qx=(iq-1)*ktot(1,ikmin)
              qy=(iq-1)*ktot(2,ikmin)
              qz=(iq-1)*ktot(3,ikmin)

              absq=sqrt(qx*qx+qy*qy+qz*qz)

c             Info file

           open(55,file='Info')
           write(55,*)'***************General***********************'
           write(55,*)''
           write(55,*)'Number of point symmetry operation is',nsim	
           write(55,88)'Wave vector (qx,qy,qz)=(',qx*Gcar,qy*Gcar,
     &     qz*Gcar,') a.u.'	                
           write(55,99)'|(qx,qy,qz)|=',absq*Gcar,'a.u.'
           if(lf.eq.1)write(55,*)'Local field effcts in z-dir'
           if(lf.eq.3)write(55,*)'Local field in all xyz-dir'
           write(55,*)'Number of local field vectors is',Nlf
           write(55,*)'Number of different K vectors in 1.B.Z. is',Ntot
           write(55,*)'Number of K vectors in I.B.Z. is',NkI
           write(55,*)'Number of bands is               ',Nband
           write(55,99)'gama dumping is ',gama*Hartree*1000.0,'meV'
           write(55,99)'Temperature is  ',T*Hartree*1000.0,'meV'
           write(55,*)''
           write(55,*)'************* Checking 1BZ integration*******'
           write(55,*)''
           write(55,12)'Number of electrons(1BZ integration)=',Nel
           write(55,*)'Number of electrons(unit cell)=',NelQE
           error=abs((NelQE-Nel)/NelQE)
           write(55,99)'Relative error=',error*100.0,'%'
           if(error.gt.0.05)then
           write(55,*)'WARRNING!!-1BZ INTEGRATION IS BAD!.'
           endif
           close(55)
88         FORMAT(A25,3F10.4,A5)     
99         FORMAT(A25,F8.4,A5)
12         FORMAT(A40,F7.4)          


           do 991 io=1,no
           do 661 iG=1,Nlf 
           do 662 jG=1,Nlf
           S0(io,iG,jG)=czero ! korelacijska funkcija (inicijalizacija)
662        continue
661        continue
991        continue
             

c              1.B.Z  LOOP STARTS HERE !!!!                   

               do 803 ik=1,Ntot ! loop over k-points in FBZ

               open(122,file='status') 
               write(122,*)'iq=',iq
               write(122,*)'ik=',ik
               close(122)                          


               kx=ktot(1,iK)
               ky=ktot(2,iK)
               kz=ktot(3,iK)

c              trazenje (kx,ky,kz) u ireducibilnoj zoni

               it=1
               if(ik.le.nkI)then
               R1=1
               K1=ik
               it=2
               else
               do 601 i=2,nsim 
               K11=RI(i,1,1)*kx+RI(i,1,2)*ky+RI(i,1,3)*kz
               K22=RI(i,2,1)*kx+RI(i,2,2)*ky+RI(i,2,3)*kz
               K33=RI(i,3,1)*kx+RI(i,3,2)*ky+RI(i,3,3)*kz
               do 602 j=1,nkI
               if(dabs(K11-KI(1,j)).le.eps)then 
               if(dabs(K22-KI(2,j)).le.eps)then 
               if(dabs(K33-KI(3,j)).le.eps)then 
               it=2
               R1=i ! pridruzena point group
               K1=j  ! trazeni vektor u IBZ
               goto 5222
               endif
               endif
               endif
602            continue
601            continue
               endif
               if(it.eq.1)then
               print*,'Can not find wave vector K=',iK,
     &         'in I.B.Z.'  
               stop 
               endif
5222           continue            
              

               it=1
               KQx=kx+qx
               KQy=ky+qy
               KQz=kz+qz

c              trazenje (KQx,KQy) prvo u 1.B.Z a onda u I.B.Z. (jer novo genrirani k+q more bit negdje vani u FBZ)

               do 701 iG=1,NG
                  do 702 jK=1,Ntot
                     if(dabs(KQx-G(1,iG)-ktot(1,jK)).le.eps)then 
                     if(dabs(KQy-G(2,iG)-ktot(2,jK)).le.eps)then 
                     if(dabs(KQz-G(3,iG)-ktot(3,jK)).le.eps)then 
                     it=2
                     iG0=iG
                     do 501 i=1,nsim 
                     K11=RI(i,1,1)*ktot(1,jK)+RI(i,1,2)*ktot(2,jK)+
     &               RI(i,1,3)*ktot(3,jK)
                     K22=RI(i,2,1)*ktot(1,jK)+RI(i,2,2)*ktot(2,jK)+
     &               RI(i,2,3)*ktot(3,jK)
                     K33=RI(i,3,1)*ktot(1,jK)+RI(i,3,2)*ktot(2,jK)+
     &               RI(i,3,3)*ktot(3,jK)
                     do 502 j=1,nkI
                     if(dabs(K11-KI(1,j)).le.eps)then 
                     if(dabs(K22-KI(2,j)).le.eps)then 
                     if(dabs(K33-KI(3,j)).le.eps)then 
                     it=3
                     R2=i  ! pridruzena point gruop
                     K2=j  ! pridruzen taj vektor pronadjen u IBZ
                     goto 2111 
                     endif
                     endif
                     endif
502            continue
501            continue
                     endif
                     endif
                     endif
702           continue
701           continue      
2111                 continue           


                if(it.eq.1)then
                print*,'Can not find wave vector K+Q=',iK,'+',iQ,
     &          'in 1.B.Z.'  
                stop 
                elseif(it.eq.2)then 
                print*,'Can not find wave vector K+Q=',iK,'+',iQ,
     &          'in I.B.Z.'  
                stop
                endif


c              R1-integer, redni broj point operacije R1 u transformaciji ''K=R1*K1''.
c              K1-integer, redni broj valnog vektora K1 u transformaciji ''K=R1*K1''.
c              iG0 i R2-integeri, redni broj vektora reciprocne restke G0 i point operacije R2 u transformaciji ''K+Q=G0+R2*K2''.
c              K2-integer, redni broj valnog vektora K2 u transformaciji  ''K+Q=G0+R2*K2''.


c              petlje po vrpcama n i m  

               do 804 n=1,9 ! loop over full bands , oprez kod metala mora ici malo iznad popunjene (n+1) 
               do 805 m=10,nband !  loop over empty bands
            




            call paths(root,K1,K2,n,m,pathK1,pathK2,bandn,bandm) ! fajl za popunit wfn. za K1 i K2 i pripadajuce vrpce n i m

               
c         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za 
c         fiksni K1,K2,n i m 

c               Otvaranje atribute za INFO
                call iotk_open_read(10,pathK1) 
                call iotk_scan_empty(10,"INFO",attr=attr) ! nalazi <info>
                call iotk_scan_attr(attr,"igwx",NG1) ! koliko imamo G-vektora /
c               Alociranje polja C1
                allocate (C1(NG1)) ! u ovo trpamo fourierove koeficijent u razvoju wfn. 
c               Ucitavanje podataka iza evc.n
                call iotk_scan_dat(10,bandn,C1) ! cita iz binranog oblika
                call iotk_close_read(10) 
C               Otvaranje atribute za INFO
                call iotk_open_read(10,pathK2)
                call iotk_scan_empty(10,"INFO",attr=attr)
                call iotk_scan_attr(attr,"igwx",NG2)
c               Alociranje polja C2
                allocate (C2(NG2))
c               Ucitavanje podataka iza evc.m
                call iotk_scan_dat(10,bandm,C2)
                call iotk_close_read(10) 

c                Konstrukcija stupca matricnih elementa MnmK1K2(G)  
          

                  if(NGd.gt.NG1)then
                  write(*,*)'NGd is bigger than NG1=',NG1
                  stop
                  elseif(NGd.gt.NG2)then
                  write(*,*)'NGd is bigger than NG2=',NG2
                  stop
                  endif              
                   

c                 matrix elements
                  iGfast=0                 
                  do 551 iG=1,Nlf
                  MnmK1K2(iG)=czero
                  do 552 iG1=1,NGd 
                  iGfast=iGfast+1
                  Gxx1=G(1,iG1)
                  Gyy1=G(2,iG1)
                  Gzz1=G(3,iG1)
                  K11=R(R1,1,1)*Gxx1+R(R1,1,2)*Gyy1+R(R1,1,3)*Gzz1
                  K22=R(R1,2,1)*Gxx1+R(R1,2,2)*Gyy1+R(R1,2,3)*Gzz1
                  K33=R(R1,3,1)*Gxx1+R(R1,3,2)*Gyy1+R(R1,3,3)*Gzz1
                  K11=K11+Glf(1,iG)
                  K22=K22+Glf(2,iG)
                  K33=K33+Glf(3,iG)
                  K11=K11+G(1,iG0)
                  K22=K22+G(2,iG0)
                  K33=K33+G(3,iG0)
                  Gxx1=RI(R2,1,1)*K11+RI(R2,1,2)*K22+RI(R2,1,3)*K33
                  Gyy1=RI(R2,2,1)*K11+RI(R2,2,2)*K22+RI(R2,2,3)*K33
                  Gzz1=RI(R2,3,1)*K11+RI(R2,3,2)*K22+RI(R2,3,3)*K33
                  if(jump.eq.1)then 
                  do 553 iG2=1,NG2
                  Gfast(iGfast)=NG2+1
                  Gxx2=G(1,iG2)
                  Gyy2=G(2,iG2)
                  Gzz2=G(3,iG2)
                  if(dabs(Gxx2-Gxx1).lt.eps)then 
                  if(dabs(Gyy2-Gyy1).lt.eps)then 
                  if(dabs(Gzz2-Gzz1).lt.eps)then 
                  Gfast(iGfast)=iG2
                  goto 1111 
                  endif
                  endif
                  endif
553               continue  
                  endif
1111              continue 
                  iG2=Gfast(iGfast)
                  if(iG2.le.NG2)then
                  MnmK1K2(iG)=MnmK1K2(iG)+
     &            conjg(C1(iG1))*C2(iG2)
                  endif
552               continue
551               continue
                  jump=2             


c                 omega loop 
                  do 802 io=1,no
                     o=(io-1)*domega
                     De=o+E(K1,n)-E(K2,m)
                      Lor=-Gama/(De*De+Gama*Gama)
                      if(dabs(Lor).ge.1.0d-3/Gama)then  
                       do 513 iG=1,Nlf                   
                       do 554 jG=1,Nlf 
                       S0(io,iG,jG)=S0(io,iG,jG)-
     &                 2.0*Lor*MnmK1K2(iG)*conjg(MnmK1K2(jG))/
     &                 (pi*Ntot*Vcell)                
554                    continue
513                    continue
                     endif   

                   
802               continue

        	deallocate(C1)
		deallocate(C2)


222             continue

                 

c                end of m do loop
805              continue

c                end of n do loop  
804              continue
                  jump=1                 
                
834              continue
c                ond of 1.B.Z do loop 
803              continue



c               Puting (qx,qy,qz) and Glf in cartezi coordinate 

                qx=Gcar*qx 
                qy=Gcar*qy
                qz=Gcar*qz

                do 703 iG=1,Nlf 
                GlfV(1,iG)=Gcar*Glf(1,iG)
                GlfV(2,iG)=Gcar*Glf(2,iG)
                GlfV(3,iG)=Gcar*Glf(3,iG)
703             continue

               

c               new sum over omega
                 do 872 io=1,no-1
c                 print*,io            
                 oi=(io-1)*domega
                 do 432 iG=1,Nlf   
                 do 433 jG=1,Nlf 
                 ReChi0=0.0
c                  static limit
                   if(io.eq.1)then
                      do 431 jo=2,no
                      oj=(jo-1)*domega
                      fact=domega/oj
                      if(jo.eq.2)fact=3.0/2.0 
                      if(jo.eq.no)fact=0.5*domega/oj
                      ReChi0=ReChi0+fact*S0(jo,iG,jG)
431                   continue
                      ReChi0=-2.0*ReChi0
                   elseif(io.eq.2)then 
                      do 434 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=1.0
                      if(jo.eq.2)fact=0.0
                      if(jo.eq.3)fact=-3.0/2.0 
                      if(jo.eq.no)fact=0.5*domega/(oi-oj)
                      ReChi0=ReChi0+fact*S0(jo,iG,jG)
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      ReChi0=ReChi0-fact*S0(jo,iG,jG)
434                   continue
                   elseif(io.eq.(no-1))then  
                      do 435 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=0.5*domega/(oi-oj)
                      if(jo.eq.(no-2))fact=3.0/2.0
                      if(jo.eq.(no-1))fact=0.0
                      if(jo.eq.no)fact=-1.0
                      ReChi0=ReChi0+fact*S0(jo,iG,jG)
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      ReChi0=ReChi0-fact*S0(jo,iG,jG)
435                   continue
                   else 
                      do 436 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=0.5*domega/(oi-oj)
                      if(jo.eq.(io-1))fact=3.0/2.0
                      if(jo.eq.io)fact=0.0
                      if(jo.eq.(io+1))fact=-3.0/2.0
                      if(jo.eq.no)fact=0.5*domega/(oi-oj)
                      ReChi0=ReChi0+fact*S0(jo,iG,jG)
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      ReChi0=ReChi0-fact*S0(jo,iG,jG)
436                   continue
                   endif 
 
                      ImChi0=-pi*S0(io,iG,jG)
                      Chi0(iG,jG)=cmplx(ReChi0,ImChi0) 


c                kraj po iG,jG              
433              continue    
432              continue 

                      

c                Calculation of the ''Chi''  by matrix invertion 


c                MATRIX V(G,G')
                
                 do 567 iG=1,Nlf
                 Gabs=sqrt((Qx+GlfV(1,iG))*(Qx+GlfV(1,iG))+
     &           (Qy+GlfV(2,iG))*(Qy+GlfV(2,iG)))
                 if(Gabs.eq.0.0)Gabs=eps
                 do 568 jG=1,Nlf 
                 V(iG,jG)=0.0
                 if(Glf(1,jG).eq.Glf(1,iG))then 
                 if(Glf(2,jG).eq.Glf(2,iG))then 
                 V(iG,jG)=4.0*pi*(1.0-exp(-Gabs*c0))/(Gabs*c0)
                 V(iG,jG)=V(iG,jG)*(Gabs*Gabs-GlfV(3,iG)*GlfV(3,jG))
                 V(iG,jG)=V(iG,jG)/(Gabs*Gabs+GlfV(3,iG)*GlfV(3,iG))
                 V(iG,jG)=V(iG,jG)/(Gabs*Gabs+GlfV(3,jG)*GlfV(3,jG))
                 V(iG,jG)=-dble(parG(iG))*dble(parG(jG))*V(iG,jG)
                 if(Glf(3,jG).eq.Glf(3,iG))then 
                 V(iG,jG)=4.0*pi/(Gabs*Gabs+GlfV(3,iG)*GlfV(3,iG))+
     &           V(iG,jG)
                 endif 
                 endif
                 endif
568              continue 
567              continue

                            

                 do 828 iG=1,Nlf
                 do 818 jG=1,Nlf
                 unit(iG,jG)=czero
818              continue
                 unit(iG,iG)=rone
828              continue

                 do 704 iG=1,Nlf
                 do 705 jG=1,Nlf
                 epsilon(iG,jG)=unit(iG,jG)
                 do 706 kG=1,Nlf
                 epsilon(iG,jG)=epsilon(iG,jG)-Chi0(iG,kG)*V(kG,jG)
706              continue
705              continue
704              continue


c                invertiranje matrice ''epsilon = 1-Chi_0*V''


                 call gjel(epsilon,Nlf,Nlfd,unit,Nlf,Nlfd)


                 do 604 iG=1,Nlf
                 do 605 jG=1,Nlf
                 Chi(iG,jG)=czero
                 do 606 kG=1,Nlf
                 Chi(iG,jG)=Chi(iG,jG)+epsilon(iG,kG)*Chi0(kG,jG)
606              continue
605              continue
604              continue

                  
c                SCREENED COULOMB INTERACTION W^T_GG'(Q,\omega) 

                 do 634 iG=1,Nlf
                 do 635 jG=1,Nlf
                 WT(io,iG,jG)=czero
                 do 636 kG1=1,Nlf
                 do 637 kG2=1,Nlf
                 WT(io,iG,jG)=WT(io,iG,jG)+
     &           V(iG,kG1)*Chi(kG1,kG2)*V(kG2,jG)
637              continue 
636              continue
                 WT(io,iG,jG)=V(iG,jG)+WT(io,iG,jG)
635              continue
634              continue

c                kraj nove petlje po omega
872              continue
 

c               ispis time ordered zasjenjene kulonske interakcije W_GG'^T(Q,\omega)                
                dato='W_Qi'
                nord=INDEX(dato,'i', BACK =.false.)
                if(iq.lt.10)then
                write(dato(nord:nord),'(i1)')iq
                elseif(iq.ge.10.and.iq.lt.100)then
                write(dato(nord:nord+1),'(i2)')iq
                else
                write(dato(nord:nord+2),'(i3)')iq
                endif

                open(74,file=dato)
                do 678 io=1,1
                o=(io-1)*domega
                write(74,*)'omega=',o,'Hartree'
                write(74,44)((WT(io,iG,jG),jG=1,Nlf),iG=1,Nlf)
678             continue  
                CLOSE(74)
44              FORMAT(10F15.5) 

                 do 665 io=1,no-1 
                 do 666 iG=1,Nlf  
                 do 667 jG=1,Nlf 
667              S0(io,iG,jG)=-(1.0/pi)*imag(WT(io,iG,jG))
666              continue 
665              continue 


                 KKS=zero
                 SKK=zero


c                new sum over omega
                 do 172 io=1,no-1
c                 print*,io            
                 oi=(io-1)*domega
                 do 132 iG=1,Nlf   
                 do 133 jG=1,Nlf 
                 W1=0.0
                 W2=0.0                
c                static limit
                   if(io.eq.1)then
                      do 131 jo=2,no
                      oj=(jo-1)*domega
                      fact=domega/oj
                      if(jo.eq.2)fact=3.0/2.0 
                      if(jo.eq.no)fact=0.5*domega/oj
                      W1=W1-fact*S0(jo,iG,jG)
131                   continue
                      W2=-W1
                   elseif(io.eq.2)then 
                      do 134 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=1.0
                      if(jo.eq.2)fact=0.0
                      if(jo.eq.3)fact=-3.0/2.0 
                      if(jo.eq.no)fact=0.5*domega/(oi-oj)
                      W1=W1+fact*S0(jo,iG,jG)
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      W2=W2+fact*S0(jo,iG,jG)
134                   continue
                   elseif(io.eq.(no-1))then  
                      do 135 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=0.5*domega/(oi-oj)
                      if(jo.eq.(no-2))fact=3.0/2.0
                      if(jo.eq.(no-1))fact=0.0
                      if(jo.eq.no)fact=-1.0
                      W1=W1+fact*S0(jo,iG,jG)
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      W2=W2+fact*S0(jo,iG,jG)
135                   continue
                   else 
                      do 136 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=0.5*domega/(oi-oj)
                      if(jo.eq.(io-1))fact=3.0/2.0
                      if(jo.eq.io)fact=0.0
                      if(jo.eq.(io+1))fact=-3.0/2.0
                      if(jo.eq.no)fact=0.5*domega/(oi-oj)
                      W1=W1+fact*S0(jo,iG,jG)
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      W2=W2+fact*S0(jo,iG,jG)
136                   continue
                   endif 
 
                      ImW=-pi*S0(io,iG,jG)
                      Gammap(iG,jG)=cmplx(W1,ImW) 
                      Gammam(iG,jG)=cmplx(-W2,0.0) 
                      if(ig.eq.1.and.jg.eq.1)W2KK=W2 
                      if(ig.eq.1.and.jg.eq.1.and.io.eq.1)G0=Gammap(1,1) 

c                kraj po iG,jG              
133              continue    
132              continue 

c                Provjera KK relacija 
                 Wind=real(WT(io,1,1)-V(1,1))
                 WindKK=real(Gammap(1,1))-W2KK
                 fact=domega
                 if(io.eq.1.or.io.eq.no-1)fact=0.5*domega
                 KKS=KKS+fact*(WindKK-Wind)*(WindKK-Wind)
                 SKK=SKK+fact*Wind*Wind


c                kraj nove petlje po omega
172              continue
                 CLOSE(74)


                 dato='Kramers-Kron_Qi'
                 nord=INDEX(dato,'i', BACK =.false.)
                 if(iq.lt.10)then
                 write(dato(nord:nord),'(i1)')iq
                 elseif(iq.ge.10.and.iq.lt.100)then
                 write(dato(nord:nord+1),'(i2)')iq
                 else
                 write(dato(nord:nord+2),'(i3)')iq
                 endif

                 open(33,file=dato)
           write(33,88)'Wave vector (qx,qy,qz)=(',qx*Gcar,qy*Gcar,
     &     qz*Gcar,') a.u.'	                
           write(33,99)'|(qx,qy,qz)|=',absq*Gcar,'a.u.'
                 write(33,*)'int(WindKK-Wind)^2 =  ',KKS
                 write(33,*)'int(Wind)^2 =  ',SKK
                 write(33,*)'****************************************'
                 write(33,*)'Kramers–Kronig relation relative error'
                 write(33,80)100.0*abs(KKS/SKK), '%'
78               format(A23,F10.5)
79               format(A16,F10.5)
80               format(5X,F7.2,A2)
                 write(33,*)'Usporedba Gamma i WT'
                 write(33,*)'real[Gamma(o=0,1,1)]=',real(G0) 
                 write(33,*)'real[WT(o=0,1,1)]/2=',
     &           real(WT(1,1,1)-V(1,1))/2.0
                 CLOSE(33)
                

c                kraj po q 
801              continue
                
                   end 

