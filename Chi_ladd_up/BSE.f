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
c        nq-number of wave vectors in W calculation     
c        Gama-Damping parameter in eV
c        Vcell-unit-cell volume in a.u.^3 
c        a0-unit cell parameter in parallel direction in a.u.  
c        c0-unit cell parameter in perpendicular direction in a.u. (z-separation between supercells)   
c        Ecut-- crystal local field cut-off for W 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        Quantum Esspresso: 
c        verbosity           = 'high' 
c        VALID JUST FOR NORM-CONSERVING PSEUDOPOTENTIALS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 


        
         use iotk_module
         implicitnone      
         character(iotk_attlenx) :: attr
         logical :: found

         INTEGER NkI,Nband,ik,i,Nk,j,jk,it,lk,Ntot,nsim,iK1,iK2,
     &   NG,io,no,iq,nq,nMPx,nMPy,nMPz,n,m,iG,R1,K1,R2,K2,iG01,iG02,
     &   Nlf,NG1,NG2,NGd,iG1,iG2,jG,Nlfd,kG,jo,jump,ni,mi,nj,mj,mod,
     &   iGfast,ikmin,NelQE,lf,kG1,kG2,nord,io1,io2,NKBZ,ref,nKd,
     &   pol,spin,nval,nvalLS,nmax,mmax 
         PARAMETER(nMPx=51,nMPy=51,nMPz=1,NkI=460,Nband=40,NelQE=26,
     &   Nk=48*NkI,NGd=4000,NG=8000,no=401,nq=121,Nlfd=50,nKd=2700)
    

c        skalars
         REAL*8 a0,c0,eps,kx,ky,kz,omin,omax,sgn,struja,
     &   qx,qy,qz,kmin,domega,o,Ef,K11,K22,K33,T,Lor,De,Gabs,kref,
     &   Eref,Ecut,Gxx1,Gyy1,Gzz1,Gxx2,Gyy2,Gzz2,fact,Vcell,oi,oj,
     &   Q,Gcar,aBohr,zero,Nel,absq,error,o1,o2,Gama,Vbare,qmax,
     &   Estar,KpointX,KpointY,KpointZ,qmin,vQ,GW,alpha,epsq,
     &   valley
         DOUBLEPRECISION pi,three,Hartree,planck,eV,eta
         PARAMETER(Hartree=2.0d0*13.6056923d0,EF=1.5703/Hartree,
     &   a0=5.9715,c0=29.8575,pi=3.141592654d0,Gcar=2.0*pi/a0,
     &   eps=1.0d-4,T=0.01/Hartree,eta=0.05/Hartree,eV=1.602176487d-19,
     &   planck=6.626196d-34,Ecut=0.0,Vcell=922.0586,
     &   aBohr=0.5291772d0,zero=0.0)
         DOUBLE COMPLEX a,ione,czero,rone,eM,G0,W,L0i,L0j,j1,j2


c        arrays
         INTEGER Gfast,Gi,parG
         REAL*8 kI,E,R,RI,k,ktot,G,Glf,KC,V,KBZ,Delta
         COMPLEX*8 MnmK1K2,WT,FockK,Xi0,Lladd,current,chi_ladd,
     &   chi_para 
         DOUBLE COMPLEX Fock,unit
         DIMENSION kI(3,NkI),E(NkI,Nband),R(48,3,3),RI(48,3,3),
     &   k(3,Nk),ktot(3,Nk),G(3,NG),Glf(3,Nlfd),MnmK1K2(2,Nlfd),
     &   Gfast(Nlfd*NGd),KC(3,3),Gi(3),WT(nq,Nlfd,Nlfd),parG(NG),
     &   KBZ(3,Nk),FockK(nKd,nKd),Delta(nKd),Fock(nKd,nKd),
     &   unit(nKd,nKd),Xi0(nKd,nKd),Lladd(nKd,nKd),current(2,nkd),
     &   chi_ladd(2),chi_para(2)
      
 
         CHARACTER*100 bandn,bandm,nis,pathK1,pathK2,dato,
     &   root,path,fajl,root1,root2
         CHARACTER*35 tag,buffer

         integer :: ngw, igwx, nbnd, nk1
         complex(kind = 8),pointer, dimension(:) :: C1,C2

         
          rone=dcmplx(1.0,0.0) 
          czero=dcmplx(0.0,0.0)
          ione=dcmplx(0.0,1.0)


c BRAVAIS LATTICE PARAMETERS 

c     bravais-lattice index     =            4
c     lattice parameter (alat)  =       5.9715  a.u.
c     unit-cell volume          =     922.0586 (a.u.)^3
c     number of atoms/cell      =            3
c     number of atomic types    =            2
c     *******************************OVO JE BROJ ELEKTRONA KAD IMAS FULL REL. PSEUDOPOTENTIAL     
c     number of electrons       =        26.00
c     ************************************************************************************** 
c     number of Kohn-Sham states=           34
c     kinetic-energy cutoff     =      50.0000  Ry
c     charge density cutoff     =     200.0000  Ry
c     convergence threshold     =      1.0E-06
c     mixing beta               =       0.7000
c     number of iterations used =            8  plain     mixing
c     Exchange-correlation      = PBE ( 1  4  3  4 0 0)
c     Non magnetic calculation with spin-orbit


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
              root='/home/vito/PROJECTS/MoS2-BSE/MoS2_51x51/'
c             Staticka zasjenjena kulonska interakcija WT je smjestena u folderu
              root1='/home/vito/PROJECTS/MoS2-BSE/W/'               



c             Crystal local field effects are included in z direction lf=1 
c             Crystal local field effects are included in x,y,z direction lf=3 
c             If mod=1 Fock kernel calculation 
c             If mod=2 Ladder ireducible current-current polarisability calculation
c             SAMO ZA MOD=2
c             For spin up put spin=1 
c             For spin down put spin=2             

              mod=2
              spin=1
              lf=1
              jump=1
              three=3.0d0
              omin=1.0d-5
              omax=4.0/Hartree
              domega=(omax-omin)/(no-1)
            

c             maximum transfer wave vector around K point..       
c             qmax=sqrt(3.0)/3.0
              qmax=0.3
              if(qmax.le.(1.0/3.0))valley=2.0
              if(qmax.gt.(1.0/3.0))valley=1.0
c             minimum transfer wave vector in WT calculation
              qmin=0.006 
c             broj valentnih vrpci bez LS vezanja  
              nval=9
c             broj valentnih vrpci sa LS vezanjem (full-rel. PP) 
              nvalLS=26
c             band gap correction (in eV)
              GW=1.0
             


c           CALL FOR POINT GROUP TRANSFORMATIONS
c           Point group transformations are in Cartesian coordinate 

            call PointR(root,nsim,R,RI)



c           Upis valnih vektora iz irreducibilne Brillouinove 
c           zone i pripadnih energijskih nivoa iz filea '****.band'. 
c           wave vectors are in Cartesian coordinate


            fajl='MoS2_LS.band'
            path=TRIM(root)//TRIM(fajl)
            open(1,FILE=path) 
            do 11 ik=1,NkI 
            if(ik.eq.1)READ(1,*)nis
            READ (1,20)kI(1,ik),kI(2,ik),kI(3,ik)
            READ (1,10)(E(ik,i),i=1,Nband)
11          continue 
            CLOSE(1)
10          FORMAT(10F8.4)
20          FORMAT(10X,F10.3,F10.3,F10.3)  


            do 333 ik=1,NkI 
            do 444 i=1,Nband 
            E(ik,i)=E(ik,i)/Hartree
            if(i.gt.nvalLS)E(ik,i)=E(ik,i)+GW/Hartree
444         continue
333         continue

            print*,
            print*,
            print*,'Be careful MoS2 gap should be GW corrected'
            print*,'Band structure is full relativistic'
            print*,




c            generator 1.B.Z. 
c            Dio programa koji pomocu operacija tockaste grupe i vektora iz 
c            I.B.Z. generira sve (MEDJUSOBNO RAZLICITE!!!) v. vektore u 1.B.Z.  
c            Ntot-Tot number of different points ''ktot'' inside 1.B.Z 

       

             jk=0
             Ntot=0
             do 13 i=1,nsim
             do 14 ik=1,NkI 
             it=1 
             jk=jk+1
             do 345 n=1,3
             k(n,jk)=zero
             do 346 m=1,3 
             k(n,jk)=k(n,jk)+R(i,n,m)*kI(m,ik)
346          continue                    
345          continue           
             if(jk.gt.1)then   
             do 17 lk=1,jk-1
             if(abs(k(1,jk)-k(1,lk)).le.eps)then 
             if(abs(k(2,jk)-k(2,lk)).le.eps)then      
             if(abs(k(3,jk)-k(3,lk)).le.eps)then    
             it=2           
             endif
             endif
             endif
17           continue
             endif 
             if(it.eq.1)then
             Ntot=Ntot+1
             ktot(1,Ntot)=k(1,jk)
             ktot(2,Ntot)=k(2,jk)
             ktot(3,Ntot)=k(3,jk)
             endif
14           continue  
13           continue


c             Checking 1BZ integration
              Nel=0
              do 883 ik=1,Ntot
              kx=ktot(1,iK)
              ky=ktot(2,iK)
              kz=ktot(3,iK)
              do 442 n=1,Nband
              if(n.eq.1)then
              it=1
              if(ik.le.nkI)then
              K1=ik
              it=2
              else
              do 631 i=2,nsim 
              K11=RI(i,1,1)*kx+RI(i,1,2)*ky+RI(i,1,3)*kz
              K22=RI(i,2,1)*kx+RI(i,2,2)*ky+RI(i,2,3)*kz
              K33=RI(i,3,1)*kx+RI(i,3,2)*ky+RI(i,3,3)*kz
              do 642 j=1,nkI
              if(dabs(K11-KI(1,j)).le.eps)then 
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
              if(E(k1,n).lt.Ef)Nel=Nel+1.0
442           continue
883           continue
              Nel=Nel/Ntot

              do 521 i=1,ntot
521           write(887,*)ktot(1,i),ktot(2,i)
              

             epsq=sqrt(4.0*pi*c0/(Vcell*Ntot)) 

    
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
            do 776  n=1,3
            G(n,iG)=zero
            do 777  m=1,3
            G(n,iG)=G(n,iG)+KC(n,m)*dble(Gi(m))
777         continue
776         continue 
            parG(iG)=Gi(3)  
25          continue
100         FORMAT(I10,I11,I11)
            CLOSE(1) 


           
c            Reciprocal vectors for crystal local field effects calculations in array ''Glf(3,Nlf)'' 


             Nlf=0
             if(lf.eq.1)then
             do 813 iG=1,NG 
             if(G(1,iG).eq.0.0.and.G(2,iG).eq.0.0)then
             Eref=Gcar*Gcar*G(3,iG)*G(3,iG)/2.0
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
             do 812 iG=1,NG 
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



c             GENERIRANJE male BZ OKO 'K' TOCKE  
c             K-point in Cartesi cordintes in  alatt 0.333333  0.577350  0.000000
              KpointX=0.333333
              KpointY=0.577350
              KpointZ=0.000000
              NKBZ=0 
              

c            generiranje 'male' BZ 

             jk=0
             NKBZ=0
             do 103 i=1,nsim
             do 104 ik=1,NkI 
             if(abs(kI(2,ik)).le.qmax)then
             it=1 
             jk=jk+1
             do 3045 n=1,3
             k(n,jk)=zero
             do 3046 m=1,3 
             k(n,jk)=k(n,jk)+R(i,n,m)*kI(m,ik)
3046         continue                    
3045         continue           
             if(jk.gt.1)then   
             do 107 lk=1,jk-1
             if(abs(k(1,jk)-k(1,lk)).le.eps)then 
             if(abs(k(2,jk)-k(2,lk)).le.eps)then      
             if(abs(k(3,jk)-k(3,lk)).le.eps)then    
             it=2           
             endif
             endif
             endif
107          continue
             endif 
             if(it.eq.1)then
             NKBZ=NKBZ+1
             KBZ(1,NKBZ)=k(1,jk)
             KBZ(2,NKBZ)=k(2,jk)
             KBZ(3,NKBZ)=k(3,jk)
             endif
             endif 
104          continue  
103          continue

c            Pomak 'male BZ' u K tocku  


              do 709 i=1,NKBZ
              KBZ(1,i)=KBZ(1,i)+KpointX
              KBZ(2,i)=KBZ(2,i)+Kpointy
709           continue

c            PROVJERA DA LI IMA ISTIH TOCAKA u OKLICI K tocke 
             do 5558 i=1,NKBZ 
             kx=KBZ(1,i)
             ky=KBZ(2,i)
             kz=KBZ(3,i) 
             do 5559 j=1,NKBZ 
             if(j.ne.i)then 
             if(abs(kx-KBZ(1,j)).le.eps)then
             if(abs(ky-KBZ(2,j)).le.eps)then
             if(abs(kz-KBZ(3,j)).le.eps)then
             Print*,'Postoje jednake k tocke u okolini K tocke'  
             stop 
             endif
             endif
             endif 
             endif
5559         continue 
5558         continue 

  


              do 5521 i=1,NKBZ
5521          write(889,*)KBZ(1,i),KBZ(2,i)
             
              
             
              


c             Info file

           open(55,file='Info')
           write(55,*)'***************General***********************'
           write(55,*)''
           write(55,*)'Number of point symmetry operation is',nsim	
           if(lf.eq.1)write(55,*)'Local field effcts in z-dir'
           if(lf.eq.3)write(55,*)'Local field in all xyz-dir'
           write(55,*)'Number of local field vectors is',Nlf
           write(55,*)'Number of different K vectors in 1.B.Z. is',Ntot
           write(55,*)'Number of K vectors in I.B.Z. is',NkI
           write(55,*)'Number of K points around K-point is',NKBZ
           write(55,*)'Number of bands is               ',Nband
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

           
              if(mod.eq.2)goto 8889



c             Upis time ordered zasjenjene kulonske interakcije W_GG'^T(Q,\omega)                
              do 5567 iq=2,nq
              Q=(iq-1)*qmin
              vQ=(2.0*pi)/Q
              dato='W_Qi'
              nord=INDEX(dato,'i', BACK =.false.)
              if(iq.lt.10)then
              write(dato(nord:nord),'(i1)')iq
              elseif(iq.ge.10.and.iq.lt.100)then
              write(dato(nord:nord+1),'(i2)')iq
              else
              write(dato(nord:nord+2),'(i3)')iq
              endif
              path=TRIM(root1)//TRIM(dato)

              OPEN(74,file=path)
              read(74,*)nis
              read(74,44)((WT(iq,iG,jG),jG=1,Nlf),iG=1,Nlf)
              CLOSE(74)
44            FORMAT(10F15.5) 

c             trazenje alpha 
              if(iq.eq.2)alpha=(c0*vQ/real(WT(2,1,1))-1.0)/qmin 
                          

5567          continue


           
                                 
                 
c             valni vektor K 
              do 801 iK1=1,NKBZ

              print*,ik1

              kx=KBZ(1,iK1)
              ky=KBZ(2,iK1)
              kz=KBZ(3,iK1)
 


c              trazenje K prvo u 1.B.Z a onda u I.B.Z.
               it=1
               do 7001 iG=1,NG
                  do 7002 jK=1,Ntot
                     if(dabs(kx-G(1,iG)-ktot(1,jK)).le.eps)then 
                     if(dabs(ky-G(2,iG)-ktot(2,jK)).le.eps)then 
                     if(dabs(kz-G(3,iG)-ktot(3,jK)).le.eps)then 
                     it=2
                     iG01=iG
                     do 5001 i=1,nsim 
                     K11=RI(i,1,1)*ktot(1,jK)+RI(i,1,2)*ktot(2,jK)+
     &               RI(i,1,3)*ktot(3,jK)
                     K22=RI(i,2,1)*ktot(1,jK)+RI(i,2,2)*ktot(2,jK)+
     &               RI(i,2,3)*ktot(3,jK)
                     K33=RI(i,3,1)*ktot(1,jK)+RI(i,3,2)*ktot(2,jK)+
     &               RI(i,3,3)*ktot(3,jK)
                     do 5002 j=1,nkI
                     if(dabs(K11-KI(1,j)).le.eps)then 
                     if(dabs(K22-KI(2,j)).le.eps)then 
                     if(dabs(K33-KI(3,j)).le.eps)then 
                     it=3
                     R1=i
                     K1=j  
                     goto 2001 
                     endif
                     endif
                     endif
5002            continue
5001            continue
                     endif
                     endif
                     endif
7002           continue
7001           continue      
2001                 continue           


                if(it.eq.1)then
                print*,'Can not find wave vector K=',iK1,
     &          'in 1.B.Z.'  
                stop 
                elseif(it.eq.2)then 
                print*,'Can not find wave vector K=',iK1,
     &          'in I.B.Z.'  
                stop
                endif
      

        
 
c             valni vektor K' 
              do 8001 iK2=1,NKBZ

                    
               kx=KBZ(1,iK2)
               ky=KBZ(2,iK2)
               kz=KBZ(3,iK2)
 

               it=1
c              trazenje K' prvo u 1.B.Z a onda u I.B.Z.

               do 701 iG=1,NG
                  do 702 jK=1,Ntot
                     if(dabs(kx-G(1,iG)-ktot(1,jK)).le.eps)then 
                     if(dabs(ky-G(2,iG)-ktot(2,jK)).le.eps)then 
                     if(dabs(kz-G(3,iG)-ktot(3,jK)).le.eps)then 
                     it=2
                     iG02=iG
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
                     R2=i
                     K2=j  
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
                print*,'Can not find wave vector K=',iK2,
     &          'in 1.B.Z.'  
                stop 
                elseif(it.eq.2)then 
                print*,'Can not find wave vector K=',iK2,
     &          'in I.B.Z.'  
                stop
                endif
      


                Qx=KBZ(1,iK2)-KBZ(1,iK1)
                Qy=KBZ(2,iK2)-KBZ(2,iK1)
                Q=Gcar*sqrt(Qx*Qx+Qy*Qy) 
                iq=Q/qmin+1 
           


                if(iq.gt.nq)goto 9999  

              
c              ''K=G01+R1*K1''.
c              ''K'=G02+R2*K2''.
c              K1 i K2 su integeri koji predstavljaju valne vektore u IBZ pridruzene 
c              valnim vektorima K i K' 
c              iG01 i iG02 su integeri rec vec. G koji translatira K i K' u 1BZ    

c             vrpca n 
              do 804  n=nval,nval+1 
       
                      if(n.eq.nval)ni=1          
                      if(n.eq.nval+1)ni=2   
c             vrpca m 
              do 8004 m=n,n

                      

             call paths(root,K1,K2,n,m,pathK1,pathK2,bandn,bandm)

    

               
c         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za 
c         fiksni K1,K2,n i m 

c               Otvaranje atribute za INFO
                call iotk_open_read(10,pathK1)
                call iotk_scan_empty(10,"INFO",attr=attr)
                call iotk_scan_attr(attr,"igwx",NG1)
c               Alociranje polja C1
                allocate (C1(NG1))
c               Ucitavanje podataka iza evc.n
                call iotk_scan_dat(10,bandn,C1)
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

c                Konstrukcija stupca matricnih elementa MnmK1K2(ni,mi,G)  
          

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
                  MnmK1K2(ni,iG)=czero
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
                  K11=K11+G(1,iG02)-G(1,iG01)
                  K22=K22+G(2,iG02)-G(2,iG01)
                  K33=K33+G(3,iG02)-G(3,iG01)
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
                  MnmK1K2(ni,iG)=MnmK1K2(ni,iG)+
     &            conjg(C1(iG1))*C2(iG2)
                  endif
552               continue
551               continue
                  jump=2             

                 
               	deallocate(C1)
 		deallocate(C2)

c               end of m loop
8004            continue


c               end of n loop  
804             continue
                 jump=1  

                
                FockK(ik1,ik2)=czero  

                if(ik1.ne.ik2)then  
                do 513 iG=1,Nlf                   
                do 554 jG=1,Nlf 
                W=WT(iq,iG,jG)
                FockK(ik1,ik2)=FockK(ik1,ik2)-
     &          W*conjg(MnmK1K2(1,iG))*MnmK1K2(2,jG)                
554             continue
513             continue

c               komponenta Q=0 
                elseif(ik1.eq.ik2)then 
                FockK(ik1,ik2)=-Vcell*Ntot*log(1.0+alpha*epsq)/alpha
                endif 
7768            continue 
7769            continue 


9999            continue             
c               end of K'  loop
8001            continue 
c               end of K loop
801             continue



c              upisivanje Fock kernela i DeltaE  
               open(234,file='FockK')
               write(234,44)((FockK(ik1,ik2),ik2=1,NKBZ),ik1=1,NKBZ)  
               close(234)
               goto 9981 



8889           continue 

c*********************************************
c              Ovdje pocinje mod=2 
c*********************************************
              

               print*,'Reading Fock kernel'

c              citanje Fock kernela i DeltaE    
               open(234,file='FockK')
               read(234,44)((FockK(ik1,ik2),ik2=1,NKBZ),ik1=1,NKBZ)  
               close(234)
     

c               GENERIRANJE STRUJNIH VRHOVA jmu                

                do 7771 pol=1,2
               
                if(pol.eq.1)print*,'Calculation of j_x'     
                if(pol.eq.2)print*,'Calculation of j_z'  

                
             
c               valni vektor K 
                do 8101 iK1=1,NKBZ
                 
                kx=KBZ(1,iK1)
                ky=KBZ(2,iK1)
                kz=KBZ(3,iK1)
 
c               trazenje K prvo u 1.B.Z a onda u I.B.Z.
                it=1
                do 7101 iG=1,NG
                  do 7102 jK=1,Ntot
                     if(dabs(kx-G(1,iG)-ktot(1,jK)).le.eps)then 
                     if(dabs(ky-G(2,iG)-ktot(2,jK)).le.eps)then 
                     if(dabs(kz-G(3,iG)-ktot(3,jK)).le.eps)then 
                     it=2
                     iG01=iG
                     do 5101 i=1,nsim 
                     K11=RI(i,1,1)*ktot(1,jK)+RI(i,1,2)*ktot(2,jK)+
     &               RI(i,1,3)*ktot(3,jK)
                     K22=RI(i,2,1)*ktot(1,jK)+RI(i,2,2)*ktot(2,jK)+
     &               RI(i,2,3)*ktot(3,jK)
                     K33=RI(i,3,1)*ktot(1,jK)+RI(i,3,2)*ktot(2,jK)+
     &               RI(i,3,3)*ktot(3,jK)
                     do 5102 j=1,nkI
                     if(dabs(K11-KI(1,j)).le.eps)then 
                     if(dabs(K22-KI(2,j)).le.eps)then 
                     if(dabs(K33-KI(3,j)).le.eps)then 
                     it=3
                     R1=i
                     K1=j  
                     goto 2101 
                     endif
                     endif
                     endif
5102            continue
5101            continue
                     endif
                     endif
                     endif
7102           continue
7101           continue      
2101                 continue           


                if(it.eq.1)then
                print*,'Can not find wave vector K=',iK1,
     &          'in 1.B.Z.'  
                stop 
                elseif(it.eq.2)then 
                print*,'Can not find wave vector K=',iK1,
     &          'in I.B.Z.'  
                stop
                endif

c             valni vektor K' 
              do 8701 iK2=iK1,ik1

               kx=KBZ(1,iK2)
               ky=KBZ(2,iK2)
               kz=KBZ(3,iK2)
 

               it=1
c              trazenje K' prvo u 1.B.Z a onda u I.B.Z.

               do 7301 iG=1,NG
                  do 7302 jK=1,Ntot
                     if(dabs(kx-G(1,iG)-ktot(1,jK)).le.eps)then 
                     if(dabs(ky-G(2,iG)-ktot(2,jK)).le.eps)then 
                     if(dabs(kz-G(3,iG)-ktot(3,jK)).le.eps)then 
                     it=2
                     iG02=iG
                     do 5301 i=1,nsim 
                     K11=RI(i,1,1)*ktot(1,jK)+RI(i,1,2)*ktot(2,jK)+
     &               RI(i,1,3)*ktot(3,jK)
                     K22=RI(i,2,1)*ktot(1,jK)+RI(i,2,2)*ktot(2,jK)+
     &               RI(i,2,3)*ktot(3,jK)
                     K33=RI(i,3,1)*ktot(1,jK)+RI(i,3,2)*ktot(2,jK)+
     &               RI(i,3,3)*ktot(3,jK)
                     do 5302 j=1,nkI
                     if(dabs(K11-KI(1,j)).le.eps)then 
                     if(dabs(K22-KI(2,j)).le.eps)then 
                     if(dabs(K33-KI(3,j)).le.eps)then 
                     it=3
                     R2=i
                     K2=j  
                     goto 2411 
                     endif
                     endif
                     endif
5302            continue
5301            continue
                     endif
                     endif
                     endif
7302           continue
7301           continue      
2411           continue           


                if(it.eq.1)then
                print*,'Can not find wave vector K=',iK2,
     &          'in 1.B.Z.'  
                stop 
                elseif(it.eq.2)then 
                print*,'Can not find wave vector K=',iK2,
     &          'in I.B.Z.'  
                stop
                endif
      
           
              
c              ''K=G01+R1*K1''.
c              ''K'=G02+R2*K2''.
c              K1 i K2 su integeri koji predstavljaju valne vektore u IBZ pridruzene 
c              valnim vektorima K i K' 
c              iG01 i iG02 su integeri rec vec. G koji translatira K i K' u 1BZ    




               do 8404 n=nval,nval        
               do 8405 m=nval+1,nval+1      
      
      
    
            call paths(root,K1,K2,n,m,pathK1,pathK2,bandn,bandm)


               
c         u ovom dijelu programa se iscitava iz binarnih fileova ''gvectors.dat'',''evc.dat'' za 
c         fiksni K1,K2,n i m 

c               Otvaranje atribute za INFO
                call iotk_open_read(10,pathK1)
                call iotk_scan_empty(10,"INFO",attr=attr)
                call iotk_scan_attr(attr,"igwx",NG1)
c               Alociranje polja C1
                allocate (C1(NG1))
c               Ucitavanje podataka iza evc.n
                call iotk_scan_dat(10,bandn,C1)
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
                   

c              matrix elements
               iGfast=0                 
               do 5551 iG=1,1
               current(pol,iK1)=czero
               do 5552 iG1=1,NGd  
               iGfast=iGfast+1
               Gxx1=G(1,iG1)
               Gyy1=G(2,iG1)
               Gzz1=G(3,iG1)
               K11=R(R1,1,1)*Gxx1+R(R1,1,2)*Gyy1+R(R1,1,3)*Gzz1
               K22=R(R1,2,1)*Gxx1+R(R1,2,2)*Gyy1+R(R1,2,3)*Gzz1
               K33=R(R1,3,1)*Gxx1+R(R1,3,2)*Gyy1+R(R1,3,3)*Gzz1
               if(pol.eq.1)then 
               struja=(2.0*kx+Glf(1,iG)+2.0*K11-2.0*G(1,iG01))*Gcar
               elseif(pol.eq.2)then
               struja=(2.0*kz+Glf(3,iG)+2.0*K33-2.0*G(3,iG01))*Gcar
               endif 
               K11=K11+Glf(1,iG)
               K22=K22+Glf(2,iG)
               K33=K33+Glf(3,iG)
               K11=K11+G(1,iG02)-G(1,iG01)
               K22=K22+G(2,iG02)-G(2,iG01)
               K33=K33+G(3,iG02)-G(3,iG01)
               Gxx1=RI(R2,1,1)*K11+RI(R2,1,2)*K22+RI(R2,1,3)*K33
               Gyy1=RI(R2,2,1)*K11+RI(R2,2,2)*K22+RI(R2,2,3)*K33
               Gzz1=RI(R2,3,1)*K11+RI(R2,3,2)*K22+RI(R2,3,3)*K33
               if(jump.eq.1)then 
               do 5553 iG2=1,NG2
               Gfast(iGfast)=NG2+1
               Gxx2=G(1,iG2)
               Gyy2=G(2,iG2)
               Gzz2=G(3,iG2)
               if(dabs(Gxx2-Gxx1).lt.eps)then 
               if(dabs(Gyy2-Gyy1).lt.eps)then 
               if(dabs(Gzz2-Gzz1).lt.eps)then 
               Gfast(iGfast)=iG2
               goto 1011 
               endif
               endif
               endif
5553           continue  
               endif
1011           continue 
               iG2=Gfast(iGfast)
               if(iG2.le.NG2)then
               current(pol,iK1)=current(pol,iK1)+
     &         0.5d0*conjg(C1(iG1))*struja*C2(iG2)
               endif
c              kraj po iG1
5552           continue
c              kraj po C.L.F. iG
5551           continue
               jump=2             
              
                 
              
               if(spin.eq.1)Delta(ik1)=E(K1,nvalLS-1)-E(K1,nvalLS+1)
               if(spin.eq.2)Delta(ik1)=E(K1,nvalLS)-E(K1,nvalLS+2)           
              

c              end of m do loop
8405           continue
               

c              end of n do loop  
8404           continue
               jump=1                 

c              ond of B.Z do loop 
8701           continue
8101           continue

c              kraj po polarizaciji 
7771           continue
           

c               Rijesavanje BSE-FOCK jednadzbe 


                if(spin.eq.1)then  
                open(334,file='Pi_ladder_up_x')
                open(335,file='Pi_ladder_up_z')
                elseif(spin.eq.2)then  
                open(334,file='Pi_ladder_down_x')
                open(335,file='Pi_ladder_down_z')
                endif  


                 do  1005  io=1,no
                 o=omin+(io-1)*domega 
                 print*,io

                 do 828 ik1=1,NKBZ
                 do 818 ik2=1,NKBZ
                 Fock(ik1,ik2)=czero
                 unit(ik1,ik2)=czero
818              continue
                 unit(ik1,ik1)=rone
828              continue


c               Konstrukcija matrice L^0*Xi^F
c               Konstrukcija matrice Xi0=L^0*Xi^FL^0

                do 1002 ik1=1,NKBZ  
                L0i=1.0/(o+Delta(ik1)+ione*eta) 
                do 1004 ik2=1,NKBZ  
                L0j=1.0/(o+Delta(ik2)+ione*eta) 
                Fock(ik1,ik2)=unit(ik1,ik2)-
     &          L0i*FockK(ik1,ik2)/(Ntot*Vcell)
                Xi0(ik1,ik2)=L0i*FockK(ik1,ik2)*L0j
1004            continue
1002            continue
 
             
                
                call gjel(Fock,NKBZ,nkd,unit,NKBZ,nkd)
                

c               konstrukcija 4-point polarizabilnosti L^ladd



                do 6002 ik1=1,NKBZ  
                do 6004 ik2=1,NKBZ  
                Lladd(ik1,ik2)=czero
                do 9001 ik=1,NKBZ
9001            Lladd(ik1,ik2)=Lladd(ik1,ik2)+Fock(ik1,ik)*Xi0(ik,ik2)
6004            continue
6002            continue



c               Konstrukcija Ladder ireduciblne current-current polarisabilnosti chi_ladd_\mu\mu

       

                do 5501 pol=1,2
                chi_ladd(pol)=czero 
                chi_para(pol)=czero 
                do 6302 ik1=1,NKBZ  
                L0i=(o/Delta(ik1))*(1.0/(o+Delta(ik1)+ione*eta)) 
                do 6304 ik2=1,NKBZ  
                j1=current(pol,iK1)
                j2=conjg(current(pol,iK2))
                chi_ladd(pol)=chi_ladd(pol)-
     &          valley*j1*Lladd(ik1,ik2)*j2/(Ntot*Vcell*Ntot*Vcell)
6304            continue
c               Paramagnetska struja-struja polarizabilost  
                chi_para(pol)= chi_para(pol)+
     &          valley*j1*L0i*conjg(j1)/(Ntot*Vcell)
6302            continue
5501            continue
               
                write(334,*)o*Hartree,chi_ladd(1)
                write(335,*)o*Hartree,chi_ladd(2)

                write(104,*)o*hartree,imag(chi_para(1)+chi_ladd(1))
                write(105,*)o*hartree,imag(chi_para(1))
                write(106,*)o*hartree,real(chi_para(1)+chi_ladd(1))
                write(107,*)o*hartree,real(chi_para(1))

                write(204,*)o*hartree,imag(chi_para(2)+chi_ladd(2))
                write(205,*)o*hartree,imag(chi_para(2))
                write(206,*)o*hartree,real(chi_para(2)+chi_ladd(2))
                write(207,*)o*hartree,real(chi_para(2))




c               omega do loop 
1005            continue                 



                close(334)
                close(335)
                close(336)



9981            continue



            

                   end 

