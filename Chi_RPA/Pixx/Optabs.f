         program Surface_current 

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
c        Vcell-unit-cell volume in a.u.^-3 
c        a0-unit cell parameter in  a.u.^-1  

        
         use iotk_module
         implicitnone      
         character(iotk_attlenx) :: attr
         logical :: found

         INTEGER NkI,Nband,ik,i,Nk,j,jk,it,lk,Ntot,iG0,nsim,
     &   NG,io,no,iq,nq,nMPx,nMPy,nMPz,n,m,iG,R1,K1,R2,K2,
     &   Nlf,NG1,NG2,NGd,iG1,iG2,jG,Nlfd,kG,jo,mod,jump,
     &   iGfast,ikmin,NelQE,pol,frac,lf,nocc
         PARAMETER(nMPx=51,nMPy=51,nMPz=1,NkI=6835,Nband=60,nocc=9,
     &   NelQE=18,Nk=48*NkI,NGd=4000,NG=8000,no=5001,nq=1,Nlfd=20)
    

c        skalars
         REAL*8 a0,eps,kx,ky,kz,KQx,KQy,KQz,qGx,qGy,qGz,omin,omax,qx,
     &   qy,qz,kmin,domega,o,Ef,K11,K22,K33,T,f1,f2,Lor,De,Gabs,
     &   kref,Eref,Ecut,Gxx1,Gyy1,Gzz1,Gxx2,Gyy2,Gzz2,fact,Vcell,
     &   oi,oj,Q,Gcar,aBohr,zero,Nel,error,Lossf,struja,Pabs,
     &   Gama_intra,Gama_inter,expo1,expo2,gama,Cond,c0,strujay,
     &   strujaz,ImChi0,ReChi0,DGW
         DOUBLEPRECISION pi,three,Hartree,planck,eV
         PARAMETER(Hartree=2.0d0*13.6056923d0,EF=0.5554/Hartree,
     &   a0=5.9715,c0=29.8575,pi=3.141592654d0,Gcar=2.0*pi/a0,
     &   eps=1.0d-4,T=0.025/Hartree,Gama_intra=0.025/Hartree,
     &   Gama_inter=0.025/Hartree,eV=1.602176487d-19,gama=1.0/137.0,
     &   planck=6.626196d-34,Ecut=0.0,Vcell=922.0586,
     &   aBohr=0.5291772d0,zero=0.0)
         DOUBLE COMPLEX ione,czero,rone,eM

c        arrays
         INTEGER Gfast,Gi
         REAL*8 kI,E,R,RI,k,ktot,G,Glf,V,KC,f,Kk
         DOUBLE COMPLEX unit,epsilon
         COMPLEX*8 MnmK1K21,Pi_dia,Pi_tot,MnmK1K22,Qeff,S0
         DIMENSION kI(3,NkI),E(NkI,Nband),R(48,3,3),RI(48,3,3),
     &   k(3,Nk),ktot(3,Nk),V(Nlfd,Nlfd),G(3,NG),
     &   Glf(3,Nlfd),MnmK1K21(Nlfd),unit(Nlfd,Nlfd),
     &   epsilon(Nlfd,Nlfd),S0(-no:no,Nlfd,Nlfd),Gfast(Nlfd*NGd),
     &   KC(3,3),Gi(3),f(48,3),MnmK1K22(Nlfd),Pi_dia(Nlfd,Nlfd),
     &   Qeff(Nlfd,Nlfd),Pi_tot(Nlfd,Nlfd),Kk(3,6) 
 
         CHARACTER*100 bandn,bandm,nis,pathK1,pathK2,dato1,
     &   root,path,fajl,dato2,dato3,root1,root2
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



c        QUANTUM ESSPRESSO IMPUTS:
         root1='/home/vito/PROJECTS/MoS2-BSE'
         root2='/MoS2_201X201/'
         root=TRIM(root1)//TRIM(root2)


        

c             For free spectral function calculation put mod=1
c             For current-current response function calculation put mod=2
c             For x polarization put pol=1
c             For y polarization put pol=2
c             For z polarization put pol=3
c             For mixed component yz put pol=4
c             Crystal local field effects are included in z direction lf=1 
c             Crystal local field effects are included in x,y,z direction lf=3 
            
              mod=2
              lf=1
              pol=1
            
c             CORRELATION FUNCTIONS, CURRENT-CURRENT RESPONSE FUNCTIONS and 
c             EFFECTIVE CHARGE CARRIERS MATRIX OUTPUTS
              
              dato1='Corrfun'
              dato2='Qeff'
              if(pol.eq.1)dato3='Pi_RPA_xx'
              if(pol.eq.2)dato3='Pi_RPA_yy'
              if(pol.eq.3)dato3='Pi_RPA_zz'

c             scissors shift
              DGW=1.0/Hartree



              jump=1
              three=3.0d0
              omin=1.0d-5
              omax=(50.0/Hartree+omin)
              domega=(omax-omin)/(no-1)
         
              
            


c      CALL FOR POINT GROUP TRANSFORMATIONS


            call PointR(root,nsim,R,RI,f)

            

c           Upis valnih vektora iz irreducibilne Brillouinove 
c           zone i pripadnih energijskih nivoa iz filea '*****.band'. 
c           wave vectors are in cart.koord.

            fajl='/MoS2.band'
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


c          

            do 333 ik=1,NkI 
            do 444 i=1,Nband 
            E(ik,i)=E(ik,i)/Hartree
            if(i.ge.nocc+1)E(ik,i)=E(ik,i)+DGW
444         continue
333         continue


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
              Nel=2.0*Nel/Ntot

c             Checking the existence of fraction translation operations  
              frac=0 
              do 120, i=1,nsim   
              if(f(i,1).ne.zero.or.f(i,2).ne.
     &        zero.or.f(i,3).ne.zero)then 
              frac=1
              endif            
120           continue          


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
             endif
812          continue
             endif              
             if(Nlf.gt.Nlfd)then
             print*,'Nlf is bigger than Nlfd'
             stop
             endif
                   
                        
 
c             q LOOP STARTS HERE!!! 

              do 801 iq=nq,nq

c             searching for the smalest 'optical' q               
              kmin=1.0
              do 913 i=1,Ntot
              kref=sqrt(ktot(1,i)*ktot(1,i)+
     &        ktot(2,i)*ktot(2,i)+ktot(3,i)*ktot(3,i)) 
              if(kref.eq.zero)goto 970
              if(kref.lt.kmin)then
              kmin=kref      
              ikmin=i
              endif
970           continue
913           continue
             
              qx=(iq-1)*ktot(1,ikmin)
              qy=(iq-1)*ktot(2,ikmin)
              qz=(iq-1)*ktot(3,ikmin)
    
c             Info file

       open(55,file='Info')
       write(55,*)'***************General***********************'
       write(55,*)' Currently we calculate         ---->',dato1
       write(55,*)' Currently we calculate         ---->',dato2
       write(55,*)''
       write(55,*)'Number of point symmetry operation is',nsim	
       if(frac.eq.0)write(55,*)'Fraction translation is not detected'	
       if(frac.eq.1)write(55,*)'Fraction translation is detected'	
       write(55,88)'Wave vector (qx,qy,qz)=(',qx*Gcar,qy*Gcar,
     & qz*Gcar,') a.u.'	                
       write(55,99)'|(qx,qy,qz)|=',sqrt(qx*qx+qy*qy+qz*qz)*Gcar,'a.u.'
       if(lf.eq.1)write(55,*)'Local field effcts in z-dir'
       if(lf.eq.3)write(55,*)'Local field in all xyz-dir'
       write(55,*)'Number of local field vectors is',Nlf
       write(55,*)'Number of different K vectors in 1.B.Z. is',Ntot
       write(55,*)'Number of K vectors in I.B.Z. is',NkI
       write(55,*)'Number of bands is               ',Nband
       write(55,99)'Gama_intra is  ',gama_intra*Hartree*1000.0,'meV'
       write(55,99)'Gama_inter is  ',gama_inter*Hartree*1000.0,'meV'  
       write(55,99)'Temperature is      ',T*Hartree*1000.0,'meV'
       write(55,*)''
       write(55,*)'-Im(Chi(io,G1,G2))/pi is in file---->',dato1
       write(55,*)' Qeff complex matrix is in file ---->',dato2
       write(55,*)' Pi_munu is in file            ---->',dato3
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
88     FORMAT(A25,3F10.4,A5)     
99     FORMAT(A25,F7.3,A5)
12     FORMAT(A40,F8.4)          


               if(mod.eq.2)goto 888



               S0=czero
               Qeff=czero

               open(74,file=dato1)
               open(75,file=dato2)               


c              1.B.Z  LOOP STARTS HERE !!!!                   

               do 803 ik=1,Ntot

               open(33,file='status')
               write(33,*)ik
               close(33)
c                print*,ik

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
               R1=i
               K1=j  
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

c              trazenje (KQx,KQy) prvo u 1.B.Z a onda u I.B.Z.

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

               

               do 804 n=1,nband          
               do 805 m=1,nband      
 


                  expo1=exp((E(K1,n)-EF)/T)
                  expo2=exp((E(K2,m)-EF)/T)
                  f1=1.0/(expo1+1.0)
                  f2=1.0/(expo2+1.0)
                  f1=f1-f2 

                  if((dabs(f1).ge.1.0d-3).or.(n.eq.m))then 

    
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
                   

c                 matrix elements
                  iGfast=0                 
                  do 551 iG=1,Nlf 
                  MnmK1K21(iG)=czero
                  MnmK1K22(iG)=czero
                  do 552 iG1=1,NGd  
                  iGfast=iGfast+1
                  Gxx1=G(1,iG1)
                  Gyy1=G(2,iG1)
                  Gzz1=G(3,iG1)
                  K11=R(R1,1,1)*Gxx1+R(R1,1,2)*Gyy1+R(R1,1,3)*Gzz1
                  K22=R(R1,2,1)*Gxx1+R(R1,2,2)*Gyy1+R(R1,2,3)*Gzz1
                  K33=R(R1,3,1)*Gxx1+R(R1,3,2)*Gyy1+R(R1,3,3)*Gzz1
                if(pol.eq.1)struja=(qx+2.0*kx+Glf(1,iG)+2.0*K11)*Gcar
                if(pol.eq.2)struja=(qy+2.0*ky+Glf(2,iG)+2.0*K22)*Gcar
                if(pol.eq.3)struja=(qz+2.0*kz+Glf(3,iG)+2.0*K33)*Gcar
                if(pol.eq.4)then 
                strujay=(qy+2.0*ky+Glf(2,iG)+2.0*K22)*Gcar
                strujaz=(qz+2.0*kz+Glf(3,iG)+2.0*K33)*Gcar
                endif
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
                  if(pol.ne.4)then  
                  MnmK1K21(iG)=MnmK1K21(iG)+
     &            0.5d0*conjg(C1(iG1))*struja*C2(iG2)
                  MnmK1K22(iG)=MnmK1K21(iG)
                  else  
                  MnmK1K21(iG)=MnmK1K21(iG)+
     &            0.5d0*conjg(C1(iG1))*strujay*C2(iG2)
                  MnmK1K22(iG)=MnmK1K22(iG)+
     &            0.5d0*conjg(C1(iG1))*strujaz*C2(iG2)
                  endif  
                  endif
c                 kraj po iG1
552               continue
c                 kraj po C.L.F. iG
551               continue
                  jump=2             
              
                    
                if(n.ne.m)then 
c                  omega loop 
                   do 802 io=-no,no
                   o=io*domega
                   De=o+E(K1,n)-E(K2,m)
                   Lor=Gama_inter/(De*De+Gama_inter*Gama_inter)
                   if(dabs(Lor).ge.1.0d-5/Gama_inter)then  
                     do 513 iG=1,Nlf                   
                     do 554 jG=1,Nlf 
c                    -1/pi*ImChi_munu-for Kramers Kronig
                     S0(io,iG,jG)=S0(io,iG,jG)-
     &               2.0*f1*Lor*MnmK1K21(iG)*conjg(MnmK1K22(jG))/
     &               (pi*Ntot*Vcell)                
554                  continue
513                  continue
                   endif   
802                continue
                elseif(n.eq.m.and.dabs(E(K1,n)-EF).le.10.0*T)then    
c                  Effective number of charge carriers (tensor)
                     fact=expo1/((expo1+1.0d0)*(expo1+1.0d0))
                     fact=-fact/T
                     do 533 iG=1,Nlf                   
                     do 574 jG=1,Nlf 
                     Qeff(iG,jG)=Qeff(iG,jG)
     &               +2.0*fact*MnmK1K21(iG)*conjg(MnmK1K22(jG))
     &               /(Ntot*Vcell)                
574                  continue                   
533                  continue 
c               end of intra/interband if loop 
                endif
    
        	deallocate(C1)
		deallocate(C2)

                endif
c               end of the occupation if loop 
      
              

c              end of m do loop
805            continue
               

c              end of n do loop  
804            continue
               jump=1                 
             
834            continue
c              ond of 1.B.Z do loop 
803            continue


c              WRITING CORFUN S0_\mu\nu 
c              WRITTING Q_eff_\mu\nu
               do 6678 io=-no,no
               o=io*domega
               write(74,*)'omega=',o,'Hartree'
               write(74,44)((S0(io,iG,jG),jG=1,Nlf),iG=1,Nlf)
6678           continue  
               write(75,44)((Qeff(iG,jG),jG=1,Nlf),iG=1,Nlf)
44             FORMAT(10F15.10) 
               CLOSE(74)
               CLOSE(75)
                 

                 if(mod.eq.1)goto 999

c               SECOND PART OF THE PROGRAM mod=2
c               Calculation of the matrix '''Pi_\mu\nu'' by using matrix ''S0_\mu\nu(G,G')'' and Kramers-Krroning relations 
                 
888             continue

                open(74,file=dato1)
                do 523 io=-no,no
                read(74,*)nis
                read(74,44)((S0(io,iG,jG),jG=1,Nlf),iG=1,Nlf)
523             continue
                CLOSE(74)

                open(75,file=dato2)
                read(75,44)((Qeff(iG,jG),jG=1,Nlf),iG=1,Nlf)
503             continue
                CLOSE(75)




c               Puting (qx,qy,qz) and Glf in cartezi coordinate 

                qx=Gcar*qx 
                qy=Gcar*qy
                qz=Gcar*qz


                do 703 iG=1,Nlf 
                Glf(1,iG)=Gcar*Glf(1,iG)
                Glf(2,iG)=Gcar*Glf(2,iG)
                Glf(3,iG)=Gcar*Glf(3,iG)
703             continue




                 open(77,file=dato3)     
c               new sum over omega
                 do 872 io=1,no-1
                print*,io            
                 oi=(io-1)*domega
                 do 432 iG=1,Nlf   
                 do 433 jG=1,Nlf 
c                Real part of the response function Re(Chi)
                 ReChi0=0.0
c                  static limit
                   if(io.eq.1)then
                      do 431 jo=2,no
                      oj=(jo-1)*domega
                      fact=domega/oj
                      if(jo.eq.2)fact=3.0/2.0 
                      if(jo.eq.no)fact=0.5*domega/oj
                      ReChi0=ReChi0+
     &                fact*(real(S0(-jo+1,iG,jG))-real(S0(jo-1,iG,jG)))
431                   continue
                   elseif(io.eq.2)then 
                      do 434 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=1.0
                      if(jo.eq.2)fact=0.0
                      if(jo.eq.3)fact=-3.0/2.0 
                      if(jo.eq.no)fact=0.5*domega/(oi-oj)
                      ReChi0=ReChi0+fact*real(S0(jo-1,iG,jG))
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      ReChi0=ReChi0+fact*real(S0(-jo+1,iG,jG))
434                   continue
                   elseif(io.eq.(no-1))then  
                      do 435 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=0.5*domega/(oi-oj)
                      if(jo.eq.(no-2))fact=3.0/2.0
                      if(jo.eq.(no-1))fact=0.0
                      if(jo.eq.no)fact=-1.0
                      ReChi0=ReChi0+fact*real(S0(jo-1,iG,jG))
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      ReChi0=ReChi0+fact*real(S0(-jo+1,iG,jG))
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
                      ReChi0=ReChi0+fact*real(S0(jo-1,iG,jG))
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      ReChi0=ReChi0+fact*real(S0(-jo+1,iG,jG))
436                   continue
                   endif 
                      ReChi0=ReChi0+pi*imag(S0(io-1,iG,jG))
                   
c                Imaginary part of the response function Im(Chi)
                   ImChi0=0.0
c                  static limit
                   if(io.eq.1)then
                      do 531 jo=2,no
                      oj=(jo-1)*domega
                      fact=domega/oj
                      if(jo.eq.2)fact=3.0/2.0 
                      if(jo.eq.no)fact=0.5*domega/oj
                      ImChi0=ImChi0+
     &                fact*(imag(S0(-jo+1,iG,jG))-imag(S0(jo-1,iG,jG)))
531                   continue
                   elseif(io.eq.2)then 
                      do 534 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=1.0
                      if(jo.eq.2)fact=0.0
                      if(jo.eq.3)fact=-3.0/2.0 
                      if(jo.eq.no)fact=0.5*domega/(oi-oj)
                      ImChi0=ImChi0+fact*imag(S0(jo-1,iG,jG))
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      ImChi0=ImChi0+fact*imag(S0(-jo+1,iG,jG))
534                   continue
                   elseif(io.eq.(no-1))then  
                      do 535 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=0.5*domega/(oi-oj)
                      if(jo.eq.(no-2))fact=3.0/2.0
                      if(jo.eq.(no-1))fact=0.0
                      if(jo.eq.no)fact=-1.0
                      ImChi0=ImChi0+fact*imag(S0(jo-1,iG,jG))
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      ImChi0=ImChi0+fact*imag(S0(-jo+1,iG,jG))
535                   continue
                   else 
                      do 536 jo=1,no
                      oj=(jo-1)*domega
                      if(jo.ne.io)fact=domega/(oi-oj)
                      if(jo.eq.1)fact=0.5*domega/(oi-oj)
                      if(jo.eq.(io-1))fact=3.0/2.0
                      if(jo.eq.io)fact=0.0
                      if(jo.eq.(io+1))fact=-3.0/2.0
                      if(jo.eq.no)fact=0.5*domega/(oi-oj)
                      ImChi0=ImChi0+fact*imag(S0(jo-1,iG,jG))
                      fact=domega/(oi+oj) 
                      if(jo.eq.1.or.jo.eq.no)fact=0.5*domega/(oi+oj) 
                      ImChi0=ImChi0+fact*imag(S0(-jo+1,iG,jG))
536                   continue
                   endif 
                     ImChi0=ImChi0-pi*real(S0(io-1,iG,jG))
            
                   if(io.eq.1)Pi_dia(iG,jG)=-cmplx(ReChi0,zero)
                   Pi_tot(iG,jG)=cmplx(ReChi0,ImChi0)
                   Pi_tot(iG,jG)=Pi_tot(iG,jG)+Pi_dia(iG,jG)

c                   dodavanje intraband clana 

                    Pi_tot(iG,jG)=Pi_tot(iG,jG)
     &              +Qeff(iG,jG)*oi/(oi+ione*Gama_intra)

                      
c                kraj po iG,jG              
433              continue    
432              continue 

c                WRITTING TOTAL RESPONSE FUNCTION Pi_xx 
               
                 write(77,*)oi*Hartree,Pi_tot(1,1)
               
         

c                end new sum over omega
872              continue
                 close(77)
              
999              continue
c                kraj po q 
801              continue
                 end 

