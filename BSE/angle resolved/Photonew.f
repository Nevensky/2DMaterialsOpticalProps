            program FOTON

c       program for the calculation of the current-current corelation function in MoS2
c       XY-plane  unit cell parameter  a0=    4.6510  a.u.  
c       Z-dir unit cell parameter c0=  23.255 a.u.    

 
      
        implicitnone

        integer n,no,iG,jG,kG,NG,Nlfd,Nlf,Nlf2,lf,nord,
     &  nMPx,nMPy,itheta,ntheta,nq


c       Nlfd je minimalno 2 zbog 2X2 blok matrice za p mod    
        parameter(nMPx=201,nMPy=201,no=2001,NG=6000,Nlfd=130)
          
        double precision zero,one,two,three,four,six,pi,eta,eV,
     &  Hartree,aBohr,planck,kb
        REAL*8 a0,c0,Eref,Gcar,gama,Ecut,domega,theta,dtheta,Vcell

        parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,
     &  four=4.0d0,six=6.0d0,pi=3.141592654d0,a0=5.9715,c0=29.8575d0,
     &  kB=1.3806503d-23,eV=1.602176487d-19,Ecut=20.0,Gcar=2.0*pi/a0,
     &  Hartree=2.0d0*13.6056923d0,aBohr=0.5291772d0,Vcell=922.0586,
     &  planck=6.626196d-34,eta=0.000001/Hartree,gama=1.0d0/137.0d0)
        double complex czero,ione,ieta,cone


c     ...scalars...
        integer i,j,k,m,iq,io,iz,jz,qk,Gx,Gy,Gz
        real*8 Ef,kf,q,o,omax,dq,omin,Ff1,Ff2,Gf1,Gf2,As,Trs,Rs,Ap,
     &  Trp,Rp,sem,rhor,z
        double complex oi,beta,Aii,Ajj,Aij,Bij,Dxx0,Dyy0,diagonal1,
     &  diagonal2
        complex*8 Trans,Refs,Absos,Dxxr,Dyyr,Dzzr,Dyzr,Dxxt,Dyyt,
     &  Dzzt,Dyzt,Dxxa,Dyya,Dzza,Dyza,Tranp,Refp,Absop,rhoc

c     ...arrays...
        REAL*8 G,Glf,KC
        INTEGER parG,Gi
        DOUBLECOMPLEX unit,TS,TP
        COMPLEX*8 unit2,Dxx,Dyy,Dzz,Dyz,Dzy,Pixx,Piyy,Piyz,Pizy,
     &  Pizz,rho
        DIMENSION Dxx(Nlfd,Nlfd),Dyy(Nlfd,Nlfd),Dzz(Nlfd,Nlfd),
     &  Dyz(Nlfd,Nlfd),Dzy(Nlfd,Nlfd),Pixx(Nlfd,Nlfd),
     &  Piyy(Nlfd,Nlfd),Piyz(Nlfd,Nlfd),Pizy(Nlfd,Nlfd),
     &  Pizz(Nlfd,Nlfd),unit2(Nlfd,Nlfd),unit(Nlfd,Nlfd),
     &  TP(Nlfd,Nlfd),TS(Nlfd,Nlfd),G(3,NG),Glf(3,Nlfd),
     &  KC(3,3),parG(NG),Gi(3),rho(Nlfd)
        character*100 nis,dato,root,path,fajl
        CHARACTER*35 tag,buffer

    
             czero=dcmplx(zero,zero)
             cone=dcmplx(one,zero)
             ione=dcmplx(zero,one)
             ieta=dcmplx(zero,eta)
          
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


c       definition of the unit reciprocal vectors rading from output datas above in units G0=2pi/a0  
             
c             QUANTUM ESSPRESSO IMPUTS:
              root='/home/vito/PROJECTS/MoS2/MoS2_201X201'

c             Crystal local field effects are included in z direction lf=1 
c             Crystal local field effects are included in x,y,z direction lf=3 
              lf=1
 
C****************************************************+
c             PLASMON-POLARITONS 
c             za slucaj skeniranja evanescent (omega<Qc) + radijativnih (onega>Qc) eigen-modova 
c             obavezno podesiti: 
c             ntheta=1               
c             nq>1

c             LIGHT SCATTERING 
c             za slucaj rasprsnja elektromagnetskog polja u ovisnosti o kutu theta 
c             obavezno podesiti:
c             ntheta>1  
c             nq=1               

              ntheta=1
              nq=21
c**********************************************************************************

             

c             trensfer wave vector step \deltaQ  
              dq=1.0d-6            
              dtheta=pi/(2.0*dble(ntheta))
              omin=1.0d-5
              omax=2.0d0

              domega=(omax-omin)/(no-1)
      


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


c            Puting Glf in cartezi coordinate 

             do 703 iG=1,Nlf 
             Glf(1,iG)=Gcar*Glf(1,iG)
             Glf(2,iG)=Gcar*Glf(2,iG)
             Glf(3,iG)=Gcar*Glf(3,iG)
703          continue

              Nlf2=2*Nlf 


            open(33,file='Pixx')
            open(34,file='Piyy')
            open(35,file='Pizz')
c           reading unscreened current-current response tensor Pi_\mu\nu
            do 939 io=1,no-1
            read(33,*)nis
            read(33,44)((Pixx(iG,jG),jG=1,Nlf),iG=1,Nlf)
            read(34,*)nis
            read(34,44)((Piyy(iG,jG),jG=1,Nlf),iG=1,Nlf)
            read(35,*)nis
            read(35,44)((Pizz(iG,jG),jG=1,Nlf),iG=1,Nlf)
939         continue
44          FORMAT(10F15.10) 
            close(33)
            close(34)
            close(35)



c*******************************************************************
c              POCETAK ITERIRANJA PO theta, Q i OMEGA
c******************************************************************
    
c             angle loop \theta starts here!!! 
              do 801 itheta=1,ntheta
              
              theta=(itheta-1)*dtheta 
                
              write(*,38)'theta=',180.0*theta/pi,'°'
38            FORMAT(A6,F5.1,A5)              
          
c              Q loop starts  here
               do 8 iq=1,nq

               
c              \omega loop starts here 
               do 9 io=1,no-1

               o=omin+(io-1)*domega  

               if(itheta.eq.1)then 
               q=(iq-1)*dq
               else 
               q=gama*o*sin(theta)
               endif


               oi=dcmplx(o,eta) 
               beta=dcmplx(gama*gama*oi*oi-q*q)
               beta=zsqrt(beta)
 
c         matrix elements D_{\mu\nu}(G,G')

           do 9001 iG=1,Nlf 
           Aii=parG(iG)*(cone-exp(ione*beta*c0))/
     &     (beta*beta-Glf(3,iG)*Glf(3,iG))
           do 9002 jG=1,Nlf 
     	   Ajj=Aii*parG(jG)/(beta*beta-Glf(3,jG)*Glf(3,jG))
           Aij=2.0d0*(beta*beta+Glf(3,iG)*Glf(3,jG))*Ajj/c0
           Bij=2.0d0*beta*(Glf(3,iG)+Glf(3,jG))*Ajj/c0
           if(jG.eq.iG)then
           diagonal1=2.0*ione*beta/(beta*beta-Glf(3,iG)*Glf(3,iG))
           diagonal2=2.0*ione*Glf(3,iG)/
     &     (beta*beta-Glf(3,iG)*Glf(3,iG))           
           else 
           diagonal1=czero
           diagonal2=czero 
           endif
c          matricni element Dxx i Dyy
           Dxx(iG,jG)=2.0d0*pi*ione*gama*gama*(Aij+diagonal1)/beta 
           Dyy(iG,jG)=2.0d0*pi*ione*beta*(Aij+diagonal1)/(oi*oi)
c          matricni elementi Dyz i Dzy
           Dyz(iG,jG)=-two*pi*ione*Q*(Bij+diagonal2)/(oi*oi)
           Dzy(iG,jG)=Dyz(iG,jG)            
c          matricni element Dzz
           Dzz(iG,jG)=two*pi*ione*Q*Q*(Aij+diagonal1)/(beta*oi*oi)
9002       continue
           Dzz(iG,iG)=Dzz(iG,iG)-four*pi/(oi*oi)
9001       continue


             
c           S-MOD

            do 118 iG=1,Nlf
            do 128 jG=1,Nlf 
            unit(iG,jG)=czero
            IF(iG.eq.jG)unit(iG,jG)=cone
            TS(iG,jG)=czero             
            do 3458 kG=1,Nlf   
            TS(iG,jG)=TS(iG,jG)+Pixx(iG,kG)*Dxx(kG,jG)
3458        continue
            TS(iG,jG)=unit(iG,jG)-TS(iG,jG)
128         continue
118         continue

c           invertiranje matrice TS
          

            call gjel(TS,Nlf,Nlfd,unit,Nlf,Nlfd)
            
                
c           KONSTRUKCIJA TP MATRICE P-MOD
            


            do 828 iG=1,Nlf2
            do 818 jG=1,Nlf2
            TP(iG,jG)=czero     
            unit2(iG,jG)=czero
818         continue
            unit2(iG,iG)=cone
828         continue
         



            do 11 iG=1,Nlf2
            do 12 jG=1,Nlf2 
            IF(iG.le.Nlf.and.jG.le.Nlf)THEN 
              do 3451 kG=1,Nlf   
              TP(iG,jG)=TP(iG,jG)+Piyy(iG,kG)*Dyy(kG,jG)+
     &        Piyz(iG,kG)*Dzy(kG,jG)
3451          continue
            ELSEIF(iG.le.Nlf.and.jG.gt.Nlf)then
              do 3452 kG=1,Nlf   
              TP(iG,jG)=TP(iG,jG)+Piyy(iG,kG)*Dyz(kG,jG-Nlf)+
     &        Piyz(iG,kG)*Dzz(kG,jG-Nlf)
3452          continue
            ELSEIF(iG.gt.Nlf.and.jG.le.Nlf)THEN
              do 3453 kG=1,Nlf   
              TP(iG,jG)=TP(iG,jG)+Pizy(iG-Nlf,kG)*Dyy(kG,jG)+
     &        Pizz(iG-Nlf,kG)*Dzy(kG,jG)
3453          continue
            ELSE 
              do 3454 kG=1,Nlf   
              TP(iG,jG)=TP(iG,jG)+Pizy(iG-Nlf,kG)*Dyz(kG,jG-Nlf)+
     &        Pizz(iG-Nlf,kG)*Dzz(kG,jG-Nlf)
3454          continue
            ENDIF
              TP(iG,jG)=unit2(iG,jG)-TP(iG,jG)
12          continue
11          continue
           

c          invertiranje matrice TP
         

         

 
           call gjel(TP,Nlf2,Nlfd,unit2,Nlf2,Nlfd)

         

c           SCREENED CURRENT-CURRENT RESPONSE MATRICES  
c           S-MOD CURRENT-CURRENT MATRIX stored in ''Dxx''           

            do 1018 iG=1,Nlf
            do 1028 jG=1,Nlf 
            Dxx(iG,jG)=czero             
            do 3045 kG=1,Nlf   
            Dxx(iG,jG)=Dxx(iG,jG)+TS(iG,kG)*Pixx(kG,jG)
3045        continue
1028        continue
1018        continue


c           P-MOD CURRENT-CURRENT MATRIX stored in ''Dyy, Dyz, Dzy, Dzz''      

          
            do 505 iG=1,Nlf
            do 605 jG=1,Nlf 
            Dyy(iG,jG)=czero             
            do 305 kG=1,Nlf   
            Dyy(iG,jG)=Dyy(iG,jG)+TP(iG,kG)*Piyy(kG,jG)+
     &      TP(iG,kG+Nlf)*Pizy(kG,jG)
305         continue
605         continue
505         continue

            do 506 iG=1,Nlf
            do 606 jG=1,Nlf 
            Dyz(iG,jG)=czero             
            do 306 kG=1,Nlf   
            Dyz(iG,jG)=Dyz(iG,jG)+TP(iG,kG)*Piyz(kG,jG)+
     &      TP(iG,kG+Nlf)*Pizz(kG,jG)
306         continue
606         continue
506         continue

            do 507 iG=1,Nlf
            do 607 jG=1,Nlf 
            Dzy(iG,jG)=czero             
            do 307 kG=1,Nlf   
            Dzy(iG,jG)=Dzy(iG,jG)+TP(iG+Nlf,kG)*Piyy(kG,jG)+
     &      TP(iG+Nlf,kG+Nlf)*Pizy(kG,jG)
307         continue
607         continue
507         continue

            do 508 iG=1,Nlf
            do 608 jG=1,Nlf 
            Dzz(iG,jG)=czero             
            do 308 kG=1,Nlf   
            Dzz(iG,jG)=Dzz(iG,jG)+TP(iG+Nlf,kG)*Piyz(kG,jG)+
     &      TP(iG+Nlf,kG+Nlf)*Pizz(kG,jG)
308         continue
608         continue
508         continue

          
c           calculation of reflected, transmited and absorbed coefficients  
            Dxx0=gama*gama*two*pi*ione/beta 
            Dyy0=gama*two*pi*ione/oi 
            Dxxr=czero
            Dyyr=czero
            Dzzr=czero
            Dyzr=czero
            Dxxt=czero
            Dyyt=czero
            Dzzt=czero
            Dyzt=czero
            Dxxa=czero
            Dyya=czero
            Dzza=czero
            Dyza=czero
            do 808 iG=1,Nlf 
            Ff1=(2.0*parG(iG)/sqrt(c0))*sin(beta*c0/2.0)
     &      /(beta+Glf(3,iG)) 
            Gf1=(2.0*parG(iG)/sqrt(c0))*sin(beta*c0/2.0)
     &      /(beta-Glf(3,iG)) 
            do 909 jG=1,Nlf
            Ff2=(2.0*parG(jG)/sqrt(c0))*sin(beta*c0/2.0)
     &      /(beta+Glf(3,jG))
            Gf2=(2.0*parG(jG)/sqrt(c0))*sin(beta*c0/2.0)
     &      /(beta-Glf(3,jG))

c            reflection
             Dxxr=Dxxr+Ff1*Dxx(iG,jG)*Gf2
             Dyyr=Dyyr+Ff1*Dyy(iG,jG)*Gf2 
             Dzzr=Dzzr+Ff1*Dzz(iG,jG)*Gf2 
             Dyzr=Dyzr+Ff1*Dyz(iG,jG)*Gf2
c            transmission
             Dxxt=Dxxt+Gf1*Dxx(iG,jG)*Ff2
             Dyyt=Dyyt+Gf1*Dyy(iG,jG)*Ff2
             Dzzt=Dzzt+Gf1*Dzz(iG,jG)*Ff2
             Dyzt=Dyzt+Gf1*Dyz(iG,jG)*Ff2
c            absorption
             Dxxa=Dxxa+Ff1*Dxx(iG,jG)*Ff2
             Dyya=Dyya+Ff1*Dyy(iG,jG)*Ff2
             Dzza=Dzza+Ff1*Dzz(iG,jG)*Ff2
             Dyza=Dyza+Ff1*Dyz(iG,jG)*Ff2
909         continue
808         continue


c************ s-mod **********************************************
              Refs=Dxx0*Dxxr
              Trans=Dxx0*Dxxt
              Absos=Dxxa
c             Absorbption
              As=(4.0d0*pi*gama/o)*Imag(Absos)
c             Reflection 
              Rs=cos(theta)*real(Refs*conjg(Refs))  
c             Transmission
              Trs=1.0-cos(theta)*(real(Trans*conjg(Trans))-
     &        2.0d0*real(Trans))   

c************ p-mod **********************************************
              Refp=Dyy0*(Dyyr*cos(theta)-
     &        Dzzr*sin(theta)*sin(theta)/cos(theta))

              Tranp=Dyy0*(Dyyt*cos(theta)-2.0*Dyzt*sin(theta)+
     &        Dzzt*sin(theta)*sin(theta)/cos(theta))

              Absop=cos(theta)*Dyya*cos(theta)+
     &        sin(theta)*Dzza*sin(theta)-
     &        2.0*cos(theta)*Dyza*sin(theta)


c             Absorbption
              Ap=(4.0d0*pi*gama/o)*Imag(Absop)
c             Reflection 
              Rp=real(Refp*conjg(Refp))  
c             Transmission
              Trp=1.0-real(Tranp*conjg(Tranp))-
     &        cos(theta)*2.0d0*real(Tranp) 

                        
              write(100+itheta,*)o*Hartree,Ap
              write(200+itheta,*)o*Hartree,1-Ap-Rp
              write(300+itheta,*)o*Hartree,Rp

              if(itheta.eq.1)then
c             vodljivost u jedinicama pi*e^2/2h   
              write(131,*)o*Hartree,real(-ione*4.0*c0*Pixx(1,1)/oi)
              write(132,*)o*Hartree,imag(-ione*4.0*c0*Pixx(1,1)/oi)
              endif 
          
c             end of omega loop
9             continue      

c             end of  q loop
8             continue         

c             end of angle loop         
801           continue


            

            

              end             





