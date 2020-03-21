            program FOTON

c       program for the calculation of the current-current reponse  function in MoS2
c       XY-plane unit cell parameter  a0  
c       Z-dir unit cell parameter c0
c       gama=e^{2}/(\hbar c) --  fine-structure constant
 
      
        implicitnone

        integer n,no,nq,noladd,nL,pol

c       Nlfd je minimalno 2 zbog 2X2 blok matrice za p mod    
        parameter(no=5001,noladd=401,nq=401,nL=30)
          
        double precision zero,one,two,three,four,six,pi,eta,eV,
     &  Hartree,aBohr,planck,kb
        REAL*8 a0,c0,gama,domega,h

        parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,
     &  four=4.0d0,six=6.0d0,pi=3.141592654d0,a0=5.9715,c0=29.8575d0,
     &  kB=1.3806503d-23,eV=1.602176487d-19,Hartree=2.0d0*13.6056923d0,
     &  aBohr=0.5291772d0,planck=6.626196d-34,
     &  eta=0.00001/Hartree,gama=1.0d0/137.0d0)
    

c     ...scalars...
        integer i,j,k,m,iq,io
        real*8 q,o,omax,dq,omin
        complex*8 Pi_ladder_down,Pi_ladder_up,Pi_RPA,oi,D0,
     &  beta,Dxx,Dyy,Pis,Pip,czero,ione,ieta,cone,Dyz,Dzz
c     ...arrays...
        COMPLEX*8 Pixx,Pizz,Pi0
        double complex eps,unit
        DIMENSION Pixx(no),Pizz(no),Pi0(nL,nL),eps(nL,nL),
     &  unit(nL,nL)        
        character*100 nis,dato,root,path,fajl,rootRPA,
     &  rootladddown,rootladdup
 
    
             czero=cmplx(zero,zero)
             cone=cmplx(one,zero)
             ione=cmplx(zero,one)
             ieta=cmplx(zero,eta)
          
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
              root='/home/vito/PROJECTS/MoS2-BSE/MoS2_51x51'

c             RPA IRREDUCIBLE POLARIZABILICI Chi_RPA FOLDER: 
              rootRPA='/home/vito/PROJECTS/MoS2-BSE/Chi_RPA'

c             LADDER IRREDUCIBLE POLARIZABILICI Chi_ladd_down FOLDER: 
              rootladddown='/home/vito/PROJECTS/MoS2-BSE/Chi_ladd_down'

c             LADDER IRREDUCIBLE POLARIZABILICI Chi_ladd_up FOLDER: 
              rootladdup='/home/vito/PROJECTS/MoS2-BSE/Chi_ladd_up'


c            udaljenost medju slojevima u Ang
             h=3.0/aBohr  


c            for s-mode pol=1  
c            for p-mode pol=2  
               
              pol=1

           
              dq=0.00001d0            
              omin=1.0d-5
              omax=(50.0/Hartree+omin)
              domega=(omax-omin)/(no-1)
      
              fajl='/Pixx/Pi_RPA_xx'
              path=TRIM(rootRPA)//TRIM(fajl)
              open(31,file=path)
              fajl='/Pi_ladder_down_x'
              path=TRIM(rootladddown)//TRIM(fajl)
              open(32,file=path)
              fajl='/Pi_ladder_up_x'
              path=TRIM(rootladdup)//TRIM(fajl)
              open(33,file=path)
              do 999 io=1,noladd
              read(31,*)o,Pi_RPA
              if(io.le.noladd)then
              read(32,*)o,Pi_ladder_down
              read(33,*)o,Pi_ladder_up
              endif
              Pixx(io)=Pi_RPA+Pi_ladder_down+Pi_ladder_up
999           continue
              close(31)
              close(32)
              close(33) 
     

              fajl='/Pizz/Pi_RPA_zz'
              path=TRIM(rootRPA)//TRIM(fajl)
              open(31,file=path)
              fajl='/Pi_ladder_down_z'
              path=TRIM(rootladddown)//TRIM(fajl)
              open(32,file=path)
              fajl='/Pi_ladder_up_z'
              path=TRIM(rootladdup)//TRIM(fajl)
              open(33,file=path)
              do 979 io=1,noladd
              read(31,*)o,Pi_RPA
              if(io.le.noladd)then
              read(32,*)o,Pi_ladder_down
              read(33,*)o,Pi_ladder_up
              endif
              Pizz(io)=Pi_RPA+Pi_ladder_down+Pi_ladder_up
979           continue
              close(31)
              close(32)
              close(33) 




c*******************************************************************
c              POCETAK PETLJI PO Q I OMEGA
c******************************************************************
                open(3,file='plot') 

c              Q loop strts here 
               do 801 iq=1,nq 
               q=(iq-1)*dq  

               print*,iq

c              omega loop starts here 
               do 9 io=3,noladd
               o=omin+(io-1)*domega  
               oi=cmplx(o,eta) 
               beta=cmplx(gama*gama*oi*oi-q*q)
               beta=sqrt(beta)
 
              

               Dxx=2.0d0*pi*ione*c0*gama*gama/beta 
               Dyy=2.0d0*pi*ione*c0*beta/(oi*oi)
               Dyz=-2.0*pi*ione*Q*c0/(oi*oi) 
               Dzz=2.0*pi*ione*Q*Q*c0/(beta*oi*oi)
               
              
              
c             S-MOD
              if(pol.eq.1)D0=Dxx              
c             P-MOD
              if(pol.eq.2)D0=Dyy        


              eps=czero
              do 777 i=1,nL
              do 888 j=1,nL 
              eps(i,j)=
     &        -Pixx(io)*D0*exp(ione*h*beta*abs(dble(i)-dble(j))) 
              unit(i,j)=czero
              Pi0(i,j)=czero
888           continue 
              eps(i,i)=cone+eps(i,i)
              unit(i,i)=cone
              Pi0(i,i)=Pixx(io)
777           continue
              
              call gjel(eps,nL,nL,unit,nL,nL)
 
              Pip=czero  
              do 909 i=1,nL 
909           Pip=Pip+eps(1,i)*Pi0(i,1)

          
c           optical conductivity in units  pi*e^2/2h s-mode
            if(iq.eq.1)write(100,*)o*Hartree,real(-ione*4.0*c0*Pip/oi)
c            if(iq.eq.1)write(200,*)o*Hartree,real(-ione*4.0*c0*Pip/oi)



         
            write(3,*)10.0*q/aBohr,
     &      o*Hartree,real(-ione*4.0*c0*Pip/oi)



c            kraj po omega
9            continue      



c            end of q loop   
801          continue

             close(3)
            

              end             





