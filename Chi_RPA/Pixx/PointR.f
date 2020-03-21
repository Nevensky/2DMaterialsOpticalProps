        SUBROUTINE PointR(root,nsim,R,RI,f)
                          
        implicitnone

        INTEGER i,nsim,is,n,m,j,k
        REAL*8 R,RI,zero,one,x,y,z,f,fx,fy,fz,fi
        double complex T,unit
        PARAMETER(zero=0.0,one=1.0)   
        DIMENSION R(48,3,3),T(3,3),RI(48,3,3),unit(3,3),f(48,3),
     &  fi(48,3)
        CHARACTER*11,buffer1,tag1
        CHARACTER*7,buffer2,tag2
        CHARACTER*100 root,fajl,path   
        

c       point group transformations R readed from 'MoS2.sc.out'
c       matrices R are in cart. coordinate system  because k 
c       from IBZ are printed  in cart.coord. 


        fajl='/MoS2.sc.out'
        path=TRIM(root)//TRIM(fajl)
        tag1='     atomic' 
        tag2=' cryst.' 
       
      

        open(1,FILE=path) 
        
        do 66 i=1,5000
        read(1,'(a)')buffer1
        if(buffer1.eq.tag1)then 
        read(1,80)
        read(1,80)
        read(1,*)nsim
        goto 999 
        endif
66      continue  
999     continue

        is=0        
        do 67 i=1,5000
        read(1,'(a)')buffer2
        if(buffer2.eq.tag2)then 
        is=is+1        
        read(1,80)
        read(1,80)
        read(1,80)
        read(1,70)x,y,z,fx
        R(is,1,1)=x
        R(is,1,2)=y
        R(is,1,3)=z
        f(is,1)=fx
        read(1,70)x,y,z,fy
        R(is,2,1)=x
        R(is,2,2)=y
        R(is,2,3)=z      
        f(is,2)=fy
        read(1,70)x,y,z,fz
        R(is,3,1)=x
        R(is,3,2)=y
        R(is,3,3)=z
        f(is,3)=fz
        if(is.eq.nsim)goto 998
        endif 
67      continue
998     continue  
        close(1)
70      FORMAT(19X,3F11.3,10X,F11.3)  
80      FORMAT(X) 


c       INVERTION
        do 68 i=1,nsim
        do 828 n=1,3
        do 818 m=1,3
        unit(n,m)=dcmplx(zero,zero)
        T(n,m)=dcmplx(R(i,n,m),zero)
818     continue
        unit(n,n)=dcmplx(one,zero)
828     continue
        call gjel(T,3,3,unit,3,3)
        do 88 n=1,3
        do 89 m=1,3
        RI(i,n,m)=real(T(n,m))
89      continue
88      continue
68      continue

c        Multiplication of f by RI
         do 130, i=1,nsim   
         do 140, j=1,3   
         fi(i,j)=zero
         do 150, k=1,3     
         fi(i,j)=fi(i,j)+RI(i,j,k)*f(i,k)     
150      continue
140      continue
130      continue          


          RETURN
           END





