        subroutine paths(root1,K1,K2,n,m,pathK1,pathK2,bandn,bandm)

               INTEGER K1,K2,n,m,nord
               CHARACTER*100 root1,root2,folder,file,pathK1,pathK2,
     &         bandn,bandm
     

     
                  root2='/MoS2.save'
                  file='/evc.dat'
        
		  folder='/K0000i'
                  nord=INDEX(folder,'i', BACK = .false.)  
		  if(K1.lt.10)then
                  write (folder(nord:nord),'(i1)')K1 
                  elseif(K1.ge.10.and.K1.lt.100)then
                  write (folder(nord-1:nord),'(i2)')K1 
                  elseif(K1.ge.100.and.K1.lt.1000)then                
                  write (folder(nord-2:nord),'(i3)')K1 
                  elseif(K1.ge.1000.and.K1.lt.10000)then                  
                  write (folder(nord-3:nord),'(i4)')K1 
                  else               
                  write (folder(nord-4:nord),'(i5)')K1 
                  endif
                  pathK1=TRIM(root1)//TRIM(root2)
     &            //TRIM(folder)//TRIM(file)

                  folder='/K0000i'
                  nord=INDEX(folder,'i', BACK = .false.)
		  if(K2.lt.10)then
                  write (folder(nord:nord),'(i1)')K2 
                  elseif(K2.ge.10.and.K2.lt.100)then
                  write (folder(nord-1:nord),'(i2)')K2
                  elseif(K2.ge.100.and.K2.lt.1000)then                
                  write (folder(nord-2:nord),'(i3)')K2 
                  elseif(K2.ge.1000.and.K2.lt.10000)then                  
                  write (folder(nord-3:nord),'(i4)')K2 
                  else               
                  write (folder(nord-4:nord),'(i5)')K2 
                  endif
                  pathK2=TRIM(root1)//TRIM(root2)
     &            //TRIM(folder)//TRIM(file)



                  bandn='evc.n'
                  nord=INDEX(bandn,'n', BACK = .false.)
		  if(n.lt.10)then
                  write (bandn(nord:nord),'(i1)')n
                  elseif(n.ge.10.and.n.lt.100)then
                  write (bandn(nord:nord+1),'(i2)')n
                  else                
                  write (bandn(nord:nord+2),'(i3)')n
                  endif
                  bandm='evc.n'
                  nord=INDEX(bandm,'n', BACK = .false.)
		  if(m.lt.10)then
                  write (bandm(nord:nord),'(i1)')m
                  elseif(m.ge.10.and.m.lt.100)then
                  write (bandm(nord:nord+1),'(i2)')m
                  else                
                  write (bandm(nord:nord+2),'(i3)')m
                  endif
                  
 

  
                   return
                      end


