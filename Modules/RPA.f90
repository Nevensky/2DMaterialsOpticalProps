module RPA
  use iso_fortran_env, only: dp => real64
  implicit none
  
  public :: genCurrentVertices
  private

contains 
  subroutine genCurrentVertices(pol, jump, eps, Gcar, qx,qy,qz, kx,ky,kz, Nlf, iG0, NG1, NG2, NGd, R1, R2, R, RI, Glf, G, Gfast, C1, C2, MnmK1K2, MnmK1K22)
    ! Konstrukcija matricnih elementa strujnih vrhova MnmK1K2(iG) i MnmK1K2(iG) 
    character(len=3), intent(in)    :: pol
    integer,          intent(in)    :: iG0, Nlf, NG1, NG2, NGd
    integer,          intent(in)    :: R1,R2
    real(kind=dp),    intent(in)    :: eps
    real(kind=dp),    intent(in)    :: Gcar
    real(kind=dp),    intent(in)    :: qx,qy,qz
    real(kind=dp),    intent(in)    :: kx,ky,kz
    real(kind=dp),    intent(in)    :: R(:,:,:)
    real(kind=dp),    intent(in)    :: RI(:,:,:)
    real(kind=dp),    intent(in)    :: Glf(:,:)
    real(kind=dp),    intent(in)    :: G(:,:)  ! polje valnih vektora G u recp. prost. za wfn.
    complex(kind=dp), intent(in)    :: C1(:)   ! dim NG
    complex(kind=dp), intent(in)    :: C2(:)   ! dim NG (bcs the input is dim iband x NG )
    integer,          intent(inout) :: jump
    integer,          intent(inout) :: Gfast(:)
    complex(kind=dp), intent(out)   :: MnmK1K2(:)
    complex(kind=dp), intent(out)   :: MnmK1K22(:)
  
    integer          :: iG, iG1, iG2
    integer          :: iGfast
    real(kind=dp)    :: K11, K22, K33
    real(kind=dp)    :: Gxx1,Gyy1,Gzz1
    real(kind=dp)    :: Gxx2,Gyy2,Gzz2
    real(kind=dp)    :: current, current_y, current_z
    
    ! neven debug
    ! print *, 'kx,ky,kz:', kx,ky,kz
    ! print *, 'qx,qy,qz:', qx,qy,qz
    ! print *, 'Glf:', Glf(1:3,1:Nlf)
    ! print *, '-----'
    ! print *, 'Gfast(1:Nlf*NGd)',Gfast(1:Nlf*NGd)
    ! print *,'iG0, Nlf, NG1, NG2, NGd'
    ! print *, iG0, Nlf, NG1, NG2, NGd
    ! print *,'R1,R2: ',R1,R2
    ! print *,'Glf(1:5,1:5)',Glf(1:5,1:5)
    ! print *, 'G(1:3,1:5,',G(1:3,1:5)
    ! print *,'R, RI,',R(1:2,1,1:3), RI(1:2,1,1:3)
    ! print *,'Gfast(1:5)',Gfast(1:5)
    ! print *,'C1(2), C2(2):',C1(2), C2(2)
  
    iGfast = 0
    MnmK1K2(1:Nlf)  = cmplx(0.0,0.0)
    MnmK1K22(1:Nlf) = cmplx(0.0,0.0)
    iG_loop: do  iG = 1,Nlf
      iG1_loop: do iG1 = 1,NGd
        iGfast = iGfast + 1
        Gxx1 = G(1,iG1)
        Gyy1 = G(2,iG1)
        Gzz1 = G(3,iG1)
        ! neven debug -> provjeri tocnost
        k11 = R(R1,1,1)*Gxx1 + R(R1,1,2)*Gyy1 + R(R1,1,3)*Gzz1
        k22 = R(R1,2,1)*Gxx1 + R(R1,2,2)*Gyy1 + R(R1,2,3)*Gzz1
        k33 = R(R1,3,1)*Gxx1 + R(R1,3,2)*Gyy1 + R(R1,3,3)*Gzz1
        ! k11 = sum( R(R1,1,1:3)*G(1:3,iG1) )
        ! k22 = sum( R(R1,2,1:3)*G(1:3,iG1) )
        ! k33 = sum( R(R1,3,1:3)*G(1:3,iG1) )
        if (pol == 'xx') then
          ! neven debug
          ! if (kx /= 0.0) then
            ! print *,'qx,kx,Glf(1,iG),k11,Gcar'
            ! print *,qx,kx,Glf(1,iG),k11,Gcar
          ! end if
          current = (qx + 2.0*kx + Glf(1,iG) + 2.0*k11)*Gcar
          ! print *,'current:',current
        else if (pol == 'yy') then
          current = (qy + 2.0*ky + Glf(2,iG) + 2.0*k22)*Gcar
        else if (pol == 'zz') then
          current = (qz + 2.0*kz + Glf(3,iG) + 2.0*k33)*Gcar
        else if (pol == 'yz') then
          current_y = (qy + 2.0*ky + Glf(2,iG) + 2.0*k22)*Gcar
          current_z = (qz + 2.0*kz + Glf(3,iG) + 2.0*k33)*Gcar
        end if
        k11 = k11 + Glf(1,iG) + G(1,iG0)
        k22 = k22 + Glf(2,iG) + G(2,iG0)
        k33 = k33 + Glf(3,iG) + G(3,iG0)
        
        Gxx1 = RI(R2,1,1)*k11 + RI(R2,1,2)*k22 + RI(R2,1,3)*k33
        Gyy1 = RI(R2,2,1)*k11 + RI(R2,2,2)*k22 + RI(R2,2,3)*k33
        Gzz1 = RI(R2,3,1)*k11 + RI(R2,3,2)*k22 + RI(R2,3,3)*k33
        if (jump == 1) then
          iG2_loop: do iG2 = 1,NG2
            Gfast(iGfast) = NG2 + 1
            Gxx2 = G(1,iG2)
            Gyy2 = G(2,iG2)
            Gzz2 = G(3,iG2)
            if (      abs(Gxx2-Gxx1) < eps &
                .and. abs(Gyy2-Gyy1) < eps &
                .and. abs(Gzz2-Gzz1) < eps ) then
              Gfast(iGfast) = iG2
              EXIT iG2_loop
            end if
          end do iG2_loop
        end if
        iG2 = Gfast(iGfast)
        if (iG2 <= NG2) then
          ! ako je polarazcija je tipa xx, yy, zz 
          if (pol == 'xx' .or. pol== 'yy' .or. pol == 'zz') then
            MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5D0*conjg(C1(iG1)) * current * C2(iG2)
            MnmK1K22(iG) = MnmK1K2(iG) ! strujni vrhovi su isti
            ! neven debug
            ! print *,'MnmK1K2(iG)',MnmK1K2(iG)
            ! print *, '0.5D0*conjg(C1(iG1))',0.5D0*conjg(C1(iG1))
            ! print *,'current',current
            ! print *,'C2',C2(iG2)
          elseif (pol =='yz' .or. pol =='zy') then ! ako  su miksani yz
            MnmK1K2(iG)  = MnmK1K2(iG)  + 0.5D0*conjg(C1(iG1)) * current_y * C2(iG2)
            MnmK1K22(iG) = MnmK1K22(iG) + 0.5D0*conjg(C1(iG1)) * current_z * C2(iG2)
          else
            print *,'WARNING Specified mixed polarization component not supported.'//adjustl(trim(pol))//' not allowed.'
            stop
          end if
        end if
      end do iG1_loop
    end do iG_loop
    jump = 2 ! za svaki valni vektor q i dani k zapamti Gfast i za svaku vrpcu preskaci taj postupak   
      
  end subroutine genCurrentVertices
end module RPA