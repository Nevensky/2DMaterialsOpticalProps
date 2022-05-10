module io_qe
  use iso_fortran_env, only: dp => real64
  implicit none

  public :: loadG, loadCs, loadCs_full
  private

contains
  subroutine loadG(savedir, KC, NG, G, parG)
    !! Reads the reciprocal vectors in crystal coordinates and transformation
    !! in Cartesian cordinates.
    character(len=*),  intent(in)   :: savedir
    real(kind=dp),     intent(in)   :: KC(3,3)
    integer,           intent(out)  :: NG                   !! No. of rec. latt. vectors
    real(kind=dp),     intent(out), allocatable  :: G(:,:)  !! rec. latt. wavevector G (for the density)
    integer, optional, intent(out), allocatable  :: parG(:) !! parity of each rec. latt. wavevector G

    integer :: iuni, ios0, ios1, ios2
    integer :: n, m, iG, Nspin
    logical :: gamma_only
    integer, allocatable :: Gi(:,:)
    character(len=218)   :: fname            

    fname = trim(savedir)//'/charge-density.dat'
    print *,'status: Reading Gvecs from file: ',adjustl(trim(fname))
    
    open(newunit=iuni,file=fname,form ='unformatted',action='read',status='old',iostat=ios0)
    if (ios0 /= 0) then
      stop "ERROR: Can\'t open G-vector file."
    end if

    read(iuni, iostat=ios1) gamma_only, NG, Nspin
    ! print *, 'Number of Gvecs (NG):', NG

    read(iuni, iostat=ios1) ! skip reading b1, b2, b3 rec.latt.vecs. 

    if (.not. allocated(Gi)) allocate(Gi(3,NG))
    read (iuni, iostat=ios2) Gi(1:3,1:NG)
    close(iuni)

    if (ios1 /= 0) then
      stop 'ERROR: Failed to read NG and Nspin from G-vector file.'
    else if (ios2 /=0 ) then
      stop 'ERROR: Failed to read Miller indices from G-vector file.'
    end if

    if (Gi(1,1) /= 0  .or.  Gi(2,1) /= 0  .or.  Gi(3,1) /= 0) then
      print *, 'WARNING: G vectors input is wrong. G(1:3,1) is not (0,0,0)!'
    end if

    ! transformation to cart.coord (now all G components are in 2pi/a0 units)
    if (.not. allocated(G)) allocate(G(3,NG))
    G(1:3,1:NG) = 0.0
    do iG=1,NG
      do n = 1,3
        do m = 1,3
          G(n,iG) = G(n,iG) + KC(n,m)*dble( Gi(m,iG) )
        end do
      end do
      if (present(parG)) then 
        if (.not. allocated(parG)) allocate(parG(NG))
        parG(iG) = Gi(3,iG) 
      end if
    end do
    close(iuni)

    deallocate(Gi)

  end subroutine loadG

  subroutine loadCs(ik, ibnd, savedir, igwx, evc)
    ! read_a_wfc(ibnd, filename, evc, ik, xk, nbnd, ispin, npol, gamma_only, ngw, igwx )
    ! read QE 6.0 and greater, wfn coefficeints
    character (len=*), intent(in)                  :: savedir
    integer,           intent(in)                  :: ik, ibnd
    integer,           intent(out)                 :: igwx
    complex(dp),       intent(inout)               :: evc(*)

    character (len=300) :: path 

    integer  :: nbnd, npol, i, ngw
    ! real(dp) :: xk(3)
    ! integer  :: ik2, ispin, dummy_int   
    ! logical  :: gamma_only 
    integer  :: ios, iuni 
    ! real(dp) :: scalef
    ! real(dp) :: b1(3), b2(3), b3(3) !, dummy_real 

    character(len=4)   :: str1 ='/wfc'
    character(len=4)   :: str3 ='.dat'
    character(len=100) :: ik_str
    write(ik_str,'(I10)') ik

    path = trim(savedir)//str1//trim(adjustl(ik_str))//str3
    ! iuni = 10 + ik*100 + ibnd*20000
    ! print *,path
    open(newunit = iuni, file = trim(adjustl(path)),action='read', form = 'unformatted', status = 'old', iostat=ios) 
    read(iuni) ! ik2, xk, ispin, gamma_only, scalef
    read(iuni) ngw, igwx, npol, nbnd
    read(iuni) ! b1, b2, b3 

    ! avoid reading miller indices of G vectors below E_cut for this kpoint 
    ! if needed allocate  an integer array of dims (1:3,1:igwx) 

    read(iuni) ! dummy_int
    ! allocate (evc(npol*igwx))
    if ( ibnd > nbnd) then 
       print '("looking for band nr. ",I7," but there are only ",I7," bands in the file")',ibnd, nbnd
       stop
    end if 
    do i = 1, nbnd 
       if ( i == ibnd ) then 
          read(iuni) evc(1:npol*igwx) 
          exit 
       else 
          read(iuni) ! dummy_real
       end if 
    end do 
    close(iuni)
!    deallocate(evc) 
  end subroutine loadCs

  subroutine loadCs_full(ik, savedir, igwx, evc)
    ! read_a_wfc(ibnd, filename, evc, ik, xk, nbnd, ispin, npol, gamma_only, ngw, igwx )
    ! read QE 6.0 and greater, wfn coefficeints
    character (len=*), intent(in)                  :: savedir
    integer,           intent(in)                  :: ik
    integer,           intent(out)                 :: igwx
    complex(dp),       intent(inout)               :: evc(:,:) ! ngk x nbnd

    character (len=300) :: path 

    integer  :: nbnd, npol, i, ngw
    ! real(dp) :: xk(3)
    ! integer  :: ik2, ispin, dummy_int   
    ! logical  :: gamma_only 
    integer  :: ios, iuni 
    ! real(dp) :: scalef
    ! real(dp) :: b1(3), b2(3), b3(3) !, dummy_real 

    character(len=4)   :: str1 ='/wfc'
    character(len=4)   :: str3 ='.dat'
    character(len=100) :: ik_str
    write(ik_str,'(I10)') ik

    path = trim(savedir)//str1//trim(adjustl(ik_str))//str3
    ! iuni = 10 + ik*100 + ibnd*20000
    ! print *,path
    open(newunit = iuni, file = trim(adjustl(path)),action='read', form = 'unformatted', status = 'old', iostat=ios) 
    read(iuni) ! ik2, xk, ispin, gamma_only, scalef
    read(iuni) ngw, igwx, npol, nbnd
    read(iuni) ! b1, b2, b3 

    ! avoid reading miller indices of G vectors below E_cut for this kpoint 
    ! if needed allocate  an integer array of dims (1:3,1:igwx) 

    read(iuni) ! dummy_int
    ! allocate (evc(npol*igwx))

    do i = 1, Nbnd 
        read(iuni) evc(1:npol*igwx,i) 
    end do 
    close(iuni)
!    deallocate(evc) 
  end subroutine loadCs_full

  subroutine loadCs_spin(ispin, ik, ibnd, savedir, igwx, evc, jspin)
    ! read_a_wfc(ibnd, filename, evc, ik, xk, nbnd, ispin, npol, gamma_only, ngw, igwx )
    ! read QE 6.0 and greater, wfn coefficeints
    character (len=*), intent(in)                  :: savedir
    integer,           intent(in)                  :: ik, ibnd
    integer,           intent(in)                  :: jspin
    integer,           intent(out)                 :: igwx
    complex(dp),       intent(inout)               :: evc(*)

    character (len=300) :: path 

    integer  :: nbnd, ispin, npol,  i, ik2, ngw
    real(dp) :: xk(3)
    ! integer  :: dummy_int   
    logical  :: gamma_only 
    integer  :: ios, iuni 
    real(dp) :: scalef
    ! real(dp) :: b1(3), b2(3), b3(3) !, dummy_real 

    character(len=4)   :: str1 ='/wfc'
    character(len=4)   :: str3 ='.dat'
    character(len=100) :: ik_str
    write(ik_str,'(I10)') ik

    path = trim(savedir)//str1//trim(adjustl(ik_str))//str3
    ! iuni = 10 + ik*100 + ibnd*20000
    ! print *,path
    open(newunit = iuni, file = trim(adjustl(path)),action='read', form = 'unformatted', status = 'old', iostat=ios) 
    read(iuni) ik2, xk, ispin, gamma_only, scalef
    read(iuni) ngw, igwx, npol, nbnd
    read(iuni) ! b1, b2, b3 

    ! avoid reading miller indices of G vectors below E_cut for this kpoint 
    ! if needed allocate  an integer array of dims (1:3,1:igwx) 

    read(iuni) ! dummy_int
    ! allocate (evc(npol*igwx))
    if ( ibnd > nbnd) then 
       print '("looking for band nr. ",I7," but there are only ",I7," bands in the file")',ibnd, nbnd
       stop
    end if 
    do i = 1, nbnd 
       if ( i == ibnd .and. ispin == jspin) then 
          read(iuni) evc(1:npol*igwx) 
          exit 
       else 
          read(iuni) ! dummy_real
       end if 
    end do 
    close(iuni)
!    deallocate(evc) 
  end subroutine loadCs_spin

  subroutine loadCs_spin_full(ik, savedir, igwx, evc, jspin)
    ! read_a_wfc(ibnd, filename, evc, ik, xk, nbnd, ispin, npol, gamma_only, ngw, igwx )
    ! read QE 6.0 and greater, wfn coefficeints
    character (len=*), intent(in)                  :: savedir
    integer,           intent(in)                  :: ik
    integer,           intent(in)                  :: jspin
    integer,           intent(out)                 :: igwx
    complex(dp),       intent(inout)               :: evc(:,:) ! ngk x nbnd

    character (len=300) :: path 

    integer  :: nbnd, ispin, npol, i, ik2, ngw
    real(dp) :: xk(3)
    ! integer  :: dummy_int   
    logical  :: gamma_only 
    integer  :: ios, iuni 
    real(dp) :: scalef
    ! real(dp) :: b1(3), b2(3), b3(3) !, dummy_real 

    character(len=4)   :: str1 ='/wfc'
    character(len=4)   :: str3 ='.dat'
    character(len=100) :: ik_str
    write(ik_str,'(I10)') ik

    path = trim(savedir)//str1//trim(adjustl(ik_str))//str3
    ! iuni = 10 + ik*100 + ibnd*20000
    ! print *,path
    open(newunit = iuni, file = trim(adjustl(path)),action='read', form = 'unformatted', status = 'old', iostat=ios) 
    read(iuni) ik2, xk, ispin, gamma_only, scalef
    read(iuni) ngw, igwx, npol, nbnd
    read(iuni) ! b1, b2, b3 

    ! avoid reading miller indices of G vectors below E_cut for this kpoint 
    ! if needed allocate  an integer array of dims (1:3,1:igwx) 

    read(iuni) ! dummy_int
    ! allocate (evc(npol*igwx))

    if (ispin==jspin) then
      do i = 1, nbnd 
          read(iuni) evc(1:npol*igwx,i) 
      end do 
    else
      evc(1:npol*igwx,1:nbnd) = dcmplx(0.0_dp,0.0_dp)
    end if

    close(iuni)
   ! deallocate(evc) 
  end subroutine loadCs_spin_full

end module io_qe
