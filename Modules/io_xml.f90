!! This modules reads the data-file-schema.xml of a Quantum Espresso
!! pwscf.x calculation and returns the requested quantities:
!! direct lattice: a
!! reciprocal lattice: b
!! Fermi energy: EFermi
!! No. of electrons: NelQE
!! No. of rotational symmetries: Nsym
!! No. of bands: Nband
!! No. k points in the IBZ: NkI
!! Monkhorst-Pack mesh dimensions: Nmp(3)
!! No. of planewaves at each k point: Npw(NkI)
!! rotational symmetries: R
!! k points in the IBZ: kI
!! band eigenvalues at each k point: eigenvals (Nband x NkI)
!! band occupations at each k point: eigenvals (Nband x NkI)
module io_xml
  use iso_fortran_env, only : iostat_end, dp => real64
  implicit none

  public :: loadXML_qe
  private 
  
contains
  subroutine loadXML_qe(path, Nband , NkI, Nrot, Nsym, Nmp, alat, a, b, R, Npw, kI, eigenvals, occupations, T, Efermi, NelQE, printOutput)
    character(len=*), intent(in) :: path
    integer,       intent(out), optional :: Nband, NkI, Nrot, Nsym, Nmp(3)
    real(kind=dp), intent(out), optional :: alat      !! lattice parameter
    real(kind=dp), intent(out), optional              :: a(3,3)           !! direct lattice
    real(kind=dp), intent(out), optional              :: b(3,3)           !! reciprocal lattice
    real(kind=dp), intent(out), optional, allocatable :: R(:,:,:)         !! array of rotational matrices (3 x 3 x Nrot)
    integer,       intent(out), optional, allocatable :: Npw(:)           !! array of G-vectors at each k-point (Nk)
    real(kind=dp), intent(out), optional, allocatable :: kI(:,:)          !! array of k-points in the ireducible Brillouin zone (3 x NkI)
    real(kind=dp), intent(out), optional, allocatable :: eigenvals(:,:)   !! eigenvalues for each band and k-point [Hartree] (Nband x Nk)
    real(kind=dp), intent(out), optional, allocatable :: occupations(:,:) !! occupations for each band and k-point (Nband x Nk)
    real(kind=dp), intent(out), optional :: T           !! electron temperature
    real(kind=dp), intent(out), optional :: Efermi      !! Fermi energy
    real(kind=dp), intent(out), optional :: NelQE       !! No. electrons
    logical,       intent(in),  optional :: printOutput !! print values/arrays read from the XML file to standard output
    
    ! integer, save, intent(out), optional :: Nband=0, NkI=0, Nrot=0, Nsym=0, Nmp(3) = 0

    integer, save :: ik_eig=0, ik_occ=0, ik_Npw=1 ! k-point counters/indices

    logical             :: printOutput_ = .false.               
    integer             :: iuni, ios_fopen, ios_fread=0 ! file i/o status vars
    integer             :: i, is, ie, ik                ! heleper vars
    integer,save        :: c_latt=0, c_reclatt=0, c_sym1=0, c_sym2=0, c_kpts=0, c_eig=0, c_occ=0 ! counters
    character(len=200)  :: buffer ! stores a line of the given file
    logical             :: noinv, force_symmorphic

    if (present(printOutput)) printOutput_ = printOutput
    if (present(Nband))  Nband = 0
    if (present(NkI))    NkI    = 0
    if (present(Nrot))   Nrot   = 0
    if (present(Nsym))   Nsym   = 0
    if (present(Nmp))    Nmp(3) = 0

    open(newunit=iuni,file=path,iostat=ios_fopen,status='old')  
    file_loop: do while (ios_fread /= iostat_end)
      read (iuni, '(A)', iostat=ios_fread, advance='NO') buffer

      call readXML_logicalval(noinv,'<noinv>',buffer)
      call readXML_logicalval(force_symmorphic,'<force_symmorphic>',buffer)

      if (present(Nrot))   call readXML_intval(Nrot,'<nrot>',buffer)
      if (present(Nsym))   call readXML_intval(Nsym,'<nsym>',buffer)
      if (present(NkI))    call readXML_intval(NkI,'<nks>',buffer)
      if (present(Nband))  call readXML_intval(Nband,'<nbnd>',buffer)
      if (present(NelQE))  call readXML_realval(NelQE,'<nelec>',buffer)
      if (present(Efermi)) call readXML_realval(Efermi,'<fermi_energy>',buffer)
      if (present(T))      call readXML_realtag(T,'<smearing',buffer)

      if (present(a))      call readXML_latt(a,'<cell>',buffer,c_latt)
      if (present(b))      call readXML_latt(b,'<reciprocal_lattice>',buffer,c_reclatt)
      
      if (present(alat) .and. present(a)) alat = a(1,1)

      if (present(Nmp)) then
        call readXML_inttag(Nmp(1),'<monkhorst_pack',buffer,'nk1')
        call readXML_inttag(Nmp(2),'<monkhorst_pack',buffer,'nk2')
        call readXML_inttag(Nmp(3),'<monkhorst_pack',buffer,'nk3')
      end if

      if (all([present(Nrot),present(R)])) then
        ! reads all rotational symmetries
        if (Nrot/=0) then
          if (.not. allocated(R)) allocate(R(3,3,Nrot))
          call readXML_rot(R,'<rotation',buffer,c_sym1, c_sym2)
          if (c_sym1==4) then
            c_sym2 = c_sym2+1
            c_sym1 = 0
          end if
        end if
      end if

      if (all([present(NkI),present(kI)])) then
        ! reads all k points in the IBZ
        if (NkI/=0) then
          if (.not. allocated(kI) ) allocate(kI(3,NkI))
          call readXML_3realval(kI,'<k_point',buffer,c_kpts)
        end if
      end if

      if (all([present(NkI),present(Npw)])) then
        ! reads number of plane waves at each k point
        if (NkI /= 0 .and. ik_Npw<=NkI) then
          if (.not. allocated(Npw)) allocate(Npw(NkI))
          call readXML_intval(Npw(ik_Npw),'<npw>',buffer,ik_Npw)
        end if
      end if

      if (all([present(Nband),present(NkI)])) then
        ! reads all band eigenvalues and their respective occupations
        if (NkI/=0 .and. Nband/=0) then
          if (present(eigenvals)) then
            if (.not. allocated(eigenvals)) allocate(eigenvals(Nband,NkI))
            call readXML_eigenvals(eigenvals,'<eigenvalues size=',buffer,c_eig,ik_eig)
          end if
          if (present(occupations)) then
            if (.not. allocated(occupations)) allocate(occupations(Nband,NkI))
            call readXML_eigenvals(occupations,'<occupations size=',buffer,c_occ,ik_occ)
          end if
        end if
      end if

    end do file_loop
    close(iuni)

    ! debug: prints all values read
    ! instead it should call a function where each 
    ! quantity is optional and if present(...) loadd it
    if (present(a) .and. present(alat) .and. printOutput_) then
      print *, 'alat: ',alat
    end if

    if (present(a) .and. printOutput_) then
      print *, 'unit_cell'
      print *, a
    end if

    if (present(b) .and. printOutput_) then
      print *, 'reciprocal_lattice'
      print *, b
    end if

    if (all([present(Nrot),present(R)]) .and. printOutput_) then
      print *, 'symmetries:'
      do is=1,Nrot
        write (*,'(A8,I2,A18)') '----- R(',is,', 1:3, 1:3)-------'
        do i=1,3
          print *, R(:,i,is)
        end do
      end do
      print *, '------------------------------'
    end if

    if (present(Nmp) .and. printOutput_) then
      print *, '--- Number of k-points in the Monkhorst-Pack mesh---'
      print *, Nmp
    end if

    if (present(Npw) .and. printOutput_) then
      print *, '--- Number of Plane Waves at each k-point ---'
      print *, Npw
    end if

    if (present(kI) .and. printOutput_) then
      print *, '----- k points -----'
      print *, kI
      print *, '------------'
    end if

    if (any([present(eigenvals),present(occupations)]) .and. printOutput_) then
      print *, 'eigenvals'
      do ik=1, NkI
        ! write(*,'(A8,I4,A5,3(F2.4),A2)'), '== ik = ',ik,' k = ',kI(:,ik), '=='
        print *,'___ ik =',ik,' ___'
        if (present(eigenvals)) then
          print *, '<eigenvalues>'
          ! write all elements before the last line
          do ie=0,Nband/5-1
            write(*,'(5f8.4)'), eigenvals(ie*5+1:ie*5+5,ik)
          enddo
          ! write elements in the last line
          do ie=1,mod(Nband,5)
            write(*,'(f8.4)',advance='no') eigenvals(Nband-mod(Nband,5)+ie,ik)
          enddo
          write(*,*) ''
        end if
        if (present(occupations)) then
          print *, '<occupations>'
          ! write all elements before the last line
          do ie=0,Nband/5-1
              write(*,'(5f8.4)'), occupations(ie*5+1:ie*5+5,ik)
          enddo
          ! write elements in the last line
          do ie=1,mod(Nband,5)
            write(*,'(f8.4)',advance='no') occupations(Nband-mod(Nband,5)+ie,ik)
          enddo
          write(*,*) ''
        end if
      end do
      print *, '------------'
    end if

    
    if (present(NkI)    .and. printOutput_) print *, 'NkI = ', NkI
    if (present(Nband)  .and. printOutput_) print *, 'Nband = ',Nband
    if (present(Nrot)   .and. printOutput_) print *, 'Nrot = ', Nrot
    if (present(Nsym)   .and. printOutput_) print *, 'Nsym = ', Nsym
    if (present(NelQE)  .and. printOutput_) print *, 'NelQE = ', NelQE
    if (present(Efermi) .and. printOutput_) print *, 'Fermi energy = ', Efermi
    if (present(T)      .and. printOutput_) print *, 'Electron smearing Temperature =',T


    if (all([present(Nsym),present(Nrot)])) then
      if (Nsym/=Nrot) then
        print *, 'WARNING: Quantum Espresso calculation must not use non-rotational symmetries. &
              & Use force_symmorphic=.true. and noinv=.true. in QE input.'
      end if
    end if

    if (.not. noinv) stop 'ERROR: Quantum espresso uses time reversal (k=-k) symmetry. &
              & Use noinv=.true. in QE input.'
    if (.not. force_symmorphic) stop 'ERROR: Quantum espresso uses fractionary translation symmetry. &
              & Use force_symmorphic=.true. =.true. in QE input.'


  end subroutine loadXML_qe

  subroutine readXML_intval(a,tag,buffer,counter)
    integer,          intent(out)   :: a
    character(len=*), intent(in)    :: buffer, tag
    integer,optional, intent(inout) :: counter
    
    character(len=200) :: dummy
    integer            :: idx_start, idx_end
    integer            :: a_tmp 

    if (index(buffer, tag) /= 0) then
      read(buffer,'(A)') dummy
      idx_start = index(dummy,'>')+1
      idx_end = index(dummy,'</')-1
      dummy = dummy(idx_start:idx_end)
      ! print *, dummy
      read(dummy,*) a_tmp
      a = a_tmp
      if (present(counter)) counter = counter+1
    end if
  end subroutine readXML_intval

  subroutine readXML_logicalval(a,tag,buffer)
    logical,          intent(out)   :: a
    character(len=*), intent(in)    :: buffer, tag
    
    character(len=200) :: dummy
    integer            :: idx_start, idx_end
    logical            :: a_tmp 

    if (index(buffer, tag) /= 0) then
      read(buffer,'(A)') dummy
      idx_start = index(dummy,'>')+1
      idx_end = index(dummy,'</')-1
      dummy = dummy(idx_start:idx_end)
      ! print *, dummy
      if (index(dummy,'true')/=0) then
        a_tmp = .true.
      else
        a_tmp = .false.
      end if
      a = a_tmp
    end if
  end subroutine readXML_logicalval

  subroutine readXML_realval(a,tag,buffer)
    real(kind=dp),     intent(out)   :: a
    character(len=*),  intent(in)    :: buffer, tag
    
    character(len=200) :: dummy
    integer            :: idx_start, idx_end
    real(kind=dp)      :: a_tmp 

    if (index(buffer, tag) /= 0) then
      read(buffer,'(A)') dummy
      idx_start = index(dummy,'>')+1
      idx_end = index(dummy,'</')-1
      dummy = dummy(idx_start:idx_end)
      ! print *, dummy
      read(dummy,*) a_tmp
      a = a_tmp
    end if
  end subroutine readXML_realval

  subroutine readXML_realtag(a,tag,buffer,tag2)
    real(kind=dp),              intent(out) :: a
    character(len=*),           intent(in)  :: buffer, tag
    character(len=*), optional, intent(in)  :: tag2
    
    character(len=200) :: dummy
    integer            :: idx_start, idx_end
    real(kind=dp)      :: a_tmp 

    if (index(buffer, tag) /= 0) then
      read(buffer,'(A)'), dummy
      if (present(tag2)) then
        idx_start = index(dummy,tag2)
        idx_start = idx_start+index(dummy(idx_start:),'"')
        idx_end = index(dummy(idx_start:),'"')+idx_start-2
        dummy = dummy(idx_start:idx_end)
      else
        read(buffer,'(A)') dummy
        idx_start = index(dummy,'"')+1
        idx_end = index(dummy(idx_start:),'"')+idx_start-2
        dummy = dummy(idx_start:idx_end)
        ! print *, dummy 
      end if
      read(dummy,*) a_tmp
      a = a_tmp
    end if
  end subroutine readXML_realtag

  subroutine readXML_inttag(a,tag,buffer,tag2)
    integer,                    intent(out) :: a
    character(len=*),           intent(in)  :: buffer, tag
    character(len=*), optional, intent(in)  :: tag2
    
    character(len=200) :: dummy
    integer            :: idx_start, idx_end,  ios
    integer            :: a_tmp

    if (index(buffer, tag) /= 0) then
      read(buffer,'(A)'), dummy
      if (present(tag2)) then
        idx_start = index(dummy,tag2)
        idx_start = idx_start+index(dummy(idx_start:),'"')
        idx_end = index(dummy(idx_start:),'"')+idx_start-2
        dummy = dummy(idx_start:idx_end)
      else
        idx_start = index(dummy,'"')
        idx_end = index(dummy(idx_start:),'"')+idx_start-2
        dummy = dummy(idx_start:idx_end)
        ! print *, dummy
      end if
      read(dummy,*,iostat=ios) a_tmp
      a = a_tmp
    end if
  end subroutine readXML_inttag

  subroutine readXML_3realval(a,tag,buffer,counter)
    implicit none
    real(kind=dp),     intent(inout) :: a(:,:)
    character(len=*),  intent(in)    :: buffer, tag
    integer,           intent(inout) :: counter
    
    character(len=200) :: dummy
    integer            :: idx_start, idx_end
    real(kind=dp)      :: a_tmp(3)

    if (index(buffer, tag) /= 0) then
      read(buffer,'(A)') dummy
      ! print *, dummy
      idx_start = index(dummy,'>')+1
      idx_end = index(dummy,'</')-1
      dummy = dummy(idx_start:idx_end)
      read(dummy,*) a_tmp
      counter = counter + 1
      a(:,counter) = a_tmp
    end if
  end subroutine readXML_3realval


  subroutine readXML_latt(a,tag,buffer,counter)
    real(kind=dp),    intent(inout) :: a(3,3)
    character(len=*), intent(in)    :: buffer, tag
    integer,          intent(inout) :: counter
    
    character(len=200) :: dummy
    integer            :: idx_start, idx_end
    real(kind=dp)      :: a_tmp(3)

    if (index(buffer, tag) /= 0) then
      ! print *, buffer
      counter = counter + 1
    else if (counter > 0 .and. counter < 4) then
      counter = counter + 1
      read(buffer,'(A)') dummy
      idx_start = index(dummy,'>')+1
      idx_end = index(dummy,'</')-1
      dummy = dummy(idx_start:idx_end)
      ! print *, dummy
      read(dummy,*) a_tmp(:)
      a(counter-1,:) = a_tmp(:)
    end if
  end subroutine readXML_latt

  subroutine readXML_rot(a,tag,buffer,counter,counter2)
    real(kind=dp),       intent(inout) :: a(:,:,:)
    character(len=200),  intent(in)    :: buffer
    character(len=*) ,   intent(in)    :: tag
    integer,             intent(inout) :: counter
    integer,             intent(in)    :: counter2

    integer       :: Ri 
    real(kind=dp) :: a_tmp(3)

    Ri = counter2 + 1

    if (index(buffer, tag) /= 0) then
      counter = counter + 1
    else if (counter > 0 .and. counter < 4) then
      read(buffer,*) a_tmp
      a(1:3,counter,Ri) = a_tmp(1:3)
      counter = counter + 1
    end if
  end subroutine readXML_rot


  subroutine readXML_eigenvals(vals,tag,buffer,counter,ik)
    real(kind=dp),    intent(inout) :: vals(:,:) ! Nband x NkI
    character(len=*), intent(in)    :: buffer, tag
    integer,          intent(inout) :: counter, ik
    
    character(len=200)   :: dummy
    integer              :: Nlines, Nband
    integer              :: ie
    real(kind=dp)        :: a_tmp(5) ! dim 5 is the No. of columns
    integer, allocatable :: Nie(:)

    Nband = size(vals,1)
    if (mod(Nband,5)==0) then
      Nlines = Nband/5
    else
      Nlines = Nband/5 + 1
      allocate(Nie(Nlines))
      Nie(1:Nlines-1) = 5
      Nie(Nlines) = mod(Nband,5)
    endif

    if (index(buffer, tag) /= 0) then
      ! print *, buffer
      counter = counter + 1
      ik = ik + 1
    else if (counter > 0 .and. counter <= Nlines) then
      read(buffer,'(A)') dummy
      dummy = trim(dummy)
      ! print *, dummy
      do ie=1,Nie(counter)
        read(dummy,*) a_tmp(1:Nie(counter))
        vals(ie+(counter-1)*5,ik) = a_tmp(ie)
      enddo
      counter = counter + 1
    else if(counter > Nlines) then
      counter = 0
    end if

  if (allocated(Nie)) deallocate(Nie)
  end subroutine readXML_eigenvals

end module io_xml