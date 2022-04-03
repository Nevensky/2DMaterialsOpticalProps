program read_xml
  use iso_fortran_env, only : dp=> real64
  use io_xml, only: loadXML_qe
  implicit none


  integer                    :: Nbands, NkI, Nrot, Nsym, Nmp(3)
  real(kind=dp)              :: alat             !! lattice parameter
  real(kind=dp)              :: a(3,3)           !! direct lattice
  real(kind=dp)              :: b(3,3)           !! reciprocal lattice
  real(kind=dp), allocatable :: R(:,:,:)         !! array of rotational matrices
  integer      , allocatable :: Npw(:)           !! array of G-vectors at each k-point
  real(kind=dp), allocatable :: kI(:,:)          !! array of k-points in the ireducible Brillouin zone
  real(kind=dp), allocatable :: eigenvals(:,:)   !! eigenvalues for each band and k-point
  real(kind=dp), allocatable :: occupations(:,:) !! occupations for each band and k-point
  real(kind=dp) :: T           !! electron temperature
  real(kind=dp) :: Efermi      !! Fermi energy
  real(kind=dp) :: NelQE       !! No. electrons
  character(len=200) :: path   !! path to data_file_schema.xml Quantum Espresso file

  call get_command_argument(1, path)
  print *,path
  if (len(trim(path)) == 0) then
    stop 'Specify path to data_file_schema.xml as first argument.'
  else
    call loadXML_qe(path, Nbands , NkI, Nrot, Nsym, Nmp, alat, a, b, R, Npw, kI, eigenvals, occupations, T, Efermi, NelQE, printOutput=.true.)
  end if

end program read_xml