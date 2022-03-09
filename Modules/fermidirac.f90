module FermiDirac
  use iso_fortran_env, only: dp => real64
  implicit none

  public :: FermiDirac
  private

contains

  pure function FermiDirac(E,Ef,T) result(Ni)
      ! Generates Fermi Dirac distribution
      ! and returns No. of fermions Ni 
      ! with energy E at temperature T

      real(kind=dp), intent(in)  :: E, Ef, T
      real(kind=dp)              :: Ni, expo

      ! real(kind=dp), parameter :: kB = 1.3806503D-23
      real(kind=dp), parameter :: kB = 1       ! Ei is in E/kB units
      real(kind=dp), parameter :: cutoff = 6.0 ! 0.99->1. 

      expo = (Ei-Ef)/(kB*T)

      if (expo < -cutoff) then
        Ni = 1.0
      else if (expo > cutoff)
        Ni = 0.0
      else
        Ni = 1/(exp( expo ) + 1)
      endif
      
  end function FermiDirac

end module FermiDirac