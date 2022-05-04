!! Functions for generating Fermi-Dirac, Bose-Einstein or Boltzmann occupation factors.

module statistics
  use iso_fortran_env, only: dp => real64
  implicit none

  public :: FermiDirac, BoseEinstein, Boltzmann
  private

contains
  pure function FermiDirac(Ei,Efermi,T,cutoff) result(Ni)
      !! Generates Fermi-Dirac distribution
      !! and returns No. of fermions Ni 
      !! with energy E at temperature T

      real(kind=dp), intent(in)            :: Ei, Efermi, T
      real(kind=dp), intent(in), optional  :: cutoff

      real(kind=dp) :: Ni, expo

      ! real(kind=dp), parameter :: kB = 1.3806503D-23
      real(kind=dp), parameter :: kB = 1.0_dp       ! Ei is in E/kB units
      real(kind=dp) :: cutoff_
      if (present(cutoff)) then
        cutoff_ = cutoff
      else
        cutoff_ = 6.0_dp ! guarantees less than 1% error 
      end if

      expo = (Ei-Efermi)/(kB*T)

      if (expo < -cutoff_) then
        Ni = 1.0_dp
      elseif (expo > cutoff_) then
        Ni = 0.0_dp
      else
        Ni = 1.0_dp/(exp(expo) + 1.0_dp)
      endif
      
  end function FermiDirac

  pure function BoseEinstein(Ei,Efermi,T,cutoff) result(Ni)
      !! Generates Bose-Einstein distribution
      !! and returns No. of bosons Ni 
      !! with energy E at temperature T

      real(kind=dp), intent(in)            :: Ei, Efermi, T
      real(kind=dp), intent(in), optional  :: cutoff

      real(kind=dp) :: Ni, expo

      ! real(kind=dp), parameter :: kB = 1.3806503D-23
      real(kind=dp), parameter :: kB = 1.0_dp       ! Ei is in E/kB units
      real(kind=dp) :: cutoff_
      if (present(cutoff)) then
        cutoff_ = cutoff
      else
        cutoff_ = 6.0_dp ! guarantees less than 1% error 
      end if

      expo = (Ei-Efermi)/(kB*T)

      if (expo < -cutoff_) then
        Ni = 1.0_dp
      elseif (expo > cutoff_) then
        Ni = 0.0_dp
      else
        Ni = 1.0_dp/(exp(expo) - 1.0_dp)
      endif
      
  end function BoseEinstein

  pure function Boltzmann(Ei,Efermi,T,cutoff) result(Ni)
      !! Generates Bose-Einstein distribution
      !! and returns No. of particles Ni 
      !! with energy E at temperature T

      real(kind=dp), intent(in)            :: Ei, Efermi, T
      real(kind=dp), intent(in), optional  :: cutoff

      real(kind=dp) :: Ni, expo

      ! real(kind=dp), parameter :: kB = 1.3806503D-23
      real(kind=dp), parameter :: kB = 1.0_dp       ! Ei is in E/kB units
      real(kind=dp) :: cutoff_
      if (present(cutoff)) then
        cutoff_ = cutoff
      else
        cutoff_ = 6.0_dp ! guarantees less than 1% error 
      end if

      expo = (Ei-Efermi)/(kB*T)

      if (expo < -cutoff_) then
        Ni = 1.0_dp
      elseif (expo > cutoff_) then
        Ni = 0.0_dp
      else
        Ni = 1.0_dp/exp(expo)
      endif
      
  end function Boltzmann

end module statistics