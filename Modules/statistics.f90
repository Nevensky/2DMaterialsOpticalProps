!! Functions for generating Fermi-Dirac, Bose-Einstein or Boltzmann occupation factors.

module statistics
  use iso_fortran_env, only: dp => real64
  implicit none

  public :: FermiDirac, BoseEinstein, Boltzmann, genOccupation
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

  ! subroutine genOccupation(Ei, Efermi, T, expo, f, eps)
  !   implicit none
  !   real(kind=dp), intent(in)  :: Ei, Efermi, T

  !   real(kind=dp), intent(out) :: expo, f
  !   real(kind=dp), intent(in), optional  :: eps

  !   real(kind=dp) :: eps_ = 20.0_dp
  !   if (present(eps)) eps_=eps

  !   expo = (Ei-Efermi)/T
  !   if (expo < -eps_) then
  !     expo = 0.0_dp
  !     f = 1.0_dp
  !   elseif(expo > eps_) then
  !     f = 0.0_dp
  !   else
  !     expo = exp(expo)
  !     f = 1.0_dp/(expo + 1.0)
  !   endif
  ! end subroutine genOccupation

  subroutine genOccupation(n, m, k1, k2, df_cut, T, Efermi, E, populated, df)
    !! Takes band indicies n, m and determines whether the transition
    !! will take place based on Fermi Dirac distribution for the electron occupation
    integer,       intent(in)  :: n         !! 1st band index
    integer,       intent(in)  :: m         !! 2nd band index
    integer,       intent(in)  :: k1        !! 1st k-point index
    integer,       intent(in)  :: k2        !! 2nd k-point index
    real(kind=dp), intent(in)  :: df_cut    !! cutoff below which a band is considered empty
    real(kind=dp), intent(in)  :: Efermi    !! Fermi energy
    real(kind=dp), intent(in)  :: T         !! electron temperature
    real(kind=dp), intent(in)  :: E(:,:)    !! energies (Nband x NkI)
    logical,       intent(out) :: populated !! band is filled? (true/false)
    real(kind=dp), intent(out), optional :: df

    real(kind=dp) :: f1, f2, df_

    if (T/=0.0_dp) then
      populated = .false.
      f1 = FermiDirac(E(n,K1),Efermi,T)
      f2 = FermiDirac(E(m,K2),Efermi,T)
      df = f1 - f2
      if ((abs(df_) >= df_cut) .or. (n == m)) populated = .true.
    else
      populated= .true.
    end if
    if (present(df)) df = df_
  end subroutine genOccupation

end module statistics