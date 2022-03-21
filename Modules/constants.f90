module constants
    use iso_fortran_env, only: dp => real64
    implicit none

    save
    real(kind=dp),    parameter :: pi      = 4.0_dp*atan(1.0_dp)
    real(kind=dp),    parameter :: eV      = 1.602176634e-19_dp            ! J
    real(kind=dp),    parameter :: kB      = 1.380649e-23_dp               ! J K^-1 
    ! real(kind=dp),    parameter :: Hartree = 2.0_dp*13.6056923_dp        ! eV
    real(kind=dp),    parameter :: Hartree = 4.3597447222071E-18_dp/eV     ! J -> eV
    real(kind=dp),    parameter :: Rydberg = 4.3597447222071E-18_dp/(2*eV) ! J -> Ry
    real(kind=dp),    parameter :: Planck  = 6.62607015e-34_DP             ! J s
    real(kind=dp),    parameter :: aBohr   = 0.529177210903e-10_dp         ! m
    real(kind=dp),    parameter :: alpha   = 0.0072973525693_dp            ! 1/137

end module constants