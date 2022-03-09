module constants
    use iso_fortran_env, only: dp => real64
    implicit none

    save
    real(kind=dp),    parameter :: pi      = 4.D0*atan(1.d0)
    real(kind=dp),    parameter :: eV      = 1.602176487D-19
    real(kind=dp),    parameter :: kB      = 1.3806503d-23
    real(kind=dp),    parameter :: Hartree = 2.0D0*13.6056923D0
    real(kind=dp),    parameter :: Planck  = 6.626196D-34
    real(kind=dp),    parameter :: aBohr   = 0.5291772d0
    real(kind=dp),    parameter :: gamma   = 1.0d0/137.0d0

end module constants