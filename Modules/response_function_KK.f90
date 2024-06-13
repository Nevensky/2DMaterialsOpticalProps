module response_function_KK
  use iso_fortran_env, only: dp => real64
  use constants, only: pi
  implicit none

  public :: genReChi0, genImChi0
  private

  ! WARNING the frequency (io) loop outside the subroutines must be from -No+1 do No-1
  ! instead of -No to No => NEEDS TO BE FIXED

  ! interface infers if correlation function S0 is real or complex
  ! and calls the appropriate _real or _cmplx subroutine
  interface genReChi0
      module procedure genReChi0_cmplx
      module procedure genReChi0_real
  end interface

  interface genImChi0
      module procedure genImChi0_cmplx
      module procedure genImChi0_real
  end interface

contains
    subroutine genReChi0_real(io, No, Nlf, iG, jG, oi, domega, S0, ReChi0)
        !! Wrapper for genChi0 returnign only the real part.
        implicit none
        integer,          intent(in)  :: io, No, Nlf
        integer,          intent(in)  :: iG, jG
        real(kind=dp),    intent(in)  :: oi, domega
        real(kind=dp),    intent(in)  :: S0(-No:No, Nlf, Nlf)
        real(kind=dp),    intent(out) :: ReChi0

        call genChi0_S0real(io, No, Nlf, iG, jG, oi, domega, S0, ReChi0, .true.)
    end subroutine genReChi0_real

    subroutine genReChi0_cmplx(io, No, Nlf, iG, jG, oi, domega, S0, ReChi0)
        !! Wrapper for genChi0 returnign only the real part.
        implicit none
        integer,          intent(in)  :: io, No, Nlf
        integer,          intent(in)  :: iG, jG
        real(kind=dp),    intent(in)  :: oi, domega
        complex(kind=dp), intent(in)  :: S0(-No:No, Nlf, Nlf)
        real(kind=dp),    intent(out) :: ReChi0

        call genChi0_S0complex(io, No, Nlf, iG, jG, oi, domega, S0, ReChi0, .true.)
    end subroutine genReChi0_cmplx

    subroutine genImChi0_real(io, No, Nlf, iG, jG, oi, domega, S0, ImChi0)
        !! Wrapper for genChi0 returnign only the real part.
        implicit none
        integer,          intent(in)  :: io, No, Nlf
        integer,          intent(in)  :: iG, jG
        real(kind=dp),    intent(in)  :: oi, domega
        real(kind=dp),    intent(in)  :: S0(-No:No, Nlf, Nlf)
        real(kind=dp),    intent(out) :: ImChi0

        call genChi0_S0real(io, No, Nlf, iG, jG, oi, domega, S0, ImChi0, .false.)
    end subroutine genImChi0_real

    subroutine genImChi0_cmplx(io, No, Nlf, iG, jG, oi, domega, S0, ImChi0)
        !! Wrapper for genChi0 returnign only the real part.
        implicit none
        integer,          intent(in)  :: io, No, Nlf
        integer,          intent(in)  :: iG, jG
        real(kind=dp),    intent(in)  :: oi, domega
        complex(kind=dp), intent(in)  :: S0(-No:No, Nlf, Nlf)
        real(kind=dp),    intent(out) :: ImChi0

        call genChi0_S0complex(io, No, Nlf, iG, jG, oi, domega, S0, ImChi0, .false.)
    end subroutine genImChi0_cmplx

    ! for optical calc
    subroutine genChi0_S0complex(io, No, Nlf, iG, jG, oi, domega, S0, Chi0, isChiReal)
        !! Computes the real or imaginary part of the response function \( \Chi^0_{G,G'}(\omega)) \)
        implicit none
        integer,          intent(in)  :: io, No, Nlf
        integer,          intent(in)  :: iG, jG
        real(kind=dp),    intent(in)  :: oi, domega
        complex(kind=dp), intent(in)  :: S0(-No:No, Nlf, Nlf)
        real(kind=dp),    intent(out) :: Chi0
        logical,          intent(in)  :: isChiReal

        integer :: jo
        real(kind=dp) :: oj, fact
        real(kind=dp) :: S0_value

        if (isChiReal) then
            Chi0 = pi * aimag(S0(io - 1, iG, jG))
        else
            Chi0 = - pi * dble(S0(io - 1, iG, jG))
        endif
        if (io == 1) then
            do jo = 2, No
                oj = (jo - 1) * domega
                fact = domega / oj
                if (jo == 2)    fact = 3.0_dp / 2.0_dp
                if (jo == No)   fact = 0.5_dp * domega / oj
                if (isChiReal) then
                    S0_value = dble(S0(-jo + 1, iG, jG)) - dble(S0(jo - 1, iG, jG))
                else
                    S0_value = aimag(S0(-jo + 1, iG, jG)) - aimag(S0(jo - 1, iG, jG))
                end if
                Chi0 = Chi0 + fact * S0_value
            end do
        elseif (io == 2) then
            do jo = 1, No
                oj   = (jo - 1) * domega
                if (jo /= io)  fact = domega / (oi - oj)
                if (jo == 1)   fact = 1.0_dp
                if (jo == 2)   fact = 0.0_dp
                if (jo == 3)   fact = -3.0_dp / 2.0_dp
                if (jo == No)  fact = 0.5_dp * domega / (oi - oj)
                if (isChiReal) then
                    S0_value = dble(S0(jo - 1, iG, jG))
                else
                    S0_value = aimag(S0(jo - 1, iG, jG))
                end if
                Chi0 = Chi0 + fact * S0_value
                fact = domega / (oi + oj)
                if (jo == 1 .or. jo == No) fact = 0.5_dp * domega / (oi + oj)
                if (isChiReal) then
                    S0_value = dble(S0(-jo + 1, iG, jG))
                else
                    S0_value = aimag(S0(-jo + 1, iG, jG))
                end if
                Chi0 = Chi0 + fact * S0_value
            end do
        elseif (io == (No - 1)) then
            do jo = 1, No
                oj   = (jo - 1) * domega
                if (jo /= io)       fact = domega / (oi - oj)
                if (jo == 1)        fact = 0.5_dp * domega / (oi - oj)
                if (jo == (No - 2)) fact = 3.0_dp / 2.0_dp
                if (jo == (No - 1)) fact = 0.0_dp
                if (jo == No)       fact = -1.0_dp
                if (isChiReal) then
                    S0_value = dble(S0(jo - 1, iG, jG))
                else
                    S0_value = aimag(S0(jo - 1, iG, jG))
                end if
                Chi0 = Chi0 + fact * S0_value
                fact = domega / (oi + oj)
                if (jo == 1 .or. jo == No) fact = 0.5_dp * domega / (oi + oj)
                if (isChiReal) then
                    S0_value = dble(S0(-jo + 1, iG, jG))
                else
                    S0_value = aimag(S0(-jo + 1, iG, jG))
                end if
                Chi0 = Chi0 + fact * S0_value
            end do
        else
            do jo = 1, No
                oj   = (jo - 1) * domega
                if (jo /= io)       fact = domega / (oi - oj)
                if (jo == 1)        fact = 0.5_dp * domega / (oi - oj)
                if (jo == (io - 1)) fact = 3.0_dp / 2.0_dp
                if (jo == io)       fact = 0.0_dp
                if (jo == (io + 1)) fact = -3.0_dp / 2.0_dp
                if (jo == No)       fact = 0.5_dp * domega / (oi - oj)
                if (isChiReal) then
                    S0_value = dble(S0(jo - 1, iG, jG))
                else
                    S0_value = aimag(S0(jo - 1, iG, jG))
                end if
                Chi0 = Chi0 + fact * S0_value
                fact = domega / (oi + oj)
                if (jo == 1 .or. jo == No) fact = 0.5_dp * domega / (oi + oj)
                if (isChiReal) then
                    S0_value = dble(S0(-jo + 1, iG, jG))
                else
                    S0_value = aimag(S0(-jo + 1, iG, jG))
                end if
                Chi0 = Chi0 + fact * S0_value
            end do
        endif

    end subroutine genChi0_S0complex

  ! for Sloss / W calc
  subroutine genChi0_S0real(io, No, Nlf, iG, jG, oi, domega, S0, Chi0, isChiReal)
    implicit none

    integer,          intent(in)  :: io, No, Nlf
    integer,          intent(in)  :: iG, jG
    real(kind=dp),    intent(in)  :: oi, domega
    real(kind=dp),    intent(in)  :: S0(-No:No, Nlf, Nlf)
    real(kind=dp),    intent(out) :: Chi0
    logical,          intent(in)  :: isChiReal
    integer       :: jo
    real(kind=dp) :: oj, fact
    real(kind=dp) :: ReChi0 !, ImChi0

    if (isChiReal) then
        ! oi = (io-1)*domega
        ! iG_loop: do iG = 1,Nlf ! moved outside subroutine
        ! jG_loop: do jG = 1,Nlf ! moved outside subroutine
        ReChi0 = 0.0_dp
        ! static limit
        if (io == 1) then
          do jo = 2,No
            oj = (jo-1)*domega
            
            fact = domega/oj
            if (jo == 2)  fact = 3.0_dp/2.0_dp
            if (jo == No) fact = 0.5_dp*domega/oj
            
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
          end do
          ReChi0 = -2.0_dp * ReChi0
        else if (io == 2) then
          do jo = 1,No
            oj = (jo-1)*domega
            if (jo /= io) then
              fact = domega/(oi-oj)
            else if (jo == 1) then
              fact = 1.0_dp
            else if (jo == 2) then
              fact = 0.0_dp
            else if (jo == 3) then
              fact = -3.0_dp/2.0_dp
            else if (jo == No) then
              fact = 0.5_dp*domega/(oi-oj)
            end if
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
            fact = domega/(oi+oj)
            if (jo == 1 .or. jo == No) then
              fact=0.5_dp*domega/(oi+oj)
            end if
            ReChi0 = ReChi0 - fact*S0(jo,iG,jG)
          end do
        else if (io == (No - 1)) then
          do  jo = 1,No
            oj = (jo - 1)*domega
            if (jo /= io)       fact = domega/(oi-oj)
            if (jo == 1)        fact = 0.5_dp*domega/(oi-oj)
            if (jo == (No - 2)) fact = 3.0_dp/2.0_dp
            if (jo == (No - 1)) fact = 0.0_dp
            if (jo == No)       fact = -1.0_dp

            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
            
            fact = domega/(oi + oj)
            if (jo == 1 .or. jo == no) fact = 0.5_dp*domega/(oi + oj)

            ReChi0 = ReChi0 - fact*S0(jo,iG,jG)
          end do
        else
          do  jo = 1,No
            oj = (jo - 1)*domega
            if (jo /= io)     fact = domega/(oi - oj)
            if (jo == 1)      fact = 0.5_dp*domega/(oi - oj)
            if (jo == (io-1)) fact = 3.0_dp/2.0_dp
            if (jo == io)     fact = 0.0_dp
            if (jo == (io+1)) fact = -3.0_dp/2.0_dp
            if (jo == No)     fact = 0.5_dp*domega/(oi - oj)
            
            ReChi0 = ReChi0 + fact*S0(jo,iG,jG)
            
            fact = domega/(oi + oj)
            if (jo == 1 .or. jo == No) fact = 0.5_dp*domega/(oi + oj)

            ReChi0 = ReChi0 - fact*S0(jo,iG,jG)
          end do
        end if
        Chi0 = ReChi0
        ! end do jG_loop
        ! end do iG_loop
    else
        Chi0 = -pi*S0(io,iG,jG) ! = ImChi0
    end if

end subroutine genChi0_S0real

end module response_function_KK