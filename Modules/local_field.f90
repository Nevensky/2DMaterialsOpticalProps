module local_field
  use iso_fortran_env, only: dp => real64
  implicit none

  public :: genGlf
  private

contains

  subroutine genGlf(lf,Ecut,Gcar,G,Nlf,Nlfd,Glf,parG)
    !! Generate Reciprocal latt. vectors for crystal local field 
    !! effects calculations in array Glf(1:3,1:Nlf)
    implicit none
    character(len=*), intent(in)    :: lf
    real(kind=dp),    intent(in)    :: Ecut
    real(kind=dp),    intent(in)    :: Gcar
    real(kind=dp),    intent(in)    :: G(:,:)
    integer,          intent(inout) :: Nlfd
    integer,          intent(out)   :: Nlf
    real(kind=dp),    intent(out)   :: Glf(:,:)
    integer,optional, intent(inout) :: parG(:)
  
    integer       :: iG
    real(kind=dp) :: Eref
    integer       :: NG = size(G,2)
  
    if (.not. allocated(Glf)) allocate(Glf(3,NG*Nlf))
    Nlf = 0
    if (lf == 'z') then
      ! local field effects included only in the perpendicular dirrection (z)
      do iG = 1, NG
        if (G(1,iG) == 0.0_dp .and. G(2,iG) == 0.0_dp) then
          Eref = Gcar**2 * G(3,iG)**2 / 2.0_dp
          if (Eref <= Ecut) then
            Nlf = Nlf + 1
            Glf(1:2,Nlf) = 0.0_dp
            Glf(3,Nlf) = G(3,iG)
            
            ! store parity of G-vectors if required
            if (present(parG)) then
              if ( (parG(iG)/2)*2 == parG(iG) ) then
                parG(Nlf) = 1
              else
                parG(Nlf) = -1
              end if
            end if

          end if
        end if
      end do
    else if (lf == 'xyz') then
      ! local fiel effects included in all directions (xyz)
      do  iG = 1, NG
        Eref = Gcar**2 * sum(G(1:3,iG)**2) / 2.0_dp
        if (Eref <= Ecut) then
          Nlf = Nlf+1
          Glf(1:3,Nlf) = G(1:3,iG)

          ! store parity of G-vectors if required
          if (present(parG)) then
            if ( (parG(iG)/2)*2 == parG(iG) ) then
              parG(Nlf) = 1
            else
              parG(Nlf) = -1
            end if
          end if

        end if
      end do
    else
      stop 'ERROR: Unsupported local-field vector orientation in genGlf()'
    end if
    if (Nlf > Nlfd) then
      stop 'Nlf is bigger than Nlfd'
    else if(Nlf<Nlfd) then
      Nlfd = Nlf
    end if
  end subroutine genGlf

end module local_field