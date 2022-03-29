module propagators
  use iso_fortran_env, only: dp => real64, sp => real32
  implicit none

  public :: genGammaPropagator
  private

contains
  subroutine genGammaPropagator(iq, Nq, No, q, domega, D, Gammap, Gammam)
  integer,          intent(in)     :: iq, Nq, No
  real(kind=dp),    intent(in)     :: q, domega
  complex(kind=dp), intent(inout)  :: D(No)
  complex(kind=dp), intent(out)    :: Gammap(No), Gammam(No)

  integer            :: io, jo
  real(kind=dp)      :: oi, oj
  real(kind=dp)      :: W1, W2, ImW, Wind, WindKK, fact
  real(kind=dp)      :: KKS, SKK
  complex(kind=dp)   :: Gammap0

  real(kind=dp)    ::  S(No)

  real(kind=dp), parameter :: pi = 3.141592654_dp

  !  Spectra of surface excitations

  S = -(1.0_dp/pi)*aimag(D)

  !  Construction of the correlation propagator \Gamma(omega)
  KKS = 0
  SKK = 0
  omega_loop: do io = 1, No-1
    oi=(io - 1)*domega
    W1 = 0.0_dp
    W2 = 0.0_dp
    ! static limit
    if(io == 1) then
      do  jo = 2,No
        oj=(jo - 1)*domega
        fact = domega/oj
        if(jo == 2) fact = 3.0_dp/2.0_dp
        if(jo == No)fact = 0.5_dp*domega/oj
        W1 = W1 - fact*S(jo)
      end do
      W2 = -W1
    else if(io == 2) then
      do  jo = 1,No
        oj=(jo - 1)*domega
        if(jo /= io) fact = domega/(oi - oj)
        if(jo == 1) fact = 1.0_dp
        if(jo == 2) fact = 0.0_dp
        if(jo == 3) fact= -3.0_dp/2.0_dp
        if(jo == No) fact = 0.5_dp*domega/(oi - oj)
        W1 = W1 + fact*S(jo)
        fact = domega/(oi + oj)
        if(jo == 1 .or. jo == No) fact = 0.5_dp*domega/(oi + oj)
        W2 = W2 + fact*S(jo)
      end do
    else if(io == (No - 1)) then
      do  jo = 1,No
        oj=(jo - 1)*domega
        if(jo /= io)fact = domega/(oi - oj)
        if(jo == 1)fact = 0.5*domega/(oi - oj)
        if(jo == (No - 2)) fact = 3.0_dp/2.0_dp
        if(jo == (No - 1)) fact = 0
        if(jo == No) fact=-1
        W1 = W1 + fact*S(jo)
        fact = domega/(oi + oj)
        if(jo == 1 .or. jo == No) fact = 0.5_dp*domega/(oi + oj)
        W2 = W2 + fact*S(jo)
      end do
    else
      do  jo = 1,No
        oj=(jo - 1)*domega
        if(jo /= io) fact = domega/(oi - oj)
        if(jo == 1) fact = 0.5_dp*domega/(oi - oj)
        if(jo == (io - 1)) fact = 3.0_dp/2.0_dp
        if(jo == io) fact = 0
        if(jo == (io + 1)) fact= -3.0_dp/2.0_dp
        if(jo == No) fact = 0.5_dp*domega/(oi - oj)
        W1 = W1 + fact*S(jo)
        fact = domega/(oi + oj)
        if(jo == 1 .or. jo == No) fact = 0.5*domega/(oi + oj)
        W2 = W2 + fact*S(jo)
      end do
    end if
    
    
    ImW = -pi*S(io)
    Gammap(io) = cmplx(W1,ImW)
    Gammam(io) = cmplx(-W2,0.0_dp)
    
    if(io == 1) Gammap0 = Gammap(1)
    
    ! Check Kramers - Kroning relations
    Wind = real(D(io))
    WindKK = real(Gammap(io)) + real(Gammam(io))
    fact = domega
    if(io == 1 .or. io == No - 1) fact = 0.5_dp*domega
    KKS = KKS + fact*(WindKK - Wind)**2
    SKK = SKK + fact*Wind**2
    
    
  end do omega_loop

  if(iq == 2 .or. iq == Nq) then
    ! write output of Kramers-Kroning relations check
    call writeKramKron_Qi(iq, q, KKS, SKK, Gammap0, D)
  end if


  end subroutine genGammaPropagator

  subroutine writeKramKron_Qi(iq, q, KKS, SKK, Gammap0, D)
    integer,          intent(in)    :: iq
    real(kind=dp),    intent(in)    :: q, KKS, SKK
    complex(kind=dp), intent(in)    :: D(:)
    complex(kind=dp), intent(in)    :: Gammap0

    integer            :: iuni, nord, ios
    character(len=100) :: dato

    dato='Kramers - Kron_Qi'
    nord = index(dato,'i', back = .false.)
    if(iq < 10) then
      write(dato(nord:nord),'(i1)')iq
    else if(iq >= 10 .and. iq < 100) then
      write(dato(nord:nord + 1),'(i2)') iq
    else
      write(dato(nord:nord + 2),'(i3)') iq
    end if
    
    
    
    open(newunit = iuni, file = dato, iostat=ios)
    if (ios /=0) then
      stop 'ERROR: Failed to write Kramers Kronging output file. (propagators.f90)'
    end if
    write(iuni,'(a25,f8.4,a5)') 'Q=', q, 'a.u.'
    write(iuni,*) 'int(DKK - D)^2 =  ', KKS
    write(iuni,*) 'int(D)^2 =  ', SKK
    write(iuni,*) '****************************************'
    write(iuni,*) 'Kramers–Kronig relation relative error'
    write(iuni,'( 5X,f7.2,a2)') 100.0_dp*abs(KKS/SKK), '%'
    write(iuni,*) 'Usporedba Gamma i D'
    write(iuni,*) 'real[Gamma(o = 0)]=', real(Gammap0)
    write(iuni,*) 'real[D(o = 0)]/2=', real(D(1))/2.0_dp
    close(iuni)

  end subroutine writeKramKron_Qi

end module propagators
