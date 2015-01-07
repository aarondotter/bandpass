module blackbody
  
  use const

  implicit none

  real(dp), parameter :: h = planck_h
  real(dp), parameter :: c = clight/AAcm

contains
  
  elemental function BBflux(lam,T)
    real(dp), intent(in) :: lam, T
    real(dp) :: l5, BBflux, l
    l = lam/AAcm
    l5 = l*l*l*l*l
    BBflux = 2d0*h*c*c/(l5*(exp(h*c/(l*kB*T))-1d0))
    BBflux = BBflux*rauch_flux_conv
  end function BBflux

end module blackbody
    
