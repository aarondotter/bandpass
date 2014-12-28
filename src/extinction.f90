      module extinction
      
      ! provides the Rv-dependent extinction curve from Rv = Av/E(B-V)

      ! Cardelli, Clayton, & Mathis (1989, ApJ, 345, 245-256)

      ! Al_div_Av_CCM provides <A(lambda)/A(V)> as a function of input 
      ! lambda and Rv, with lambda in Angstroms.

      implicit none

      private

      integer, parameter :: sp = selected_real_kind(p=5)
      integer, parameter :: dp = selected_real_kind(p=15)

      public :: Al_div_Av_CCM

      contains

      !elemental: can be called on a scalar or an array in the same way
      elemental function Al_div_Av_CCM(lambda,Rv) result(Alam)
         real(dp), intent(in) :: lambda !in angstrom
         real(dp), intent(in) :: Rv ! Av/E(B-V)
         real(dp) :: a, b, x, y, Alam
         real(dp), parameter :: small = 1d-8, min_x = 0d0 ! 2d-1
        
         ! x is wavenumber in 1/micron : x = 10^4/lambda
         ! impose an upper bound in the UV
         x = min(1d1,1d4/lambda)

         ! CCM break up the extinction curve into 4 regions in x
         if (x>=min_x .and. x<= 1.1d0) then ! Infrared
         ! strictly speaking, CCM imposes a limit on x >= 0.2 but we ignore this here
         ! in order to allow for mid-IR wavelengths.
            a = 0.574d0*x**1.61d0
            b = -0.527d0*x**1.61d0
         else if ( x> 1.1d0 .and. x<=3.3d0) then !Optical
            y = x - 1.82d0            
            a = 1d0 + 1.7699d-1*y - 5.0447d-1*y**2 - 2.427d-2*y**3 &
            + 7.2085d-1*y**4 + 1.979d-2*y**5 - 7.7530d-1*y**6 + 3.2999d-1*y**7
            b = 1.41338d0*y + 2.28305d0*y**2 + 1.07233d0*y**3 &
            - 5.38434d0*y**4 - 6.2251d-1*y**5 + 5.30260d0*y**6 - 2.09002d0*y**7 
         else if (x> 3.3d0 .and. x<8d0) then !Ultraviolet
            a = 1.752d0 - 3.16d-1*x - 1.04d-1/((x-4.67d0)**2 + 3.41d-1) + Fa(x)
            b = -3.090d0 + 1.825d0*x + 1.206d0/( (x-4.62d0)**2 + 2.63d-1) + Fb(x)
         else if (x> 8d0 .and. x <= 1d1) then !Far-Ultraviolet
            y = x - 8d0
            a = -1.073d0 - 6.28d-1*y + 1.37d-1 * y**2 - 7d-2*y**3
            b = 13.67d0 + 4.257d0*y - 4.2d-1*y**2 + 3.74d-1*y**3
         else !return zero if outside defined range
            a = 0d0
            b = 0d0
         endif
         Alam = a + b/Rv   
      end function Al_div_Av_CCM
         
      !special functions Fa and Fb for UV region
      elemental function Fa(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         if(x<5.9d0)then
            y = 0d0
         else
            y = -4.473d-2*(x-5.9d0)**2 - 9.779d-3*(x-5.9d0)**3
         endif
      end function Fa

      elemental function Fb(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         if(x<5.9d0)then
            y = 0d0
         else
            y = 2.13d-1*(x-5.9)**2 + 1.207d-1*(x-5.9d0)**3
         endif
      end function Fb

      end module extinction
