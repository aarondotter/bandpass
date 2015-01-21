program test

  implicit none

  integer :: i
  integer, parameter :: numT = 106
  double precision :: Teff(numT)

  Teff(1) = 1.0d3; Teff(2) = 1.5d3; Teff(3) = 2.0d3; Teff(4)=2.5d3
  Teff(5) = 2.8d3; Teff(6) = 3.0d3; Teff(7) = 3.2d3; Teff(8)=3.5d3
       
  !  3,750 to    13,000
  do i=9,46
     Teff(i) = Teff(i-1) + 2.5d2
  enddo

  ! 14,000 to    50,000
  do i=47,83
     Teff(i) = Teff(i-1) + 1.0d3
  enddo
  
  ! 60,000 to   200,000
  do i=84,98
     Teff(i) = Teff(i-1) + 1.0d4
  enddo

  !200,000 to 1,000,000
  do i=99,numT
     Teff(i) = Teff(i-1) + 1.0d5
  enddo
  
  do i=1,numT
     write(*,*) i, Teff(i)
  enddo

end program test
