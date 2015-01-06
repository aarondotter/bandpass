  module const

    implicit none
    
    integer, parameter :: sp = selected_real_kind(p=5)
    integer, parameter :: dp = selected_real_kind(p=15)

    ! important files: for Vega, 005 is newer than 004
    ! also note that Vega/alpha_lyr*.dat contain Vega flux in 1st column, AB in 2nd, ST in 3rd
    character(len=256), parameter :: vega_filename = 'Vega/alpha_lyr_stis_005.dat'
    !character(len=256), parameter :: vega_filename = 'Vega/alpha_lyr_stis_004.dat'
    
    ! all constants' units are cgs
    real(dp), parameter :: Msolbol = 4.75d0, Lsun = 3.8418d33, Msun=1.9891d33
    real(dp), parameter :: G = 6.67384d-8, sigma = 5.670373d-5 ! 2010 CODATA
    real(dp), parameter :: clight=2.99792458d18      ! speed of light in Angstrom/s
    real(dp), parameter :: planck_h = 6.62606957d-27 ! Planck's constant (erg s)
    real(dp), parameter :: pc = 3.0856770322224d18, pc10=1d1*pc ! 10 parsec
    real(dp), parameter :: pi = 3.141592653589793d0, pi4=4d0*pi
    real(dp), parameter :: AAcm = 1d8  ! AA per cm
    real(dp), parameter :: phx_flux_conv = 1d0/AAcm   ! PHOENIX unit conversion
    real(dp), parameter :: rauch_flux_conv = pi/AAcm !Rauch unit conversion
    real(dp), parameter :: ST_flux_const = 3.631d-09 ! erg/s/cm^2/AA
    real(dp), parameter :: AB_flux_const = 3.631d-20 ! erg/s/cm^2/Hz
    real(dp), parameter :: solar_const = Lsun / (pi4 * pc10 * pc10) ! for BC

  end module const
