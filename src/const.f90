  module const

    implicit none
    
    integer, parameter :: sp = selected_real_kind(p=5)
    integer, parameter :: dp = selected_real_kind(p=15)

     !index for each data type
     integer, parameter :: VEGAZP=0
     integer, parameter :: PHOENIX=1
     integer, parameter :: CK2003=2
     integer, parameter :: ATLAS_spec=3
     integer, parameter :: ATLAS_SED=4
     integer, parameter :: RAUCH=5
     integer, parameter :: BB=6

     !index for each photometric system
     integer, parameter :: HST_WFC3 = 1
     integer, parameter :: HST_ACS_WFC = 2
     integer, parameter :: HST_ACS_HRC = 3
     integer, parameter :: HST_WFPC2 = 4
     integer, parameter :: SDSS = 5
     integer, parameter :: CFHT = 6
     integer, parameter :: UBVRIJHKsKp = 7
     integer, parameter :: WashDDOuvby = 8
     integer, parameter :: UKIDSS = 9 
     integer, parameter :: WISE = 10
     integer, parameter :: PanSTARRS = 11
     integer, parameter :: SPITZER = 12
     integer, parameter :: SkyMapper = 13
     integer, parameter :: LSST = 14
     integer, parameter :: Swift = 15
     integer, parameter :: FSPS  = 16

     !index for types of zeropoints
     integer, parameter :: zero_point_Vega = 1
     integer, parameter :: zero_point_AB = 2
     integer, parameter :: zero_point_ST = 3

    ! important files: for Vega, 005 is newer than 004
    ! also note that Vega/alpha_lyr*.dat contain Vega flux in 1st column, AB in 2nd, ST in 3rd
    character(len=256), parameter :: vega_filename = 'Vega/alpha_lyr_stis_005.dat'
    !character(len=256), parameter :: vega_filename = 'Vega/alpha_lyr_stis_004.dat'
    
    ! all constants' units are cgs
    real(dp), parameter :: Msolbol = 4.75d0, Lsun = 3.8418d33, Msun=1.9891d33
    real(dp), parameter :: G = 6.67384d-8 ! 2010 CODATA
    real(dp), parameter :: sigma = 5.670373d-5 ! 2010 CODATA
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
    real(dp), parameter :: kB=1.3806488d-16 !erg/K, 2010 CODATA

  end module const
