  module const

    implicit none
    
    integer, parameter :: sp = selected_real_kind(p=5)
    integer, parameter :: dp = selected_real_kind(p=15)

     !index for each data type
     integer, parameter :: VEGAZP=0
     integer, parameter :: BB=1
     integer, parameter :: CK2003=2
     integer, parameter :: PHOENIX=3
     integer, parameter :: ATLAS_spec=4
     integer, parameter :: ATLAS_SED=5
     integer, parameter :: ATLAS_flux=9
     integer, parameter :: RAUCH=6
     integer, parameter :: KOESTER=7     
     integer, parameter :: PHOENIX_ACES=8

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
     integer, parameter :: DECam = 17
     integer, parameter :: GALEX = 18
     integer, parameter :: JHKLM = 19
     integer, parameter :: JWST  = 20
     integer, parameter :: WFIRST = 21
     integer, parameter :: RoboAO = 22
     integer, parameter :: Phill = 23

     !index for types of zeropoints
     integer, parameter :: zero_point_unknown = 0
     integer, parameter :: zero_point_Vega    = 1   
     integer, parameter :: zero_point_AB      = 2
     integer, parameter :: zero_point_ST      = 3
     character(len=2), parameter :: AB_string = 'AB'
     character(len=4), parameter :: Vega_string = 'Vega'
     character(len=2), parameter :: ST_string = 'ST'

    ! important files: for Vega, 005 is newer than 004
    ! also note that Vega/alpha_lyr*.dat contain Vega flux in 1st column, AB in 2nd, ST in 3rd
    character(len=256), parameter :: vega_filename = 'Vega/alpha_lyr_stis_005.dat'
    !character(len=256), parameter :: vega_filename = 'Vega/alpha_lyr_stis_004.dat'
    
    ! all constants' units are cgs
    real(dp), parameter :: Lsun = 3.8418d33, Msun=1.9891d33
    real(dp), parameter :: G = 6.67384d-8 ! 2010 CODATA
    real(dp), parameter :: sigma = 5.670373d-5 ! 2010 CODATA
    real(dp), parameter :: clight=2.99792458d18      ! speed of light in Angstrom/s
    real(dp), parameter :: kB=1.3806488d-16 !erg/K, 2010 CODATA
    real(dp), parameter :: planck_h = 6.62606957d-27 ! Planck's constant (erg s)
    real(dp), parameter :: pc = 3.0856770322224d18, pc10=1d1*pc ! 10 parsec
    real(dp), parameter :: pi = 3.141592653589793d0, pi4=4d0*pi
    real(dp), parameter :: AAcm = 1d8  ! AA per cm
    real(dp), parameter :: phx_flux_conv = 1d0/AAcm   ! PHOENIX unit conversion
    real(dp), parameter :: WD_flux_conv = pi/AAcm !Rauch/Koester unit conversion
    real(dp), parameter :: ST_flux_const = 3.631d-09 ! erg/s/cm^2/AA
    real(dp), parameter :: AB_flux_const = 3.631d-20 ! erg/s/cm^2/Hz
    !from G. Torres 2015, IAU resolution to define flux,Mbol standards
    real(dp), parameter :: Mbol_ref = 0d0
    real(dp), parameter :: flux_ref = 3.012d35 !erg/s, note that this implies Lsun=3.828d33
    !!the old version, with no reference, difference between these two is 0.0142 mag
    !!real(dp), parameter :: Mbol_ref = 4.75d0
    !!real(dp), parameter :: flux_ref=3.8418d33 !Lsun
    real(dp), parameter :: flux_ref_const = flux_ref / (pi4 * pc10 * pc10) ! for BC


  end module const
