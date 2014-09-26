      module bandpass

      use extinction
      
      implicit none

      ! important files: for Vega, 005 is newer than 004
      ! also note that Vega/alpha_lyr*.dat contain Vega flux in 1st column, AB in 2nd, ST in 3rd
      character(len=256), parameter :: vega_filename = 'Vega/alpha_lyr_stis_005.dat'
      !character(len=256), parameter :: vega_filename = 'Vega/alpha_lyr_stis_004.dat'
      character(len=256) :: data_dir
      logical :: data_dir_set = .false.

      integer, parameter :: VEGAZP=0, PHOENIX=1, CK2003=2, ATLAS_spec=3, ATLAS_SED=4

      ! all constants' units are cgs
      double precision, parameter :: Msolbol = 4.75d0, Lsun = 3.8418d33, Msun=1.9891d33
      double precision, parameter :: G = 6.67428d-8, sigma = 5.6704d-5 ! 2006 CODATA
      double precision, parameter :: clight=2.99792458d18      ! speed of light in Angstrom/s
      double precision, parameter :: planck_h = 6.62606896d-27 ! Planck's constant (erg s)
      double precision, parameter :: pc = 3.0856770322224d18, pc10=10d0*pc ! 10 parsec
      double precision, parameter :: pi4 = 1.6d1*atan(1d0)
      double precision, parameter :: AAcm = 1d8  ! AA per cm
      double precision, parameter :: phx_flux_conv = AAcm   ! PHOENIX unit conversion
      double precision, parameter :: ST_flux_const = 3.631d-09 ! erg/s/cm^2/AA
      double precision, parameter :: AB_flux_const = 3.631d-20 ! erg/s/cm^2/Hz
      double precision, parameter :: const = Lsun / (pi4 * pc10 * pc10) ! for BC
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
      integer, parameter :: zero_point_Vega = 1, zero_point_AB = 2, zero_point_ST = 3
      integer :: zero_point_type
      logical, parameter :: debug = .false., read_on_the_fly = .true., do_check_total_flux=.false.

      type spectrum
         character(len=256) :: filename
         integer :: npts
         double precision :: FeH, Teff, logg, R, M, Fbol, alpha_Fe ! R in cm, M in Msun
         double precision, allocatable :: wave(:), flux(:), extinction(:)
      end type spectrum                              

      contains

      subroutine set_data_dir(my_data_dir)
      character(len=256), intent(in) :: my_data_dir
      data_dir = my_data_dir
      data_dir_set = .true.
      end subroutine set_data_dir

      subroutine write_bin_file( s, binfile,ierr)
      type(spectrum), intent(in) :: s
      character(len=256), intent(in) :: binfile
      integer, intent(out) :: ierr
      integer :: io=69
      ierr=0
      open(io,file=trim(binfile),form='unformatted',iostat=ierr)
      if(ierr/=0) then
         ierr=-1
         return
      endif
      write(io) s% filename, s% npts, s% feh, s% Teff, s% logg, s% R, s% M, s% Fbol
      write(io) s% wave, s% flux
      close(io)
      write(0,*) '   wrote binary file ' , trim(binfile)
      end subroutine write_bin_file

      subroutine read_vega(vega,AB,ST,ierr)
      type(spectrum), intent(inout) :: vega, AB, ST
      integer, intent(out) :: ierr
      character(len=256) :: filename
      integer :: w

      ierr=0
      filename = trim(data_dir) // '/' // trim(vega% filename)
      write(0,*)"reading: ",trim(filename)
      open(1,file=trim(filename),iostat=ierr)
      if(ierr/=0) return
      read(1,'(i6)',iostat=ierr) vega% npts
      if(ierr/=0) return
      AB% npts = vega% npts
      ST% npts = vega% npts
      allocate(vega% wave(vega% npts),vega% flux(vega% npts), vega% extinction(vega% npts))
      allocate(AB% wave(vega% npts), AB% flux(vega% npts), AB% extinction(AB% npts))
      allocate(ST% wave(vega% npts), ST% flux(vega% npts), ST% extinction(ST% npts))
      do w=1,vega% npts
         read(1,*) vega% wave(w), vega% flux(w), AB% flux(w), ST% flux(w)
      enddo
      close(1)
      AB% wave = vega% wave
      ST% wave = vega% wave
      vega% extinction = 1d0
      AB% extinction = 1d0
      ST% extinction = 1d0
      end subroutine read_vega

      double precision function STflux(wave)
      double precision, intent(in) :: wave
      double precision :: Junk
      junk=wave !so compiler doesn't complain
      STflux = ST_flux_const
      end function STflux

      double precision function ABflux(wave)
      double precision, intent(in) :: wave
      ABflux = AB_flux_const * clight / (wave*wave)
      end function ABflux

      subroutine read_filters(list,set,num,ierr)
      character(len=256), intent(in) :: list
      character(len=256) :: filename, filter_filename
      type(spectrum), allocatable, intent(out) :: set(:)
      integer, intent(out) :: num, ierr
      integer :: f, w, pts, skip, num_cols, col_wave, col_flux
      double precision, pointer :: line(:)
      
      ierr=0
      filename=trim(data_dir)//'/'//trim(list)
      open(1,file=trim(filename),iostat=ierr)
      if(ierr/=0) return
      read(1,*) num
      allocate(set(num))
      read(1,*) num_cols, col_wave, col_flux
      allocate(line(num_cols))
      do f = 1, num
         read(1,'(a40,i5,1x,i5)') set(f)% filename, pts, skip
         set(f)% npts = pts - skip
         filter_filename = ''
         filter_filename = trim(data_dir)//'/'//'filters/'//set(f)% filename
         open(2,file=trim(filter_filename),iostat=ierr)
         if(ierr/=0) return
         allocate(set(f)% wave(set(f)% npts),set(f)% flux(set(f)% npts))
         if(skip > 0)then
            do w = 1,skip       !skip the header, if any
               read(2,*)
            enddo
         endif
         do w = 1,set(f)% npts
            read(2,*) line
            set(f)% wave(w) = line(col_wave)
            set(f)% flux(w) = max(line(col_flux), 0d0) !no negative thru-put allowed
         enddo
         close(2)
      enddo
      close(1)
      deallocate(line)
      end subroutine read_filters

      subroutine read_CK2003(file,set,num,ierr)
      !CK2003 flux in ergs/cm**2/s/hz/ster
      character(len=256), intent(in) :: file
      type(spectrum), pointer, intent(out) :: set(:)
      integer, intent(out) :: num, ierr
      integer, parameter :: nwav=1221, natm=476
      integer :: i
      double precision :: wl(nwav),fl(nwav),bb(nwav),feh

      ierr=0
      num= natm
      !read header and get [Fe/H]
      open(1,file=trim(file),iostat=ierr)
      if(ierr/=0) return
      read(1,'(14x,f4.1)')feh
      do i=1,21
         read(1,*)
      enddo

      !wavelength in nm -> convert to AA
      read(1,'(8f10.2)') wl
      wl=wl*1.0d1

      nullify(set)
      allocate(set(natm))

      !loop over all the atms in the file
      do i=1,natm
         allocate(set(i)% wave(nwav), set(i)% flux(nwav), set(i)% extinction(nwav))
         set(i)% npts = nwav
         set(i)% wave = wl
         set(i)% feh = feh
         set(i)% extinction = 1d0
      !read in Teff and Log g
         read(1,'(6x,f6.0,10x,f7.5,50x)',iostat=ierr) set(i)% Teff, set(i)% logg
         if(ierr/=0) exit
      !read in fluxes and convert to proper units
         read(1,'(8e10.4)')fl
         read(1,'(8e10.4)')bb
      !convert flux from /hz/ster to /AA -> 4*pi*c/lambda^2 
         set(i)% flux = pi4*clight*fl/(wl*wl) 
         set(i)% Fbol = sigma * set(i)% Teff**4
         if(debug) write(*,*) i, set(i)% npts, set(i)% Teff, set(i)% logg
      enddo
      close(1)
      end subroutine read_CK2003

      subroutine read_ATLAS_spec(list,set,num,ierr)
      character(len=256), intent(in) :: list
      type(spectrum), pointer, intent(out) :: set(:)
      integer, intent(out) :: num, ierr
      integer :: s
      ierr=0
      open(1,file=trim(list),iostat=ierr)
      if(ierr/=0) return
      read(1,*) num
      nullify(set)
      allocate(set(num))
      do s=1,num
         read(1,'(a)') set(s)% filename
         if(.not.read_on_the_fly)then
            call load_ATLAS_spec(set(s),ierr)
            if(set(s)% npts == 0) then
               write(0,*) '  WARNING: empty spectrum file: ', trim(set(s)% filename)
               ierr=-1
            endif   
         endif
      enddo
      close(1)
      end subroutine read_ATLAS_spec


      subroutine read_ATLAS_sed(list,set,num,ierr)
      character(len=256), intent(in) :: list
      type(spectrum), pointer, intent(out) :: set(:)
      integer, intent(out) :: num, ierr
      integer :: s
      ierr=0
      open(1,file=trim(list),iostat=ierr)
      if(ierr/=0) return
      read(1,*) num
      nullify(set)
      allocate(set(num))
      do s=1,num
         read(1,'(a)') set(s)% filename
         if(.not.read_on_the_fly)then
            call load_ATLAS_sed(set(s),ierr)
            if(set(s)% npts == 0) then
               write(0,*) '  WARNING: empty spectrum file: ', trim(set(s)% filename)
               ierr=-1
            endif   
         endif
      enddo
      close(1)
      end subroutine read_ATLAS_sed


      subroutine load_ATLAS_spec(s,ierr)
      type(spectrum), intent(inout) :: s
      integer, intent(out) :: ierr
      character(len=256) :: binfile, filename
      integer, parameter :: nwav=1700600 !nwav=30000
      integer :: i
      ierr=0
      filename = trim(data_dir) // '/' // trim(s% filename)
      binfile=trim(filename)//'.bin'
      open(2,file=trim(binfile),iostat=ierr,form='unformatted',status='old')
      if(ierr/=0) then  !'no binary file; open ascii file and write binfile
         close(2)
         open(2,file=trim(filename),iostat=ierr,status='old')
         if(ierr/=0) then
            write(*,*) trim( filename)
            return
         endif
         s% filename = s% filename
         call read_teff_logg_from_spec_file(s)
         s% feh = 0d0
         s% npts = nwav
         allocate(s% wave(s% npts), s% flux(s% npts), s% extinction(s% npts))
         do i=1,nwav
            read(2,*) s% wave(i), s% flux(i)
         enddo
         s% M = 1d0
         s% R = 1d0
      !convert flux from /hz/ster to /AA -> 4*pi*c/lambda^2 
         s% flux = pi4*clight*s% flux/(s% wave * s% wave) 
         s% Fbol = sigma * s% Teff**4
         call write_bin_file(s, binfile,ierr)
      else
         read(2) s% filename, s% npts, s% feh, s% Teff, &
         s% logg, s% R, s% M, s% Fbol
         allocate(s% wave(s% npts), s% flux(s% npts), s% extinction(s% npts))
         read(2) s% wave, s% flux
      endif
      close(2)
      if(do_check_total_flux) call check_total_flux(s)
      end subroutine load_ATLAS_spec

      subroutine read_teff_logg_from_spec_file(s)
      type(spectrum), intent(inout) :: s
      integer :: tloc, gloc, teff
      character(len=4) :: tchar, gchar
      tloc = index(s% filename,"_t",back=.true.)
      tchar= s% filename(tloc+2:tloc+5)
      gloc = index(s% filename,"g",back=.true.)
      gchar= s% filename(gloc+1:gloc+4)
      read(gchar,'(f4.2)') s% logg
      read(tchar,'(i4)') teff
      s% Teff = dble(teff)
      end subroutine read_teff_logg_from_spec_file

      subroutine load_ATLAS_sed(s,ierr)
      type(spectrum), intent(inout) :: s
      integer, intent(out) :: ierr
      character(len=256) :: binfile, filename
      integer, parameter :: nwav=26500
      integer :: i
      ierr=0
      filename = trim(data_dir) // '/' // trim(s% filename)
      binfile=trim(filename)//'.bin'
      open(2,file=trim(binfile),iostat=ierr,form='unformatted',status='old')
      if(ierr/=0) then  !'no binary file; open ascii file and write binfile
         close(2)
         open(2,file=trim(filename),iostat=ierr,status='old')
         if(ierr/=0) then
            write(*,*) trim( filename)
            return
         endif
         s% filename = s% filename
         call read_teff_logg_from_sed_file(s)
         s% feh = 0d0
         s% npts = nwav
         allocate(s% wave(s% npts), s% flux(s% npts), s% extinction(s% npts))
         do i=1,nwav
            read(2,*) s% wave(i), s% flux(i)
         enddo
         s% M = 1d0
         s% R = 1d0
        !convert flux from /hz/ster to /AA -> 4*pi*c/lambda^2 
         s% flux = pi4*clight*s% flux/(s% wave * s% wave) 
         s% Fbol = sigma * s% Teff**4
         call write_bin_file(s, binfile,ierr)
      else
         read(2) s% filename, s% npts, s% feh, s% Teff, &
         s% logg, s% R, s% M, s% Fbol
         allocate(s% wave(s% npts), s% flux(s% npts), s% extinction(s% npts))
         read(2) s% wave, s% flux
      endif
      close(2)
      if(do_check_total_flux) call check_total_flux(s)
      end subroutine load_ATLAS_sed
      
      subroutine read_teff_logg_from_sed_file(s)
      !ATLAS/ckc/at12_feh+0.0_afe+0.0_t5500g5.00.sed
      type(spectrum), intent(inout) :: s
      integer :: tloc, gloc, mloc, aloc, teff
      character(len=4) :: gchar, mchar, achar
      character(len=5) :: tchar
      mloc = index(s% filename,"feh",back=.true.)
      mchar= s% filename(mloc+3:mloc+6)
      aloc = index(s% filename,"afe",back=.true.)
      achar = s% filename(aloc+3:aloc+6)
      gloc = index(s% filename,"g",back=.true.)
      gchar= s% filename(gloc+1:gloc+4)
      tloc = index(s% filename,"_t",back=.true.)
      tchar = s% filename(tloc+2:gloc-1)
      read(gchar,'(f4.2)') s% logg
      read(tchar,'(i5)') teff
      read(mchar,'(f4.2)') s% FeH
      read(mchar,'(f4.2)') s% alpha_Fe
      s% Teff = dble(teff)
      !write(*,*) s% Teff, s% logg, s% FeH, s% alpha_Fe
      end subroutine read_teff_logg_from_sed_file

      subroutine read_phoenix(list,set,num,ierr)
      character(len=256), intent(in) :: list
      type(spectrum), pointer, intent(out) :: set(:)
      integer, intent(out) :: num, ierr
      integer :: s
      ierr=0
      open(1,file=trim(list),iostat=ierr)
      if(ierr/=0) return
      read(1,*) num
      nullify(set)
      allocate(set(num))
      do s=1,num
         read(1,'(a)') set(s)% filename
         if(.not.read_on_the_fly)then
            call load_phoenix_spec(set(s),ierr)
            if(set(s)% npts == 0) then
               write(0,*) '  WARNING: empty spectrum file: ', trim(set(s)% filename)
               ierr=-1
            endif   
         endif
      enddo
      close(1)
      end subroutine read_phoenix

      subroutine load_spec(choice,s,ierr)
      integer, intent(in) :: choice
      type(spectrum), intent(inout) :: s
      integer, intent(out) :: ierr
      ierr=0
      select case(choice)
      case(PHOENIX) 
         call load_phoenix_spec(s,ierr)
      case(ATLAS_spec)
         call load_ATLAS_spec(s,ierr)
      case(ATLAS_SED)
         call load_ATLAS_sed(s,ierr)
      end select
      end subroutine load_spec

      subroutine load_phoenix_spec(s,ierr)
      type(spectrum), intent(inout) :: s
      integer, intent(out) :: ierr
      character(len=256) :: binfile, filename
      ierr=0
      filename = trim(data_dir) // '/' // trim(s% filename)
      binfile=trim(filename)//'.bin'
      open(2,file=trim(binfile),iostat=ierr,form='unformatted',status='old')
      if(ierr/=0) then  !'no binary file; open ascii file and write binfile
         close(2)
         open(2,file=trim(filename),iostat=ierr,status='old')
         if(ierr/=0) then
            write(*,*) trim( filename)
            return
         endif
         read(2,*) s% Teff, s% logg, s% feh
         read(2,*) s% npts
         allocate(s% wave(s% npts), s% flux(s% npts), s% extinction(s% npts))
         read(2,*) s% wave
         read(2,*) s% flux
         s% M = 1d0
         s% R = 1d0
         s% flux = s% flux/phx_flux_conv
         s% Fbol = sigma * s% Teff**4
         call write_bin_file(s, binfile,ierr)
      else
         read(2) s% filename, s% npts, s% feh, s% Teff, &
         s% logg, s% R, s% M, s% Fbol
         allocate(s% wave(s% npts), s% flux(s% npts), s% extinction(s% npts))
         read(2) s% wave, s% flux
      endif
      close(2)
      if(do_check_total_flux) call check_total_flux(s)
      end subroutine load_phoenix_spec

      subroutine check_total_flux(s)
      !check integrated flux vs. expected from sigma*Teff^4
      type(spectrum), intent(in) :: s
      double precision :: flux1, flux2
      integer :: i
      flux1 = s% Fbol
      flux2 = 0d0
      do i=2,s% npts
         flux2=flux2 + s% flux(i) * (s% wave(i) - s%wave(i-1))
      enddo
      write(*,*) s% Teff, flux1, flux2, flux1/flux2
      end subroutine check_total_flux

      subroutine unload_spec(s,ierr)
      type(spectrum), intent(inout) :: s
      integer, intent(out) :: ierr
      ierr=0
      deallocate(s% wave, s% flux, s% extinction)
      end subroutine unload_spec

      subroutine extinction_for_spectrum(s,Av,Rv)
      type(spectrum), intent(inout) :: s
      double precision, intent(in) :: Av, Rv
      if(Av > 0d0 .and. Rv > 0d0)then
         s% extinction = 1d1**(-0.4d0*Av*Al_div_Av_CCM(s% wave,Rv))
      else
         s% extinction = 1d0
      endif
      end subroutine extinction_for_spectrum

      subroutine smooth_spectrum(s,n)
      type(spectrum), intent(inout) :: s
      integer, intent(in) :: n !number of smoothing passes
      integer :: i, j
      do j=1,n
         do i=1,s% npts
            if(i==1 .or. i==s% npts) then
               cycle
            else if(i==2 .or. i==s% npts-1) then
               s% flux(i) = (s% flux(i-1) + s% flux(i) + s% flux(i+1))/3d0
            else
               s% flux(i) = (s% flux(i-2) + s% flux(i-1) + s% flux(i) + s% flux(i+1) + s% flux(i+2))/5d0
            endif
         enddo
      enddo
      end subroutine smooth_spectrum

      subroutine write_spectrum(s,out,smooth,ierr)
      type(spectrum), intent(in) :: s
      logical, intent(in) :: smooth
      character(len=256), intent(in) :: out
      integer, intent(out) :: ierr
      integer :: i, pass
      type(spectrum), pointer :: s1
      logical, parameter :: long_format = .false.

      ierr=0

      write(*,*) s% filename

      if(long_Format)then
         pass = 10
         nullify(s1)
         allocate(s1)
         s1 = s
      
         if(smooth) call smooth_spectrum(s1,pass)

         open(1,file=trim(out),iostat=ierr)
         if(ierr/=0) return
         write(*,'(a1,3f12.3,i10)') '#', s1% Teff, s1% logg, s1% FeH, s1% npts
         do i=1,s1% npts
            write(1,'(1p3e20.10)') s1% wave(i), s1% flux(i), s1% extinction(i)
         enddo
         close(1)
         deallocate(s1)
      endif
      end subroutine write_spectrum

      function integrate_bandpass(s,f,ierr) result(integral)
      type(spectrum), intent(in) :: s, f
      integer, intent(out) :: ierr
      integer :: i
      double precision :: integral
      ierr = 0
      integral = 0d0
      !check that filter and spectrum are compatible
      if(s% npts < 2) return
      if( f% wave(1) > s% wave(s% npts) .or. f% wave(f% npts) < s% wave(1)) then
         write(0,*) "*WARNING: No overlap between spectrum and filter*"
         return
      endif
      do i=2,s% npts !integrate over s: integral = flux* filter *wave * dwave
         integral=integral+filter_interp(f,s% wave(i))*s% extinction(i)* &
                            s% flux(i)*s% wave(i)*(s% wave(i)-s% wave(i-1))
      enddo
      end function integrate_bandpass

      subroutine R_from_Teff_and_logg(s)
      type(spectrum), intent(inout) :: s
      s% R = sqrt(G * s% M * Msun/10d0**s% logg)
      end subroutine R_from_Teff_and_logg

      double precision function filter_interp(f,w)
      type(spectrum), intent(in) :: f
      double precision, intent(in) :: w
      integer :: loc, i
      double precision :: q(4)
      if(w < f% wave(1) .or. w > f% wave(f% npts))then
         filter_interp = 0d0
         return
      endif
      loc = 1
      loc = binary_search(f% npts, f% wave, loc, w)
      loc = max(loc,2)
      loc = min(loc, f% npts - 2)
      call interp(f% wave(loc-1:loc+2),q,w,4)
      forall(i=1:4) q(i) = q(i) * f% flux(loc-2+i)
      filter_interp = sum(q)
      end function filter_interp

      include 'num_binary_search.dek'

      subroutine interp(a,b,x,n)
! {a} are the tabulated values for use in interpolation
! {b} are coefficients of the interpolating polynomial
!  x  is the abscissa to be interpolated
!  n  is the number of points to be used, interpolating polynomial
!     has order n-1 
      integer, intent(in) :: n
      double precision, intent(in) :: a(n), x
      double precision, intent(out) :: b(n)
      integer :: i,j
      do i=1,n
         b(i)=1.0d0
         do j=1,n
            if(j/=i) b(i)=b(i)*(x-a(j))/(a(i)-a(j))
         enddo
      enddo
      end subroutine interp

      subroutine interp_4pt_pm(x, y, a) 
         ! returns coefficients for monotonic cubic interpolation from x(2) to x
         double precision, intent(in)    :: x(4)    ! junction points, strictly monotoni
         double precision, intent(in)    :: y(4)    ! data values at x's
         double precision, intent(out)   :: a(3)    ! coefficients
         double precision :: h1, h2, h3, s1, s2, s3, p2, p3, as2, ss2, yp2, yp3
         ! for x(2) <= x <= x(3) and dx = x-x(2), 
         ! y(x) = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
         h1 = x(2)-x(1)
         h2 = x(3)-x(2)
         h3 = x(4)-x(3)
         s1 = (y(2)-y(1))/h1
         s2 = (y(3)-y(2))/h2
         s3 = (y(4)-y(3))/h3
         p2 = (s1*h2+s2*h1)/(h1+h2)
         p3 = (s2*h3+s3*h2)/(h2+h3)
         as2 = abs(s2)
         ss2 = sign(1d0, s2)
         yp2 = (sign(1d0, s1)+ss2)*min(abs(s1), as2, 0.5d0*abs(p2))
         yp3 = (ss2+sign(1d0, s3))*min(as2, abs(s3), 0.5d0*abs(p3))
         a(1) = yp2
         a(2) = (3*s2-2*yp2-yp3)/h2
         a(3) = (yp2+yp3-2*s2)/h2**2
      end subroutine interp_4pt_pm

      end module bandpass
