   module bandpass
     
     use const
     use blackbody
     use extinction
      
     implicit none

     character(len=256) :: data_dir
     logical :: data_dir_set = .false., read_on_the_fly = .false.
     logical, parameter :: debug = .false., do_check_total_flux=.false.

     integer :: zero_point_type
     
     type spectrum
        character(len=256) :: filename
        integer :: npts
        real(dp) :: FeH, Teff, logg, R, M, Fbol, alpha_Fe ! R in cm, M in Msun
        real(dp), allocatable :: wave(:), flux(:), extinction(:)
     end type spectrum

   contains

     subroutine set_data_dir(my_data_dir)
       character(len=256), intent(in) :: my_data_dir
       data_dir = my_data_dir
       data_dir_set = .true.
     end subroutine set_data_dir

     subroutine create_BBs(set,num,ierr)
       type(spectrum), pointer, intent(out) :: set(:)
       integer, intent(out) :: num, ierr
       integer, parameter :: nwav=20000
       integer :: i
       integer, parameter :: numT = 106
       real(dp) :: wave(nwav), dw
       real(dp) :: Teff(numT)

       ierr = 0
       num = numT
       
       nullify(set)
       allocate(set(num))

       !set the temperature scale for MIST bolometric correction tables
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
  
       if(debug) then
          do i=1,numT
             write(*,*) i, Teff(i)
          enddo
       endif

       !set the wavelength scale
       wave(1) = 1d-1
       wave(nwav) = 1d10
       dw = log(wave(nwav)/wave(1))/real(nwav-1,kind=dp)
       do i=2,nwav-1
          wave(i)=exp(log(wave(i-1))+dw)
       enddo

       set(:)% logg = 5d0

       do i=1,num
          set(i)% Teff = Teff(i)
          allocate(set(i)% wave(nwav), set(i)% flux(nwav), set(i)% extinction(nwav))
          set(i)% npts = nwav
          set(i)% wave = wave
          set(i)% feh = 0d0
          set(i)% extinction = 1d0
          set(i)% Fbol = sigma * set(i)% Teff * set(i)% Teff * set(i)% Teff * set(i)% Teff
          set(i)% flux = BBflux(wave,set(i)% Teff)
       enddo

     end subroutine create_BBs
     
     subroutine write_bin_file(s, binfile,ierr)
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

     subroutine read_bin_file(io,s)
       integer, intent(in) :: io
       type(spectrum), intent(out) :: s
       read(io) s% filename, s% npts, s% feh, s% Teff, s% logg, s% R, s% M, s% Fbol
       allocate(s% wave(s% npts), s% flux(s% npts), s% extinction(s% npts))
       read(io) s% wave, s% flux
     end subroutine read_bin_file
     
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

     !not used
     real(dp) function STflux(wave)
       real(dp), intent(in) :: wave
       real(dp) :: Junk
       junk=wave !so compiler doesn't complain
       STflux = ST_flux_const
     end function STflux

     !not used
     real(dp) function ABflux(wave)
       real(dp), intent(in) :: wave
       ABflux = AB_flux_const * clight / (wave*wave)
     end function ABflux
     
     subroutine read_filters(list,set,num,ierr)
       character(len=256), intent(in) :: list
       character(len=256) :: filename, filter_filename
       type(spectrum), allocatable, intent(out) :: set(:)
       integer, intent(out) :: num, ierr
       integer :: f, w, pts, skip, num_cols, col_wave, col_flux
       real(dp), pointer :: line(:)
       
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
       real(dp) :: wl(nwav),fl(nwav),bb(nwav),feh
       
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

     subroutine read_spec(choice,list,set,num,ierr)
       integer, intent(in) :: choice
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
             call load_spec(choice,set(s),ierr)
             if(set(s)% npts == 0) then
                write(0,*) '  WARNING: empty spectrum file: ', trim(set(s)% filename)
                ierr=-1
             endif
          endif
       enddo
       close(1)
     end subroutine read_spec     

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
       case(RAUCH)
          call load_rauch_spec(s,ierr)
       case(KOESTER)
          call load_koester_spec(s,ierr)
       end select
     end subroutine load_spec

     subroutine load_ATLAS_spec(s,ierr)
       type(spectrum), intent(inout) :: s
       integer, intent(out) :: ierr
       character(len=256) :: binfile, filename
       integer, parameter :: nwav=1700600 !nwav=30000
       integer :: i
       ierr=0
       filename = trim(s% filename)
       binfile=trim(filename)//'.bin'
       open(2,file=trim(binfile),iostat=ierr,form='unformatted',status='old')
       if(ierr/=0) then  !no binary file; open ascii file and write binfile
          close(2)
          open(2,file=trim(filename),iostat=ierr,status='old')
          if(ierr/=0) then
             write(*,*) trim( filename)
             return
          endif
          s% filename = s% filename
          call teff_logg_from_spec_file(s)
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
          call read_bin_file(2,s)
       endif
       close(2)
       if(do_check_total_flux) call check_total_flux(s)
     end subroutine load_ATLAS_spec

     subroutine teff_logg_from_spec_file(s)
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
     end subroutine teff_logg_from_spec_file

     subroutine load_ATLAS_sed(s,ierr)
       type(spectrum), intent(inout) :: s
       integer, intent(out) :: ierr
       character(len=256) :: binfile, filename
       !integer, parameter :: nwav=26500
       integer, parameter :: nwav=47378
       integer :: i
       ierr=0
       filename = trim(s% filename)
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
          call teff_logg_from_sed_file(s)
          s% feh = 0d0
          s% npts = nwav
          allocate(s% wave(s% npts), s% flux(s% npts), s% extinction(s% npts))
          do i=1,nwav
             read(2,'(f14.3,e12.4)',iostat=ierr) s% wave(i), s% flux(i)
             !check for NaNs, replace with zero
             if(s% flux(i) - s% flux(i) /= 0) s% flux(i) = 0d0
          enddo
          !check for less than zero:
          s% flux = abs(s% flux)

          s% M = 1d0
          s% R = 1d0
          !convert flux from /hz/ster to /AA -> 4*pi*c/lambda^2 
          s% flux = pi4*clight*s% flux/(s% wave * s% wave) 
          s% Fbol = sigma * s% Teff**4
          call write_bin_file(s, binfile,ierr)
       else
          call read_bin_file(2,s)
       endif
       close(2)
       if(do_check_total_flux) call check_total_flux(s)
     end subroutine load_ATLAS_sed
      
     subroutine teff_logg_from_sed_file(s)
       type(spectrum), intent(inout) :: s
       integer :: tloc, gloc, mloc, aloc, teff
       character(len=4) :: gchar, mchar, achar
       character(len=5) :: tchar
       mloc = index(s% filename,"feh",back=.true.)
       if(mloc>0) then
          mchar= s% filename(mloc+3:mloc+6)
          read(mchar,'(f4.2)') s% FeH
       else
          s% FeH = 0
       endif
       
       aloc = index(s% filename,"afe",back=.true.)
       if(aloc>0)then
          achar = s% filename(aloc+3:aloc+6)
          read(achar,'(f4.2)') s% alpha_Fe
       else
          s% alpha_Fe = 0
       endif

       gloc = index(s% filename,"g",back=.true.)
       gchar= s% filename(gloc+1:gloc+4)
       read(gchar,'(f4.2)') s% logg

       tloc = index(s% filename,"_t",back=.true.)
       tchar = s% filename(tloc+2:gloc-1)
       read(tchar,'(i5)') teff
       s% Teff = dble(teff)

       write(*,*) s% Teff, s% logg, s% FeH, s% alpha_Fe
     end subroutine teff_logg_from_sed_file

     subroutine load_rauch_spec(s,ierr)
       type(spectrum), intent(inout) :: s
       integer, intent(out) :: ierr
       character(len=256) :: binfile, filename
       character(len=60) :: line
       integer :: i
       ierr=0
       filename=trim(s% filename)
       binfile=trim(filename) // '.bin'
       open(2,file=trim(binfile),iostat=ierr,form='unformatted',status='old')
       if(ierr/=0)then
          close(2)
          open(2,file=trim(filename),iostat=ierr,status='old')
          if(ierr/=0)then
             write(*,*) trim(filename)
             return
          endif
          
          !read header
          do i=1,35
             read(2,'(a)') line
             if(i==6) then       !get Teff
                read(line(10:18),*) s% Teff
             elseif(i==7)then
                read(line(10:18),*) s% logg
             elseif(i==35)then
                read(line(3:7),'(i5)') s% npts
             endif
          enddo
          
          s% FeH = 0d0
          
          allocate(s% wave(s% npts), s% flux(s% npts), s% extinction(s% npts))
          
          do i=1,s% npts
             read(2,'(e18.12,e12.4)') s% wave(i), s% flux(i)
          enddo
          close(2)
          
          s% flux = WD_flux_conv * s% flux
          s% Fbol = sigma * s% Teff**4
          s% R    = 1d0
          s% M    = 1d0

          call write_bin_file(s, binfile, ierr)
          
       else
          call read_bin_file(2,s)
       endif
       close(2)
       if(do_check_total_flux) call check_total_flux(s)
     end subroutine load_rauch_spec

     subroutine load_koester_spec(s,ierr)
       type(spectrum), intent(inout) :: s
       integer, intent(out) :: ierr
       character(len=256) :: binfile, filename, line
       integer :: i, lo, hi
       integer, parameter :: num_extra = 100
       real(dp) :: wave_extra(num_extra), flux_extra(num_extra), dwave
       ierr=0
       filename=trim(s% filename)
       binfile=trim(filename) // '.bin'
       open(2,file=trim(binfile),iostat=ierr,form='unformatted',status='old')
       if(ierr/=0) then  !'no binary file; open ascii file and write binfile
          close(2)
          open(2,file=trim(filename),iostat=ierr,status='old')
          if(ierr/=0) then
             write(*,*) trim( filename)
             return
          endif
          !read the header and get some useful info
          do i=1,33
             read(2,'(a)') line
             !write(*,*) i, line(1:7)
             if(trim(line(1:7))==  'TEFF') read(line(12:30),*) s% Teff
             if(trim(line(1:7))== 'LOG_G') read(line(12:30),*) s% logg
             if(trim(line(1:7))=='NAXIS2') read(line(12:30),*) s% npts
          enddo
          s% npts = s% npts + num_extra !fudge to add a blackbody tail beyond 3 microns
          allocate(s% wave(s% npts), s% flux(s% npts), s% extinction(s% npts))
          !now read wavelength and flux
          s% flux = 0d0
          do i=1,s% npts-num_extra
             read(2,'(f11.3,e13.5)') s% wave(i), s% flux(i)
          enddo
          s% flux = WD_flux_conv * s% flux !then convert flux to match sigma*T^4
          s% FeH = 0d0
          s% M = 1d0
          s% R = 1d0
          s% Fbol = sigma * s% Teff * s% Teff * s% Teff * s% Teff

          !graft on the extra blackbody flux
          if(num_extra > 0)then
             wave_extra(1) = s% wave(s% npts - num_extra) + 5d1
             wave_extra(num_extra) = 1d7
             dwave = log(wave_extra(num_extra)/wave_extra(1))/real(num_extra-1,kind=dp)
             do i=2,num_extra-1
                wave_extra(i) = exp(log(wave_extra(i-1))+dwave)
             enddo
             lo=s% npts - num_extra + 1
             hi=s% npts 
             flux_extra = BBflux(wave_extra, s% Teff)
             s% flux(lo:hi) = flux_extra
             s% wave(lo:hi) = wave_extra
          endif

          !resume normal operations....
          call write_bin_file(s, binfile,ierr)
       else
          call read_bin_file(2,s)
       endif
       close(2)
       if(do_check_total_flux) call check_total_flux(s)
     end subroutine load_koester_spec

     subroutine load_phoenix_spec(s,ierr)
       type(spectrum), intent(inout) :: s
       integer, intent(out) :: ierr
       character(len=256) :: binfile, filename
       ierr=0
       filename = trim(s% filename)
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
          s% flux = phx_flux_conv * s% flux
          s% Fbol = sigma * s% Teff**4
          call write_bin_file(s, binfile,ierr)
       else
          call read_bin_file(2,s)
       endif
       close(2)
       if(do_check_total_flux) call check_total_flux(s)
     end subroutine load_phoenix_spec

     subroutine check_total_flux(s)
       !check integrated flux vs. expected from sigma*Teff^4
       type(spectrum), intent(in) :: s
       real(dp) :: flux1, flux2
       flux1 = s% Fbol
       flux2 = tsum(s% wave, s% flux)
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
       real(dp), intent(in) :: Av, Rv
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
       real(dp) :: integral
       real(dp), allocatable :: filter(:)
       ierr = 0
       integral = 0d0
       !check that filter and spectrum are compatible
       if(s% npts < 2) return
       if( f% wave(1) > s% wave(s% npts) .or. f% wave(f% npts) < s% wave(1)) then
          write(0,*) "*WARNING: No overlap between spectrum and filter*"
          write(*,*) f% wave(1), f% wave(f% npts)
          write(*,*) s% wave(1), s% wave(s% npts)
          write(*,*) s% filename
          stop
          return
       endif
       allocate(filter(s% npts))
       do i=1,s% npts
          filter(i) = filter_interp(f,s% wave(i))
       enddo
       integral = tsum(s% wave, filter*s% flux*s% extinction*s% wave)
       deallocate(filter)
     end function integrate_bandpass

     subroutine R_from_Teff_and_logg(s)
       type(spectrum), intent(inout) :: s
       s% R = sqrt(G * s% M * Msun/10d0**s% logg)
     end subroutine R_from_Teff_and_logg
     
     real(dp) function filter_interp(f,w)
       type(spectrum), intent(in) :: f
       real(dp), intent(in) :: w
       integer :: loc, i
       real(dp) :: q(4)
       if(w < f% wave(1) .or. w > f% wave(f% npts))then
          filter_interp = 0d0
          return
       endif
       loc = 1
       loc = binary_search(f% npts, f% wave, loc, w)
       loc = max(loc,2)
       loc = min(loc, f% npts - 2)
       call interp(f% wave(loc-1:loc+2),q,w,4)
       do i=1,4
          q(i) = q(i)*f% flux(loc-2+i)
       enddo
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
       real(dp), intent(in) :: a(n), x
       real(dp), intent(out) :: b(n)
       integer :: i,j
       do i=1,n
          b(i)=1.0d0
          do j=1,n
             if(j/=i) b(i)=b(i)*(x-a(j))/(a(i)-a(j))
          enddo
       enddo
     end subroutine interp
     
     subroutine interp_4pt_pm(x, y, a) !from MESA/interp_1d by Bill Paxton
       ! returns coefficients for monotonic cubic interpolation from x(2) to x
       real(dp), intent(in)    :: x(4)    ! junction points, strictly monotoni
       real(dp), intent(in)    :: y(4)    ! data values at x's
       real(dp), intent(out)   :: a(3)    ! coefficients
       real(dp) :: h1, h2, h3, s1, s2, s3, p2, p3, as2, ss2, yp2, yp3
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

     !from FSPS by Charlie Conroy, simple trapezoidal integration
     !                                               \int y(x) dx
     function tsum(x,y)
       real(dp), intent(in) :: x(:), y(:)
       real(dp) :: tsum
       integer :: n
       n=size(x)
       tsum = 0.5d0*sum( abs(x(2:n) - x(1:n-1)) * (y(2:n) + y(1:n-1)) )
     end function tsum

   end module bandpass
