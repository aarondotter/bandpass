      program wash_fix

      implicit none
      integer :: i, ierr
      double precision :: wave, thru
      character(len=128) :: infile

      call getarg(1,infile)

      open(1,file=trim(infile))
      i=1
      ierr=0
      do while(.true.)
         read(1,*,iostat=ierr) wave, thru
         if(ierr/=0) exit
         write(*,'(1p2e20.10)') wave*1d1, thru/wave
      enddo
      close(1)
      
      end program wash_fix
