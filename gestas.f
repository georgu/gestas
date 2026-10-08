
!------------------------------------------------------------
!
! box model gestas
!
! needs files boxes*.txt from shyfem simulation
!
! todo:
!
! consider rain and evaporation for volume error
!
!------------------------------------------------------------

!============================================================
	module boxes
!============================================================

	implicit none

	integer, save :: nvers_file = 0	!version of files

	integer, save :: nbox = 0	!total number of boxes
	integer, save :: nsect = 0	!total number of sections
	integer, save :: nosect = 0	!total number of OB sections
	integer, save :: nbound = -1	!total number of OB boundaries
	double precision, save :: dt = 0.

	integer, save :: imax = 0	!maximum number of time steps
	integer, save :: irec = 0

	!integer, parameter :: ifreq = 1    !frequency of output (every time)
	!integer, parameter :: ifreq = 24   !frequency of output (daily)
	!integer, parameter :: ifreq = 168  !frequency of output (weekly)
	integer, parameter :: ifreq = 730  !frequency of output (monthly)
	!integer, parameter :: ifreq = 8766 !frequency of output (yearly)

	integer, save :: date0,time0
	double precision, save :: atime0

	double precision, save :: dtime
	character*20, save :: stime
	character*20, save :: stime_start, stime_end

	integer, save, allocatable :: nlayer(:)
	double precision, save, allocatable :: area(:)
	double precision, save, allocatable :: volume0(:)
	double precision, save, allocatable :: depth(:)
	integer, save, allocatable :: isect(:,:)	!section index
	integer, save, allocatable :: iosect(:,:)	!OB section index

	double precision, save, allocatable :: zold(:)
	double precision, save, allocatable :: znew(:)
	double precision, save, allocatable :: vold(:)
	double precision, save, allocatable :: vnew(:)

	double precision, save, allocatable :: temp(:)
	double precision, save, allocatable :: salt(:)
	double precision, save, allocatable :: cspeed(:)

	double precision, save, allocatable :: meteo(:,:)

	double precision, save, allocatable :: flux(:,:)  !section fluxes
	double precision, save, allocatable :: oflux(:,:) !OB section fluxes

	integer :: iheader = 30*max(1,ifreq)
	character*80, save :: header = ' iterations                dtime'
     &			// '             date and time             verr'

!============================================================
	end module boxes
!============================================================

	program gestas

	use boxes

	implicit none

	logical bout
	integer ierr
	double precision verr
	double precision time_start,time_end

	bout = ifreq == 0

	call cpu_time(time_start)

	call read_stats
	call initialize_time
	call get_open_boundary_conditions
	call alloc_arrays
	call read_zinit
	call scalar_init

	write(6,*) '==========================='
	write(6,*) 'start of time loop'
	write(6,*) '==========================='
	write(6,*) trim(header)

	do
	  call copy_new2old
	  call read_time_record(ierr)
	  call read_meteo_record(ierr)
	  if( ierr /= 0 ) exit
	  irec = irec + 1
	  call check_volume(verr)
	  if( ifreq /= 0 ) bout = ( mod(irec,ifreq) == 0 )
	  if( mod(irec,iheader) == 0 ) write(6,*) trim(header)
	  if( bout ) write(6,*) irec,dtime,stime,sngl(verr)
	  call hydro_out
	  call meteo_out
	  call scalar_run
	  if( imax > 0 .and. irec > imax ) exit
	end do

	write(6,*) '==========================='
	write(6,*) 'end of time loop'
	write(6,*) '==========================='

	write(6,*) irec,' time records read'
	write(6,*) 'start time of simulation: ',stime_start
	write(6,*) 'end time of simulation:   ',stime_end

	call scalar_finalize

	call cpu_time(time_end)
	write(6,*) time_end-time_start,' seconds cpu time used'

	end

!************************************************************

	subroutine check_volume(verr)

	use boxes

	implicit none

	double precision verr

	integer is,ib
	double precision vdiff,vaver,vaaver,vadiff
	double precision volume(nbox)

	volume = vold

	do is=1,nsect
	  !write(6,*) isect(:,is)
	end do

	do is=1,nsect
	  ib = isect(1,is)	! box from
	  if( ib > 0 ) then
	    volume(ib) = volume(ib) + dt*flux(3,is)
	    volume(ib) = volume(ib) - dt*flux(2,is)
	  end if
	  ib = isect(2,is)	! box to
	  if( ib > 0 ) then
	    volume(ib) = volume(ib) + dt*flux(2,is)
	    volume(ib) = volume(ib) - dt*flux(3,is)
	  end if
	end do

	if( nvers_file <= 5 ) then
	  do is=1,nosect
	    ib = iosect(1,is)	! box from
	    if( ib > 0 ) then
	      volume(ib) = volume(ib) + dt*oflux(3,is)
	      volume(ib) = volume(ib) - dt*oflux(2,is)
	    end if
	    ib = iosect(2,is)	! box to
	    if( ib > 0 ) then
	      volume(ib) = volume(ib) + dt*oflux(2,is)
	      volume(ib) = volume(ib) - dt*oflux(3,is)
	    end if
	  end do
	end if

	vaaver = 0.
	vadiff = 0.
	do ib=1,nbox
	  vaver = 0.5*(vnew(ib)+volume(ib))
	  vdiff = abs(vnew(ib)-volume(ib))
	  vaaver = vaaver + vaver
	  vadiff = vadiff + vdiff
	  !write(6,*) ib,vnew(ib),volume(ib),vdiff,vdiff/vaver
	end do

	verr = vadiff / vaaver

	end

!************************************************************

	subroutine copy_new2old

	use boxes

	implicit none

	zold = znew
	vold = vnew

	end

!************************************************************

	subroutine alloc_arrays

	use boxes

	implicit none

	allocate(zold(nbox),znew(nbox))
	allocate(vold(nbox),vnew(nbox))
	allocate(temp(nbox),salt(nbox))
	allocate(cspeed(nbox))

	allocate(flux(3,nsect))

	allocate(meteo(8,nbox))

	write(6,*) 'finished allocating arrays'

	end

!************************************************************

	subroutine get_open_boundary_conditions

! we have to read first part of 2d file to find out about nbound

	use boxes

	implicit none

	integer ierr

	integer, save :: iu = 0
	character*80 :: file

	integer nb,ns,nvars,nob,nvers
	integer i,is,from

	iu = 11
	file='boxes_2d.txt'
	call check_file_exists(file)
	open(iu,file=file,status='old',form='formatted')
	call check_file(iu,2,nvers,ierr)
	if( ierr < 0 ) then
	  stop 'error stop get_open_boundary_conditions: ' //
     &			'reading boxes_2d.txt'
	end if
	if( nvers_file == 0 ) nvers_file = nvers
	if( nvers /= nvers_file ) then
	  stop 'error stop get_open_boundary_conditions: ' //
     &			'inconsistent nvers'
	end if

	call skip_lines(iu,3)		!skip time record
	read(iu,*) nb,nvars
	call skip_lines(iu,nb+2)	!skip boxes
	read(iu,*) ns
	if( nvers <= 5 ) then
	  call skip_lines(iu,ns+2)	!skip sections
	  read(iu,*) nob
	else
	  call skip_lines(iu,1)		!skip sections
	  nob = 0
	  do i=1,ns
	    read(iu,*) is,from
	    if( i /= is ) then
	      write(6,*) 'i,is: ',i,is
	      stop 'error stop get_open_boundary_conditions: i/=is'
	    end if
	    if( from < 0 ) nob = nob + 1	!boundary
	  end do
	end if
	close(iu)

	nbound = nob
	nosect = nob
	write(6,*) 'open boundary conditions: ',nbound

	end

!************************************************************

	subroutine skip_lines(iu,n)

	implicit none

	integer iu,n

	integer i

	do i=1,n
	  read(iu,*)
	end do

	end

!************************************************************

	subroutine read_time_record(ierr)

	use boxes

	implicit none

	integer ierr

	integer, save :: iu = 0
	character*80 :: file
	character*80 :: string

	logical, save :: bfirst = .true.

	integer ib,is
	integer i,ii,i1,i2
	integer nb,ns,no,nob,nvars
	integer, save :: nvers = 0
	double precision a

	nob = 0		!number of open boundary conditions

	if( bfirst ) then
	  iu = 2
	  file='boxes_2d.txt'
	  call check_file_exists(file)
	  open(iu,file=file,status='old',form='formatted')
	  call check_file(iu,2,nvers,ierr)
	  if( ierr < 0 ) return
	  if( nvers_file == 0 ) nvers_file = nvers
	  if( nvers /= nvers_file ) then
	    stop 'error stop read_time_record: inconsistent nvers'
	  end if
	end if

	read(iu,'(a)',end=1)
	read(iu,*) dtime,stime

	call check_time
	if( bfirst ) stime_start = stime
	stime_end = stime

	read(iu,*)
	read(iu,*) nb,nvars
	if( nb /= nbox ) stop 'error stop read_time_record: nb/=nbox'
	if( nvars /= 7 ) stop 'error stop read_time_record: nvars/=7'

	read(iu,*)
	do ib=1,nbox
	  read(iu,*) i,temp(ib),salt(ib),cspeed(ib),znew(ib),a,a,vnew(ib)
	  if( ib /= i ) stop 'error stop read_time_record: ib/=i'
	end do

	read(iu,*)
	read(iu,*) ns
	if( ns > nsect ) stop 'error stop read_time_record: ns/=nsect'

	flux = 0.
	read(iu,*)
	do is=1,ns
	  if( nvers <= 5 ) then
	    i = is
	    read(iu,*) i1,i2,(flux(ii,is),ii=1,3)
	  else
	    read(iu,*) i,i1,i2,(flux(ii,is),ii=1,3)
	  end if
	  if( is /= i ) stop 'error stop read_time_record: is/=i'
	  nob = min(nob,i1,i2)
	end do

	if( nvers <= 5 ) then

	read(iu,*)
	read(iu,*) no

	if( bfirst ) then
	  nosect = no
	  allocate(oflux(3,nosect),iosect(2,nosect))
	  !write(6,*) 'open boundary sections: ',nosect
	end if

	if( no /= nosect ) then
	  write(6,*) nosect,no
	  write(6,*) 'number of OB sections has changed'
	  stop 'error stop read_time_record: nosect /= no'
	end if

	read(iu,*)
	do is=1,no
	  read(iu,*) (iosect(ii,is),ii=1,2),(oflux(ii,is),ii=1,3)
	  nob = min(nob,iosect(1,is),iosect(2,is))
	end do

	end if

	nob = -nob
	if( nbound == 0 ) then
	  nbound = nob
	  !write(6,*) 'open boundary conditions: ',nbound
	end if

	if( nbound /= nob ) then
	  write(6,*) nbound,nob
	  write(6,*) 'number of OBC have changed'
	  stop 'error stop read_time_record: nbound /= nob'
	end if

	bfirst = .false.

	ierr = 0
	return
    1	continue
	ierr = -1
	return
	end

!************************************************************

	subroutine read_meteo_record(ierr)

	use boxes

	implicit none

	integer ierr

	integer, save :: iu = 0
	character*80 :: file

	integer ib,is
	integer i,ii,i1,i2
	integer nb,ns,no,nob,nvars,nvers
	double precision a

	nob = 0		!number of open boundary conditions

	if( iu == 0 ) then
	  iu = 1
	  file='boxes_meteo.txt'
	  call check_file_exists(file)
	  open(iu,file=file,status='old',form='formatted')
	  call check_file(iu,1,nvers,ierr)
	  if( ierr < 0 ) return
	  if( nvers_file == 0 ) nvers_file = nvers
	  if( nvers /= nvers_file ) then
	    stop 'error stop read_meteo_record: inconsistent nvers'
	  end if
	end if

	read(iu,*,end=1)
	read(iu,*) dtime,stime

	call check_time

	read(iu,*)
	read(iu,*) nb,nvars
	if( nb /= nbox ) stop 'error stop read_time_record: nb/=nbox'
	if( nvars /= 8 ) stop 'error stop read_time_record: nvars/=8'

	read(iu,*)
	do ib=1,nbox
	  read(iu,*) i,meteo(1:8,ib)
	  if( ib /= i ) stop 'error stop read_time_record: ib/=i'
	end do

	ierr = 0
	return
    1	continue
	ierr = -1
	return
	end

!************************************************************

	subroutine check_file(iu,ft,nvers,ierr)

	implicit none

	integer iu
	integer ft
	integer nvers
	integer ierr

	integer idbox,ftype
	integer, parameter :: idbox_expect = 473226
	integer, parameter :: nvers_expect = 7
	character*80 string,file
	character*2 aux

	ierr = 0

	read(iu,*,end=97)
	read(iu,'(a)') string
	string(1:) = string(3:)
	!write(6,*) '|',trim(string),'|'
	!write(6,*) '|','123456789012345678901234567890','|'
	!write(6,*) 'gguu: ',iu,'  ',trim(string)
	read(string,'(3i10)') idbox,ftype,nvers

	if( idbox /= idbox_expect ) goto 99
	if( ft /= ftype ) goto 99
	if( nvers == 3 ) then
	  !write(6,*) 'old version... using compatibility mode: ',nvers
	else if( nvers /= nvers_expect ) then
	  goto 99
	end if

	return
   97	continue
	ierr = -1
	return
   99	continue
	write(6,*) 'error reading header of file'
	call get_file_name(ft,file)
	write(6,*) 'file name: ',trim(file)
	write(6,*) '                    id       ftype       nvers'
	write(6,*) 'read:     ',idbox,ftype,nvers
	write(6,*) 'expected: ',idbox_expect,ft,nvers_expect
	stop 'error stop check_file: error in header'
	end

!************************************************************

	subroutine check_file_exists(file)

	implicit none

	character*(*) file

	logical bex

	inquire(file=file,exist=bex)

	if( .not. bex ) then
	  write(6,*) 'file not existing: ',trim(file)
	  stop 'error stop check_file_exists: no such file'
	end if

	end

!************************************************************

	subroutine get_file_name(ftype,file)

	implicit none

	integer ftype
	character*(*) file

	if( ftype == 0 ) then
	  file = 'boxes_init.txt'
	else if( ftype == 1 ) then
	  file = 'boxes_meteo.txt'
	else if( ftype == 2 ) then
	  file = 'boxes_2d.txt'
	else if( ftype == 3 ) then
	  file = 'boxes_3d.txt'
	else if( ftype == 7 ) then
	  file = 'boxes_geom.txt'
	else if( ftype == 8 ) then
	  file = 'boxes_stats.txt'
	else
	  write(6,*) 'ftype = ',ftype
	  stop 'error stop get_file_name: unknown file type'
	end if

	end

!************************************************************

	subroutine read_stats

	use boxes

	implicit none

	integer :: iu = 1
	character*80 :: file

	logical bvers3
	integer ib,is,nvers,ierr
	integer i,ii,is_aux

	file='boxes_stats.txt'
	call check_file_exists(file)
	open(iu,file=file,status='old',form='formatted',err=99)
	call check_file(iu,8,nvers,ierr)
	if( ierr < 0 ) stop 'error stop read_stats: ' //
     &			'reading boxes_stats.txt'
	if( nvers_file == 0 ) nvers_file = nvers
	if( nvers /= nvers_file ) then
	  stop 'error stop read_stats: inconsistent nvers'
	end if

	bvers3 = ( nvers == 3 )

	read(iu,*)
	read(iu,*) date0,time0
	read(iu,*) dt

	read(iu,*)
	read(iu,*)

	read(iu,*)
	read(iu,*) nbox
	read(iu,*)

	allocate(nlayer(nbox),area(nbox),volume0(nbox),depth(nbox))

	do ib=1,nbox
	  read(iu,*) i,nlayer(ib),area(ib),volume0(ib),depth(ib)
	  if( i /= ib ) then
	    write(6,*) i,ib
	    stop 'error stop read_stats: i /= ib'
	  end if
	end do

	read(iu,*)
	read(iu,*) nsect
	read(iu,*)

	allocate(isect(5,nsect))

	do is=1,nsect
	  if( bvers3 ) then
	    read(iu,*) isect(:,is)
	    is_aux = is
	  else
	    read(iu,*) is_aux,isect(:,is)
	  end if
	  if( is /= is_aux ) goto 98
	  if( isect(1,is) > nbox ) goto 97
	  if( isect(2,is) > nbox ) goto 97
	end do
	
	close(iu)

	write(6,*) 'nbox,nsect: ',nbox,nsect
	write(6,*) 'finished reading ',trim(file)

	return
   97	continue
	write(6,*) 'error in box number in section'
	write(6,*) isect(:,is)
	stop 'error stop read_stats: ib > nbox'
   98	continue
	write(6,*) 'inconsistency reading section information'
	write(6,*) is,is_aux
	stop 'error stop read_stats: is /= is_aux'
   99	continue
	write(6,*) 'error opening file: ',trim(file)
	stop 'error stop read_stats: opening file'
	end

!************************************************************

	subroutine read_zinit

	use boxes

	implicit none

	integer :: iu = 1
	character*80 :: file

	integer ib,nb,nvers
	integer i,ierr

	file='boxes_init.txt'
	call check_file_exists(file)
	open(iu,file=file,status='old',form='formatted')
	call check_file(iu,0,nvers,ierr)
	if( ierr < 0 ) stop 'error stop read_zinit: ' //
     &			'reading boxes_init.txt'
	if( nvers_file == 0 ) nvers_file = nvers
	if( nvers /= nvers_file ) then
	  stop 'error stop read_zinit: inconsistent nvers'
	end if

	read(iu,*)
	read(iu,*) nb
	read(iu,*)

	if( nb /= nbox ) then
	  write(6,*) nb,nbox
	  stop 'error stop read_zinit: nb /= nbox'
	end if

	do ib=1,nbox
	  read(iu,*) i,zold(ib)
	  if( i /= ib ) then
	    write(6,*) i,ib
	    stop 'error stop read_zinit: i /= ib'
	  end if
	end do

	vold = volume0 + area*zold

	znew = zold
	vnew = vold

	close(iu)

	write(6,*) 'finished reading ',trim(file)

	end

!************************************************************
!************************************************************
!************************************************************

	subroutine initialize_time

	use boxes
	use mod_date
	use mod_lintp

	implicit none

	integer year0
	character*20 aline

	year0 = date0
	call date_from_year(year0,date0)
	call date_to_atime(date0,time0,atime0)
	call date_format(atime0,aline)

	write(6,*) 'reference time: ',aline
	write(6,*) date0,time0,atime0

	call lintp_set_reference(atime0)

	end

!************************************************************

	subroutine check_time

	use boxes
	use mod_date

	implicit none

	double precision atime
	character*20 aline

	atime = atime0 + dtime
	call date_format(atime,aline)

	!write(6,*) aline,'  ',stime

	if( aline /= stime ) then
	  write(6,*) 'in time record ',irec+1
	  write(6,*) 'expected: ',aline
	  write(6,*) 'read    : ',stime
	  stop 'error stop check_time: time stamp is wrong'
	end if

	end

!************************************************************

	subroutine hydro_out

	use boxes

	implicit none

	double precision atime
	character*20 aline

	integer, save :: iubase = 10
	integer, save :: iuaver = 999
	integer, save :: icall = 0
	integer, save :: nhydro = 4
	double precision :: taver,saver,zaver,vaver
	double precision, allocatable, save :: hydro(:,:)
	character*80 format

	if( icall == 0 ) then
	  allocate(hydro(nbox,nhydro))
	  call scalar_out_header(iubase+1,'water temperature [C]')
	  call scalar_out_header(iubase+2,'salinity [psu]')
	  call scalar_out_header(iubase+3,'current speed [m/s]')
	  call scalar_out_header(iubase+4,'water level [m]')

	  write(iuaver,'(a)')  '# averages of hydro data over all boxes'
          write(iuaver,'(4a)') '#      date_and_time'
     &		// '   temperature'
     &		// '      salinity'
     &		// '   water_level'
     &		// '        volume'

	  icall = 1
	end if

	hydro(:,1) = temp(:)
	hydro(:,2) = salt(:)
	hydro(:,3) = cspeed(:)
	hydro(:,4) = znew(:)

	call scalar_out(iubase,nhydro,hydro)

	call box_average(temp,taver)
	call box_average(salt,saver)
	call box_average(znew,zaver)
	call box_average(vnew,vaver)

	format='(a,3f14.3,e14.6)'
	write(iuaver,format) stime,sngl(taver),sngl(saver)
     &				,sngl(zaver),sngl(vaver)

	end

!************************************************************

	subroutine meteo_out

	use boxes

	implicit none

	double precision atime
	character*20 aline

	integer, save :: iubase = 14
	integer, save :: iuaver = 998
	integer, save :: icall = 0
	integer, save :: nmeteo = 5
	double precision :: taver,saver,haver,waver,caver
	double precision, allocatable, save :: lmeteo(:,:)
	character*80 format

	if( icall == 0 ) then
	  allocate(lmeteo(nbox,nmeteo))
	  call scalar_out_header(iubase+1,'solar radiation [W/m**2]')
	  call scalar_out_header(iubase+2,'air temperature [C]')
	  call scalar_out_header(iubase+3,'realtive humidity [%]')
	  call scalar_out_header(iubase+4,'wind speed [m/s]')
	  call scalar_out_header(iubase+5,'cloud cover [0-1]')

	  write(iuaver,'(a)')  '# averages of meteo data over all boxes'
          write(iuaver,'(4a)') '#      date_and_time'
     &		// '       srad'
     &		// '      atemp'
     &		// '       rhum'
     &		// '     wspeed'
     &		// '     ccover'

	  icall = 1
	end if

	lmeteo(:,1) = meteo(1,:)
	lmeteo(:,2) = meteo(2,:)
	lmeteo(:,3) = meteo(3,:)
	lmeteo(:,4) = meteo(4,:)
	lmeteo(:,5) = meteo(5,:)

	call scalar_out(iubase,nmeteo,lmeteo)

	call box_average(lmeteo(:,1),saver)
	call box_average(lmeteo(:,2),taver)
	call box_average(lmeteo(:,3),haver)
	call box_average(lmeteo(:,4),waver)
	call box_average(lmeteo(:,5),caver)

	format='(a,5f11.3)'
	write(iuaver,format) stime,sngl(saver),sngl(taver)
     &			,sngl(haver),sngl(waver),sngl(caver)

	end

!************************************************************

	subroutine box_average(val,valaver)

	use boxes

	implicit none

	double precision val(nbox)
	double precision valaver

	integer ib
	double precision aver,vol

	aver = 0.
	vol = 0.

	do ib=1,nbox
	  aver = aver + val(ib)*vnew(ib)
	  vol = vol + vnew(ib)
	end do

	valaver = aver / vol
	
	end

!************************************************************

