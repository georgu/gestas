
!------------------------------------------------------------
!
! box model estas2
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

	integer, save :: nbox = 0	!total number of boxes
	integer, save :: nsect = 0	!total number of sections
	integer, save :: nosect = 0	!total number of OB sections
	integer, save :: nbound = -1	!total number of OB boundaries
	real, save :: dt = 0.

	integer, save :: irec = 0

	integer, save :: date0,time0
	double precision, save :: atime0

	double precision, save :: dtime
	character*20, save :: stime

	integer, save, allocatable :: nlayer(:)
	real, save, allocatable :: area(:)
	real, save, allocatable :: volume0(:)
	real, save, allocatable :: depth(:)
	integer, save, allocatable :: isect(:,:)	!section index
	integer, save, allocatable :: iosect(:,:)	!OB section index

	real, save, allocatable :: zold(:)
	real, save, allocatable :: znew(:)
	real, save, allocatable :: vold(:)
	real, save, allocatable :: vnew(:)

	real, save, allocatable :: temp(:)
	real, save, allocatable :: salt(:)
	real, save, allocatable :: cspeed(:)

	real, save, allocatable :: meteo(:,:)

	real, save, allocatable :: flux(:,:)		!section fluxes
	real, save, allocatable :: oflux(:,:)		!OB section fluxes

!============================================================
	end module boxes
!============================================================

	program estas2

	use boxes

	implicit none

	integer ierr
	real verr
	real time_start,time_end

	call cpu_time(time_start)

	call read_stats
	call initialize_time
	call get_open_boundary_conditions
	call alloc_arrays
	call read_zinit
	call scalar_init

	do
	  call copy_new2old
	  call read_time_record(ierr)
	  call read_meteo_record(ierr)
	  if( ierr /= 0 ) exit
	  irec = irec + 1
	  call check_volume(verr)
	  write(6,*) irec,dtime,stime,verr
	  call hydro_out
	  call scalar_run
	  !if( irec > 3 ) exit
	end do

	write(6,*) irec,' time records read'

	call scalar_finalize

	call cpu_time(time_end)
	write(6,*) time_end-time_start,' seconds cpu time used'

	end

!************************************************************

	subroutine check_volume(verr)

	use boxes

	implicit none

	real verr

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

	integer nb,ns,nvars,nob

	iu = 11
	file='boxes_2d.txt'
	open(iu,file=file,status='old',form='formatted')
	call check_file(iu,2,ierr)
	if( ierr < 0 ) return

	call skip_lines(iu,3)		!skip time record
	read(iu,*) nb,nvars
	call skip_lines(iu,nb+2)	!skip boxes
	read(iu,*) ns
	call skip_lines(iu,ns+2)	!skip sections
	read(iu,*) nob
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
	real a

	nob = 0		!number of open boundary conditions

	if( bfirst ) then
	  iu = 2
	  file='boxes_2d.txt'
	  open(iu,file=file,status='old',form='formatted')
	  call check_file(iu,2,ierr)
	  if( ierr < 0 ) return
	end if
	!call check_file(iu,2,ierr)
	!if( ierr < 0 ) return

	read(iu,'(a)',end=1)
	read(iu,*) dtime,stime

	call check_time

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
	  read(iu,*) i1,i2,(flux(ii,is),ii=1,3)
	  nob = min(nob,i1,i2)
	end do

	read(iu,*)
	read(iu,*) no

	if( bfirst ) then
	  nosect = no
	  allocate(oflux(3,nosect),iosect(2,nosect))
	  write(6,*) 'open boundary sections: ',nosect
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

	nob = -nob
	if( nbound == 0 ) then
	  nbound = nob
	  write(6,*) 'open boundary conditions: ',nbound
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
	integer nb,ns,no,nob,nvars
	real a

	nob = 0		!number of open boundary conditions

	if( iu == 0 ) then
	  iu = 1
	  file='boxes_meteo.txt'
	  open(iu,file=file,status='old',form='formatted')
	  call check_file(iu,1,ierr)
	  if( ierr < 0 ) return
	end if
	!call check_file(iu,1,ierr)
	!if( ierr < 0 ) return

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

	subroutine check_file(iu,ft,ierr)

	implicit none

	integer iu
	integer ft
	integer ierr

	integer idbox,ftype,nvers
	integer, parameter :: idbox_expect = 473226
	integer, parameter :: nvers_expect = 5
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
	if( nvers /= nvers_expect ) goto 99

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

	integer ib,is,ierr
	integer i,ii,is_aux

	file='boxes_stats.txt'
	open(iu,file=file,status='old',form='formatted',err=99)

	call check_file(iu,8,ierr)

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
	  read(iu,*) is_aux,isect(:,is)
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

	integer ib,nb
	integer i,ierr

	file='boxes_init.txt'
	open(iu,file=file,status='old',form='formatted')

	call check_file(iu,0,ierr)
	if( ierr < 0 ) return

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

	implicit none

	integer year0
	character*20 aline

	year0 = date0
	call date_from_year(year0,date0)
	call date_to_atime(date0,time0,atime0)
	call date_format(atime0,aline)

	write(6,*) 'reference time: ',aline
	write(6,*) date0,time0,atime0

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
	integer, save :: icall = 0
	integer, save :: nhydro = 4
	double precision, allocatable, save :: hydro(:,:)

	if( icall == 0 ) then
	  allocate(hydro(nbox,nhydro))
	  call scalar_out_header(iubase+1,'temperature')
	  call scalar_out_header(iubase+2,'salinity')
	  call scalar_out_header(iubase+3,'current speed')
	  call scalar_out_header(iubase+4,'water level')
	  icall = 1
	end if

	hydro(:,1) = temp(:)
	hydro(:,2) = salt(:)
	hydro(:,3) = cspeed(:)
	hydro(:,4) = znew(:)

	call scalar_out(iubase,nhydro,hydro)

	end

!************************************************************

	subroutine box_average
	end
!************************************************************

