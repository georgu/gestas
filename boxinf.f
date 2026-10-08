
!***************************************************************
!
! boxinf - routine to check one and compare two box files
!
!***************************************************************

	module mod_boxinf

	logical, save :: bdebug = .false.

	logical, save :: bsilent = .false.
	logical, save :: bquiet = .false.
	logical, save :: bnodiff = .false.
	logical, save :: bverbose = .false.
	logical, save :: bnostop = .false.
	logical, save :: bstop = .true.
	logical, save :: b2dwrite = .false.
	double precision, save :: maxdiff = 0.

	double precision, save :: dstart,dend,dt
	character*20, save :: astart,aend

	end module mod_boxinf

!***************************************************************

	program boxinf

	use clo
	use mod_boxinf

	implicit none

	integer nc,ierr

	call boxinf_init

        nc = clo_number_of_files()

        if( nc == 0 ) then
          call clo_usage
        else if( nc == 1 ) then
          call read_box_file
        else if( nc == 2 ) then
          call compare_box_files(ierr)
        else
          write(6,*) 'nc = ',nc
          stop 'error stop boxinf: wrong number of files'
        end if

        if( ierr > 0 ) then
          if( ierr == 99 ) ierr = 100   !terrible hack - FIXME
          call exit(ierr)
        else
          call exit(99)
        end if

	end

!***************************************************************

	subroutine read_box_file

	use mod_boxinf
	use mod_boxutil
	use clo

	implicit none

	logical bwrite				!write single time records
	logical b2d,b3d,bmeteo,bhydro,bobs
	double precision dtime
	double precision rain,evap,period,dyear
	integer nbox,ibox,ftype,i,nvar,nsect,lmax,nvar3d,nvar2d
	integer iu,iuw
	integer ios,ierr,nrec,nvals,ival,nvers,nobs,ib,nout,nobc
	character*20 aline
	character*80 sh,format,format77,formatb
	character*80 name_one

	double precision, allocatable :: acum(:,:),aver(:)
	double precision, allocatable :: boxvals(:,:)
	double precision, allocatable :: boxvals3d(:,:,:)
	double precision, allocatable :: vals2d(:,:)
	double precision, allocatable :: sectvals(:,:)
	double precision, allocatable :: obcvals(:,:)
	double precision, allocatable :: sectvals3d(:,:,:)
	double precision, allocatable :: obcvals3d(:,:,:)
	double precision, allocatable :: sectout(:,:)
	integer, allocatable :: lbmax(:)
	integer, allocatable :: isvals(:,:)
	integer, allocatable :: iovals(:,:)
	integer, allocatable :: isout(:,:)

	bwrite = bverbose
	nrec = 0
	nbox = 0
	nvar = 0
	nsect = 0
	nobs = 0
	nobc = 0
	ierr = 0
	dt = 0.

!--------------------------------------------------------------
! open file
!--------------------------------------------------------------

        call clo_get_file(1,name_one)

	iu = 1
	open(iu,file=name_one,status='old',form='formatted',iostat=ios)

        if( ios /= 0 ) then
	  write(6,*) 'cannot open file: ',trim(name_one)
	  stop 'error opening file'
	end if

!--------------------------------------------------------------
! check file type
!--------------------------------------------------------------

	call read_file_header(iu,ftype,nvers)
	b2d = ( ftype == 2 )
	b3d = ( ftype == 3 )
	bhydro = b2d .or. b3d
	bmeteo = ( ftype == 1 )
	bobs = ( nvers <= 5 )
	if( .not. (b2d .or. b3d .or. bmeteo ) ) then
	  stop 'error stop: ftype/=[123]'
	end if
	call make_formats(ftype,sh,format,format77,formatb)

	call get_box_dimensions(iu,ftype,nbox,nvar
     &				,nsect,nobc,lmax,ierr)
	nobs = abs(nobc)
	if( ierr /= 0 ) stop 'error stop: cannot determine dimensions'

	if( bverbose ) then
	  write(6,*) 'ftype,nvers:  ',ftype,nvers
	  write(6,*) 'nbox,nvar:    ',nbox,nvar
	  write(6,*) 'nsect,nobs:   ',nsect,nobs
	  write(6,*) 'lmax:         ',lmax
	end if

!--------------------------------------------------------------
! allocate arrays
!--------------------------------------------------------------

	if( b3d ) then
	  nvar = 7
	  allocate(boxvals(nvar,nbox))
	  allocate(sectvals(3,nsect))
	end if
	if( b2dwrite ) then
	  nout = nsect + nobs
	  allocate(sectout(3,nout))
	  allocate(isout(3,nout))
	end if
	if( nvers > 5 ) then		!avoid compiler error
	  allocate(obcvals(3,nobs))
	  allocate(iovals(3,nobs))
	end if

	allocate(acum(nvar,nbox),aver(nvar))
	acum = 0.

!--------------------------------------------------------------
! open output file
!--------------------------------------------------------------

	if( b2dwrite ) then
	  iuw = 3
	  open(iuw,file='new_2d.txt',status='unknown',form='formatted')
	  call write_file_header(iuw,2)
	  write(6,*) 'writing 2D file'
	end if

!--------------------------------------------------------------
! time loop
!--------------------------------------------------------------

	do
	  call read_time_record(iu,dtime,aline,ierr)
	  if( ierr /= 0 ) exit

	  call handle_time(nrec,dtime,aline)

	  if( bmeteo ) then
	    call read_box_2d(iu,nbox,nvar,boxvals)
	  else if( b2d ) then
	    call read_box_2d(iu,nbox,nvar,boxvals)
	    call read_flux_2d(iu,nsect,isvals,sectvals)
	    if( bobs ) call read_flux_2d(iu,nobs,iovals,obcvals)
	  else if( b3d ) then
	    call read_box_3d(iu,nbox,nvar3d,lmax,lbmax,boxvals3d,vals2d)
	    call read_flux_3d(iu,nsect,lmax,isvals,sectvals3d)
	    if( bobs ) call read_flux_3d(iu,nobs,lmax,iovals,obcvals3d)
	  end if

	  if( b3d ) then	!average of layers
            call aver3dbox(nbox,nvar,nvar3d,lmax,lbmax
     &				,boxvals3d,vals2d,boxvals)
            call aver3dflux(nsect,lmax,isvals
     &				,sectvals3d,sectvals)
	  end if

	  acum(:,:) = acum(:,:) + boxvals(:,:)
	  do ival=1,nvar	!average over boxes
	    aver(ival) = sum(boxvals(ival,:)) / nbox
	  end do

	  if( bwrite ) write(6,*) 'time = ',aline,dtime
	  write(77,format77) aline,aver
	  !call debug_boxvals(88,aline,nbox,nvar,boxvals)
	  !if( .not. bmeteo ) call debug_isvals(89,nsect,isvals)
	  !if( bobs ) call debug_isvals(89,nobs,iovals)
	  !stop


	  if( b2dwrite ) then
	    nvar2d = 7
	    call write_time_record(iuw,dtime,aline)
	    call write_box_records_2d(iuw,nbox,nvar2d,boxvals)
	    call combine_2d_fluxes(iu,nbox,nsect,nobs,isvals,iovals
     &					,sectvals,obcvals
     &					,nout,isout,sectout)
	    call write_flux_records_2d(iuw,nout,isout,sectout)
	  end if
	end do

!--------------------------------------------------------------
! end of time loop
!--------------------------------------------------------------

	if( ierr > 0 ) stop 'error stop: read error'

	if(  bsilent ) stop

	write(6,*) 'final message:'
	write(6,*) 'type of file:      ',ftype
	write(6,*) 'version of file:   ',nvers
	write(6,*) 'boxes:             ',nbox
	write(6,*) 'sections:          ',nsect
	write(6,*) 'boundary sections: ',nobs
	write(6,*) 'open boundaries:   ',abs(nobc)
	write(6,*) 'box values:        ',nvar
	write(6,*) 'levels read:       ',lmax
	write(6,*) 'records read:      ',nrec
	write(6,*) 'start date of simulation:  ',astart
	write(6,*) 'end date of simulation:    ',aend
	write(6,*) 'time step:              ',dt

	if( .not. bquiet ) then
	do ival=1,nvar
	  aver(ival) = sum(acum(ival,:)) / (nbox*nrec)
	end do
	write(6,*) 'time average over box values:'
	write(6,'(a)') trim(sh)
	write(6,format) aver

	write(6,*) 'time average of box values:'
	write(6,'(2a)') '  box',trim(sh)
	acum = acum / nrec
	do ib=1,nbox
	  aver(:) = acum(:,ib)
	  write(6,formatb) ib,aver
	end do

	if( bmeteo ) call compute_rain(aver)
	end if

	if( b2dwrite ) write(6,*) '2D file has been written to new_2d.txt'

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!***************************************************************

	subroutine compare_box_files(ierr)

	use mod_boxinf
	use mod_boxutil
	use clo

	implicit none

	integer ierr

	logical bwrite				!write single time records
	logical b2d,b3d,bmeteo,bhydro,bobs
	double precision dtime,dtime1,dtime2
	double precision rerr,rerr1,rerr2,rerrmax
	integer ibox,i,iu1,iu2
	integer i1,i2
	integer ftype,ftype1,ftype2
	integer nvers,nvers1,nvers2
	integer nbox,nbox2
	integer nvar,nvar2
	integer nsect,nsect2
	integer lmax,lmax2
	integer nobs,nobs2
	integer nobc,nobc2
	integer ios,nrec,ival
	integer ierr1,ierr2
	character*20 aline
	character*80 sh,format,format77,formatb
	character*80 name_one,name_two

	double precision, allocatable :: boxvals1(:,:),boxvals2(:,:)
	double precision, allocatable :: sectvals1(:,:),sectvals2(:,:)
	integer, allocatable :: isvals1(:,:),isvals2(:,:)

	bwrite = bverbose
	nrec = 0
	ierr = 0
	dt = 0.

!--------------------------------------------------------------
! open files
!--------------------------------------------------------------

        call clo_get_file(1,name_one)
        call clo_get_file(2,name_two)

	iu1 = 1
	open(iu1,file=name_one,status='old',form='formatted',iostat=ios)

        if( ios /= 0 ) then
	  write(6,*) 'cannot open file: ',trim(name_one)
	  stop 'error opening file'
	end if

	iu2 = 2
	open(iu2,file=name_two,status='old',form='formatted',iostat=ios)

        if( ios /= 0 ) then
	  write(6,*) 'cannot open file: ',trim(name_two)
	  stop 'error opening file'
	end if

!--------------------------------------------------------------
! check file types
!--------------------------------------------------------------

	call read_file_header(iu1,ftype,nvers)
	call read_file_header(iu2,ftype2,nvers2)
	call compare_params('ftype',ftype,ftype2)
	call compare_params('nvers',nvers,nvers2)
	b2d = ( ftype == 2 )
	b3d = ( ftype == 3 )
	bhydro = b2d .or. b3d
	bmeteo = ( ftype == 1 )
	bobs = ( nvers <= 5 )
	if( .not. (b2d .or. b3d .or. bmeteo ) ) then
	  stop 'error stop: ftype/=[123]'
	end if
	call make_formats(ftype,sh,format,format77,formatb)

	call get_box_dimensions(iu1,ftype,nbox,nvar
     &				,nsect,nobc,lmax,ierr1)
	nobs = abs(nobc)
	call get_box_dimensions(iu2,ftype,nbox2,nvar2
     &				,nsect2,nobc2,lmax2,ierr2)
	nobs2 = abs(nobc2)
	if( ierr1 /= 0 .or. ierr2 /= 0 ) then
	  write(6,*) ierr1,ierr2
	  stop 'error stop: cannot determine dimensions'
	end if

	call compare_params('nbox',nbox,nbox2)
	call compare_params('nvar',nvar,nvar2)
	call compare_params('nsect',nsect,nsect2)
	call compare_params('nobs',nobs,nobs2)
	call compare_params('lmax',lmax,lmax2)

	if( bverbose ) then
	  write(6,*) 'ftype,nvers:  ',ftype,nvers
	  write(6,*) 'nbox,nvar:    ',nbox,nvar
	  write(6,*) 'nsect,nobs:   ',nsect,nobs
	  write(6,*) 'lmax:         ',lmax
	end if

!--------------------------------------------------------------
! time loop
!--------------------------------------------------------------

	do
	  call read_time_record(iu1,dtime,aline,ierr)
	  if( ierr /= 0 ) exit
	  call read_time_record(iu2,dtime2,aline,ierr)
	  if( ierr /= 0 ) exit
	  if( dtime /= dtime2 ) stop 'error stop: dtime1/=dtime2'

	  call handle_time(nrec,dtime,aline)

	  call read_box_2d(iu1,nbox,nvar,boxvals1)
	  call read_box_2d(iu2,nbox,nvar,boxvals2)

	  call compare_values('box',nvar,nbox,boxvals1,boxvals2,rerr)
	  rerrmax = max(rerrmax,rerr)

	  rerr1 = 0
	  rerr2 = 0
	  if( b2d ) then
	    call read_flux_2d(iu1,nsect,isvals1,sectvals1)
	    call read_flux_2d(iu2,nsect,isvals2,sectvals2)
	    call compare_values('flux',3,nsect,sectvals1,sectvals2,rerr1)
	    rerrmax = max(rerrmax,rerr1)
	    if( bobs ) then
	      call read_flux_2d(iu1,nsect,isvals1,sectvals1)
	      call read_flux_2d(iu2,nsect,isvals2,sectvals2)
	      call compare_values('flux',3,nsect,sectvals1,sectvals2
     &						,rerr2)
	      rerrmax = max(rerrmax,rerr2)
	    end if
	  end if

	  if( bwrite ) write(6,1000) aline,rerr,rerr1,rerr2
 1000	  format(a,3g16.8)
	end do

!--------------------------------------------------------------
! end of time loop
!--------------------------------------------------------------

	if( ierr > 0 ) stop 'error stop: read error'

	write(6,*) 'final message:'
	write(6,*) 'boxes:        ',nbox
	write(6,*) 'records read: ',nrec
	write(6,*) 'values read:  ',nvar
	write(6,*) 'start date of simulation:  ',astart
	write(6,*) 'end date of simulation:    ',aend
	write(6,*) 'time step:              ',dt
	write(6,*) 'maximum relative difference:    ',rerrmax

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*****************************************************

	subroutine compare_values(type,n1,n2,vals1,vals2,rerr)

	use mod_boxinf

	implicit none

	character*(*) type
	integer n1,n2
	double precision vals1(n1,n2)
	double precision vals2(n1,n2)
	double precision rerr

	integer i1,i2
	double precision atot,diff,aux

	rerr = 0.

	if( type == 'box' ) then
	  do i2=1,n2
	    do i1=1,n1
	      diff = abs( vals1(i1,i2) - vals2(i1,i2) )
	      atot = abs(vals1(i1,i2)) + abs(vals2(i1,i2))
	      if( diff > 0 .and. atot > 0. ) then
	        aux = 2.*diff/atot
	        rerr = max(rerr,aux)
		if( bdebug ) write(6,*) i1,i2,diff,aux
	      end if
	    end do
	  end do
	else if( type == 'flux' ) then
	  do i2=1,n2
	    do i1=1,n1
	      diff = abs( vals1(i1,i2) - vals2(i1,i2) )
	      atot = abs(vals1(i1,i2)) + abs(vals2(i1,i2))
	      if( diff > 0 .and. atot > 0. ) then
	        aux = 2.*diff/atot
	        rerr = max(rerr,aux)
		if( bdebug ) write(6,*) i1,i2,diff,aux
	      end if
	    end do
	  end do
	else
	  stop 'error stop compare_values: unknown section type'
	end if
	
	end

!*****************************************************

	subroutine handle_time(nrec,dtime,aline)

	use mod_boxinf

	implicit none

	integer nrec
	double precision dtime
	character*(*) aline

	if( nrec == 0 ) then
	  dstart = dtime
	  astart = aline
	else if( nrec == 1 ) then
	  dt = dtime - dstart
	else
	  if( dt /= dtime - dend ) then
	    write(6,*) 'changing time step: ',dt,dtime-dend
	    stop 'error stop: time step'
	  end if
	end if

	dend = dtime
	aend = aline

	nrec = nrec + 1

	end

!*****************************************************

	subroutine compute_rain(vals)

	implicit none

	double precision vals(*)	!averaged values

	double precision, parameter :: dyear = 365.25
	double precision :: rain,evap

	rain = vals(6) * dyear		!converts from mm/day to mm/year
	evap = vals(7) * dyear
	write(6,*) 'average rain per year [mm]: ',rain
	write(6,*) 'average evap per year [mm]: ',evap

	end
	
!*****************************************************

	subroutine make_formats(ftype,sh,format,format77,formatb)

	implicit none

	integer ftype
	character*80 sh,format,format77,formatb

	if( ftype == 0 ) then
          sh = '  init_eta'
	  format   = '(    f10.3)'
	  format77 = '(a20,f10.3)'
	else if( ftype == 1 ) then
          sh = '    srad    tair    rhum  wspeed'
     +		 // '     cc   rain   evap    ice'
	  format   = '(    4f8.3,4f7.3)'
	  formatb  = '(    i5,4f8.3,4f7.3)'
	  format77 = '(a20,4f8.3,4f7.3)'
	else if( ftype == 2 ) then
          sh = '     temp     salt  c_speed  act_eta'
     +		 // ' aver_eta  bstress      volume'
	  format   = '(    6f9.3,1e12.4)'
	  formatb  = '(    i5,6f9.3,1e12.4)'
	  format77 = '(a20,6f8.3,1e12.4)'
	else if( ftype == 3 ) then
          sh = '     temp     salt  c_speed  act_eta'
     +		 // ' aver_eta  bstress      volume'
	  format   = '(    6f9.3,1e12.4)'
	  formatb  = '(    i5,6f9.3,1e12.4)'
	  format77 = '(a20,6f8.3,1e12.4)'
	else if( ftype == 7 ) then
          sh = '        area   depth'
	  format   = '(    e12.4,f8.3)'
	  format77 = '(a20,e12.4,f8.3)'
	else if( ftype == 8 ) then
          sh = '   layers        area      volume    depth'
	  format   = '(    1f9.3,2e12.4,1f9.3)'
	  format77 = '(a20,1f9.3,2e12.4,1f9.3)'
	else
	  write(6,*) 'cannot handle this file type: ',ftype
	  !stop 'error stop make_formats: ftype not supported'
	end if

	end

!*****************************************************

        subroutine aver3dbox(nbox,nval,nval3d,lmax,lbmax
     &				,boxvals3d,vals2d,boxvals)

	implicit none

	integer nbox,nval,nval3d,lmax
	integer lbmax(nbox)
	double precision boxvals3d(lmax,nval3d,nbox)
	double precision vals2d(3,nbox)
	double precision boxvals(nval,nbox)

	integer ib,l,lb
	double precision vol,volacum
	double precision acum(nval3d)

	if( nval3d /= 4 ) stop 'error stop aver3d: nval3d/=4'

	do ib=1,nbox
	  lb = lbmax(ib)
	  volacum = 0.
	  acum = 0.
	  do l=1,lb
	    vol = boxvals3d(l,4,ib)
	    volacum = volacum + vol
	    acum(1:3) = acum(1:3) + boxvals3d(l,1:3,ib) * vol
	  end do
	  acum(1:3) = acum(1:3) / volacum	!do not average volume

	  boxvals(1:3,ib) = acum(1:3)
	  boxvals(4,ib) = vals2d(2,ib)
	  boxvals(5,ib) = vals2d(1,ib)
	  boxvals(6,ib) = vals2d(3,ib)
	  boxvals(7,ib) = volacum
	end do

	end

!*****************************************************

        subroutine aver3dflux(nsect,lmax,isvals
     &				,sectvals3d,sectvals)

	implicit none

	integer nsect,lmax
	integer isvals(3,nsect)
	double precision sectvals3d(lmax,3,nsect)
	double precision sectvals(3,nsect)

	integer is,l,ls
	double precision acum(3)

	do is=1,nsect
	  ls = isvals(3,is)
	  acum = 0.
	  do l=1,ls
	    acum(:) = acum(:) + sectvals3d(l,:,is)
	  end do
	  sectvals(:,is) = acum(:)
	end do

	end

!*****************************************************

	subroutine combine_2d_fluxes(iu,nbox,nsect,nobs,isvals,iovals
     &					,sectvals,obcvals
     &					,nout,isout,sectout)

	implicit none

	integer iu
	integer nbox,nsect,nobs
	integer isvals(3,nsect)
	integer iovals(3,nobs)
	double precision sectvals(3,nsect)
	double precision obcvals(3,nobs)
	integer nout
	integer isout(3,nout)
	double precision sectout(3,nout)

	integer is,it

	do is=1,nsect
	  isout(:,is) = isvals(:,is)
	  sectout(:,is) = sectvals(:,is)
	end do

	it = nsect
	do is=1,nobs
	  if( any( iovals(1:2,is) == 0 ) ) cycle
	  it = it + 1
	  isout(:,it) = iovals(:,is)
	  sectout(:,it) = obcvals(:,is)
	end do

	nout = it

	end

!*****************************************************

	subroutine compare_params(text,n1,n2)

	implicit none

	character*(*) text
	integer n1,n2

	if( n1 /= n2 ) then
	  stop 'error stop: ' // trim(text) // ' are different'
	end if

	end

!*****************************************************

	subroutine debug_isvals(iu,nsect,isvals)

	implicit none

	integer iu
	integer nsect
	integer isvals(3,nsect)

	integer is

	write(iu,*) nsect
	do is=1,nsect
	  write(iu,*) is,isvals(:,is)
	end do

	end

!*****************************************************

	subroutine debug_boxvals(iu,aline,nbox,nvar,boxvals)

	implicit none

	integer iu
	character*20 aline
	integer nbox,nvar
	double precision boxvals(nvar,nbox)

	integer ib

	write(iu,'(a20,i10)') aline,nbox
	do ib=1,nbox
	  write(iu,*) boxvals(1,ib),boxvals(2,ib),boxvals(7,ib)
	end do

	end

!*****************************************************

        subroutine boxinf_init

        use mod_boxinf
        use clo

        implicit none

        logical baux
        character*80 version

        version = '1.0'

        call clo_init('boxinf','box-file(s)',trim(version))

        call clo_add_info('checks box files')

        call clo_add_sep('general options:')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('nodiff',.false.,'do not show differences')
        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('nostop',.false.,'do not stop at error')
        call clo_add_option('write2d',.false.,'write 2D file')
        call clo_add_option('maxdiff',0.,'maximum tolerated difference')

        call clo_parse_options

        call clo_get_option('silent',bsilent)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('nodiff',bnodiff)
        call clo_get_option('verbose',bverbose)
        call clo_get_option('nostop',baux)
        call clo_get_option('write2d',b2dwrite)
        call clo_get_option('maxdiff',maxdiff)

        if( baux ) bstop = .false.
        if( bsilent ) bquiet = .true.
        if( bquiet ) bverbose = .false.

        end

!*****************************************************

