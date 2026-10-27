
!--------------------------------------------------------------------------
!
!    Copyright (C) 2024-2026  Georg Umgiesser
!
!    This file is part of GESTAS.
!
!    GESTAS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    GESTAS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

!------------------------------------------------------------
!
! box model gestas - administration of box model
!
! info :
!
! needs files boxes*.txt from shyfem simulation
!
! todo :
!
! consider rain and evaporation for volume error
! implement 3D
!
! revision log :
!
! 01.01.2024    ggu     written from scratch
! 13.10.2026    ggu     some new functionalities
! 27.10.2026    ggu     code cleaned
!
!------------------------------------------------------------

!============================================================
	module boxes
!============================================================

	implicit none

	!integer, parameter :: ifreq = 1    !frequency of output (every time)
	!integer, parameter :: ifreq = 24   !frequency of output (daily)
	!integer, parameter :: ifreq = 168  !frequency of output (weekly)
	!integer, parameter :: ifreq = 730  !frequency of output (monthly)
	!integer, parameter :: ifreq = 8766 !frequency of output (yearly)

	double precision, save :: ptime0 = 0
	double precision, save :: pextra = 0	!if not full year available
	integer, save :: pyear = 1	!how many years to simulate
	character*20, save :: pline		!date for perpetual year

!------------------------------------------------------------
! do not change anything beyond this point
!------------------------------------------------------------

	integer, save :: nvers_file = 0	!version of files

	integer, save :: nbox = 0	!total number of boxes
	integer, save :: nsect = 0	!total number of sections
	integer, save :: nosect = 0	!total number of OB sections
	integer, save :: nbound = -1	!total number of OB boundaries
	double precision, save :: dt = 0.

	integer, save :: irec = 0

	integer, save :: date0,time0
	double precision, save :: atime0

	double precision, save :: dtime
	double precision, save :: dlast
	character*20, save :: stime
	character*20, save :: stime_start, stime_end

	integer, save, allocatable :: nlayer(:)
	double precision, save, allocatable :: area(:)
	double precision, save, allocatable :: volume0(:)
	double precision, save, allocatable :: depth(:)
	double precision, save, allocatable :: vinit(:)
	double precision, save, allocatable :: zinit(:)
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

!------------------------------------------------------------

	logical, save :: bsilent = .false.
	logical, save :: bquiet = .false.
	logical, save :: bverbose = .false.
	integer, save :: ifreq = -1		!frequency of output
	integer, save :: imax = 0		!maximum number of time steps

        character*80, save :: boxfiles = ' '

	logical, save :: bconz = .false.
	logical, save :: bwrt = .false.
	logical, save :: bsteady = .false.
	logical, save :: bweutro = .false.
	logical, save :: baquabc = .false.
	logical, save :: bbfm = .false.

	character*80, save :: weutro_input
	integer, save :: wrt_freq

	logical, save :: bpe = .false.	  !use rain/evap for mass error
        logical, save :: blimit = .false. !limit concentration
        logical, save :: badvect = .true. !advect concentration

!============================================================
	end module boxes
!============================================================

	program gestas

	use clo
	use boxes

	implicit none

	logical bout
	integer ierr
	integer iheader
	integer iyear
	double precision verr
	double precision time_start,time_end

	call cpu_time(time_start)

	call gestas_init	!handle command line
	call read_stats
	call initialize_time
	call get_open_boundary_conditions
	call check_bound(nbound,nosect)
	call alloc_arrays
	call read_zinit
	call scalar_init

	if( .not. bquiet ) then
	  write(6,*) 'running with time step ',dt

	  write(6,*) '==========================='
	  write(6,*) 'start of time loop'
	  write(6,*) '==========================='
	end if

	do iyear=1,pyear

	if( pyear > 1 .and. .not. bsilent ) then
	  write(6,*) 'running perpetual year ',iyear
	end if

	do
	  call copy_new2old
	  call read_2d_record(ierr)
	  call read_meteo_record(ierr)
	  if( ierr /= 0 ) exit
	  irec = irec + 1
	  call check_volume(verr)
	  call hydro_write_to_terminal(verr)
	  call hydro_out
	  call meteo_out
	  call fluxes_out
	  call scalar_run
	  if( imax > 0 .and. irec > imax ) exit
	  if( ierr /= 0 ) exit			!finished perpetual
	end do

	call handle_perpetual(dlast,ierr)
	if( ierr /= 0 ) exit

	end do

	if( .not. bquiet ) then
	  write(6,*) '==========================='
	  write(6,*) 'end of time loop'
	  write(6,*) '==========================='

	  write(6,*) irec,' time records read'
	  write(6,*) 'start time of simulation: ',stime_start
	  write(6,*) 'end time of simulation:   ',stime_end
	end if

	call hydro_finalize
	call scalar_finalize

	if( .not. bquiet ) then
	  call cpu_time(time_end)
	  write(6,*) time_end-time_start,' seconds cpu time used'
	end if

	end

!************************************************************

	subroutine hydro_write_to_terminal(verr)

	use boxes

	implicit none

	logical bout
	integer iheader
	double precision verr

        character*80, save :: header =                                  &
     &             ' iterations                dtime'                   &
     &          // '             date and time             verr'

	iheader = 30*max(1,ifreq)

	if( ifreq == 0 ) return

	if( irec == 1 .or. mod(irec,iheader) == 0 ) then
	  write(6,*) trim(header)
	end if
	bout = ( mod(irec,ifreq) == 0 )
	if( bout ) write(6,*) irec,dtime,pline,sngl(verr)

	end

!************************************************************

	subroutine hydro_finalize

	use boxes

	implicit none

	integer n
	character*80 text

	if( bquiet ) return

	call make_max_file('fort.',700,nsect,text)

	write(6,*) 'hydro module has been run'
	write(6,*) '  results for water variables are in fort.101-fort.104'
	write(6,*) '  results for meteo variables are in fort.105-fort.109'
	write(6,*) '  fluxes between boxes are in fort.701-' // trim(text)
	write(6,*) '  averages over boxes are in fort.911 and fort.912'

	end

!************************************************************

	subroutine check_volume(verr)

	use boxes

	implicit none

	double precision verr

	integer is,ib
	double precision aarea,rain,evap
	double precision vdiff,vaver,vaaver,vadiff
	double precision volume(nbox)

	volume = vold

!-------------------------------------------------------
! some debug
!-------------------------------------------------------

	do is=1,nsect
	  !write(6,*) is,isect(:,is)
	end do

!-------------------------------------------------------
! loop over sections to compute new volume
!-------------------------------------------------------

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

!-------------------------------------------------------
! for old version we had extra section for the boundary conditions
!-------------------------------------------------------

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

!-------------------------------------------------------
! consider rain and evaporation (not working yet)
!-------------------------------------------------------

	do ib=1,nbox
	  aarea = area(ib)
	  rain = meteo(6,ib)
	  evap = meteo(7,ib)
	  if( bpe ) volume(ib) = volume(ib) + aarea * dt * ( rain - evap )
	end do

!-------------------------------------------------------
! compute volume error
!-------------------------------------------------------

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

!-------------------------------------------------------
! end of routine
!-------------------------------------------------------

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

	if( .not. bquiet ) write(6,*) 'finished allocating arrays'

	end

!************************************************************

	subroutine check_version(text,nvers)

	use boxes

	implicit none

	character*(*) text
	integer nvers

	if( nvers_file == 0 ) nvers_file = nvers

	if( nvers /= nvers_file ) then
	  write(6,*) trim(text)
	  stop 'error stop check_version: inconsistent nvers'
	end if

	end

!************************************************************

	subroutine check_bound(nbound,nosect)

	implicit none

	integer nbound
	integer nosect

        if( nbound < 0 ) then
          write(6,*) 'nbound = ',nbound
          write(6,*) 'nbound not initialized'
          stop 'error stop check_bound: nbound<0'
        else if( nosect /= nbound ) then
          write(6,*) 'nbound = ',nbound
          write(6,*) 'nosect = ',nosect
          write(6,*) 'nbound and nosect are different'
          stop 'error stop check_bound: nbound/=nosect'
        end if

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
	call date_set_reference(atime0)

	if( .not. bquiet ) then
	  write(6,*) 'reference time: ',aline
	  write(6,*) date0,time0,atime0
	end if

	call lintp_set_reference(atime0)

	end

!************************************************************

	subroutine gestas_copyright(routine)

	implicit none

	character*(*) routine

	character*50 copyright,acronym

	acronym = ' '
	copyright = 'Copyright (C) Georg Umgiesser 2024-2026'

        write(6,*) '---------------------------------------'
        write(6,*) 'GESTAS' ! // ' - ' // trim(acronym)
        write(6,*) copyright
        write(6,*) routine
        write(6,*) '---------------------------------------'

	end

!************************************************************
!************************************************************
!************************************************************

	subroutine gestas_init

        use clo
	use boxes

        implicit none

	logical bnoadvect
        integer nfiles

        call clo_init('gestas','[box-files]','2.0')

        call clo_add_info('box model for shyfem')

        call clo_add_sep('output options')

        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('silent',.false.,'be silent')

        call clo_add_option('freq ifreq',ifreq &
     &				,'frequency of output [time steps]')

        call clo_add_sep('modules')

        call clo_add_option('weutro input-files',' ','run weutro module')
        call clo_add_option('aquabc',.false.,'run aquabc module')
        call clo_add_option('bfm',.false.,'run bfm module')
        call clo_add_option('wrt freq',-1,'run wrt module ' &
     &				// 'with reset frequency freq')
        call clo_add_option('conz',.false.,'run conz module')
        call clo_add_option('steady',.false.,'run steady state module')

        call clo_add_sep('modifiers')

        call clo_add_option('noadvect',.false.,'do not advect scalars')
        call clo_add_option('limit',.false.,'limit scalars between min/max')
        call clo_add_option('imax imax',0,'do only imax time steps')
        call clo_add_option('years iyear',0,'run iyear perpetual years')

        call clo_hide_next_options
        call clo_add_sep('extra options')
        call clo_add_option('pextra pextra',0,'use pextra for perpetual years')
        call clo_show_next_options

        call clo_add_sep('additional information')
        call clo_add_com('  if no box-files is given,' &
     &			// ' standard boxes_*.txt are used')
 
        call clo_parse_options

        call clo_get_option('verbose',bverbose)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)

        call clo_get_option('freq',ifreq)
	if( ifreq == -1 ) then
	  if( bsilent ) then
	    ifreq = 0
	  else if( bverbose ) then
	    ifreq = 1
	    ifreq = 730
	    ifreq = 168
	  else
	    ifreq = 730
	  end if
	end if

        if( bsilent ) bquiet = .true.
        if( bquiet ) bverbose = .false.

        call clo_get_option('weutro',weutro_input)
	if( weutro_input /= ' ' ) bweutro = .true.
        call clo_get_option('aquabc',baquabc)
        call clo_get_option('bfm',bbfm)
        call clo_get_option('wrt',wrt_freq)
	if( wrt_freq >= 0 ) bwrt = .true.
        call clo_get_option('conz',bconz)
        call clo_get_option('steady',bsteady)

        call clo_get_option('noadvect',bnoadvect)
	badvect = .not. bnoadvect
        call clo_get_option('limit',blimit)
        call clo_get_option('imax',imax)
        call clo_get_option('years',pyear)
	if( pyear < 1 ) pyear = 1
        call clo_get_option('pextra',pextra)

	nfiles = clo_number_of_files()
	if( nfiles > 0 ) then
          call clo_get_file(1,boxfiles)
	  write(6,*) 'box-files: ',trim(boxfiles)
	  stop 'error stop gestas_init: box-files not yet ready'
	end if

        if( .not. bsilent ) then
          call gestas_copyright('gestas - box model for shyfem')
        end if

	end

!************************************************************
!************************************************************
!************************************************************

