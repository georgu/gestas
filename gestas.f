
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
!
! revision log :
!
! 01.01.2024    ggu     written from scratch
! 13.10.2026    ggu     some new functionalities
!
!------------------------------------------------------------

!============================================================
	module boxes
!============================================================

	implicit none

	logical, save :: bpe = .false.	  !use rain/evap for mass error
        logical, save :: blimit = .false. !limit concentration
        logical, save :: badvect = .true. !advect concentration

	!integer, save :: imax = 50	!maximum number of time steps
	integer, save :: imax = 0	!maximum number of time steps

	!integer, parameter :: ifreq = 1    !frequency of output (every time)
	!integer, parameter :: ifreq = 24   !frequency of output (daily)
	!integer, parameter :: ifreq = 168  !frequency of output (weekly)
	integer, parameter :: ifreq = 730  !frequency of output (monthly)
	!integer, parameter :: ifreq = 8766 !frequency of output (yearly)

	double precision, save :: ptime0 = 0
	double precision, save :: pextra = 3600	!if not full year available
	double precision, save :: pyear = 1	!how many years to simulate
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
	call check_bound(nbound,nosect)
	call alloc_arrays
	call read_zinit
	call scalar_init

	write(6,*) 'running with time step ',dt

	write(6,*) '==========================='
	write(6,*) 'start of time loop'
	write(6,*) '==========================='
	write(6,*) trim(header)

	do

	do
	  call copy_new2old
	  call read_2d_record(ierr)
	  call read_meteo_record(ierr)
	  if( ierr /= 0 ) exit
	  irec = irec + 1
	  call check_volume(verr)
	  if( ifreq /= 0 ) bout = ( mod(irec,ifreq) == 0 )
	  if( mod(irec,iheader) == 0 ) write(6,*) trim(header)
	  if( bout ) write(6,*) irec,dtime,pline,sngl(verr)
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

	write(6,*) '==========================='
	write(6,*) 'end of time loop'
	write(6,*) '==========================='

	write(6,*) irec,' time records read'
	write(6,*) 'start time of simulation: ',stime_start
	write(6,*) 'end time of simulation:   ',stime_end

	call hydro_finalize
	call scalar_finalize

	call cpu_time(time_end)
	write(6,*) time_end-time_start,' seconds cpu time used'

	end

!************************************************************

	subroutine hydro_finalize

	use boxes

	implicit none

	integer n
	character*80 text

	call make_max_file('fort.',700,nsect,text)

	write(6,*) 'results for water variables are in fort.101-fort.104'
	write(6,*) 'results for meteo variables are in fort.105-fort.109'
	write(6,*) 'fluxes between boxes are in fort.701-' // trim(text)
	write(6,*) 'averages over boxes are in fort.911 and fort.912'

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

	write(6,*) 'finished allocating arrays'

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

	write(6,*) 'reference time: ',aline
	write(6,*) date0,time0,atime0

	call lintp_set_reference(atime0)

	end

!************************************************************
!************************************************************
!************************************************************

