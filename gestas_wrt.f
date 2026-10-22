
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
! box model gestas - WRT (water residence time) computations
!
! this simulates the water residence time and age of the water
!
!	scal(:,1)	concentration for WRT
!	scal(:,2)	concentration for confinement
!	scal(:,3)	age
!
!	fort.201	concentration for WRT
!	fort.202	concentration for confinement
!	fort.203	age
!	fort.204	WRT computed through integration
!	fort.205	confinement computed through integration
!	fort.206	nothing usefull
!	fort.207	WRT computed through logarithmic regression
!	fort.208	confinement computed through logarithmic regression
!	fort.209	nothing usefull
!
!------------------------------------------------------------

!============================================================
        module boxes_wrt
!============================================================

        implicit none

        integer, parameter :: nscal = 3		!total number of scalars

        integer, save :: iubase = 200		!unit for output
        integer, save :: iumass = 920		!unit for mass output

	logical, parameter :: bgoro = .false.	!goro application

	!double precision, save :: dtwrt = 90	!reset time for WRT (in days)
	!double precision, save :: dtwrt = 20	!reset time for WRT (in days)
	!double precision, save :: dtwrt = 180	!reset time for WRT (in days)
	double precision, save :: dtwrt = 90	!reset time for WRT (in days)
	double precision, save :: wrtmax = 999	!maximum for WRT

	integer, save :: iadv(nscal) = 1	!advection of scalars

!	-----------------------------------------------------

	integer, save :: ifrqwrt = 0
	integer, save :: iacc = 0

	double precision, save :: dtime_reset = 0
	double precision, save :: tacu = 0
	double precision, save, allocatable :: wrtacc(:,:)
	double precision, save, allocatable :: wrtlog(:,:)

	double precision, parameter :: spd = 86400	!seconds per day

	double precision, save, allocatable :: scal(:,:) !scalar variables
	double precision, save, allocatable :: scal_bound(:,:)	!scalar BC

!============================================================
        end module boxes_wrt
!============================================================

!************************************************************

	subroutine wrt_run

! runs wrt routines

	use boxes
	use boxes_wrt

	implicit none

	call wrt_boundary_condition
	call wrt_react
	call scalar_advect(nscal,iadv,scal,scal_bound)
	call wrt_custom					!this is for goro
	call scalar_mass(iumass,nscal,scal)
	call scalar_out(iubase,nscal,scal)

	end
	
!************************************************************

	subroutine wrt_finalize

! finalizes WRT module

	use boxes_wrt

	implicit none

	character*80 text

        call make_max_file('fort.201-fort.',200,9,text)

	write(6,*) 'water residence time (WRT) module has been run'
	if( bgoro ) then
	  write(6,*) 'special setting for Goro has been used'
	end if
	write(6,*) 'results for wrt module are in '//trim(text)
	write(6,*) 'total mass for wrt module is in ',iumass

	end

!************************************************************

	subroutine wrt_init

! initializes scalar for WRT computations

	use boxes
	use boxes_wrt

	implicit none

	integer iu

	if( nscal /= 3 ) then
	  stop 'error stop wrt_init: nscal /= 3'
	end if

	if( nbound < 0 ) then
	  write(6,*) 'nbound = ',nbound
	  write(6,*) 'number of open boundary conditions not initialized'
	  stop 'error stop wrt_init: nbound'
	end if

	allocate(scal(nbox,nscal))
	allocate(scal_bound(-nbound:0,nscal))
	scal = 0.
	scal_bound = 0.

	allocate(wrtacc(nbox,nscal))
	allocate(wrtlog(nbox,nscal))

	call wrt_reset

	ifrqwrt = nint( spd * dtwrt / dt )	!frequency of reset

	iu = iumass
	call scalar_mass_header(iu,nscal,'total mass for WRT')

	iu = iubase
	call scalar_out_header(iu+1,'concentration for WRT')
	call scalar_out_header(iu+2,'concentration for confinement')
	call scalar_out_header(iu+3,'water age [days]')

	iu = iubase + 3
	call scalar_out_header(iu+1,'WRT [days] (integration)')
	call scalar_out_header(iu+2,'confinement [days] (integration)')
	call scalar_out_header(iu+3,'water age [days] (snap shot)')

	iu = iubase + 6
	call scalar_out_header(iu+1,'WRT [days] (log)')
	call scalar_out_header(iu+2,'confinement [days] (log)')
	call scalar_out_header(iu+3,'water age [days] (snap shot)')

	write(6,*) 'finished initializing wrt: ',nscal,ifrqwrt

	end
	
!************************************************************

	subroutine wrt_reset

! resets scalar for WRT computations

	use boxes
	use boxes_wrt

	implicit none

	dtime_reset = dtime
	tacu = 0.

	iacc = 0
	wrtacc = 0.
	wrtlog = 0.
	scal(:,1) = 1.
	scal(:,2) = 0.

	end

!************************************************************

	subroutine wrt_boundary_condition

! sets boundary conditions for WRT - open boundary is 1

	use boxes
	use boxes_wrt

	implicit none

	scal_bound = 0.
	scal_bound(-nbound:0,2) = 1.		!OBC for confinement

	end
	
!************************************************************

	subroutine wrt_react

! does computations for WRT

	use boxes
	use boxes_wrt

	implicit none

	integer iu,ib
	double precision time
	double precision rl(nbox)
	character*80 format

!--------------------------------------------------------
! compute parameters
!--------------------------------------------------------

	iacc = iacc + 1

	time = dtime - dtime_reset
	tacu = tacu + time*time

!--------------------------------------------------------
! compute parameters
!--------------------------------------------------------

	rl = -log(scal(:,1))
	where( scal(:,1) <= 0. ) rl = 0.
	where( scal(:,1) >= 1. ) rl = 1.
	wrtlog(:,1) = wrtlog(:,1) + rl*time
	where( scal(:,2) < 1. ) rl = -log(1.-scal(:,2))
	where( scal(:,2) <= 0. ) rl = 0.
	where( scal(:,2) >= 1. ) rl = 1.
	wrtlog(:,2) = wrtlog(:,2) + rl*time

	where( scal(:,1) < 0. ) scal(:,1) = 0.
	where( scal(:,1) > 1. ) scal(:,1) = 1.
	where( scal(:,2) < 0. ) scal(:,2) = 0.
	where( scal(:,2) > 1. ) scal(:,2) = 1.

	wrtacc(:,1) = wrtacc(:,1) + dt*scal(:,1)		!WRT
	wrtacc(:,2) = wrtacc(:,2) + dt*(1.-scal(:,2))		!confinement
	scal(:,3) = scal(:,3) + dt / spd			!age [days]
	!scal(:,3) = scal(:,3) + dt				!age [s]

!--------------------------------------------------------
! output results
!--------------------------------------------------------

	if( iacc == ifrqwrt ) then
	  write(format,'(a,i4,a)') '(a20,',nbox,'f8.2)'
	  scal(:,1) = wrtacc(:,1) / spd
	  scal(:,2) = wrtacc(:,2) / spd
	  iu = iubase + 3
	  write(iu+1,format) stime,scal(:,1)
	  write(iu+2,format) stime,scal(:,2)
	  write(iu+3,format) stime,scal(:,3)

	  scal(:,1) = (tacu/wrtlog(:,1)) / spd
	  scal(:,2) = (tacu/wrtlog(:,2)) / spd
	  where( wrtlog(:,1) == 0. ) scal(:,1) = 0.
	  where( wrtlog(:,2) == 0. ) scal(:,2) = 0.
	  where( scal(:,1) > wrtmax ) scal(:,1) = 0.
	  where( scal(:,2) > wrtmax ) scal(:,2) = 0.
	  iu = iubase + 6
	  write(iu+1,format) stime,scal(:,1)
	  write(iu+2,format) stime,scal(:,2)
	  write(iu+3,format) stime,scal(:,3)
	  call wrt_reset	!resets all vars, also iacc
	end if

!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	end

!************************************************************

	subroutine wrt_custom

	use boxes
	use boxes_wrt

	implicit none

	if( .not. bgoro ) return

	scal(6:nbox,1) = 0.
	scal(6:nbox,2) = 1.
	scal(6:nbox,3) = 0.

	end

!************************************************************

