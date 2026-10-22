
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
! box model gestas - bfm ecological model
!
!------------------------------------------------------------

!============================================================
        module boxes_bfm
!============================================================

        implicit none

        integer, parameter :: nscal = 49	!total number of scalars

        integer, save :: iubase = 100		!unit for output

        integer, save :: iadv(nscal) = 1        !advection of scalars

!	-----------------------------------------------------

	double precision, save, allocatable :: scal(:,:)	!scalar vars
	double precision, save, allocatable :: scal_bound(:,:)	!scalar BC

!============================================================
        end module boxes_bfm
!============================================================

!************************************************************

	subroutine bfm_run

! runs bfm routines

	use boxes
	use boxes_bfm

	implicit none

	call bfm_boundary_condition
	call scalar_advect(nscal,iadv,scal,scal_bound)
	call bfm_react
	call scalar_mass(iubase,nscal,scal)
	call scalar_out(iubase,nscal,scal)

	end
	
!************************************************************

        subroutine bfm_finalize

! finalizes BFM module

        write(6,*) 'BFM module has been run'

        end

!************************************************************

	subroutine bfm_init

! initializes scalars for BFM model

	use boxes
	use boxes_bfm

	implicit none

	if( nbound < 0 ) then
	  write(6,*) 'nbound = ',nbound
	  write(6,*) 'number of open boundary conditions not initialized'
	  stop 'error stop bfm_init: nbound'
	end if

	allocate(scal(nbox,nscal))
	allocate(scal_bound(-nbound:0,nscal))
	scal = 0.
	scal_bound = 0.

	! here also adjust iadv for benthic variables
	! iadv(5) = 0

	! internal initialization of BFM model

	write(6,*) 'finished initializing bfm: ',nscal

	end
	
!************************************************************

	subroutine bfm_boundary_condition

! sets boundary conditions for bfm model

	use boxes
	use boxes_bfm

	implicit none

	scal_bound = 0.

	end

!************************************************************

	subroutine bfm_react

	use boxes
	use boxes_bfm

	implicit none

	integer ib
	double precision s,t,cs,sr,dn,do

! does reactions with bfm model

	do ib=1,nbox
	  s = salt(ib)			!salinity [psu]
	  t = temp(ib)			!temperature [C]
	  cs = cspeed(ib)		!current speed [m/s]
	  sr = meteo(1,ib)		!solar radiation [W/m**2]
	  !call BFM reactor
	end do

	end

!************************************************************

	subroutine bfm_get_forcing(ib,ts,ss,sw,ws,cs)

	use boxes

	implicit none

	integer ib
	double precision ts,ss,sw,ws,cs

        ts = temp(ib)                 !temperature [C]
        ss = salt(ib)                 !salinity [psu]
        sw = meteo(1,ib)              !solar radiation [W/m**2]
        ws = meteo(4,ib)              !wind speed [m/s]
        cs = cspeed(ib)               !current speed [m/s]

	end

!************************************************************

