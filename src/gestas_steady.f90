
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
! box model gestas - steady state computations
!
! this runs scalar simulations to see how long it takes to reach steady state
!
! boundary conditions are set to 100
! initial conditions are set to 0 and 200
!
!------------------------------------------------------------

!============================================================
        module boxes_steady
!============================================================

        implicit none

        integer, parameter :: nscal = 2		!total number of scalars

	integer, save :: iubase = 150		!unit for output
	integer, save :: iumass = 915		!unit for output of mass

	double precision, save :: sinit1 = 0.		!initial for scalar 1
	double precision, save :: sinit2 = 200.		!initial for scalar 2
	double precision, save :: sbound = 100.		!boundary condition

        integer, save :: iadv(nscal) = 1                !advection of scalars

!	-----------------------------------------------------

	double precision, save, allocatable :: scal(:,:) !scalar variables
	double precision, save, allocatable :: scal_bound(:,:)	!scalar BC

!============================================================
        end module boxes_steady
!============================================================

!************************************************************

	subroutine steady_run

! runs concentration routines

	use boxes
	use boxes_steady

	implicit none

	call steady_boundary_condition
	call scalar_advect(nscal,iadv,scal,scal_bound)
	call steady_react
	call scalar_mass(iumass,nscal,scal)	!writes to iubase
	call scalar_out(iubase,nscal,scal)	!writes to iubase+iscal

	end
	
!************************************************************

        subroutine steady_finalize

	use boxes
	use boxes_steady

! finalizes steady module

        character*80 text

	if( bquiet ) return

        call make_max_file('fort.151-fort.',150,3,text)

        write(6,*) 'steady state module has been run'
        write(6,*) '  results for steady module are in '//trim(text)
        write(6,*) '  total mass for wrt module is in ',iumass


        end

!************************************************************

	subroutine steady_init

! initializes scalar for concentration computations

	use boxes
	use boxes_steady

	implicit none

	integer iu

	if( nbound < 0 ) then
	  write(6,*) 'nbound = ',nbound
	  write(6,*) 'number of open boundary conditions not initialized'
	  stop 'error stop steady_init: nbound'
	end if

	allocate(scal(nbox,nscal))
	allocate(scal_bound(-nbound:0,nscal))
	scal = 0.
	scal_bound = 0.

	scal(:,1) = sinit1
	scal(:,2) = sinit2

        iu = iumass
        call scalar_mass_header(iu,nscal,'total mass for steady')

        iu = iubase
        call scalar_out_header(iu+1,'steady scalar low init')
        call scalar_out_header(iu+2,'steady scalar high init')
        call scalar_out_header(iu+3,'steady scalar difference')

	if( .not. bquiet ) write(6,*) 'finished initializing steady: ',nscal

	end
	
!************************************************************

	subroutine steady_boundary_condition

! sets boundary conditions for concentration

	use boxes
	use boxes_steady

	implicit none

	integer i

	scal_bound = sbound

	end

!************************************************************

	subroutine steady_react

! does reactions on concentrations

	use boxes
	use boxes_steady

	implicit none

	integer iu
	double precision sdiff(nbox)
	character*80 format

	sdiff = ( scal(:,2) - scal(:,1) ) / sbound
	
	iu = iubase + 3
        write(format,'(a,i4,a)') '(a20,',nbox,'f8.2)'
	write(iu,format) stime,sdiff(:)

	end

!************************************************************

