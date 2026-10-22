
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
! box model gestas - concentration computations
!
! general tracer computation routines
!
! concentration inputs are from every boundary separately
! concentration nbound+1 uses input from all boundaries
!
!------------------------------------------------------------

!============================================================
        module boxes_conz
!============================================================

        implicit none

        integer, save :: nscal = 0			!total number of scalars

	integer, save :: iubase = 400			!unit for output

        integer, save, allocatable :: iadv(:)           !advection of scalars

!	-----------------------------------------------------

	double precision, save, allocatable :: scal(:,:) !scalar variables
	double precision, save, allocatable :: scal_bound(:,:)	!scalar BC

!============================================================
        end module boxes_conz
!============================================================

!************************************************************

	subroutine conz_run

! runs concentration routines

	use boxes
	use boxes_conz

	implicit none

	call conz_boundary_condition
	call scalar_advect(nscal,iadv,scal,scal_bound)
	call conz_react
	call scalar_mass(iubase,nscal,scal)
	call scalar_out(iubase,nscal,scal)

	end
	
!************************************************************

        subroutine conz_finalize

! finalizes concentration module

        write(6,*) 'concentration module has been run'

        end

!************************************************************

	subroutine conz_init

! initializes scalar for concentration computations

	use boxes
	use boxes_conz

	implicit none

	integer iu,is
	character*80 what

	if( nbound < 0 ) then
	  write(6,*) 'nbound = ',nbound
	  write(6,*) 'number of open boundary conditions not initialized'
	  stop 'error stop conz_init: nbound'
	end if

	nscal = nbound + 1

	allocate(scal(nbox,nscal))
	allocate(scal_bound(-nbound:0,nscal))
	allocate(iadv(nscal))
	scal = 0.
	scal_bound = 0.
	iadv = 1

        iu = iubase
        call scalar_mass_header(iu,nscal,'total mass for concentration')

	do is=1,nscal
          iu = iubase+is
	  write(what,'(a,i3)') '# concentration ',is
          call scalar_out_header(iu,what)
	end do

	write(6,*) 'finished initializing conz: ',nscal

	end
	
!************************************************************

	subroutine conz_boundary_condition

! sets boundary conditions for concentration

	use boxes
	use boxes_conz

	implicit none

	integer i

	scal_bound = 0.
	do i=1,nbound
	  scal_bound(-i,i) = 1.			!single conz per boundary
	end do
	scal_bound(-nbound:-1,nscal) = 1.	!all bounds with this OBC

	end

!************************************************************

	subroutine conz_react

! does reactions on concentrations

	end

!************************************************************

