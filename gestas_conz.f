
!------------------------------------------------------------
!
! box model gestas - concentration computations
!
! general tracer computation routines
! boundary conditions must be initialized in conz_boundary_condition()
!
!------------------------------------------------------------

!============================================================
        module boxes_conz
!============================================================

        implicit none

        integer, parameter :: nscal = 4			!total number of scalars

	integer, save :: iubase = 40			!unit for output

        integer, save :: iadv(nscal) = 1                !advection of scalars

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

	if( nscal /= nbound ) then
	  write(6,*) 'nscal,nbound: ',nscal,nbound
	  stop 'error stop conz_init: nscal /= nbound'
	end if

	if( nbound < 0 ) then
	  write(6,*) 'nbound = ',nbound
	  write(6,*) 'number of open boundary conditions not initialized'
	  stop 'error stop conz_init: nbound'
	end if

	allocate(scal(nbox,nscal))
	allocate(scal_bound(-nbound:0,nscal))
	scal = 0.
	scal_bound = 0.

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
	  scal_bound(-i,i) = 1.		!OBC for concentration
	end do

	end

!************************************************************

	subroutine conz_react

! does reactions on concentrations

	end

!************************************************************

