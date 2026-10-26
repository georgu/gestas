
!------------------------------------------------------------
! fake aquabc module to be able to compile
!------------------------------------------------------------

!************************************************************

	module mod_aquabc_fake

	integer, save :: nbox_fake = 0
	integer, save :: nstate_fake = 0
	integer, save :: nforce_fake = 0

	end module mod_aquabc_fake

!************************************************************

        subroutine aquabc_init(nbox,nstate,n_driving_functions,surf_box)

	use mod_aquabc_fake

	implicit none

	integer nbox
	integer nstate
	integer n_driving_functions
	integer surf_box(nbox)

	nbox_fake = nbox
	nstate_fake = nstate
	nforce_fake = nstate

	write(6,*) 'aquabc files have not been linked... cannot run'
	stop 'error stop aquabc_init: no aquabc linked'

	end

!************************************************************

	subroutine  aquabc_run(time,time_step,scal,ph,forcing,sed_flux)

	use mod_aquabc_fake

	implicit none

	double precision time
	double precision time_step
	double precision scal(nbox_fake,nstate_fake)
	double precision ph(nbox_fake)
	double precision forcing(nbox_fake,nforce_fake)
	double precision sed_flux(nbox_fake,nstate_fake)

	end

!************************************************************

	subroutine aquabc_get_flags(nflags,flags)

	implicit none

	integer nflags
	integer flags(nflags)

	end

!************************************************************

	subroutine aquabc_set_flags(nflags,flags)

	implicit none

	integer nflags
	integer flags(nflags)

	end

!************************************************************

	subroutine aquabc_read_constants(file)

	implicit none

	character*(*) file

	end

!************************************************************

	subroutine aquabc_write_constants(file)

	implicit none

	character*(*) file

	end

!************************************************************

