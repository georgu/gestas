
!------------------------------------------------------------
!
! box model gestas - weutro variables
!
! administers runs with various weutro variables
! provides common routines advect, mass, output
!
! to add a new module do the following:
!	new logical variable (e.g., bwrt)
!	new calls in weutro_init, weutro_run, weutro_finalize
!
!------------------------------------------------------------

!============================================================
        module boxes_weutro
!============================================================

        implicit none

        integer, parameter :: nscal = 9        !total number of scalars
        integer, parameter :: nforce = 10       !total number of forcings

	integer, parameter :: nstate = nscal

        integer, save :: iubase = 300                   !unit for output

        integer, save :: iadv(nscal) = 1                !advection of scalars

        logical, save :: bdebug = .false.

!       -----------------------------------------------------

        integer, parameter :: time_steps_per_day = 24
        integer, parameter :: seconds_per_day = 86400
        double precision, save :: time = 1
        double precision, save :: time_step = 1./time_steps_per_day

        double precision, save, allocatable :: scal(:,:)
        double precision, save, allocatable :: scal_old(:,:)
        double precision, save, allocatable :: forcing(:,:)
        integer, save, allocatable :: surf_box(:)

        double precision, save, allocatable :: volumes(:,:)
        double precision, save, allocatable :: flow_rates(:,:)
        double precision, save, allocatable :: scal_bound(:,:)

        integer, save, allocatable :: ids(:)

!============================================================
        end module boxes_weutro
!============================================================

	subroutine weutro_init

! initializes weutro

	use boxes
	use boxes_weutro

	implicit none

	character*80, save :: param_file = 'weutro.nml'

        if( nbound < 0 ) then
          write(6,*) 'nbound = ',nbound
          write(6,*) 'nbound not initialized'
          stop 'error stop weutro_init: nbound'
        end if

        write(6,*) 'initializing weutro: ',nscal

        !allocate(vals(nscal))

        allocate(scal(nbox,nscal))
        allocate(scal_old(nbox,nscal))
        !allocate(forcing(nbox,nforce))
        !allocate(ph(nbox))
        !allocate(sed_flux(nbox,nscal))
        allocate(surf_box(nbox))

        !allocate(volumes(nbox,nscal))
        !allocate(flow_rates(nbox,nscal))
        allocate(scal_bound(-nosect:0,nscal))

        !allocate(description(nscal))

        allocate(ids(nosect))

        scal = 0.
        scal_old = 0.
        scal_bound = 0.

	call handle_param_init(param_file)
	call eutroini
        call weutro_initial_conditions

	end
	
!************************************************************

	subroutine weutro_run

! runs weutro

	use boxes
	use boxes_weutro

	implicit none

        call weutro_boundary_conditions
        call scalar_advect(nscal,iadv,scal,scal_bound)
        call weutro_react
        call scalar_mass(iubase,nscal,scal)
        call scalar_out(iubase,nscal,scal)

	end
	
!************************************************************

	subroutine weutro_finalize

! finalizes weutro

	use boxes
	use boxes_weutro

	implicit none

	write(6,*) 'weutro module has finished to run'

	end
	
!************************************************************
!************************************************************
!************************************************************

        subroutine weutro_initial_conditions

! sets boundary conditions for weutro

        use boxes
        use boxes_weutro

        implicit none

	real, save :: bcs(nscal) = 
     &		(/0.02,0.05,0.03,0.06,1.0,7.8,0.03,0.005,0.012/)
        integer i

	if( nbound /= nosect ) then
	  write(6,*) 'applying OBCs: ',nbound,nosect
	  stop 'error stop weutro_boundary_conditions: nbound /= nosect'
	end if

        scal = 0.
        do i=1,nbox
          scal(i,:) = bcs(:)         !initial conditions for state variables
        end do

        end

!************************************************************

        subroutine weutro_boundary_conditions

! sets boundary conditions for weutro

        use boxes
        use boxes_weutro

        implicit none

	real, save :: bcs(nscal) = 
     &         (/0.2560,0.5516,0.0181,0.0065
     &		,1.5000,6.0168,9.0168,0.0118,0.0100/)
        integer i

	if( nbound /= nosect ) then
	  write(6,*) 'applying OBCs: ',nbound,nosect
	  stop 'error stop weutro_boundary_conditions: nbound /= nosect'
	end if

        scal_bound = 0.
        do i=1,nbound
          scal_bound(-i,:) = bcs(:)         !OBC for state variables
        end do

        end

!************************************************************

        subroutine weutro_react

! does reactions on weutro

        use boxes
        use boxes_weutro

	implicit none

	integer ib,id
	real td,dtd
	real vol,dep,vel,uws,stp,sal,qss
	real c(nstate)
	real, save :: loads(nstate) = 0

!# box     srad    tair    rhum  wspeed      cc    rain    evap     ice

	td = dtime / 86400.	!td is time in days
	dtd = dt / 86400.	!dtd is time step in days

        do ib=1,nbox

	  id = ib*1000 + 1		!prepare to recover node and level

	  c(:) = scal(ib,:)

	  vol = vnew(ib)		!volume of box
	  dep = depth(ib) + znew(ib)	!depth of box
	  vel = cspeed(ib)		!current speed of box
	  uws = meteo(8,ib)		!wind speed of box
	  stp = temp(ib)		!water temperature of box
	  sal = salt(ib)		!water salinity of box
	  qss = meteo(1,ib)		!solar radiation of box

	  call eutro0d(id,td,dtd,vol,dep,vel,uws,stp,sal,qss,c,loads)

	end do

        end

!************************************************************
!************************************************************
!************************************************************
! utility routines
!************************************************************
!************************************************************
!************************************************************

	subroutine get_act_timeline(aline)

	use boxes
	use mod_date

	implicit none

        double precision atime
        character*20 aline

        atime = atime0 + dtime
        call date_format(atime,aline)

	end

	function ipint(k)

	implicit none

	integer ipint
	integer k

	ipint = k

	end

!************************************************************

