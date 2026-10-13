
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

        integer, save :: iadv(nscal) = 1       !advection of scalars (pelagic)
        logical, save :: bdebug = .false.
	real, save :: bfact = 1.0
	!real, save :: bfact = 0.3		!factor for boundary values
	!real, save :: bfact = 0.0		!factor for boundary values
	logical, save :: badvect = .true.	!do advection of scalars?
	!logical, save :: badvect = .false.	!if .false. simulates lake

!       -----------------------------------------------------
!	please do not change anything below this point
!       -----------------------------------------------------

	character*10, save :: what = 'weutro'
	character*13, save :: descrp(nscal) = (/ 
     &		 '   nh3 [mg/L]','   no3 [mg/L]','  opo4 [mg/L]'
     &		,' phyto [mg/L]','  cbod [mg/L]','    do [mg/L]'
     &		,'    on [mg/L]','    op [mg/L]','   zoo [mg/L]'
     &						/)

!       -----------------------------------------------------

	integer, parameter :: nstate = nscal   !number of state variables
        integer, save :: iubase = 300          !base unit for output
        integer, save :: iumass = 930          !base unit for mass output

        integer, parameter :: time_steps_per_day = 24
        integer, parameter :: seconds_per_day = 86400
        double precision, save :: time = 1
        double precision, save :: time_step = 1./time_steps_per_day

        double precision, save, allocatable :: scal(:,:)
        integer, save, allocatable :: surf_box(:)
        double precision, save, allocatable :: scal_bound(:,:)

        integer, save, allocatable :: ids(:)

	character*80, save :: param_file = ' '
	character*80, save :: init_file = ' '
	character*80, save, allocatable :: bound_files(:)

!============================================================
        end module boxes_weutro
!============================================================

	subroutine weutro_init

! initializes weutro

	use boxes
	use boxes_weutro

	implicit none

	integer iu,is

	character*80 text
	character*4 svar

        write(6,*) 'initializing weutro: ',nscal

        if( nbound < 0 ) then
          write(6,*) 'nbound = ',nbound
          write(6,*) 'nbound not initialized'
          stop 'error stop weutro_init: nbound'
	else if( nosect /= nbound ) then
          write(6,*) 'nbound = ',nbound
          write(6,*) 'nosect = ',nosect
          write(6,*) 'nbound and nosect are different'
          stop 'error stop weutro_init: nbound/=nosect'
        end if

        allocate(scal(nbox,nscal))
        allocate(scal_bound(-nosect:0,nscal))
        allocate(surf_box(nbox))

        allocate(ids(nosect))
	allocate(bound_files(nbound))

        scal = 0.
        scal_bound = 0.
	surf_box = 1
	ids = 0
	bound_files(:) = ' '

	call read_input_file_names(nscal,nbound
     &			,param_file,init_file,bound_files)

	call handle_param_init(param_file)
        call weutro_initial_conditions
        call weutro_boundary_conditions

	call eutroini

	text = 'total mass for ' // trim(what)
	call scalar_mass_header(iumass,nscal,text)

	do is=1,nscal
	  iu = iubase + is
	  write(svar,'(i2)') is
	  text = trim(what) // ' state variable ' // svar // descrp(is)
	  call scalar_out_header(iu,text)
	end do

	end
	
!************************************************************

	subroutine weutro_run

! runs weutro

	use boxes
	use boxes_weutro

	implicit none

        call weutro_boundary_conditions
        if( badvect ) call scalar_advect(nscal,iadv,scal,scal_bound)
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

	write(6,*) 'weutro module has successfully finished'

	end
	
!************************************************************
!************************************************************
!************************************************************

        subroutine weutro_initial_conditions

! sets boundary conditions for weutro

        use boxes
        use boxes_weutro

        implicit none

!	real, save :: bcs(nscal) = 
!     &		(/0.02,0.05,0.03,0.06,1.0,7.8,0.03,0.005,0.012/)
!
	if( nbound /= nosect ) then
	  write(6,*) 'applying OBCs: ',nbound,nosect
	  stop 'error stop weutro_boundary_conditions: nbound /= nosect'
	end if

        scal = 0.
	call read_init_box(init_file,nbox,nscal,scal)

        end

!************************************************************

        subroutine weutro_boundary_conditions

! sets boundary conditions for weutro

        use boxes
        use boxes_weutro

	use mod_lintp

        implicit none

!	real, save :: bcs(nscal) = 
!     &         (/0.2560,0.5516,0.0181,0.0065
!     &		,1.5000,6.0168,9.0168,0.0118,0.0100/)

        integer i,id,ierr
	double precision bcs(nscal)
	integer, save :: icall = 0
	character*80 file
	character*20 aline

!------------------------------------------------------------
! first call - initialize
!------------------------------------------------------------

	if( icall == 0 ) then
	  if( nbound /= nosect ) then
	    write(6,*) 'applying OBCs: ',nbound,nosect
	    stop 'error stop weutro_boundary_conditions: nbound /= nosect'
	  end if

          do i=1,nbound
	    file = bound_files(i)
	    call lintp_init(file,id)
	    ids(i) = id
	  end do

	  icall = 1
        end if

!------------------------------------------------------------
! interpolate
!------------------------------------------------------------

        scal_bound = 0.
        do i=1,nbound
          id = ids(i)
          !call lintp_info(id,'in boundary')
          call lintp_interpolate(id,dtime,bcs,ierr)
	  if( ierr /= 0 ) goto 99
          scal_bound(-i,:) = bfact * bcs(:)         !OBC for state variables
	end do

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	return
   99	continue
	file = bound_files(i)
	write(6,*) 'error interpolating boundary file ',trim(file)
	call get_act_timeline(aline)
	write(6,*) 'dtime = ',dtime,'  date = ',aline
	stop 'error stop weutro_boundary_conditions: interpolation error'
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

	td = dtime / 86400.	!td is time in days
	dtd = dt / 86400.	!dtd is time step in days

        do ib=1,nbox

	  id = ib*1000 + 1		!prepare to recover node and level

	  c(:) = scal(ib,:)

	  vol = vnew(ib)		!volume of box
	  dep = depth(ib) + znew(ib)	!depth of box
	  vel = cspeed(ib)		!current speed of box
	  uws = meteo(4,ib)		!wind speed of box
	  stp = temp(ib)		!water temperature of box
	  sal = salt(ib)		!water salinity of box
	  qss = meteo(1,ib)		!solar radiation of box

	  call eutro0d(id,td,dtd,vol,dep,vel,uws,stp,sal,qss,c,loads)

	  scal(ib,:) = c(:)
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

! atime0 and dtime are in module boxes

	use boxes
	use mod_date

	implicit none

        double precision atime
        character*20 aline

        atime = atime0 + dtime
        call date_format(atime,aline)

	end

!************************************************************

	function ipint(k)

	use boxes

	implicit none

	integer ipint
	integer k

	if( k > nbox ) then
	  ipint = 0
	else
	  ipint = k
	end if

	end

!************************************************************

