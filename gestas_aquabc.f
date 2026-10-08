
!------------------------------------------------------------
!
! box model gestas - aqua_bc model
!
!------------------------------------------------------------

!============================================================
        module boxes_aqua_bc
!============================================================

        implicit none

        integer, parameter :: nscal = 32	!total number of scalars
        integer, parameter :: nforce = 10	!total number of forcings

        integer, save :: iubase = 200			!unit for output

        integer, save :: iadv(nscal) = 1                !advection of scalars

	logical, save :: bdebug = .false.

!	-----------------------------------------------------

	integer, parameter :: time_steps_per_day = 24
	integer, parameter :: seconds_per_day = 86400
	double precision, save :: time = 1
	double precision, save :: time_step = 1./time_steps_per_day

	double precision, save, allocatable :: scal(:,:)
	double precision, save, allocatable :: scal_old(:,:)
	double precision, save, allocatable :: forcing(:,:)
	double precision, save, allocatable :: ph(:)
	double precision, save, allocatable :: sed_flux(:,:)
	integer, save, allocatable :: surf_box(:)

	double precision, save, allocatable :: volumes(:,:)
	double precision, save, allocatable :: flow_rates(:,:)
	double precision, save, allocatable :: scal_bound(:,:)

	integer, save, allocatable :: ids(:)

!============================================================
        end module boxes_aqua_bc
!============================================================

!************************************************************

	subroutine aqua_bc_run

! runs aqua_bc routines

	use boxes
	use boxes_aqua_bc

	implicit none

	call aqua_bc_boundary_condition
	call scalar_advect(nscal,iadv,scal,scal_bound)
	call aqua_bc_react
	call scalar_mass(iubase,nscal,scal)
	call scalar_out(iubase,nscal,scal)
	call scalar_out(iubase+88,1,ph)

	end
	
!************************************************************

        subroutine aqua_bc_finalize

! finalizes AquaBC module

        write(6,*) 'AquaBC module has been run'

        end

!************************************************************

	subroutine aqua_bc_init

! initializes scalars for AquaBC model

	use boxes
	use boxes_aqua_bc
	use mod_lintp

	implicit none

	integer nstate,n_driving_functions,iscal
	integer ibound
	integer iu,is,id,n
	integer zeros
	integer, parameter :: nflags = 12
	integer, save :: flags(nflags)
	double precision, parameter :: zero = 0.
	double precision, allocatable :: vals(:)
	character*80, allocatable :: description(:)
	character*128 what
	character*80 initfile,basename,filename,constfile
	character*80 descrptfile,extension
	character*80 firstline

!-------------------------------------------------
	zeros = 0
	extension = '.txt'
	initfile  = 'data/NEW_INIT.txt'
	basename  = 'data/NEW_FORC_TS_'
	constfile = 'data/const_CL.txt'
	constfile = 'data/WCONST_04.txt'
	descrptfile = 'data/description.txt'	! not strictly necessary
!-------------------------------------------------

	if( nbound < 0 ) then
	  write(6,*) 'nbound = ',nbound
	  write(6,*) 'number of open boundary conditions not initialized'
	  stop 'error stop aqua_bc_init: nbound'
	end if

	write(6,*) 'initializing aquabc: ',nscal

	allocate(vals(nscal))

	allocate(scal(nbox,nscal))
	allocate(scal_old(nbox,nscal))
	allocate(forcing(nbox,nforce))
	allocate(ph(nbox))
	allocate(sed_flux(nbox,nscal))
	allocate(surf_box(nbox))

	allocate(volumes(nbox,nscal))
	allocate(flow_rates(nbox,nscal))
	allocate(scal_bound(-nosect:0,nscal))

	allocate(description(nscal))

	allocate(ids(nosect))

	scal = 0.
	scal_old = 0.
	forcing = 0.
	ph = 0.
	sed_flux = 0.
	surf_box = 1.

	volumes = 0.
	flow_rates = 0.
	scal_bound = 0.

	nstate = nscal
	n_driving_functions = nforce
	time_step = dt / seconds_per_day

	call aquabc_init(nbox,nstate,n_driving_functions,surf_box)

	call aquabc_get_flags(nflags,flags)
	flags(7) = 0				!ADVANCED_REDOX_SIMULATION
	call aquabc_set_flags(nflags,flags)

	!call read_init_var('data/NEW_INIT.txt',nscal,vals)
	!do iscal=1,nscal
	!  scal(:,iscal) = vals(iscal)
	!end do

	!call read_init_var('data/states_bound_sea.txt',nscal,vals)
	!scal_bound(-1,:) = vals(:)
	
	!call read_init_var('data/states_bound_river.txt',nscal,vals)
	!do ibound=2,nosect
	!  scal_bound(-ibound,:) = vals(:)
	!end do
	
	write(6,*) 'opening file for initial conditions: ',trim(initfile)
	call get_first_line(initfile,firstline)
	write(6,*) 'init file: ',trim(firstline)
	call read_init_box(initfile,nbox,nscal,scal)

	write(6,*) 'opening files for boundary conditions: ',nosect
	do ibound=1,nosect
	  n = ibound
	  call make_file_name(basename,extension,n,zeros,filename)
	  call get_first_line(filename,firstline)
	  write(6,'(a,i2,2a)') 'boundary ',n,'  ',trim(firstline)
	  call lintp_init(filename,id)
	  write(6,*) ibound,id,trim(filename)
	  ids(ibound) = id
	  !scal_bound(-ibound,:) = vals(:)
	end do

	call aquabc_read_constants(constfile)
	call aquabc_write_constants('const_out.txt')
	call aqua_bc_read_description(descrptfile,description)

	!call initialize_output('OUTPUT.csv')
	!call write_output(time,nbox,nstate,scal)

        iu = iubase
        call scalar_mass_header(iu,nscal,'total mass for aquabc')

        do is=1,nscal
          iu = iubase+is
          write(what,'(a,i3)') 'aquabc state variable ',is
	  what = trim(what) // trim(description(is))
          call scalar_out_header(iu,what)
        end do
	iu = iubase+88+1
        write(what,'(a,i3)') 'aquabc ph'
	call scalar_out_header(iu,what)

	write(6,*) 'finished initializing aqua_bc: ',nscal,time_step

	end
	
!************************************************************

	subroutine aqua_bc_boundary_condition

! sets boundary conditions for aqua_bc model

	use boxes
	use boxes_aqua_bc
	use mod_lintp

	implicit none

	integer ibound,id,ierr
	double precision vals(nscal)

	!write(6,*) 'setting boundary conditions for time: ',dtime

	do ibound=1,nosect
	  id = ids(ibound)
	  !call lintp_info(id,'in boundary')
	  call lintp_interpolate(id,dtime,vals,ierr)
	  if( ierr /= 0 ) goto 99
	  !write(6,*) ibound,vals
	  scal_bound(-ibound,:) = vals(:)
	end do

	!scal_bound = 0.

	return
   99	continue
	write(6,*) '*** error reading boundary condition'
	write(6,*) id,ierr
	call lintp_info(id)
	stop 'error stop aqua_bc_boundary_condition: boundary'
	end

!************************************************************

	subroutine aqua_bc_react

! does reactions with aqua_bc model

	use boxes
	use boxes_aqua_bc

	implicit none

	integer ib
	integer ibox,iscal
	integer nstate

	call aqua_bc_set_forcing

	time = time + time_step
	scal_old = scal

        call aquabc_run(time,time_step,scal,ph,forcing,sed_flux)

	end

!************************************************************

	subroutine aqua_bc_set_forcing

!	sets forcing variables from boxes file

!    DRIVING_FUNCTIONS(:,1 ) =   14.0D0 ! Water temperature [C]
!    DRIVING_FUNCTIONS(:,2 ) =    2.0D0 ! Salinity [psu]
!    DRIVING_FUNCTIONS(:,3 ) =  200.0D0 ! Available light [W/m**2]
!    DRIVING_FUNCTIONS(:,4 ) =    0.5D0 ! Fraction of day [0-1]
!    DRIVING_FUNCTIONS(:,5 ) =   12.0D0 ! Air temperature [C]
!    DRIVING_FUNCTIONS(:,6 ) =    3.0D0 ! Wind speed [m/s]
!    DRIVING_FUNCTIONS(:,7 ) =    2.0D0 ! Elevation [m]
!    DRIVING_FUNCTIONS(:,8 ) =    3.6D0 ! Depth [m]
!    DRIVING_FUNCTIONS(:,9 ) =    0.1D0 ! Background light extinction [1/m]
!    DRIVING_FUNCTIONS(:,10) =    0.0D0 ! Ice cover [0-1]

	use boxes
	use boxes_aqua_bc

	implicit none

	integer ib

	do ib=1,nbox
	  forcing(ib,1)  = temp(ib)
	  forcing(ib,2)  = salt(ib)
	  forcing(ib,3)  = meteo(1,ib)
	  forcing(ib,4)  = 1.
	  forcing(ib,5)  = meteo(2,ib)
	  forcing(ib,6)  = meteo(4,ib)
	  forcing(ib,7)  = znew(ib)
	  forcing(ib,8)  = depth(ib)
	  forcing(ib,9)  = 0.1
	  forcing(ib,10) = meteo(8,ib)
	end do

	end

!************************************************************

	subroutine aqua_bc_read_description(descrptfile,description)

	use boxes_aqua_bc

	implicit none

	character*80 descrptfile
	character*80 description(nscal)
	
	integer is,ios,j
	character*80 acron,text

	description = ' '

	write(6,*) 'reading description file'
	open(1,file=descrptfile,status='old',form='formatted',iostat=ios)
	if( ios /= 0 ) return

	do is=1,nscal
	  read(1,'(i12,1x,a15,a60)') j,acron,text
	  if( is /= j ) stop 'error stop description: is/=j'
	  if( text(1:3) == '!!!' ) text = '(no description)'
	  if( bdebug ) write(6,*) j,acron(1:15),'  ',trim(text)
	  description(is) = ' ('//trim(acron)//' - '//trim(text)//')'
	end do

	end

!************************************************************
!************************************************************
!************************************************************
! routines for debug
!************************************************************
!************************************************************
!************************************************************

        subroutine aquabc_debug(time,time_step 
     &          ,nbox,nstate,n_driving_functions 
     &          ,STATE_VARIABLES,PH 
     &          ,DRIVING_FUNCTIONS,SEDIMENT_FLUXES,SURFACE_BOXES)

        double precision time,time_step
        integer nbox,nstate,n_driving_functions
        double precision STATE_VARIABLES(nbox,nstate)
        double precision PH(nbox)
        double precision DRIVING_FUNCTIONS(nbox,n_driving_functions)
        double precision SEDIMENT_FLUXES(nbox,nstate)
        integer SURFACE_BOXES(nbox)

        integer ib,is,id,iu
	double precision daux
	real vstate

	return

	iu = 188

        write(iu,*) '=============================='
        write(iu,*) nbox,nstate,n_driving_functions
        write(iu,*) time,time_step
        do ib=1,nbox
          write(iu,*) ib,ph(ib),SURFACE_BOXES(ib)
        end do
        do ib=1,nbox
          do id=1,n_driving_functions
            write(iu,*) ib,id,DRIVING_FUNCTIONS(ib,id)
          end do
          do is=1,nstate
	    vstate = STATE_VARIABLES(ib,is)
	    daux = SEDIMENT_FLUXES(ib,is)
            write(iu,*) ib,is,vstate,daux
          end do
        end do
        write(iu,*) '=============================='

        end

!************************************************************

	subroutine set_fake_forcings(nbox,nforce,forcing)

	implicit none

	integer nbox,nforce
	double precision forcing(nbox,nforce)

	forcing(:,1 ) =   14.0D0 ! Water temperature [C]
	forcing(:,2 ) =    2.0D0 ! Salinity [psu]
	forcing(:,3 ) =  200.0D0 ! Available light [W/m**2]
	forcing(:,4 ) =    0.5D0 ! Fraction of day [0-1]
	forcing(:,5 ) =   12.0D0 ! Air temperature [C]
	forcing(:,6 ) =    3.0D0 ! Wind speed [m/s]
	forcing(:,7 ) =    2.0D0 ! Elevation [m]
	forcing(:,8 ) =    3.6D0 ! Depth [m]
	forcing(:,9 ) =    0.1D0 ! Background light extinction [1/m]
	forcing(:,10) =    0.0D0 ! Ice cover [0-1]

	end

!************************************************************

