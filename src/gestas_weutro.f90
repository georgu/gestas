
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
! box model gestas - weutro ecological model
!
! administers runs with weutro variables
!
!	1	fort.301	   nh3 [mg/L]		ammonia
!	2	fort.302	   no3 [mg/L]		nitrate
!	3	fort.303	  opo4 [mg/L]		orto phosphate
!	4	fort.304	 phyto [mg/L]		phytoplankton
!	5	fort.305	  cbod [mg/L]		carboneus BOD
!	6	fort.306	    do [mg/L]		dissolved oxigen
!	7	fort.307	    on [mg/L]		organic nitrogen
!	8	fort.308	    op [mg/L]		organic phosphorus
!	9	fort.309	   zoo [mg/L]		zooplankton
!
!------------------------------------------------------------

! revision log :
!
! 10.10.2026    ggu     written from scratch
! 27.10.2026    ggu     use bw to decide what to write

!============================================================
        module boxes_weutro
!============================================================

        implicit none

        integer, parameter :: nscal = 9        !total number of scalars

        integer, save :: iadv(nscal) = 1       !advection of scalars (pelagic)
        logical, save :: bdebug = .false.
	real, save :: bfact = 1.0		!factor for boundary values

!       -----------------------------------------------------
!	please do not change anything below this point
!       -----------------------------------------------------

	character*10, save :: what = 'weutro'
        character*13, save :: descrp(nscal) = (/                        &
     &           '   nh3 [mg/L]','   no3 [mg/L]','  opo4 [mg/L]'        &
     &          ,' phyto [mg/L]','  cbod [mg/L]','    do [mg/L]'        &
     &          ,'    on [mg/L]','    op [mg/L]','   zoo [mg/L]'        &
     &                                          /)

!       -----------------------------------------------------

	integer, parameter :: nstate = nscal   !number of state variables
        integer, save :: iubase = 300          !base unit for output
        integer, save :: iumass = 930          !base unit for mass output

        integer, parameter :: spd = 86400

        double precision, save, allocatable :: scal(:,:)
        double precision, save, allocatable :: scal_bound(:,:)
        integer, save, allocatable :: surf_box(:)	!not yet used
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

	logical bw
	integer iu,is
	character*80 text
	character*4 svar

	bw = .not. bquiet
        if( bw ) write(6,*) 'initializing weutro: ',nscal

!       -----------------------------------------------------
!	allocate and initialize arrays
!       -----------------------------------------------------

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

!       -----------------------------------------------------
!	get file names and initialize weutro
!       -----------------------------------------------------

        call read_input_file_names(weutro_input &
     &			,bw                         &
     &			,nbox,nscal,nbound                         &
     &                  ,param_file,init_file,bound_files)

	call handle_weutro_param_init(param_file,bw)
        call weutro_initial_conditions
        call weutro_boundary_conditions

	call eutroini

!       -----------------------------------------------------
!	prepare for output
!       -----------------------------------------------------

	text = 'total mass for ' // trim(what)
	call scalar_mass_header(iumass,nscal,text)

	do is=1,nscal
	  iu = iubase + is
	  write(svar,'(i2)') is
	  text = trim(what) // ' state variable ' // svar // descrp(is)
	  call scalar_out_header(iu,text)
	end do

!       -----------------------------------------------------
!	end of routine
!       -----------------------------------------------------

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
        call scalar_mass(iumass,nscal,scal)
        call scalar_out(iubase,nscal,scal)

	end
	
!************************************************************

	subroutine weutro_finalize

! finalizes weutro

	use boxes
	use boxes_weutro

	implicit none

        character*80 text

	if( bquiet ) return

        call make_max_file('fort.301-fort.',300,nscal,text)

	write(6,*) 'weutro module has been run'
	write(6,*) '  results for weutro module are in '//trim(text)
        write(6,*) '  total mass for weutro module is in ',iumass

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

	td = dtime / spd	!td is time in days
	dtd = dt / spd		!dtd is time step in days

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

