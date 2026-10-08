
!------------------------------------------------------------
!
! box model gestas - scalar variables
!
! administers runs with various scalar variables
! provides common routines advect, mass, output
!
! to add a new module do the following:
!	new logical variable (e.g., bwrt)
!	new calls in scalar_init, scalar_run, scalar_finalize
!
!------------------------------------------------------------

!============================================================
        module boxes_scalar
!============================================================

        implicit none

        integer, save :: irun = 0		!how often called

	integer, save :: ibstabtot = 0			!stability index box
	double precision, save :: timestabtot = 0	!stability index time
	double precision, save :: stabtot = 1.e+30	!stability index value

	logical, save :: bwrt = .false.		!water residence time
	logical, save :: bconz = .false.	!concentrations
	logical, save :: bbfm = .false.		!bfm model
	logical, save :: bsteady = .false.	!steady state
	logical, save :: baqua_bc = .false.	!aquabc model
	logical, save :: bweutro = .true.	!weutro model

!============================================================
        end module boxes_scalar
!============================================================

	subroutine scalar_init

! initializes scalar

	use boxes
	use boxes_scalar

	implicit none

	write(6,*) 'initializing scalars'

	if( bwrt ) call wrt_init
	if( bconz ) call conz_init
	!if( bbfm ) call bfm_init
	if( bsteady ) call steady_init
	!if( baqua_bc ) call aqua_bc_init
	if( bweutro ) call weutro_init

	write(6,*) 'finished initializing scalars'

	end
	
!************************************************************

	subroutine scalar_run

! runs scalar

	use boxes
	use boxes_scalar

	implicit none

	irun = irun + 1

	if( bwrt ) call wrt_run
	if( bconz ) call conz_run
	!if( bbfm ) call bfm_run
	if( bsteady ) call steady_run
	!if( baqua_bc ) call aqua_bc_run
	if( bweutro ) call weutro_run

	end
	
!************************************************************

	subroutine scalar_finalize

! finalizes scalar

	use boxes
	use boxes_scalar
	use mod_date

	implicit none

	double precision atime
	character*20 aline

	write(6,*) 'the following modules have been run:'

	if( bwrt ) call wrt_finalize
	if( bconz ) call conz_finalize
	!if( bbfm ) call bfm_finalize
	if( bsteady ) call steady_finalize
	!if( baqua_bc ) call aqua_bc_finalize
	if( bweutro ) call weutro_finalize

        atime = atime0 + timestabtot
        call date_format(atime,aline)

	write(6,*) 'the minimum stability index is ',stabtot
	write(6,*) 'box ',ibstabtot,'   time ',trim(aline)

	end
	
!************************************************************
!************************************************************
!************************************************************

	subroutine scalar_stability(stabmin)

! advects all scalar variables

	use boxes
	use boxes_scalar

	implicit none

	double precision stabmin

	integer is,ib
	integer if,it
	integer ibstabmin
	real vaux,faux,stab
	double precision pos,neg
	double precision fluxbox(-nosect:nbox)

	fluxbox = 0

        do is=1,nsect
          if = isect(1,is)      ! box from
          it = isect(2,is)      ! box to
	  pos = dt*flux(2,is)
	  neg = dt*flux(3,is)
          fluxbox(if) = fluxbox(if) + neg - pos
          fluxbox(it) = fluxbox(it) + pos - neg
	end do

	stabmin = 1.e+30
	ibstabmin = 0
	do ib = 1,nbox
	  vaux = min(vold(ib),vnew(ib))
	  faux = abs(fluxbox(ib))
	  stab = vaux/faux
	  if( stab < stabmin ) then
	    stabmin = stab
	    ibstabmin = ib
	  end if
	  !write(6,*) 'stability: ',ib,vaux,faux,stab
	end do

	if( stabmin < stabtot ) then
	  stabtot = stabmin
	  ibstabtot = ibstabmin
	  timestabtot = dtime
	end if

	end

!************************************************************

	subroutine scalar_advect(nscal,iadv,scal,scal_bound)

! advects all scalar variables

	use boxes

	implicit none

	integer nscal
	integer iadv(nscal)
	double precision scal(nbox,nscal)
	double precision scal_bound(-nosect:0,nscal)

	integer is,ib
	integer if,it
	double precision pos,neg,stab
	double precision time_aux
	double precision s(-nosect:nbox)

	call scalar_stability(stab)
	!write(6,*) 'stability index for boxes: ',stab

	if( stab < 1. ) then
	  write(6,*) 'stability index: ',stab
	  write(6,*) 'not yet ready for substepping...'
	  stop 'error stop scalar_advect: stab < 1'
	end if

	do is=1,nscal
	  if( iadv(is) == 0 ) cycle			!no advection
	  s(1:nbox) = scal(1:nbox,is)
	  s(-nosect:0) = scal_bound(-nosect:0,is)
	  !call scalar_print('advect before:',is,s)
	  call scalar_advect_single(is,s)
	  !call scalar_print('advect after:',is,s)
	  scal(1:nbox,is) = s(1:nbox)
	end do

	end

!************************************************************

	subroutine scalar_advect_single(iscal,s)

! advects single scalar variable

	use boxes

	implicit none

	integer iscal
	double precision s(-nosect:nbox)

	logical bdebug
	integer if,it,is,ib
	double precision pos,neg,vaux
	double precision mass(-nosect:nbox)

	bdebug = .false.

	mass = 0.
	mass(1:nbox) = vold(1:nbox) * s(1:nbox)

        do is=1,nsect
          if = isect(1,is)      ! box from
          it = isect(2,is)      ! box to
	  !bdebug = ( if == -1 )
	  pos = dt*flux(2,is)
	  neg = dt*flux(3,is)
	  if( bdebug ) then
	    write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~'
	    write(6,*) iscal
	    write(6,*) if,it
	    write(6,*) 'pos,neg: ',pos,neg
	    vaux = 0
	    if( if > 0 ) vaux = vold(if)
	    write(6,*) 'vol: ',vaux,vold(it)
	    write(6,*) 's before: ',s(if),s(it)
	    write(6,*) 'mass flux: ',neg*s(it) - pos*s(if)
	    write(6,*) 'mass before: ',mass(if),mass(it)
	  end if
          mass(if) = mass(if) + neg*s(it) - pos*s(if)
          mass(it) = mass(it) + pos*s(if) - neg*s(it)
	  if( bdebug ) then
	    vaux = 1.
	    write(6,*) 'mass after: ',mass(if),mass(it)
	    write(6,*) 's after: ',mass(if)/vaux,mass(it)/vnew(it)
	    write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~'
	  end if
        end do

	if( nvers_file <= 5 ) then
          do is=1,nosect
            if = iosect(1,is)      ! box from
            it = iosect(2,is)      ! box to
	    pos = dt*oflux(2,is)
	    neg = dt*oflux(3,is)
            mass(if) = mass(if) + neg*s(it) - pos*s(if)
            mass(it) = mass(it) + pos*s(if) - neg*s(it)
          end do
	end if

	!if( bpe ) then	!consider precipitation and evaporation
	!  do ib=1,nbox
	!    mass(ib) = mass(ib) 
	!  end do
	!end if

	s(1:nbox) = mass(1:nbox) / vnew(1:nbox)

	end
	
!************************************************************
!************************************************************
!************************************************************

	subroutine scalar_mass_header(iu,nscal,what)

! writes header for mass scalar output

	use boxes

	implicit none

	integer iu,nscal
	character*(*) what

	integer i
	character*80 formhead

	write(iu,'(a,a)') '# ',trim(what)

	formhead = ' '
	write(formhead,'(a,i2,a)') '(a,',nscal,'(a,i2.2))'

	write(iu,formhead) '#        dtime'
     +		,('        scalar',i,i=1,nscal)

	end

!************************************************************

	subroutine scalar_mass(iu,nscal,scal)

! computes total mass of scalars for check

	use boxes

	implicit none

	integer iu		!unit to write to
	integer nscal
	double precision scal(nbox,nscal)

	integer is,ib
	double precision mass(nscal)
	double precision m
	character*80 format

	do is=1,nscal
	  m = 0.
	  do ib=1,nbox
	    m = m + vnew(ib) * scal(ib,is)
	  end do
	  mass(is) = m
	end do

	write(format,'(a,i4,a)') '(f14.3,',nscal,'e16.6)'
	write(iu,format) dtime,mass

	end
	
!************************************************************

	subroutine scalar_out_header(iu,what)

! writes header for scalar output

	use boxes

	implicit none

	integer iu
	character*(*) what

	integer i
	character*80 formhead

	write(iu,'(a,a)') '# ',trim(what)

	formhead = ' '
	write(formhead,'(a,i2,a)') '(a,',nbox,'(a,i2.2))'

	write(iu,formhead) '#      date_and_time'
     +		,('    box',i,i=1,nbox)

	end

!************************************************************

	subroutine scalar_out(iubase,nscal,scal)

! outputs scalar

	use boxes

	implicit none

	integer iubase		!base unit to write to
	integer nscal
	double precision scal(nbox,nscal)

	integer is,i,iu
	character*80 format

	write(format,'(a,i4,a)') '(a20,',nbox,'f9.3)'

	do is=1,nscal
	  iu = iubase + is
	  write(iu,format) stime,scal(:,is)
	end do

	end
	
!************************************************************

	subroutine scalar_print(text,is,s)

	use boxes

	implicit none

	character*(*) text
	integer is
	double precision s(-nosect:nbox)

	character*80 format

	format = '(a,10f7.2)'

	write(6,*) trim(text),is,nosect,nbox
	write(6,format) 'boxes:  ',s(1:nbox)
	write(6,format) 'zero:   ',s(0)
	write(6,format) 'bounds: ',s(-nosect:-1)

	end

!************************************************************

