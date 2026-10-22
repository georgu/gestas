
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
! box model gestas - handle scalar variables
!
! administers runs with various scalar variables
! provides common routines advect, mass, output
!
! to add a new module do the following:
!	new logical variable (e.g., bwrt)
!	new calls in scalar_init, scalar_run, scalar_finalize
!
! more than one module can be run together
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

	logical, save :: bwrt = .true.		!water residence time
	logical, save :: bconz = .true.		!concentrations
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

! computes stabilty index for boxes

	use boxes
	use boxes_scalar

	implicit none

	double precision stabmin

	integer is,ib
	integer if,it
	integer ibstabmin
	real vaux,faux,stab,alpha
	double precision pos,neg
	double precision fluxbox(2,-nosect:nbox)
	double precision qin(-nosect:nbox)

!	1 is into box, 2 is out of box

	fluxbox = 0
	qin = 0

        do is=1,nsect
          if = isect(1,is)      ! box from
          it = isect(2,is)      ! box to
	  pos = flux(2,is)
	  neg = flux(3,is)
          qin(if) = qin(if) + neg
          qin(it) = qin(it) + pos
          fluxbox(1,if) = fluxbox(1,if) + dt*neg
          fluxbox(1,it) = fluxbox(1,it) + dt*pos
          fluxbox(2,if) = fluxbox(2,if) + dt*pos
          fluxbox(2,it) = fluxbox(2,it) + dt*neg
	end do

	stabmin = 1.e+30
	ibstabmin = 0
	do ib = 1,nbox
	  alpha = dt * qin(ib) / vnew(ib)
	  vaux = min(vold(ib),vnew(ib))
	  faux = max(fluxbox(1,ib),fluxbox(2,ib))
	  faux = abs(faux)
	  stab = vaux/faux
	  if( stab < stabmin ) then
	    stabmin = stab
	    ibstabmin = ib
	  end if
	  if( alpha > 1 ) then
	    write(6,*) 'stability error: alpha > 1 ',alpha
	    write(6,*) stime
	    write(6,*) fluxbox(1,ib),vnew(ib)
	    stop 'error stop scalar_stability: alpha > 1'
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

	if( stab < 1. ) then
	  write(6,*) 'stability index: ',stab
	  write(6,*) 'not yet ready for substepping...'
	  stop 'error stop scalar_advect: stab < 1'
	end if

	do is=1,nscal
	  if( iadv(is) == 0 ) cycle			!no advection
	  s(1:nbox) = scal(1:nbox,is)
	  s(-nosect:0) = scal_bound(-nosect:0,is)
	  call scalar_advect_single(is,s)
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
	double precision pos,neg,smax,smin
	double precision qpos,qneg,qdiff
	double precision mold(-nosect:nbox)
	double precision mnew(-nosect:nbox)
	double precision macc(-nosect:nbox)
	double precision sold(-nosect:nbox)

	bdebug = .false.

	sold = s
	mold = 0.
	mnew = 0.
	macc = 0.
	mold(1:nbox) = vold(1:nbox) * s(1:nbox)

        do is=1,nsect
          if = isect(1,is)      ! box from
          it = isect(2,is)      ! box to
	  qpos = dt*flux(2,is)*s(if)
	  qneg = dt*flux(3,is)*s(it)
	  qdiff = qpos - qneg
          macc(if) = macc(if) + qneg - qpos
          macc(it) = macc(it) + qpos - qneg
        end do


	if( nvers_file <= 5 ) then
          do is=1,nosect
            if = iosect(1,is)      ! box from
            it = iosect(2,is)      ! box to
	    pos = dt*oflux(2,is)
	    neg = dt*oflux(3,is)
            macc(if) = macc(if) + neg*s(it) - pos*s(if)
            macc(it) = macc(it) + pos*s(if) - neg*s(it)
          end do
	end if

	!if( bpe ) then	!consider precipitation and evaporation
	!  do ib=1,nbox
	!    mass(ib) = mass(ib) 
	!  end do
	!end if

	mnew(:) = mold(:) + macc(:)
	s(1:nbox) = mnew(1:nbox) / vnew(1:nbox)

	if( blimit ) then
	  smax = maxval(sold)
	  smin = minval(sold)
	  where( s > smax ) s = smax
	  where( s < smin ) s = smin
	end if

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

	formhead = ' '
	write(formhead,'(a,i2,a)') '(a,',nscal,'(a,i2.2))'

	write(iu,'(a,a)') '# ',trim(what)
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

	write(format,'(a,i4,a)') '(a,',nscal,'e16.6)'
	write(iu,format) pline,mass

	end
	
!************************************************************

	subroutine scalar_out_header(iu,what)

! writes header for scalar output

	use boxes

	implicit none

	integer iu
	character*(*) what
	integer nscal

	integer i,is
	character*80 formhead

	formhead = ' '
	write(formhead,'(a,i2,a)') '(a,',nbox,'(a,i2.2))'

	write(iu,'(a,a)') '# ',trim(what)
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

	write(format,'(a,i4,a)') '(a20,',nbox,'f9.4)'

	do is=1,nscal
	  iu = iubase + is
	  write(iu,format) pline,scal(:,is)
	end do

	end
	
!************************************************************

	subroutine scalar_print(text,is,s)

! print scalars

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
!************************************************************
!************************************************************

