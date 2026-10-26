
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
! box model gestas - reading routines
!
!------------------------------------------------------------

!************************************************************

	subroutine get_open_boundary_conditions

! we have to read first part of 2d file to find out about nbound

	use boxes

	implicit none

	integer ierr

	integer, save :: iu = 0
	character*80 :: file

	integer nb,ns,nvars,nob,nvers
	integer i,is,from

	iu = 11
	file='boxes_2d.txt'
	call check_file_exists(file)
	open(iu,file=file,status='old',form='formatted')
	call check_file_header(iu,2,nvers,ierr)
	if( ierr < 0 ) then
          stop 'error stop get_open_boundary_conditions: ' //           &
     &                  'reading boxes_2d.txt'
	end if
	call check_version('get_open_boundary_conditions',nvers)

	call skip_lines(iu,3)		!skip time record
	read(iu,*) nb,nvars
	call skip_lines(iu,nb+2)	!skip boxes
	read(iu,*) ns
	if( nvers <= 5 ) then
	  call skip_lines(iu,ns+2)	!skip sections
	  read(iu,*) nob
	else
	  call skip_lines(iu,1)		!skip sections
	  nob = 0
	  do i=1,ns
	    read(iu,*) is,from
	    if( i /= is ) then
	      write(6,*) 'i,is: ',i,is
	      stop 'error stop get_open_boundary_conditions: i/=is'
	    end if
	    if( from < 0 ) nob = nob + 1	!boundary
	  end do
	end if
	close(iu)

	nbound = nob
	nosect = nob
	if( .not. bsilent ) write(6,*) 'open boundary conditions: ',nbound

	end

!************************************************************

	subroutine skip_lines(iu,n)

	implicit none

	integer iu,n

	integer i

	do i=1,n
	  read(iu,*)
	end do

	end

!************************************************************

	subroutine read_2d_record(ierr)

	use boxes

	implicit none

	integer ierr

	integer, save :: iu = 0
	character*80 :: file
	character*80 :: string

	integer ib,is
	integer i,ii,i1,i2
	integer nb,ns,no,nob,nvars
	integer, save :: nvers = 0
	double precision a

	nob = 0		!number of open boundary conditions

!	--------------------------------------------------
!	first time called
!	--------------------------------------------------

	if( iu == 0 ) then
	  iu = 2
	  file='boxes_2d.txt'
	  call check_file_exists(file)
	  open(iu,file=file,status='old',form='formatted')
	  call check_file_header(iu,2,nvers,ierr)
	  if( ierr < 0 ) return
	  call check_version('read_2d_record',nvers)
	end if

!	--------------------------------------------------
!	read time record
!	--------------------------------------------------

	read(iu,'(a)',end=1)
	read(iu,*) dtime,stime

	call check_time

!	--------------------------------------------------
!	read boxes
!	--------------------------------------------------

	read(iu,*)
	read(iu,*) nb,nvars
	if( nb /= nbox ) stop 'error stop read_2d_record: nb/=nbox'
	if( nvars /= 7 ) stop 'error stop read_2d_record: nvars/=7'

	read(iu,*)
	do ib=1,nbox
	  read(iu,*) i,temp(ib),salt(ib),cspeed(ib),znew(ib),a,a,vnew(ib)
	  if( ib /= i ) stop 'error stop read_2d_record: ib/=i'
	end do

!	--------------------------------------------------
!	read sections
!	--------------------------------------------------

	read(iu,*)
	read(iu,*) ns
	if( ns > nsect ) stop 'error stop read_2d_record: ns/=nsect'

	flux = 0.
	read(iu,*)
	do is=1,ns
	  if( nvers <= 5 ) then
	    i = is
	    read(iu,*) i1,i2,(flux(ii,is),ii=1,3)
	  else
	    read(iu,*) i,i1,i2,(flux(ii,is),ii=1,3)
	  end if
	  if( is /= i ) stop 'error stop read_2d_record: is/=i'
	  nob = min(nob,i1,i2)
	end do

!	--------------------------------------------------
!	next only for old version
!	--------------------------------------------------

	if( nvers <= 5 ) then

	read(iu,*)
	read(iu,*) no

	if( .not. allocated(oflux) ) then
	  allocate(oflux(3,nosect),iosect(2,nosect))
	  !write(6,*) 'open boundary sections: ',nosect
	end if

	if( no /= nosect ) then
	  write(6,*) nosect,no
	  write(6,*) 'number of OB sections has changed'
	  stop 'error stop read_2d_record: nosect /= no'
	end if

	read(iu,*)
	do is=1,no
	  read(iu,*) (iosect(ii,is),ii=1,2),(oflux(ii,is),ii=1,3)
	  nob = min(nob,iosect(1,is),iosect(2,is))
	end do

	end if

!	--------------------------------------------------
!	check number of open boundaries
!	--------------------------------------------------

	nob = -nob
	if( nbound == 0 ) then
	  nbound = nob
	  !write(6,*) 'open boundary conditions: ',nbound
	end if

	if( nbound /= nob ) then
	  write(6,*) nbound,nob
	  write(6,*) 'number of OBC have changed'
	  stop 'error stop read_2d_record: nbound /= nob'
	end if

	ierr = 0

!	--------------------------------------------------
!	end of routine
!	--------------------------------------------------

	return
    1	continue
	ierr = -1
	return
	end

!************************************************************

	subroutine read_meteo_record(ierr)

	use boxes

	implicit none

	integer ierr

	integer, save :: iu = 0
	character*80 :: file

	integer ib,is
	integer i,ii,i1,i2
	integer nb,ns,no,nob,nvars,nvers
	double precision a

	nob = 0		!number of open boundary conditions

	if( iu == 0 ) then
	  iu = 1
	  file='boxes_meteo.txt'
	  call check_file_exists(file)
	  open(iu,file=file,status='old',form='formatted')
	  call check_file_header(iu,1,nvers,ierr)
	  if( ierr < 0 ) return
	  call check_version('read_meteo_record',nvers)
	end if

	read(iu,*,end=1)
	read(iu,*) dtime,stime

	call check_time

	read(iu,*)
	read(iu,*) nb,nvars
	if( nb /= nbox ) stop 'error stop read_meteo_record: nb/=nbox'
	if( nvars /= 8 ) stop 'error stop read_meteo_record: nvars/=8'

	read(iu,*)
	do ib=1,nbox
	  read(iu,*) i,meteo(1:8,ib)
	  if( ib /= i ) stop 'error stop read_meteo_record: ib/=i'
	end do

	ierr = 0
	return
    1	continue
	ierr = -1
	return
	end

!************************************************************

	subroutine check_file_header(iu,ft,nvers,ierr)

	implicit none

	integer iu
	integer ft
	integer nvers
	integer ierr

	integer idbox,ftype
	integer, parameter :: idbox_expect = 473226
	integer, parameter :: nvers_expect = 7
	character*80 string,file
	character*2 aux

	ierr = 0

	read(iu,*,end=97)
	read(iu,'(a)') string
	string(1:) = string(3:)
	!write(6,*) '|',trim(string),'|'
	!write(6,*) '|','123456789012345678901234567890','|'
	!write(6,*) 'gguu: ',iu,'  ',trim(string)
	read(string,'(3i10)') idbox,ftype,nvers

	if( idbox /= idbox_expect ) goto 99
	if( ft /= ftype ) goto 99
	if( nvers == 3 ) then
	  !write(6,*) 'old version... using compatibility mode: ',nvers
	else if( nvers /= nvers_expect ) then
	  goto 99
	end if

	return
   97	continue
	ierr = -1
	return
   99	continue
	write(6,*) 'error reading header of file'
	call get_file_name(ft,file)
	write(6,*) 'file name: ',trim(file)
	write(6,*) '                    id       ftype       nvers'
	write(6,*) 'read:     ',idbox,ftype,nvers
	write(6,*) 'expected: ',idbox_expect,ft,nvers_expect
	stop 'error stop check_file_header: error in header'
	end

!************************************************************

	subroutine check_file_exists(file)

	implicit none

	character*(*) file

	logical bex

	inquire(file=file,exist=bex)

	if( .not. bex ) then
	  write(6,*) 'file not existing: ',trim(file)
	  stop 'error stop check_file_exists: no such file'
	end if

	end

!************************************************************

        subroutine check_time

        use boxes
        use mod_date

        implicit none

	logical, save :: bfirst = .true.
        double precision atime
        double precision ptime
        character*20 aline

        atime = atime0 + dtime
        call date_format(atime,aline)

        ptime = ptime0 + dtime          !print time
        call date_format(ptime,pline)

        dlast = dtime                   !last good time read

        !write(6,*) aline,'  ',stime

        if( aline /= stime ) then
          write(6,*) 'in time record ',irec+1
          write(6,*) 'expected: ',aline
          write(6,*) 'read    : ',stime
          stop 'error stop check_time: time stamp is wrong'
        end if

	if( bfirst ) stime_start = stime
	stime_end = stime

	bfirst = .false.

        end

!************************************************************

	subroutine get_file_name(ftype,file)

	implicit none

	integer ftype
	character*(*) file

	if( ftype == 0 ) then
	  file = 'boxes_init.txt'
	else if( ftype == 1 ) then
	  file = 'boxes_meteo.txt'
	else if( ftype == 2 ) then
	  file = 'boxes_2d.txt'
	else if( ftype == 3 ) then
	  file = 'boxes_3d.txt'
	else if( ftype == 7 ) then
	  file = 'boxes_geom.txt'
	else if( ftype == 8 ) then
	  file = 'boxes_stats.txt'
	else
	  write(6,*) 'ftype = ',ftype
	  stop 'error stop get_file_name: unknown file type'
	end if

	end

!************************************************************

	subroutine read_stats

	use boxes

	implicit none

	integer :: iu = 1
	character*80 :: file

	logical bvers3
	integer ib,is,nvers,ierr
	integer i,ii,is_aux

	file='boxes_stats.txt'
	call check_file_exists(file)
	open(iu,file=file,status='old',form='formatted',err=99)
	call check_file_header(iu,8,nvers,ierr)
        if( ierr < 0 ) stop 'error stop read_stats: ' //                &
     &                  'reading boxes_stats.txt'
	call check_version('read_stats',nvers)

	bvers3 = ( nvers == 3 )

	read(iu,*)
	read(iu,*) date0,time0
	read(iu,*) dt

	read(iu,*)
	read(iu,*)

	read(iu,*)
	read(iu,*) nbox
	read(iu,*)

	allocate(nlayer(nbox),area(nbox),volume0(nbox),depth(nbox))
	allocate(vinit(nbox),zinit(nbox))

	do ib=1,nbox
	  read(iu,*) i,nlayer(ib),area(ib),volume0(ib),depth(ib)
	  if( i /= ib ) then
	    write(6,*) i,ib
	    stop 'error stop read_stats: i /= ib'
	  end if
	end do

	read(iu,*)
	read(iu,*) nsect
	read(iu,*)

	allocate(isect(5,nsect))

	do is=1,nsect
	  if( bvers3 ) then
	    read(iu,*) isect(:,is)
	    is_aux = is
	  else
	    read(iu,*) is_aux,isect(:,is)
	  end if
	  if( is /= is_aux ) goto 98
	  if( isect(1,is) > nbox ) goto 97
	  if( isect(2,is) > nbox ) goto 97
	end do
	
	close(iu)

	if( .not. bsilent ) write(6,*) 'nbox,nsect: ',nbox,nsect
	if( .not. bquiet ) write(6,*) 'finished reading ',trim(file)

	return
   97	continue
	write(6,*) 'error in box number in section'
	write(6,*) isect(:,is)
	stop 'error stop read_stats: ib > nbox'
   98	continue
	write(6,*) 'inconsistency reading section information'
	write(6,*) is,is_aux
	stop 'error stop read_stats: is /= is_aux'
   99	continue
	write(6,*) 'error opening file: ',trim(file)
	stop 'error stop read_stats: opening file'
	end

!************************************************************

	subroutine read_zinit

	use boxes

	implicit none

	integer :: iu = 1
	character*80 :: file

	integer ib,nb,nvers
	integer i,ierr

	file='boxes_init.txt'
	call check_file_exists(file)
	open(iu,file=file,status='old',form='formatted')
	call check_file_header(iu,0,nvers,ierr)
        if( ierr < 0 ) stop 'error stop read_zinit: ' //                &
     &                  'reading boxes_init.txt'
	call check_version('read_zinit',nvers)

	read(iu,*)
	read(iu,*) nb
	read(iu,*)

	if( nb /= nbox ) then
	  write(6,*) nb,nbox
	  stop 'error stop read_zinit: nb /= nbox'
	end if

	do ib=1,nbox
	  read(iu,*) i,zold(ib)
	  if( i /= ib ) then
	    write(6,*) i,ib
	    stop 'error stop read_zinit: i /= ib'
	  end if
	end do

	vold = volume0 + area*zold

	vinit = vold
	zinit = zold

	znew = zold
	vnew = vold

	close(iu)

	if( .not. bquiet ) write(6,*) 'finished reading ',trim(file)

	end

!************************************************************
!************************************************************
!************************************************************

