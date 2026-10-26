
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
! box model gestas - utilities for reading and writing box files
!
!------------------------------------------------------------

!---------------------------------------------------------
	module mod_boxutil
!---------------------------------------------------------

	implicit none

	integer, save :: nvers_box = 0
	integer, save :: nvers_max = 6
	integer, save :: id_box = 473226

!---------------------------------------------------------
	contains
!---------------------------------------------------------

!*********************************************************
!*********************************************************
!*********************************************************
! deal with file header
!*********************************************************
!*********************************************************
!*********************************************************

	subroutine read_file_header(iu,ftype,nvers)

	implicit none

	integer iu
	integer ftype
	integer nvers

	integer i
	integer idbox
	character*80 string

        character (len=8),dimension(0:8) :: strings = (/                &
     &          'init    ','meteo   ','2d      ','3d      ',            &
     &          'vertical','        ','        ','geom    ','stats   '  &
     &          /)

	read(iu,'(a)') string
	i = index(string,'idbox')

	if( i <= 0 ) then
	  write(6,*) 'old box format... cannot handle'
	  write(6,*) trim(string),i
	  stop 'error stop read_header: old format'
	end if

	read(iu,'(a2,3i10)') string,idbox,ftype,nvers
	if( idbox /= 473226 ) then
	  write(6,*) 'wrong idbox'
	  write(6,*) 'idbox=',idbox,'   expected=473226'
	  stop 'error stop read_header: idbox'
	end if

	if( nvers < 3 .or. nvers > 6 ) then
	  write(6,*) 'wrong nvers - must be between 3 and 6'
	  write(6,*) 'nvers=',nvers
	  stop 'error stop read_header: nvers'
	end if

	if( ftype < 0 .or. ftype > 8 ) stop 'error stop: ftype'

	string = strings(ftype)	
	if( string == ' ' ) stop 'no string found'
	!write(6,*) 'box file is of type ',ftype,' (',trim(string),')'

	nvers_box = nvers

	end

!*********************************************************
!*********************************************************
!*********************************************************
! deal with time record
!*********************************************************
!*********************************************************
!*********************************************************

	subroutine read_time_record(iu,dtime,aline,ierr)

	implicit none

	integer iu
	double precision dtime
	character*(*) aline
	integer ierr

	character*80 string

	read(iu,'(a)',end=88,err=99) string
	if( trim(string) /= '# time record' ) then
	  write(6,*) 'error reading time string...'
	  write(6,*) trim(string)
	  stop 'error stop read_time_record: wrong string'
	end if

	read(iu,'(d21.8,a)') dtime,string
	aline = adjustl(string)
	ierr = 0

	return
   99	continue
	write(6,*) 'error reading time string...'
	stop 'error stop read_time_record: generic error'
   88	continue
	ierr = -1
	return
	end

!*********************************************************
!*********************************************************
!*********************************************************
! get dimensions for file
!*********************************************************
!*********************************************************
!*********************************************************

        subroutine get_box_dimensions(iu,ftype                          &
     &                          ,nbox,nvar                              &
     &                          ,nsect,nobs                             &
     &                          ,lmax                                   &
     &                          ,ierr)

	implicit none

	integer iu,ftype
	integer nbox,nvar
	integer nsect,nobs
	integer lmax
	integer ierr

	nsect = 0
	nobs = 0
	lmax = 0

	if( ftype == 1 ) then
	  call get_meteo_dimensions(iu,nbox,nvar,ierr)
	else if( ftype == 2 ) then
	  call get_2d_dimensions(iu,nbox,nvar,nsect,nobs,ierr)
	else if( ftype == 3 ) then
	  call get_3d_dimensions(iu,nbox,nvar,lmax,nsect,nobs,ierr)
	else
	end if

	end

!*********************************************************

	subroutine get_meteo_dimensions(iu,nbox,nvar,ierr)

	implicit none

	integer iu
	integer nbox,nvar
	integer ierr

	integer ftype,nvers,naux
	double precision dtime
	character*80 aline

	call read_time_record(iu,dtime,aline,ierr)
	if( ierr /= 0 ) return

	call read_box_header_2d(iu,nbox,nvar)
	call skip_block_2d(iu,nbox,.false.,naux)

	rewind(iu)
	call read_file_header(iu,ftype,nvers)

	end

!*********************************************************

	subroutine get_2d_dimensions(iu,nbox,nvar,nsect,nobs,ierr)

	implicit none

	integer iu
	integer nbox,nvar
	integer nsect,nobs
	integer ierr

	integer ftype,nvers,nobc,naux
	double precision dtime
	character*80 aline

	call read_time_record(iu,dtime,aline,ierr)
	if( ierr /= 0 ) return

	call read_box_header_2d(iu,nbox,nvar)
	call skip_block_2d(iu,nbox,.false.,naux)

	call read_flux_header_2d(iu,nsect)
	call skip_block_2d(iu,nsect,.true.,nobc)
	nobs = nobc

	if( nvers_box <= 5 ) then
	  call read_flux_header_2d(iu,nobs)
	  call skip_block_2d(iu,nobs,.true.,nobc)
	end if

	rewind(iu)
	call read_file_header(iu,ftype,nvers)

	end

!*********************************************************

	subroutine get_3d_dimensions(iu,nbox,nvar,lmax,nsect,nobs,ierr)

	implicit none

	integer iu
	integer nbox,nvar,lmax
	integer nsect,nobs
	integer ierr

	integer ftype,nvers,nobc,naux
	double precision dtime
	character*80 aline

	call read_time_record(iu,dtime,aline,ierr)
	if( ierr /= 0 ) return

	call read_box_header_3d(iu,nbox,nvar,lmax)
	call skip_block_3d(iu,nbox,.false.,naux)

	call read_flux_header_3d(iu,nsect,lmax)
	call skip_block_3d(iu,nsect,.true.,nobc)
	nobs = nobc

	if( nvers_box <= 5 ) then
	  call read_flux_header_3d(iu,nobs,lmax)
	  call skip_block_3d(iu,nobs,.true.,nobc)
	end if

	rewind(iu)
	call read_file_header(iu,ftype,nvers)

	end

!*********************************************************

	subroutine skip_block_2d(iu,nskip,bflux,nobc)

	implicit none

	integer iu
	integer nskip
	logical bflux
	integer nobc

	logical bv6
	integer i,j,from

	bv6 = ( nvers_box >= 6 )

	nobc = 0

	read(iu,*)
	do i=1,nskip
	  if( bflux ) then
	    read(iu,*) j,from
	  else
	    read(iu,*) j
	  end if
	  nobc = min(nobc,from)
	  if( bv6 .and. i/=j ) stop 'error stop skip_block_2d: i/=j'
	end do
	
	end

!*********************************************************

	subroutine skip_block_3d(iu,nskip,bflux,nobc)

	implicit none

	integer iu
	integer nskip
	logical bflux
	integer nobc

	logical bv6
	integer i,j,l,lmax,from,iaux

	bv6 = ( nvers_box >= 6 )

	nobc = 0

	do i=1,nskip
	  read(iu,*)
	  if( bflux ) then
	    if( bv6 ) then
	      read(iu,*) j,from,iaux,lmax
	    else
	      j = 0
	      read(iu,*) from,iaux,lmax
	    end if
	  else
	    read(iu,*) j,lmax
	  end if
	  nobc = min(nobc,from)
	  if( j>0 .and. i/=j ) stop 'error stop skip_block_3d: i/=j'

	  read(iu,*)
	  do l=1,lmax
	    read(iu,*) j
	    if( j>0 .and. l/=j ) stop 'error stop skip_block_3d: l/=j'
	  end do
	end do
	
	end

!*********************************************************
!*********************************************************
!*********************************************************
! deal with 2d box data
!*********************************************************
!*********************************************************
!*********************************************************

	subroutine read_box_2d(iu,nbox,nvar,vals)

	implicit none

	integer iu
	integer nbox,nvar
	double precision, allocatable :: vals(:,:)

	call read_box_header_2d(iu,nbox,nvar)
	call alloc_box_arrays_2d(nbox,nvar,vals)
	call read_box_records_2d(iu,nbox,nvar,vals)

	end

!*********************************************************

	subroutine alloc_box_arrays_2d(nbox,nvar,vals)

	implicit none

	integer nbox,nvar
	double precision, allocatable :: vals(:,:)

	logical balloc
	integer n1,n2

	balloc = .false.
	if( .not. allocated(vals) ) then
	  balloc=.true.
	else
	  n1 = size(vals,dim=1)
	  n2 = size(vals,dim=2)
	  if( n1 /= nvar .or. n2 /= nbox ) then
	    deallocate(vals)
	    balloc=.true.
	  end if
	end if

	if( balloc ) then
	  allocate(vals(nvar,nbox))
	end if

	end

!*********************************************************

	subroutine read_box_header_2d(iu,nbox,nvar)

	implicit none

	integer iu
	integer nbox,nvar

	character*80 string
	character*80 string_error

	string_error = 'read_box_header_2d'

	read(iu,'(a)') string
	call check_string(string,'#      boxes       nvars',string_error)

	read(iu,*) nbox,nvar

	end

!*********************************************************

	subroutine read_box_records_2d(iu,nbox,nvar,vals)

	implicit none

	integer iu
	integer nbox,nvar
	double precision, allocatable :: vals(:,:)

	logical balloc
	integer i,ib
	character*80 string
	character*80 string_error

	string_error = 'read_box_records_2d'

	read(iu,'(a)') string
	call check_string(string,'# box',string_error)

	do ib=1,nbox
	  read(iu,*) i,vals(:,ib)
	  if( i /= ib ) call error_stop(string_error,'i/=ib')
	end do

	end

!*********************************************************
!*********************************************************
!*********************************************************
! deal with 2d flux data
!*********************************************************
!*********************************************************
!*********************************************************

	subroutine read_flux_2d(iu,nsect,isvals,vals)

	implicit none

	integer iu
	integer nsect
	integer, allocatable :: isvals(:,:)
	double precision, allocatable :: vals(:,:)

	call read_flux_header_2d(iu,nsect)
	call alloc_flux_arrays_2d(nsect,isvals,vals)
	call read_flux_records_2d(iu,nsect,isvals,vals)

	end

!*********************************************************

	subroutine alloc_flux_arrays_2d(nsect,isvals,vals)

	implicit none

	integer nsect
	integer, allocatable :: isvals(:,:)
	double precision, allocatable :: vals(:,:)

	logical balloc
	integer n1,n2

	balloc = .false.
	if( .not. allocated(vals) ) then
	  balloc=.true.
	else
	  n2 = size(vals,dim=2)
	  if( nsect > n2 ) then
	    deallocate(vals)
	    deallocate(isvals)
	    balloc=.true.
	  end if
	end if

	if( balloc ) then
	  allocate(vals(3,nsect))
	  allocate(isvals(3,nsect))
	end if

	end

!*********************************************************

	subroutine read_flux_header_2d(iu,nsect)

	implicit none

	integer iu
	integer nsect

	character*80 string

	read(iu,'(a)') string
	if( is_string_equal(string,'#   sections' ) ) then
	else if( is_string_equal(string,'#        OBC' ) ) then
	else
	  write(6,*) 'error reading flux string...'
	  write(6,*) trim(string)
	  stop 'error stop read_flux_header_2d: wrong string'
	end if

	read(iu,*) nsect

	end

!*********************************************************

	subroutine read_flux_records_2d(iu,nsect,isvals,vals)

	implicit none

	integer iu
	integer nsect
	integer, allocatable :: isvals(:,:)
	double precision, allocatable :: vals(:,:)

	logical balloc
	integer i,is,isc
	character*80 string
	character*80 string_error

	string_error = 'read_flux_records_2d'

	read(iu,'(a)') string
	if( nvers_box <= 5 ) then
	  call check_string(string,'# from    to',string_error)
	else
	  call check_string(string,'# section  from    to',string_error)
	end if

	if( nvers_box <= 5 ) then
	  do is=1,nsect
	    read(iu,*) isvals(1:2,is),vals(:,is)
	  end do
	else
	  do is=1,nsect
	    read(iu,*) i,isvals(1:2,is),vals(:,is)
	    if( i /= is ) call error_stop(string_error,'i/=is')
	  end do
	end if

	isvals(3,:) = 1		!2d has 1 layer

	end

!*********************************************************
!*********************************************************
!*********************************************************
! deal with 3d box data
!*********************************************************
!*********************************************************
!*********************************************************

	subroutine read_box_3d(iu,nbox,nvar,lmax,lbmax,vals,vals2d)

	implicit none

	integer iu
	integer nbox,nvar,lmax
	integer, allocatable :: lbmax(:)
	double precision, allocatable :: vals(:,:,:)
	double precision, allocatable :: vals2d(:,:)

	call read_box_header_3d(iu,nbox,nvar,lmax)
	call alloc_box_arrays_3d(nbox,nvar,lmax,lbmax,vals,vals2d)
	call read_box_records_3d(iu,nbox,nvar,lmax,lbmax,vals,vals2d)

	end

!*********************************************************

	subroutine alloc_box_arrays_3d(nbox,nvar,lmax,lbmax,vals,vals2d)

	implicit none

	integer nbox,nvar,lmax
	integer, allocatable :: lbmax(:)
	double precision, allocatable :: vals(:,:,:)
	double precision, allocatable :: vals2d(:,:)

	logical balloc
	integer n1,n2,n3

	balloc = .false.
	if( .not. allocated(vals) ) then
	  balloc=.true.
	else
	  n1 = size(vals,dim=1)
	  n2 = size(vals,dim=2)
	  n3 = size(vals,dim=3)
	  if( n1 /= lmax .or. n2 /= nvar .or. n3 /= nbox ) then
	    deallocate(lbmax)
	    deallocate(vals)
	    deallocate(vals2d)
	    balloc=.true.
	  end if
	end if

	if( balloc ) then
	  allocate(lbmax(nbox))
	  allocate(vals(lmax,nvar,nbox))
	  allocate(vals2d(3,nbox))
	end if

	end

!*********************************************************

	subroutine read_box_header_3d(iu,nbox,nvar,lmax)

	implicit none

	integer iu
	integer nbox,nvar,lmax

	character*80 string
	character*80 string_error

	string_error = 'read_box_header_3d'

	read(iu,'(a)') string
	call check_string(string,'#      boxes  max_layers',string_error)

	read(iu,*) nbox,lmax,nvar

	end

!*********************************************************

        subroutine read_box_records_3d(iu,nbox,nvar,lmax                &
     &                                  ,lbmax,vals,vals2d)

	implicit none

	integer iu
	integer nbox,nvar,lmax
	integer, allocatable :: lbmax(:)
	double precision, allocatable :: vals(:,:,:)
	double precision, allocatable :: vals2d(:,:)

	logical balloc
	integer i,ib,ibb,l,il
	integer llmax
	character*80 string
	character*80 string_error

	string_error = 'read_box_records_3d'

	do ib=1,nbox

	  read(iu,'(a)') string
	  call check_string(string,'#   box layers',string_error)

	  if( nvers_box <= 5 ) then
	    read(iu,*) i,llmax,vals2d(1,ib)
	    vals2d(2,ib) = vals2d(1,ib)
	    vals2d(3,ib) = 0.
	  else
	    read(iu,*) i,llmax,vals2d(:,ib)
	  end if
	  if( i /= ib ) call error_stop(string_error,'i/=ib')
	  lbmax(ib) = llmax

	  read(iu,'(a)') string
	  call check_string(string,'# layer      temp',string_error)

	  do l=1,llmax
	    read(iu,*) i,vals(l,:,ib)
	    if( i /= l ) call error_stop(string_error,'i/=l')
	  end do
	end do

	end

!*********************************************************
!*********************************************************
!*********************************************************
! deal with 3d flux data
!*********************************************************
!*********************************************************
!*********************************************************

	subroutine read_flux_3d(iu,nsect,lmax,isvals,vals)

	implicit none

	integer iu
	integer nsect
	integer lmax
	integer, allocatable :: isvals(:,:)
	double precision, allocatable :: vals(:,:,:)

	call read_flux_header_3d(iu,nsect,lmax)
	call alloc_flux_arrays_3d(nsect,lmax,isvals,vals)
	call read_flux_records_3d(iu,nsect,lmax,isvals,vals)

	end

!*********************************************************

	subroutine alloc_flux_arrays_3d(nsect,lmax,isvals,vals)

	implicit none

	integer nsect,lmax
	integer, allocatable :: isvals(:,:)
	double precision, allocatable :: vals(:,:,:)

	logical balloc
	integer n1,n2,n3

	balloc = .false.
	if( .not. allocated(vals) ) then
	  balloc=.true.
	else
	  n1 = size(vals,dim=1)
	  n2 = size(vals,dim=2)
	  n3 = size(vals,dim=3)
	  if( n1 /= lmax .or. n2 /= 3 .or. n3 /= nsect ) then
	    deallocate(isvals)
	    deallocate(vals)
	    balloc=.true.
	  end if
	end if

	if( balloc ) then
	  allocate(isvals(3,nsect))
	  allocate(vals(lmax,3,nsect))
	end if

	end

!*********************************************************

	subroutine read_flux_header_3d(iu,nsect,lmax)

	implicit none

	integer iu
	integer nsect,lmax

	character*80 string

	read(iu,'(a)') string
	if( is_string_equal(string,'#   sections' ) ) then
	else if( is_string_equal(string,'#        OBC' ) ) then
	else
	  write(6,*) 'error reading flux string...'
	  write(6,*) trim(string)
	  stop 'error stop read_flux_header_3d: wrong string'
	end if

	read(iu,*) nsect,lmax

	end

!*********************************************************

        subroutine read_flux_records_3d(iu,nsect,lmax                   &
     &                                  ,isvals,vals)

	implicit none

	integer iu
	integer nsect,lmax
	integer, allocatable :: isvals(:,:)
	double precision, allocatable :: vals(:,:,:)

	logical balloc
	integer i,is,l
	integer llmax
	character*80 string
	character*80 string_error,string_test

	string_error = 'read_flux_records_3d'

	do is=1,nsect

	  read(iu,'(a)') string

	  if( nvers_box <= 5 ) then
	    string_test = '#       from          to'
	    call check_string(string,trim(string_test),string_error)
	    read(iu,*) isvals(:,is)
	    i = is
	  else
	    string_test = '#    section        from'
	    call check_string(string,trim(string_test),string_error)
	    read(iu,*) i,isvals(:,is)
	  end if
	  if( i /= is ) call error_stop(string_error,'i/=is')

	  read(iu,'(a)') string
	  string_test = '# layer         total'
	  call check_string(string,trim(string_test),string_error)

	  llmax = isvals(3,is)
	  do l=1,llmax
	    read(iu,*) i,vals(l,:,is)
	    if( i == 0 .and. llmax == 1 ) cycle		!old style
	    if( i /= l ) call error_stop(string_error,'i/=l')
	  end do
	end do

	end

!*********************************************************
!*********************************************************
!*********************************************************
! writing 2d files
!*********************************************************
!*********************************************************
!*********************************************************

	subroutine write_file_header(iu,ftype)

	implicit none

	integer iu
	integer ftype

	write(iu,'(a)') '#      idbox file_type     nvers'
	write(iu,'(a,3i10)') '# ',id_box,ftype,nvers_max

	end

!*********************************************************

	subroutine write_time_record(iu,dtime,aline)

	implicit none

	integer iu
	double precision dtime
	character*(*) aline

	write(iu,'(a)') '# time record'
	write(iu,*) dtime,trim(aline)

	end

!*********************************************************

	subroutine write_box_records_2d(iu,nbox,nvar,vals)

	implicit none

	integer iu
	integer nbox,nvar
	double precision, allocatable :: vals(:,:)

	integer ib
	character*80 s1,s2

	write(iu,'(a)') '#      boxes       nvars'
	write(iu,*) nbox,nvar

        s1 = '# box      temp      salt cur_speed'
        s2 = '   act_eta  aver_eta   bstress        volume'
        write(iu,'(a)') trim(s1)//trim(s2)

	do ib=1,nbox
	  write(iu,1000) ib,vals(:,ib)
	end do

	return
 1000   format(i5,6f10.5,e14.6)
	end

!*********************************************************

	subroutine write_flux_records_2d(iu,nsect,isvals,vals)

	implicit none

	integer iu
	integer nsect
	integer, allocatable :: isvals(:,:)
	double precision, allocatable :: vals(:,:)

	integer is
	character*80 s1,s2

	write(iu,'(a)') '#   sections'
	write(iu,*) nsect

        s1 = '# section  from    to'
        s2 = '         total      positive      negative'
        write(iu,'(a)') trim(s1)//trim(s2)

	do is=1,nsect
	  write(iu,2000) is,isvals(1:2,is),vals(:,is)
	end do

	return
 2000   format(i9,2i6,3e14.6)
	end

!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
! utility routines
!*********************************************************
!*********************************************************
!*********************************************************

	function is_string_equal(string,string_compare)

	implicit none

	logical is_string_equal
	character*(*) string,string_compare

	integer isc

	isc = len_trim(string_compare)
	if( string(1:isc) == string_compare(1:isc) ) then
	  is_string_equal = .true.
	else
	  is_string_equal = .false.
	end if

	end

!*********************************************************

	subroutine check_string(string,string_compare,string_error)

	implicit none

	character*(*) string,string_compare,string_error

	integer isc

	isc = len_trim(string_compare)
	if( string(1:isc) == string_compare(1:isc) ) return

	write(6,*) 'error ' // trim(string_error)
	write(6,*) 'string read:     ',string(1:isc)
	write(6,*) 'string expected: ',string_compare(1:isc)
	write(6,*) 'box file version: ',nvers_box
	stop 'error stop ' // trim(string_error) // ': wrong string'

	end

!*********************************************************

	subroutine error_stop(routine,what)

	implicit none

	character*(*) routine,what

	character*80 string

	string = 'error stop ' // trim(routine)

	stop trim(string) // ': ' // trim(what)

	end

!*********************************************************
!*********************************************************
!*********************************************************
! end of module
!*********************************************************
!*********************************************************
!*********************************************************

!---------------------------------------------------------
	end module mod_boxutil
!---------------------------------------------------------

        subroutine read_record_2d(iu,aline,dtime,nbox,nvar,ns,ns_obc    &
     &                  ,vals                                           &
     &                  ,isvals,fluxes                                  &
     &                  ,isvals_obc,fluxes_obc                          &
     &                  ,ierr)

	use mod_boxutil

	implicit none

	integer iu
	character*20 aline
	double precision dtime
	integer nbox,nvar,ns,ns_obc
	double precision, allocatable :: vals(:,:)
	integer, allocatable :: isvals(:,:)
	double precision, allocatable :: fluxes(:,:)
	integer, allocatable :: isvals_obc(:,:)
	double precision, allocatable :: fluxes_obc(:,:)
	integer ierr

	call read_time_record(iu,dtime,aline,ierr)
	if( ierr < 0 ) return
	if( ierr > 0 ) goto 99

	call read_box_2d(iu,nbox,nvar,vals)
	call read_flux_2d(iu,ns,isvals,fluxes)
	call read_flux_2d(iu,ns_obc,isvals_obc,fluxes_obc)

	return
   99	continue
	write(6,*) 'error reading time record: ',ierr
	stop 'error stop read_record_2d: time record'
	end

!*********************************************************

