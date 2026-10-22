
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
! box model gestas - handling perpetual year
!
!------------------------------------------------------------

!	module boxes
!	double precision, save :: dtime
!	double precision, save :: atime0 = 0
!	double precision, save :: ptime0 = 0
!	double precision, save :: pyear = 3
!	double precision, save :: dt = 3600
!	end module boxes

!*************************************************************

	subroutine handle_perpetual(dtime_final,ierr)

	use boxes

	implicit none

	integer iu
	double precision :: dtime_final

	double precision, save :: yearsec = 86400*365

	integer ierr,itype
	double precision :: dstart
	character*80 text
	integer, save :: iyear = 0
	integer newyear

	ierr = -1
	if( pyear == 0 ) return		!no pertual year needed

	iyear = iyear + 1
	if( iyear >= pyear ) then
	  !write(6,*) 'iyear = ',iyear,pyear
	  write(6,*) 'finished running perpetual years ',iyear
	  return
	end if
	newyear = iyear + 1

	dstart = dtime_final + dt + pextra - yearsec
	!write(6,*) 'looking for time ',dstart,newyear
	write(6,*) 'simulating perpetual year ',newyear,dstart

	iu = 2
	itype = 2
	call skip_to_start(iu,dstart,itype)

	iu = 1
	itype = 1
	call skip_to_start(iu,dstart,itype)

	ptime0 = ptime0 + yearsec
	ierr = 0

	vold = vinit
	vnew = vinit
	zold = zinit
	znew = zinit

	end
	
!*************************************************************
!*************************************************************
!*************************************************************

	subroutine get_time_record(iu,dtime,ierr)

	implicit none

	integer iu
	double precision dtime
	integer ierr

	dtime = -1
	read(iu,*,iostat=ierr)
	if( ierr /= 0 ) return
	read(iu,*,iostat=ierr) dtime

	end

!*************************************************************

	subroutine skip_header(iu)

	implicit none

	integer iu
	!integer itype

	!integer idbox,ftype,nvers
	!character*80 string

        !read(iu,*,end=97)
        !read(iu,'(a)') string
        !string(1:) = string(3:)
        !read(string,'(3i10)') idbox,ftype,nvers

	!itype = ftype

	call skip_line(iu,2)

	end

!*************************************************************

	subroutine skip_block(iu)

	implicit none

	integer iu

	integer n

	call skip_line(iu,1)		!read extra before nbox or nsect info
	read(iu,*) n			!read nbox or nsect
	call skip_line(iu,n+1)		!read extra and block

	end

!*************************************************************

	subroutine skip_back_time(iu)

	implicit none

	integer iu

	backspace(iu)
	backspace(iu)

	end

!*************************************************************

	subroutine skip_line(iu,n)

	implicit none

	integer iu
	integer n

	integer i

        do i=1,n
          read(iu,*)
        end do

	end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine skip_2d_record(iu,dtime,ierr)

	implicit none

	integer iu
	double precision dtime
	integer ierr

	call get_time_record(iu,dtime,ierr)
	if( ierr /= 0 ) return
	call skip_block(iu)
	call skip_block(iu)

	end

!*************************************************************

	subroutine skip_meteo_record(iu,dtime,ierr)

	implicit none

	integer iu
	double precision dtime
	integer ierr

	call get_time_record(iu,dtime,ierr)
	if( ierr /= 0 ) return
	call skip_block(iu)

	end

!*************************************************************

	subroutine skip_to_start(iu,dstart,itype)

	implicit none

	integer iu
	double precision dstart
	integer itype			!file type: 1: meteo  2: 2d

	integer ierr
	double precision dtime
	
	rewind(iu)
	call skip_header(iu)		!skip header

	do
	  call get_time_record(iu,dtime,ierr)
	  if( ierr /= 0 ) goto 99
	  !write(6,*) 'searching: ',dtime
	  if( dtime > dstart ) goto 98
	  if( dtime == dstart ) exit
	  call skip_block(iu)
	  if( itype == 2 ) call skip_block(iu)
	end do

	call skip_back_time(iu)		!restore time record

	return
   98	continue
	write(6,*) 'cannot find time record ',dstart
	write(6,*) 'first record found is ',dtime
	stop 'error stop skip_to_start: cannot find time record'
   99	continue
	stop 'error stop skip_to_start: error reading 2d_record'
	end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine test_perpetual

	use boxes
	use mod_date

	implicit none

	integer iu,ierr
	!double precision dtime
	!double precision atime0
	double precision dtime_last
	double precision ptime
	character*80 file
	character*20 aline
	character*20, save :: aline0 = '2015-01-01::00:00:00'

	call date_unformat(aline0,atime0,ierr)
	call date_set_reference(atime0)
	if( ierr /= 0 ) then
	  stop 'error stop test_perpetual: cannot parse date'
	end if

	iu = 2
	file = 'boxes_2d.txt'

	open(iu,file=file,status='old',form='formatted')

	call skip_header(iu)

	do

	do
	  call skip_2d_record(iu,dtime,ierr)
	  if( ierr /= 0 ) exit
	  dtime_last = dtime
	  ptime = ptime0 + dtime
	  call date_format(ptime,aline)
	  write(6,*) dtime,ptime,'  ',aline
	end do

	call handle_perpetual(dtime_last,ierr)
	if( ierr /= 0 ) exit

	end do

	close(iu)

	end

!*************************************************************
!	program main_perpetual
!	call test_perpetual
!	end
!*************************************************************

