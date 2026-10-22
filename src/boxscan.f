
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

!-----------------------------------------------------------
!
! box model gestas - scans boxe files and rewrites them
!
!-----------------------------------------------------------

!***********************************************************

	program boxscan

! scans through box file and re-writes it if necessary

	implicit none

	integer ftype

	call read_header(ftype)
	!write(6,*) 'ftype = ',ftype

	do
	  call read_time
	  call read_boxes
	  if( ftype == 2 ) then
	    call read_sections
	    call read_sections
	  end if
	  call skip_header
	end do

	end

!***********************************************************

	subroutine read_header(ftype)

	implicit none

	integer ftype

	integer id,nvers
	character*80 string

	call scan(1)
	read(5,'(a)') string
	write(6,'(a)') trim(string)

	read(string,'(2x,3i10)') id,ftype,nvers
	!write(6,*) id,ftype,nvers

	end

!***********************************************************

	subroutine skip_header

	implicit none

	call skip(2)

	end

!***********************************************************

	subroutine read_time

	implicit none

	call scan(2)

	end

!***********************************************************

	subroutine read_boxes

	implicit none

	integer nbox,nvars

	call scan(1)
	read(5,*) nbox,nvars
	write(6,*) nbox,nvars
	call scan(nbox+1)

	end

!***********************************************************

	subroutine read_sections

	implicit none

	integer nsect

	call scan(1)
	read(5,*) nsect
	write(6,*) nsect
	call scan(nsect+1)

	end

!***********************************************************

	subroutine scan(n)

	implicit none

	integer n

	integer i
	character*80 string

	do i=1,n
	  read(5,'(a)',end=1) string
	  write(6,'(a)') trim(string)
	end do

	return
    1	continue
	stop 'end of file'
	end

!***********************************************************

	subroutine skip(n)

	implicit none

	integer n

	integer i

	do i=1,n
	  read(5,*,end=1)
	end do

	return
    1	continue
	stop 'end of file'
	end

!***********************************************************

