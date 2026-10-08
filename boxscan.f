
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

