
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
! box model gestas - utility routines
!
!------------------------------------------------------------

!************************************************************

	subroutine read_init_box(file,nbox,nscal,val)

! reads file with initial conditions for scalar variables
!
! file format is:
!
! 1    scal_1 scal_2 scal_3 ... state_nscal
! 2    scal_1 scal_2 scal_3 ... state_nscal
! ...
! nbox scal_1 scal_2 scal_3 ... state_nscal
!
! or
!
! 0    scal_1 scal_2 scal_3 ... state_nscal
!
! if the first value is 0 it will set the initial values for all boxes the same
! in this case no other lines are allowed in the file
!
! comments lines start with #

	implicit none

	character*(*) file
	integer nbox
	integer nscal
	double precision val(nbox,nscal)

	double precision vals(0:nscal+1)
	character*1024 line
	integer ibox,n,ib

	integer iscand

	open(1,file=file,status='old',form='formatted',err=97)

	!write(6,*) 'nscal = ',nscal

	ibox = 0
	do
	  read(1,'(a)') line
	  if( line(1:1) == '#' ) cycle	!comment
	  if( line == ' ' ) cycle	!empty line
	  ibox = ibox + 1
	  n = iscand(line,vals,nscal+2)	!reads line - max nscal+2 values
	  n = n - 1
	  if( n /= nscal ) goto 99
	  ib = nint(vals(0))
	  if( ib == 0 ) then		!set all boxes to the same values
	    do ib=1,nbox
	      val(ib,:) = vals(1:nscal)
	    end do
	    exit			!we are done - do not read more lines
	  end if
	  if( ib /= ibox ) goto 98
	  val(ibox,:) = vals(1:nscal)
	  if( ibox == nbox ) exit
	end do

	close(1)

	return
   97	continue
	write(6,*) 'error reading file ',trim(file)
	write(6,*) 'no such file...'
	stop 'error stop read_init_box: no such file'
   98	continue
	write(6,*) 'error reading file ',trim(file)
	write(6,*) trim(line)
	write(6,*) ibox,ib
	write(6,*) 'first column is not box number'
	stop 'error stop read_init_box: error in box number'
   99	continue
	write(6,*) 'error reading file ',trim(file)
	write(6,*) trim(line)
	write(6,*) n,nscal
	write(6,*) 'wrong number of values read'
	stop 'error stop read_init_box: nscal /= n'
	end

!************************************************************

	subroutine read_init_var(file,nscal,val)

! reads file with initial conditions for scalar variables
!
! scalar variable is for all boxes
!
! file format is:
!
! 1	scal_1
! 2	scal_2
! ...
! nscal	scal_nscal
!
! comments lines start with #

	implicit none

	character*(*) file
	integer nscal
	double precision val(nscal)

	double precision vals(3)
	character*1024 line
	integer iscal,n,is,ios

	integer iscand

	open(1,file=file,status='old',form='formatted',err=97)

	iscal = 0
	do
	  read(1,'(a)') line
	  if( line(1:1) == '#' ) cycle	!comment
	  if( line == ' ' ) cycle	!empty line
	  iscal = iscal + 1
	  n = iscand(line,vals,3)	!reads line - max 3 values
	  if( n < 0 ) n = -n - 1
	  if( n /= 2 ) goto 99
	  is = nint(vals(1))
	  if( is /= iscal ) goto 98
	  val(iscal) = vals(2)
	  if( iscal == nscal ) exit
	end do

	close(1)

	return
   97	continue
	write(6,*) 'error reading file ',trim(file)
	write(6,*) 'no such file...'
	stop 'error stop read_init_var: no such file'
   98	continue
	write(6,*) 'error reading file ',trim(file)
	write(6,*) trim(line)
	write(6,*) iscal,is
	write(6,*) 'first column is not scalar number'
	stop 'error stop read_init_var: error in scalar number'
   99	continue
	write(6,*) 'error reading file ',trim(file)
	write(6,*) trim(line)
	write(6,*) n,2
	write(6,*) 'wrong number of values read'
	stop 'error stop read_init_var: n /= 2'
	end

!************************************************************
!************************************************************
!************************************************************

	subroutine write_state_box(time,nbox,nstate,state)

! writes state variables 12-14 of some boxes to unit 99

	implicit none

	double precision time
	integer nbox,nstate
	double precision state(nbox,nstate)

	integer iu,ib,is

	iu = 99

	write(iu,*) '======================================'
	write(iu,*) 'box ',time
	write(iu,*) '======================================'

	do is=1,nstate
	  if( is < 11 .or. is > 15 ) cycle
	  write(iu,*) 'state = ',is
	  do ib=1,nbox,10
	    write(iu,*) is,ib,state(ib,is)
	  end do
	end do

	end

!************************************************************

	subroutine write_state_bound(time,nbound,nstate,state)

! writes state variables 12-14 of some boundaries to unit 99

	implicit none

	double precision time
	integer nbound,nstate
	double precision state(-nbound:0,nstate)

	integer iu,ib,is

	iu = 99

	write(iu,*) '======================================'
	write(iu,*) 'boundary ',time
	write(iu,*) '======================================'

	do is=1,nstate
	  if( is < 11 .or. is > 15 ) cycle
	  write(iu,*) 'state = ',is
	  do ib=0,nbound,10
	    write(iu,*) is,-ib,state(-ib,is)
	  end do
	end do

	end

!************************************************************

	subroutine write_state(time,what,nn,nstate,is,state)

! writes state variable is of all boxes to unit 99

	implicit none

	double precision time
	character*(*) what
	integer nn,nstate,is
	double precision state(nn,nstate)

	integer iu,ib

	iu = 99

	write(iu,*) '======================================'
	write(iu,*) 'state ',trim(what),' ',time,is
	write(iu,*) '======================================'

	  do ib=1,nn
	    write(iu,*) ib,state(ib,is)
	  end do

	end

!************************************************************
!************************************************************
!************************************************************

	subroutine make_file_name(basename,extension,n,zeros,filename)

! makes file name with given number n

	implicit none

	character*(*) basename
	character*(*) extension
	integer n			!number to be used in file name
	integer zeros			!total length of number string
	character*(*) filename

	integer i
	character*10 string

	filename = adjustl(basename)
	write(string,'(i10)') n

	do i=10-zeros+1,10
	  if( string(i:i) /= ' ' ) exit
	  string(i:i) = '0'
	end do

	filename = trim(adjustl(basename))
	filename = trim(filename) // trim(adjustl(string))
	filename = trim(filename) // trim(adjustl(extension))

	end

!************************************************************

	subroutine get_first_line(filename,first)

! get first line of file

	implicit none

	character*(*) filename,first

	open(1,file=filename,status='old',form='formatted')
	read(1,'(a)') first
	close(1)

	end

!************************************************************

	subroutine make_max_file(filename,iubase,n,outfile)

! make filename of max file

	implicit none

	character*(*) filename
	integer iubase			!base unit
	integer n			!max number of files
	character*(*) outfile

	integer iu
        character*80 text

        iu = iubase + n
        write(text,'(i3)') iu
        outfile = trim(filename) // adjustl(text)

	end

!************************************************************
!************************************************************
!************************************************************

	subroutine read_input_file_names(nscal,nbound
     &			,param_file,init_file,bound_files)

! reads file input_file_names.txt and looks for input files

	implicit none

	integer nscal,nbound
	character*80 param_file
	character*80 init_file
	character*80 bound_files(nbound)

	logical bex,bdebug
	integer iline,i,in,ios,ibound
	double precision d(1)
	character*80 file
	character*80 string

	integer iscand

	file = 'input_files.txt'

	bdebug = .true.
	iline = 0
	param_file = ' '
	init_file = ' '
	bound_files = ' '

	inquire(file=file,exist=bex)
	if( .not. bex ) goto 99

	open(1,file=file,status='old',form='formatted',iostat=ios)
	if( ios /= 0 ) goto 98

	do
	  iline = iline + 1
	  read(1,'(a)',iostat=ios) string
	  if( ios > 0 ) goto 98
	  if( ios < 0 ) exit
	  if( string == ' ' ) cycle		!empty line
	  if( string(1:1) == '#' ) cycle	!comment
	  i = index(string,':')
	  if( i == 0 ) goto 97
	  if( string(1:i-1) == 'param' ) then
	    param_file = adjustl(string(i+1:))
	  else if( string(1:i-1) == 'init' ) then
	    init_file = adjustl(string(i+1:))
	  else if( string(1:5) == 'bound' ) then
	    in = iscand(string(6:i-1),d,1)
	    if( in /= 1 ) goto 96
	    ibound = nint(d(1))
	    if( ibound < 1 .or. ibound > nbound ) goto 95
	    bound_files(ibound) = adjustl(string(i+1:))
	  end if
	end do

	close(1)

	if( .not. bdebug ) return

	write(6,*) 'param: ',trim(param_file)
	write(6,*) 'init: ',trim(init_file)
	write(6,*) 'bound: '
	do i=1,nbound
	  write(6,*) i,trim(bound_files(i))
	end do

	file = param_file
	inquire(file=file,exist=bex)
	if( .not. bex ) goto 99
	file = init_file
	inquire(file=file,exist=bex)
	if( .not. bex ) goto 99
	do i=1,nbound
	  file = bound_files(i)
	  inquire(file=file,exist=bex)
	  if( .not. bex ) goto 99
	end do
	
	return
   95	continue
	write(6,*) 'boundary number out of range: ',trim(file)
	write(6,*) 'boundary number: ',ibound
	write(6,*) 'possible min/max boundary number: ',1,nbound
	write(6,*) 'line number: ',iline
	write(6,*) 'line: ',trim(string)
	stop 'error stop read_input_file_names: boundary number'
   96	continue
	write(6,*) 'error reading number of boundary: ',trim(file)
	write(6,*) 'line number: ',iline
	write(6,*) 'line: ',trim(string)
	stop 'error stop read_input_file_names: parsing boundary'
   97	continue
	write(6,*) 'error parsing line: ',trim(file)
	write(6,*) 'line number: ',iline
	write(6,*) 'line: ',trim(string)
	stop 'error stop read_input_file_names: parsing error'
   98	continue
	write(6,*) 'error reading file: ',trim(file)
	write(6,*) 'line number: ',iline
	stop 'error stop read_input_file_names: read error'
   99	continue
	write(6,*) 'file is not existing: ',trim(file)
	stop 'error stop read_input_file_names: no such file'
	end

!************************************************************
!************************************************************
!************************************************************
! test routines
!************************************************************
!************************************************************
!************************************************************

	subroutine test_make_file_name

	implicit none

	character*80 basename,filename,extension
	integer n,zeros

	basename = 'init_bound_'
	extension = '.txt'

	do zeros=0,3
	  write(6,*) 'zeros = ',zeros
	  do n=1,120,15
	    call make_file_name(basename,extension,n,zeros,filename)
	    write(6,*) trim(filename)
	  end do
	end do

	end

!************************************************************
!	program test_main
!	call test_make_file_name
!	end
!************************************************************

