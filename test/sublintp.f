
!--------------------------------------------------------------
!
! Usage:
!
! call lintp_init(file,iu,nvar,id)
! call lintp_init(file,nvar,id)
! call lintp_init(file,id)
!
! call lintp_reference(id,aline0)
! call lintp_reference(id,atime0)
!
! call lintp_interpolate(id,aline,values,ierr)
! call lintp_interpolate(id,atime,values,ierr)
! call lintp_interpolate(id,aline,value,ierr)
! call lintp_interpolate(id,atime,value,ierr)
!
! file		file name
! iu		initial unit number (default 20)
! nvar		number of columns of values (return)
! id		id of opened file 
! aline		date string
! atime		absolute time
! aline0	reference date string
! atime0	reference absolute time
! values	returned values (double precision :: values(nvar))
! value		returned value  (double precision :: value)
! ierr		error code (0:ok, <0:eof, >0:error)
!
! id is returned in lintp_init, given in lintp_interpolate
! atime can also be relative time if reference date has been set
!
!--------------------------------------------------------------

!==============================================================
        module mod_lintp
!==============================================================

        implicit none

        private

        logical, save :: bdebug = .true.
        integer, parameter :: ndim = 300
	integer, save :: idfill = 0

        type, private :: info
          integer :: iunit = 0
          integer :: nvar = 0
          logical :: bfirst = .true.
          logical :: btimestring = .true.
          logical :: bone_record = .false.
	  double precision :: atime1
	  double precision :: atime2
	  double precision, allocatable :: values1(:)
	  double precision, allocatable :: values2(:)
          double precision :: atime0               !reference time
	  character*80 file
	end type info

	!type(info), save, target, allocatable :: pinfo(:)
	type(info), save :: pinfo(ndim)

        INTERFACE lintp_init
        MODULE PROCEDURE         lintp_init_iu_nvar
     +                          ,lintp_init_nvar
     +                          ,lintp_init_id
        END INTERFACE

        INTERFACE lintp_reference
        MODULE PROCEDURE         lintp_reference_s
     +                          ,lintp_reference_d
        END INTERFACE

        INTERFACE lintp_interpolate
        MODULE PROCEDURE         lintp_interpolate_d
     +                          ,lintp_interpolate_s
     +                          ,lintp_interpolate_d1
     +                          ,lintp_interpolate_s1
        END INTERFACE

	public lintp_init, lintp_interpolate, lintp_reference

!==============================================================
        contains
!==============================================================

	subroutine lintp_init_iu_nvar(file,iu,nvar,id)

	implicit none

	character*(*) file
	integer iu
	integer nvar
	integer id

	integer n
	integer ierr,ios
	double precision atime
	character*80 s(1)
	character*1024 string

	integer iscans

	id = 0
        call lintp_find_unit(iu)
	if( iu <= 0 ) return

	idfill = idfill + 1
	if( idfill > ndim ) stop 'error stop lintp_init: ndim'
	id = idfill

	open(iu,file=file,status='old',form='formatted')

	do			!look for first data line skipping comments
	  read(iu,'(a)',iostat=ios) string
	  if( ios /= 0 ) goto 99
	  string = adjustl(string)
	  if( string == ' ' ) cycle			!empty line
	  if( string(1:1) == '#' ) cycle		!comment
	  exit
	end do

	call  set_values(string,id,nvar)

	pinfo(id)%file = file
	pinfo(id)%iunit = iu
	pinfo(id)%nvar = nvar

	allocate(pinfo(id)%values1(nvar))
	allocate(pinfo(id)%values2(nvar))

	call lintp_info(id,'after set_values')

	call parse_values(string,id,atime,nvar,pinfo(id)%values1,ierr)
	if( ierr /= 0 ) goto 98
	pinfo(id)%atime1 = atime

	read(iu,'(a)',iostat=ios) string
	if( ios > 0 ) goto 99
	if( ios < 0 ) then
	  pinfo(id)%bone_record = .true.
	else
	  call parse_values(string,id,atime,nvar,pinfo(id)%values2,ierr)
	  if( ierr /= 0 ) goto 98
	  pinfo(id)%atime2 = atime
	end if

	call lintp_info(id,'end of init')

	return
   98	continue
	write(6,*) 'file: ',trim(file)
	write(6,*) 'line: ',trim(string)
	stop 'error stop lintp_init: read error'
   99	continue
	write(6,*) 'file: ',trim(file)
	stop 'error stop lintp_init: no data or read error'
	end

!**************************************************************

	subroutine lintp_init_nvar(file,nvar,id)

	implicit none

	character*(*) file
	integer nvar
	integer id

	integer iu

	iu = 20
	call lintp_init_iu_nvar(file,iu,id,nvar)

	end

!**************************************************************

	subroutine lintp_init_id(file,id)

	implicit none

	character*(*) file
	integer id

	integer iu
	integer nvar

	iu = 20
	call lintp_init_iu_nvar(file,iu,id,nvar)

	end

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine lintp_interpolate_s1(id,aline,value,ierr)

	use mod_date

	implicit none

	integer id
	character*(*) aline
	double precision value
	integer ierr

	double precision values(1)

	call lintp_interpolate_s(id,aline,values,ierr)
	value = values(1)

	end

!**************************************************************

	subroutine lintp_interpolate_d1(id,atime,value,ierr)

	use mod_date

	implicit none

	integer id
	double precision atime
	double precision value
	integer ierr

	double precision values(1)

	call lintp_interpolate_d(id,atime,values,ierr)
	value = values(1)

	end

!**************************************************************

	subroutine lintp_interpolate_s(id,aline,values,ierr)

	use mod_date

	implicit none

	integer id
	character*(*) aline
	double precision values(:)
	integer ierr

	double precision atime

	call date_unformat(aline,atime,ierr)
	if( ierr /= 0 ) return

	call lintp_interpolate_d(id,atime,values,ierr)

	end

!**************************************************************

	subroutine lintp_interpolate_d(id,atime,values,ierr)

	implicit none

	integer id
	double precision atime
	double precision values(:)
	integer ierr

	double precision da,ra1,ra2

	if( pinfo(id)%bone_record ) then
	  values = pinfo(id)%values1
	  return
	end if

	if( atime < pinfo(id)%atime1 ) then	!cannot go backward
	  ierr = 3
	  return
	end if

	do
	  if( atime <= pinfo(id)%atime2 ) exit
	  pinfo(id)%atime1 = pinfo(id)%atime2
	  pinfo(id)%values1 = pinfo(id)%values2
	  call read_values(id,pinfo(id)%atime2,pinfo(id)%values2,ierr)
	  if( ierr /= 0 ) return
	end do

	da = pinfo(id)%atime2 - pinfo(id)%atime1
	ra1 = ( atime - pinfo(id)%atime1 ) / da
	ra2 = ( pinfo(id)%atime2 - atime ) / da
	
	values = ra2 * pinfo(id)%values1 + ra1 * pinfo(id)%values2

	end

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine lintp_reference_s(id,aline0)

	use mod_date

	implicit none

	integer id
	character*(*) aline0

	integer ierr
	double precision atime0

	call date_unformat(aline0,atime0,ierr)
	if( ierr /= 0 ) stop 'error stop: cannot parse date'
	call lintp_reference_d(id,atime0)

	end

!**************************************************************

	subroutine lintp_reference_d(id,atime0)

	implicit none

	integer id
	double precision atime0

	pinfo(id)%atime0 = atime0

	end

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine set_values(string,id,nvar)

	use mod_date

	implicit none

	character*(*) string
	integer id
	integer nvar

	integer ioff,iaux,ierr
	integer i,n
	logical btimestring
	character*80 stime
	character*80 s(1)

	integer istot,istod,iscans
	double precision atime

	n = iscans(string,s,0)	!count number of tokens on line
	nvar = n - 1

	ioff = 1
	i = istot(string,stime,ioff)	!get first token
	if( i <= 0 ) goto 99

	btimestring = .false.
	iaux = 1
	i = istod(stime,atime,iaux)
	if( i < 0 ) btimestring = .true.
	pinfo(id)%btimestring = btimestring
	pinfo(id)%bfirst = .false.
	
	!call lintp_info(id,'in set_values')

	if( btimestring ) then
	  call date_unformat(stime,atime,ierr)
	  if( ierr /= 0 ) goto 98
	else
	  iaux = 1
	  i = istod(stime,atime,iaux)
	  if( i <= 0 ) goto 98
	end if

	return
   98	continue
	write(6,*) stime
	stop 'error stop lintp_init: error reading first time stamp'
   99	continue
	write(6,*) string
	stop 'error stop lintp_init: error reading first line'
	end

!**************************************************************

	subroutine read_values(id,atime,values,ierr)

	use mod_date

	implicit none

	integer id
	double precision atime
	integer nvar
	double precision values(:)
	integer ierr

	integer ios
	integer iu
	character*1024 string

	iu = pinfo(id)%iunit
	nvar = pinfo(id)%nvar

	read(iu,'(a)',iostat=ierr) string
	if( ierr > 0 ) goto 99
	if( ierr < 0 ) return

	call parse_values(string,id,atime,nvar,values,ierr)

	return
   99	continue
	write(6,*) trim(string)
	stop 'error stop lintp_init: error reading line'
	end

!**************************************************************

	subroutine parse_values(string,id,atime,nvar,values,ierr)

	use mod_date

	implicit none

	character*(*) string
	integer id
	double precision atime
	integer nvar
	double precision values(nvar)
	integer ierr

	integer ioff,iaux
	integer i
	character*80 stime

	integer istot,istod,iscans

	ierr = 1

	ioff = 1
	i = istot(string,stime,ioff)	!get first token
	if( i <= 0 ) goto 99

	if( pinfo(id)%btimestring ) then
	  call date_unformat(stime,atime,ierr)
	  if( ierr /= 0 ) goto 98
	else
	  iaux = 1
	  i = istod(stime,atime,iaux)
	  if( i <= 0 ) goto 98
	end if
	
	read(string(ioff:),*,iostat=ierr) values

	return
   98	continue
	write(6,*) stime
	stop 'error stop lintp_init: error reading time stamp'
   99	continue
	write(6,*) string
	stop 'error stop lintp_init: error reading line'
	end

!**************************************************************

	subroutine lintp_info(id,text)

	implicit none

	integer id
	character*(*), optional :: text

	character*80 :: string

	string = ' '
	if( present(text) ) string = text

	write(6,*) '-------------------------------------'
	write(6,*) 'lintp_info: ',trim(string)
	write(6,*) '-------------------------------------'

	write(6,*) 'id: ',id
	write(6,*) 'file: ',trim(pinfo(id)%file)
	write(6,*) 'unit: ',pinfo(id)%iunit
	write(6,*) 'nvar: ',pinfo(id)%nvar
	write(6,*) 'atime0: ',pinfo(id)%atime0
	write(6,*) 'atime1: ',pinfo(id)%atime1
	write(6,*) 'atime2: ',pinfo(id)%atime2
	write(6,*) 'bfirst: ',pinfo(id)%bfirst
	write(6,*) 'btimestring: ',pinfo(id)%btimestring
	write(6,*) 'bone_record: ',pinfo(id)%bone_record

	end

!**************************************************************

        subroutine lintp_find_unit(iunit)

        implicit none

        integer iunit

        logical opened
        integer, save :: iumax = 1000
        integer iu

        iu = iunit
        if( iu .le. 0 ) iu = 20                 !set standard unit

        do
                if( iu .eq. 5 ) iu = 7          !safeguard units 5 and 6
                if( iu > iumax ) then
                        write(6,*) 'no unit available to open file'
                        write(6,*) 'maximum unit tried: ',iumax
                        iu = 0
                        exit
                end if
                inquire(iu,opened=opened)
                if( .not. opened ) exit
                iu=iu+1
        end do

        iunit = iu

        end

!==============================================================
        end module mod_lintp
!==============================================================

	subroutine test_main

	use mod_lintp

	integer id,ierr
	double precision val
	character*80 file

	file='river_1.txt'
	call lintp_init(file,id)

	call lintp_reference(id,'2017-01-01')

	call lintp_interpolate(id,'2017-01-11::12:00:00',val,ierr)
	write(6,*) 'val = ',val
	call lintp_interpolate(id,'2017-01-12::12:00:00',val,ierr)
	write(6,*) 'val = ',val
	call lintp_interpolate(id,'2017-01-18::12:00:00',val,ierr)
	write(6,*) 'val = ',val

	end

	program test_lintp
	call test_main
	end

