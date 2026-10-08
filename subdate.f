
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016,2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
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

! date routines (dealing with absolute time)
!
! revision log :
!
! 23.03.2021    ggu     newly created from subdts.f
! 10.04.2021	ggu	adapted to relative time
!
!**************************************************************

!--------------------------------------------------------------------
!
! Usage:
!
! date_to_atime(ys,atime)			!converts date to atime
! date_to_atime(date,time,atime)
! date_to_atime(date,atime)
!
! date_from_atime(atime,ys)			!converts atime to date
! date_from_atime(atime,date,time,msecs)
! date_from_atime(atime,date,time)
! date_from_atime(atime,date)
!
! date_format(ys,aline)				!converts date to string
! date_format(date,time,aline)
! date_format(date,aline)
! date_format(atime,aline)
!
! date_unformat(aline,ys,ierr)			!converts string to date
! date_unformat(aline,date,time,ierr)
! date_unformat(aline,date,ierr)
! date_unformat(aline,atime,ierr)
!
! date_to_jd(date,jd)				!converts date to julian days
! date_from_jd(jd,year,date)			!converts julian days to date
!
! routines dealing with relative time (time relative to a reference time)
!
! date_set_reference(ys)			!set reference date
! date_set_reference(date,time)
! date_set_reference(date)
! date_set_reference(year)
! date_set_reference(atime)
!
! date_from_year(year,date)			!convert year to date
!
! the following routines can be used with rtime instead atime
!
! date_from_atime(rtime,..)			!converts rtime to date
! date_format(rtime,aline)			!converts rtime to string
!
! date		YYYYMMDD
! time		hhmmss
! atime		absolute time - seconds after 1.1.1 (double precision)
! rtime		relative time - seconds after reference time (double precision)
! aline		date string (character*20, YYYY-MM-DD::hh:mm:ss)
! ys		year,month,day,hour,minute,second,msecs,tzone (integer(8))
! msecs		milliseconds
! year		year
! jd		julian day in year
! ierr		error code (0 if no error)
!
! all variables (except atime, rtime, and aline) are integer
!
!--------------------------------------------------------------------

!==============================================================
	module mod_date
!==============================================================

	implicit none

	private

	logical, save :: bdebug = .true.

	integer, save :: errors = 0
	logical, save :: bignore_error = .false.

        double precision, save :: atime0 = 0.		!reference time

        double precision, parameter :: secs_in_day = 86400.
        double precision, parameter :: ttime = 1000 * 365 * secs_in_day
	integer, parameter          :: tdate = 10000

        INTERFACE date_to_atime
        MODULE PROCEDURE         date_to_atime_ys
     +                          ,date_to_atime_dt
     +                          ,date_to_atime_d
        END INTERFACE

        INTERFACE date_from_atime
        MODULE PROCEDURE         date_from_atime_ys
     +                          ,date_from_atime_dtms
     +                          ,date_from_atime_dt
     +                          ,date_from_atime_d
        END INTERFACE

        INTERFACE date_format
        MODULE PROCEDURE         date_format_ys
     +                          ,date_format_dt
     +                          ,date_format_d
     +                          ,date_format_a
        END INTERFACE

        INTERFACE date_unformat
        MODULE PROCEDURE         date_unformat_ys
     +                          ,date_unformat_dt
     +                          ,date_unformat_d
     +                          ,date_unformat_a
        END INTERFACE

        INTERFACE date_set_reference
        MODULE PROCEDURE         date_set_reference_ys
     +                          ,date_set_reference_dt
     +                          ,date_set_reference_d
     +                          ,date_set_reference_a
        END INTERFACE

	public :: 		 date_to_atime
     +				,date_from_atime
     +				,date_format
     +				,date_unformat
     +				,date_set_reference
     +				,date_to_jd
     +				,date_from_jd
     +				,date_from_year
     +				,date_ignore_error
     +				,date_has_error

!==============================================================
	contains
!==============================================================

!**************************************************************
!**************************************************************
!**************************************************************
! conversion routines to atime
!**************************************************************
!**************************************************************
!**************************************************************

	subroutine date_to_atime_ys(ys,atime)

	implicit none

	integer ys(8)
	double precision atime

	integer date,time

	date = 10000*ys(1) + 100*ys(2) + ys(3)
	time = 10000*ys(4) + 100*ys(5) + ys(6)
	call date_to_atime_dt(date,time,atime)
	atime = atime + ys(7)

	end subroutine date_to_atime_ys

!**************************************************************

	subroutine date_to_atime_dt(date,time,atime)

	implicit none

	integer date,time
	double precision atime

	integer adays,secs

	call date_check_date(date)
	call date_check_time(time)

        call date2adays(date,adays)
        call time2secs(time,secs)

	atime = secs_in_day * adays + secs

	end subroutine date_to_atime_dt

!**************************************************************

	subroutine date_to_atime_d(date,atime)

	implicit none

	integer date
	double precision atime

	integer time

	time = 0
	call date_to_atime_dt(date,time,atime)

	end subroutine date_to_atime_d

!**************************************************************
!**************************************************************
!**************************************************************
! internal routines
!**************************************************************
!**************************************************************
!**************************************************************

        subroutine date2adays(date,adays)

! given date returns total days from 1/1/1

        implicit none

        integer date
        integer adays

        integer year,month,day
        integer iy,jd

	call date2ymd(date,year,month,day)

        iy=year-1
        adays=365*iy+iy/4-iy/100+iy/400

        call ymd2jd(year,month,day,jd)

        adays = adays + jd - 1

        end

!**************************************************************

	subroutine date2ymd(date,year,month,day)

! splits date [YYYYMMDD] into year, month, day

	implicit none

	integer date
	integer year,month,day

        month = date / 100
        day = date - month * 100
        year = month / 100
        month = month - year * 100

	end subroutine date2ymd

!**************************************************************

        subroutine time2hms(time,hour,min,sec)

! splits time [HHMMSS] into hour, min, sec

        implicit none

        integer time
        integer hour,min,sec

        min = time / 100
        sec = time - min * 100
        hour = min / 100
        min = min - hour * 100

	end subroutine time2hms

!**************************************************************

        subroutine time2secs(time,secs)

! converts time to seconds

        implicit none

        integer time,secs
        integer hour,min,sec

        call time2hms(time,hour,min,sec)

	secs = hour*3600 + min*60 + sec

	end subroutine time2secs

!**************************************************************

        subroutine adays2ymd(days,year,month,day)

! splits absolute days from 1/1/1 into year, month, day

        implicit none

        integer days
        integer year,month,day

        integer iy,auxdays

!       estimate year (using average days in 400 years)

        iy = (days*400) / (365*400+97) + 2

!       correct year

        auxdays = days + 1
        do while( auxdays .gt. days )
          iy = iy - 1
          auxdays = 365*iy + iy/4 - iy/100 + iy/400
        end do

!       year found, auxdays days left

        year = iy + 1
        auxdays = days - auxdays

!       find month

        month = 0
        do while( auxdays .ge. 0 )
          month = month + 1
          auxdays = auxdays - idmonth(year,month)
        end do

!       find day

        day = auxdays + idmonth(year,month) + 1

        end

!**************************************************************

        subroutine secs2hms(secs,hour,min,sec)

! splits secs into hour, min, sec

        implicit none

        integer secs
        integer hour,min,sec

        min = secs / 60
        sec = secs - min * 60

        hour = min / 60
        min = min - hour * 60

        end

!**************************************************************
!**************************************************************
!**************************************************************
! conversion routines from atime
!**************************************************************
!**************************************************************
!**************************************************************

        subroutine date_from_atime_ys(artime,ys)

! given absolute time converts to ys structure

        implicit none

        double precision artime
        integer ys(8)

        integer days,secs
        double precision atime
        double precision dsecs

	atime = artime
	if( atime < ttime ) then
	  if( atime0 == 0 ) then
	    write(6,*) 'relative time given but no reference date'
	    stop 'error stop date_from_atime_ys: no reference date'
	  end if
	  atime = atime + atime0
	end if

        days = atime/secs_in_day
        dsecs = atime - secs_in_day*days
        secs = dsecs

        call adays2ymd(days,ys(1),ys(2),ys(3))
        call secs2hms(secs,ys(4),ys(5),ys(6))
        ys(7) = 1000.*(dsecs-secs)
        ys(8) = 0

        end

!**************************************************************

        subroutine date_from_atime_dtms(atime,date,time,msecs)

! given absolute time converts to date,time,msecs

        implicit none

        double precision atime
	integer date,time,msecs

        integer ys(8)

        call date_from_atime_ys(atime,ys)
	date = 10000*ys(1) + 100*ys(2) + ys(3)
	time = 10000*ys(4) + 100*ys(5) + ys(6)
	msecs = ys(7)

        end

!**************************************************************

        subroutine date_from_atime_dt(atime,date,time)

! given absolute time converts to date,time

        implicit none

        double precision atime
	integer date,time

        integer ys(8)

        call date_from_atime_ys(atime,ys)
	date = 10000*ys(1) + 100*ys(2) + ys(3)
	time = 10000*ys(4) + 100*ys(5) + ys(6)

	end

!**************************************************************

        subroutine date_from_atime_d(atime,date)

! given absolute time converts to date

        implicit none

        double precision atime
	integer date

        integer ys(8)

        call date_from_atime_ys(atime,ys)
	date = 10000*ys(1) + 100*ys(2) + ys(3)

        end

!**************************************************************
!**************************************************************
!**************************************************************
! julian date routines
!**************************************************************
!**************************************************************
!**************************************************************

        subroutine date_to_jd(date,jd)

	implicit none

	integer date,jd

	integer year,month,day

	call date2ymd(date,year,month,day)
        call ymd2jd(year,month,day,jd)

	end

!**************************************************************

        subroutine ymd2jd(year,month,day,jd)

c converts date to julian days of year

        implicit none

        integer year,month,day
        integer jd

        integer jdmon,idmon

        if( year .eq. 0 ) then
          jd = day
          return
        end if

	call date_check_ymd(year,month,day)
        jd = day + jdmonth(year,month-1)

        end

!**************************************************************

        subroutine date_from_jd(jd,year,date)

! given julian days and year returns month and day

        implicit none

        integer jd,year
	integer date

        integer month,day

        call jd2ymd(jd,year,month,day)
	date = 10000 * year + 100 * month + day

	end

!**************************************************************

        subroutine jd2ymd(jd,year,month,day)

! given julian days and year returns month and day

        implicit none

        integer jd
        integer year,month,day

        integer maxday

	maxday = 365
	if( bises(year) ) maxday = maxday + 1

        if( jd .lt. 1 .or. jd .gt. maxday ) then
	  errors = errors + 1
	  month = 1
	  day = 1
	  if( bignore_error ) return
          write(6,*) '*** year,jd,maxday: ',year,jd,maxday
          stop 'error stop date_from_jd: jd'
        end if

        month = jd / 30
        if( jd .le. jdmonth(year,month) ) month = month - 1
        day = jd - jdmonth(year,month)
        month = month + 1

        end

!**************************************************************

        function jdmonth(year,month)

c returns total days from January to end of month

        implicit none

        integer jdmonth
        integer year,month

        integer, save :: jmm(0:12) = 
     +		(/0,31,59,90,120,151,181,212,243,273,304,334,365/)

        jdmonth=jmm(month)

        if( month .ge. 2 .and. bises(year) ) jdmonth = jdmonth + 1

	end

!**************************************************************

        function idmonth(year,month)

c returns total days from January to end of month

        implicit none

        integer idmonth
        integer year,month

        integer, save :: imm(0:12) = 
     +		(/0,31,28,31,30,31,30,31,31,30,31,30,31/)

        idmonth=imm(month)

        if( month .eq. 2 .and. bises(year) ) idmonth = idmonth + 1

        end

!**************************************************************

        function bises(year)

! true if year is bisestial

        logical bises
        integer year

	bises = .false.
	if( year <= 0 ) return

        if(  ( mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0 )
     +                  .or. mod(year,400) .eq. 0 ) then
                bises = .true.
        end if

        end

!**************************************************************
!**************************************************************
!**************************************************************
! format routines
!**************************************************************
!**************************************************************
!**************************************************************

        subroutine date_format_ys(ys,aline)

! formats date

        implicit none

        integer ys(8)
        character*(*) aline

        character*20 alineaux
        integer i

        alineaux = ' '

	call date_check_ys(ys)

	!write(6,*) ys

        if( ys(1) .gt. 0 ) then		!has year
          write(alineaux,1000) ys(1),ys(2),ys(3),ys(4),ys(5),ys(6)
          do i=1,20
            if( alineaux(i:i) .eq. ' ' ) alineaux(i:i) = '0'
          end do
        else				!start with days
          write(alineaux,2000) ys(3),ys(4),ys(5),ys(6)
          do i=11,20
            if( alineaux(i:i) .eq. ' ' ) alineaux(i:i) = '0'
          end do
        end if

	aline = alineaux

        return
 1000   format(i4,1h-,i2,1h-,i2,2h::,i2,1h:,i2,1h:,i2)
 2000   format(i10,2h::,i2,1h:,i2,1h:,i2)
	end

!**************************************************************

        subroutine date_format_dt(date,time,aline)

! formats date

        implicit none

	integer date,time
        character*(*) aline

        integer ys(8)

	ys = 0
	call date2ymd(date,ys(1),ys(2),ys(3))
	call time2hms(time,ys(4),ys(5),ys(6))
        call date_format_ys(ys,aline)

	end

!**************************************************************

        subroutine date_format_d(date,aline)

! formats date

        implicit none

	integer date
        character*(*) aline

	integer time

	time = 0
        call date_format_dt(date,time,aline)

	end

!**************************************************************

        subroutine date_format_a(artime,aline)

! formats date

        implicit none

	double precision artime
        character*(*) aline

	integer ys(8)
	double precision atime

	atime = artime
	if( atime < ttime ) then
	  if( atime0 == 0 ) then
	    write(6,*) 'relative time given but no reference date'
	    stop 'error stop date_format_a: no reference date'
	  end if
	  atime = atime + atime0
	end if

	call date_from_atime_ys(atime,ys)
        call date_format_ys(ys,aline)

	end

!**************************************************************
!**************************************************************
!**************************************************************

        subroutine date_unformat_ys(aline,ys,ierr)

c unformats date and time -> from string to date

        implicit none

        character*(*) aline
	integer ys(8)
        integer ierr

        integer year,month,day,hour,min,sec
        character*20 ll
        integer i,n
        logical bdebug

        bdebug = .false.

	ys = 0

        ll = aline
        n = len(aline)
        if( n .gt. 20 ) n = 20

        do i=n,1,-1
          if( ll(i:i) .ne. ' ' ) exit
        end do
        n = i

        ierr = 99
        if( n .lt. 10 ) goto 9  !we insist having at least the full date

        ierr = 1
        if( ll(5:5) .ne. '-' .or. ll(8:8) .ne. '-' ) goto 9

        ierr = 3
        read(ll(1:4),'(i4)',err=9) year
        read(ll(6:7),'(i2)',err=9) month
        read(ll(9:10),'(i2)',err=9) day

        hour = 0
        min = 0
        sec = 0

	ys(1) = year
	ys(2) = month
	ys(3) = day

        ierr = 0
        if( n .le. 10 ) return

        ierr = 5
        if( ll(11:12) .ne. '::' ) goto 9
        if( n .ge. 15 .and. ll(15:15) .ne. ':' ) goto 9
        if( n .ge. 18 .and. ll(18:18) .ne. ':' ) goto 9

        ierr = 7
        if( n .eq. 13 .or. n .eq. 16 .or. n .eq. 19 ) goto 9

        ierr = 9
        if( n .ge. 14 ) read(ll(13:14),'(i2)',err=9) hour
        if( n .ge. 17 ) read(ll(16:17),'(i2)',err=9) min
        if( n .ge. 20 ) read(ll(19:20),'(i2)',err=9) sec

	ys(4) = hour
	ys(5) = min
	ys(6) = sec

        ierr = 0

        return
    9   continue
        if( bdebug ) then
          write(6,*) '*** cannot parse date: ',ierr,aline(1:20)
          write(6,*) '    format should be YYYY-MM-DD::hh:mm:ss'
          write(6,*) '    possible also YYYY-MM-DD[::hh[:mm[:ss]]]'
        end if
        return
        end

!**************************************************************

        subroutine date_unformat_dt(aline,date,time,ierr)

c unformats date and time -> from string to date

        implicit none

        character*(*) aline
	integer date,time
        integer ierr

	integer ys(8)

        call date_unformat_ys(aline,ys,ierr)
	date = 10000 * ys(1) + 100 * ys(2) + ys(3)
	time = 10000 * ys(4) + 100 * ys(5) + ys(6)

	end

!**************************************************************

        subroutine date_unformat_d(aline,date,ierr)

c unformats date and time -> from string to date

        implicit none

        character*(*) aline
	integer date
        integer ierr

	integer time

        call date_unformat_dt(aline,date,time,ierr)

	end

!**************************************************************

        subroutine date_unformat_a(aline,atime,ierr)

c unformats date and time -> from string to date

        implicit none

        character*(*) aline
	double precision atime
        integer ierr

	integer ys(8)

	atime = 0.
        call date_unformat_ys(aline,ys,ierr)
	if( ierr /= 0 ) return
	call date_to_atime_ys(ys,atime)

	end

!**************************************************************
!**************************************************************
!**************************************************************
! reference date routines
!**************************************************************
!**************************************************************
!**************************************************************

	subroutine date_set_reference_ys(ys)

	implicit none

	integer ys(8)

	call date_to_atime_ys(ys,atime0)

	end

!**************************************************************

	subroutine date_set_reference_dt(date,time)

	implicit none

	integer date,time

	call date_to_atime_dt(date,time,atime0)

	end

!**************************************************************

	subroutine date_set_reference_d(ydate)

	implicit none

	integer ydate

	integer date

	date = ydate
	if( date < tdate ) date = 10000 * date + 101 	!is year

	call date_to_atime_d(date,atime0)

	end

!**************************************************************

	subroutine date_set_reference_a(atime)

	implicit none

	double precision atime

	atime0 = atime

	end

!**************************************************************
!**************************************************************
!**************************************************************
! conversion routines
!**************************************************************
!**************************************************************
!**************************************************************

	subroutine date_from_year(year,date)

	implicit none

	integer year,date

	date = year
	if( date < tdate ) date = 10000 * date + 101 	!is year

	end

!**************************************************************
!**************************************************************
!**************************************************************
! check routines
!**************************************************************
!**************************************************************
!**************************************************************

	subroutine date_check_ys(ys)

	implicit none

	integer ys(8)

	call date_check_ymd(ys(1),ys(2),ys(3))
	call date_check_hms(ys(4),ys(5),ys(6))

	end

!**************************************************************

	subroutine date_check_date(date)

	implicit none

	integer date

	integer year,month,day

	call date2ymd(date,year,month,day)
	call date_check_ymd(year,month,day)

	end

!**************************************************************

	subroutine date_check_ymd(year,month,day)

	implicit none

	integer year,month,day

	integer error
	integer maxday

	error = 0
        if( month .lt. 1 .or. month .gt. 12 ) error = error + 1
	maxday = idmonth(year,month)
        if( day .lt. 1 .or. day .gt. maxday ) error = error + 1

	if( error > 0 ) then
	  errors = errors + error
	  if( bignore_error ) return
          write(6,*) '*** year,month,day: ',year,month,day
          !call divide_by_zero(jd)
          stop 'error stop date_check_ymd: errors found'
	end if

	end subroutine date_check_ymd

!**************************************************************

	subroutine date_check_time(time)

	implicit none

	integer time

        integer hour,min,sec

	call time2hms(time,hour,min,sec)
	call date_check_hms(hour,min,sec)

	end

!**************************************************************

	subroutine date_check_hms(hour,min,sec)

	implicit none

        integer hour,min,sec

	integer error

	error = 0
	if( hour < 0 .or. hour >= 24 ) error = error + 1
	if( min < 0 .or. min >= 60 ) error = error + 1
	if( sec < 0 .or. sec >= 60 ) error = error + 1

	if( error > 0 ) then
	  errors = errors + error
	  if( bignore_error ) return
          write(6,*) '*** hour,min,sec: ',hour,min,sec
          stop 'error stop date_check_hms: errors found'
	end if

	end subroutine date_check_hms

!**************************************************************
!**************************************************************
!**************************************************************
! error handling routines
!**************************************************************
!**************************************************************
!**************************************************************

	subroutine date_ignore_error(bignore)

	implicit none

	logical bignore

	bignore_error = bignore
	errors = 0

	end

!**************************************************************

	function date_has_error()

	implicit none

	logical date_has_error

	date_has_error = errors > 0

	end

!==============================================================
	end module mod_date
!==============================================================

!**************************************************************
!**************************************************************
!**************************************************************
! test routines
!**************************************************************
!**************************************************************
!**************************************************************

	subroutine atime_test

	use mod_date

	implicit none

	integer, save :: nmax = 1000
	integer i
	double precision atime
	double precision r
	double precision amax
        double precision, parameter :: secs_in_day = 86400.

	amax = 3000 * 365 * secs_in_day

	do i=1,nmax
	  call random_number(r)
	  atime = 1. + r * amax
	  atime = aint(atime)			! no fractions of seconds
	  call atime_test_single(atime)
	  call format_test_single(atime)
	end do

	write(6,*) 'successful testing atime_test: ',nmax

	end

!***********************************

	subroutine atime_test_single(atime)

	use mod_date

	implicit none

	double precision atime

	integer date,time
	double precision atime1

        call date_from_atime(atime,date,time)
	call date_to_atime(date,time,atime1)

	!write(6,*) atime,date,time

	if( atime /= atime1 ) then
	  write(6,*) 'date,time: ',date,time
	  write(6,*) 'atimes: ',atime,atime1
	  stop 'error stop atime_test_single'
	end if

	end

!***********************************

	subroutine format_test_single(atime)

	use mod_date

	implicit none

	double precision atime

	integer date,time
	double precision atime1
	integer ierr
	character*20 aline

	call date_format(atime,aline)
	call date_unformat(aline,atime1,ierr)

	call date_from_atime(atime,date,time)
	!write(6,*) date,time,aline

	if( atime /= atime1 ) then
	  write(6,*) 'aline: ',aline
	  !write(6,*) 'date,time: ',date,time
	  write(6,*) 'atimes: ',atime,atime1
	  stop 'error stop format_test_single'
	end if

	end

!***********************************

	subroutine format_test_special

	implicit none

	call format_test_special_one('2020=12-01::00:00:00',.false.)
	call format_test_special_one('2020-13-01::00:00:00',.false.)
	call format_test_special_one('2020-11-01::25:00:00',.false.)
	call format_test_special_one('2020-11-01::24:00:00',.false.)
	call format_test_special_one('2020-11-21::23:00:00',.true.)
	call format_test_special_one('2020-11-41::23:00:00',.false.)

	write(6,*) 'successful testing date format'

	end

!***********************************

	subroutine format_test_special_one(aline,bvalid)

	use mod_date

	implicit none

	character*20 aline
	logical bvalid

	logical berror
	double precision atime
	integer ierr
	character*20 aline1

	!write(6,*) 'testing aline: ',aline

	call date_ignore_error(.true.)
	call date_unformat(aline,atime,ierr)
	if( ierr /= 0 .and. .not. bvalid ) return
	call date_format(atime,aline1)
	berror = date_has_error()
	if( berror .and. .not. bvalid ) return
	call date_ignore_error(.false.)

	if( .not. bvalid ) then
	  write(6,*) 'aline should be invalid: ',aline
	  stop 'error stop format_test_special_one: format error'
	else if( bvalid .and. berror ) then
	  write(6,*) 'aline should be valid: ',aline
	  stop 'error stop format_test_special_one: format error'
	end if

	end

!**************************************************************

	subroutine julian_day_test

	implicit none

	integer, save :: nmax = 1000
	integer i
	integer jd,year
	real r

	do i=1,nmax
	  call random_number(r)
	  jd = 1 + r*365
	  call random_number(r)
	  year = 1 + r*2000
	  call julian_day_test_day(jd,year)
	end do

	call julian_day_test_day(366,2004)

	write(6,*) 'successful testing julian_day_test: ',nmax

	end

!***********************************

	subroutine julian_day_test_day(jd,year)

	use mod_date

	implicit none

	integer jd,year
	integer date
	integer jd1

        call date_from_jd(jd,year,date)
        call date_to_jd(date,jd1)

	!write(6,*) year,month,day,jd,jd1

	if( jd /= jd1 ) then
	  write(6,*) date,jd,jd1
	  stop 'error stop julian_day_test_day'
	end if

	end

!***********************************

	subroutine date_test
	call julian_day_test
	call atime_test
	call format_test_special
	end

!***********************************
!	program main_date_test
!	call date_test
!	end
!**************************************************************

