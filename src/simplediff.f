
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
! box model gestas - simple routine to compute differences
!
!------------------------------------------------------------

!************************************************************

	program diff

	implicit none

	integer, parameter :: nval = 3
	integer, parameter :: nboxmax = 100
	double precision val1(nval)
	double precision val2(nval)

	logical bwrite
	integer iu1,iu2
	integer nbox1,nbox2,nbox
	integer ib,nc,irec
	character*20 string1,string2,string
	character*80 file1,file2
	double precision dval(nval)
	double precision rval(nval)
	double precision rmax,rmaxval
	double precision, parameter :: eps = 1.e-5

	bwrite = .true.
	bwrite = .false.

        nc = command_argument_count()
	if( nc /= 2 ) stop 'usage: diff file1 file2'

	iu1 = 1
	iu2 = 2

        call get_command_argument(1,file1)
        call get_command_argument(2,file2)

	open(iu1,file=file1,status='old',form='formatted')
	open(iu2,file=file2,status='old',form='formatted')

	write(6,*) 'file1: ',trim(file1)
	write(6,*) 'file2: ',trim(file2)

	irec = 0
	rmaxval = 0.

	do
	  read(iu1,'(a20,i10)',end=99) string1,nbox1
	  read(iu2,'(a20,i10)',end=99) string2,nbox2
	  if( string1 /= string2 ) stop 'string'
	  if( nbox1 /= nbox2 ) stop 'nbox'
	  string = string1
	  nbox = nbox1
	  if( nbox > nboxmax ) stop 'nboxmax'
	  irec = irec + 1

	  if( bwrite ) write(6,*) string
	  do ib=1,nbox
	    read(iu1,*) val1(:)
	    read(iu2,*) val2(:)
	    dval = abs(val1-val2)
	    rval = abs(val1)+abs(val2)
	    rval = 0.5*dval/rval
	    rmax = maxval(rval)
	    if( rmax > rmaxval ) rmaxval = rmax
	    if( rmax > eps ) goto 98
	    if( bwrite ) write(6,*) rval
	  end do

	end do

   99	continue

	write(6,*) 'total number of time records read: ',irec
	write(6,*) 'number of boxes: ',nbox
	write(6,*) 'maximum relative error: ',rmaxval

	return
   98	continue
	write(6,*) '*** max relative error too high'
	write(6,*) 'time: ',string
	write(6,*) 'irec: ',irec
	write(6,*) 'box:  ',ib
	write(6,*) 'eps:  ',eps
	write(6,*) 'rmax: ',rmax
	stop 'error stop: rmax>eps'
	end

!************************************************************

