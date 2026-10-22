
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
! box model gestas - program to analyse file boxes.txt
!
!------------------------------------------------------------

!************************************************************

	program boxindex

! analyses boxes.txt file and writes info to terminal

	implicit none

	integer nel,nbox,nmax
	integer is,nodes,from,to
	integer iaux,i,ns
	integer ie,iaa,ib,ietot
	integer, allocatable :: ia(:),ic(:)

	read(5,*) nel,nbox,nmax			! reading header
	write(6,*) 'elements: ',nel
	write(6,*) 'nbox,nmax: ',nbox,nmax

	allocate(ia(nel),ic(nmax))

	read(5,*) ia				! reading area code

	ic = 0
	do ie=1,nel
	  iaa = ia(ie)
	  if( iaa < 1 .or. iaa > nmax ) goto 97
	  ic(iaa) = ic(iaa) + 1
	end do

	write(6,*) 'boxes:'
	write(6,*) '        box    elements'
	ietot = 0
	do ib=1,nmax
	  write(6,*) ib,ic(ib)
	  ietot = ietot + ic(ib)
	end do
	write(6,*) 'total number of boxes:    ',nbox
	write(6,*) 'total number of elements: ',ietot
	if( nel /= ietot ) goto 98

	write(6,*) 'sections:'
	write(6,*) '    section       nodes    from box      to box'
	ns = 0
	do
	  read(5,*) is,nodes,from,to		! reading info on section
	  if( is == 0 ) exit
	  ns = ns + 1
	  if( is /= ns ) goto 99
	  write(6,*) is,nodes,from,to
	  read(5,*) (iaux,i=1,nodes)		! skipping nodes
	end do
	write(6,*) 'total number of sections: ',ns

	stop
   97	continue
	stop 'error stop: internal error (1)'
   98	continue
	stop 'error stop: internal error (2)'
   99	continue
	stop 'error stop: internal error (3)'
	end

!************************************************************

