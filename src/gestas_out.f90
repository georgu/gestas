
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
! box model gestas - output routines
!
!------------------------------------------------------------

!************************************************************

	subroutine hydro_out

	use boxes

	implicit none

	character*20 aline

	integer, save :: iubase = 100
	integer, save :: iuaver = 911
	integer, save :: icall = 0
	integer, save :: nhydro = 4
	double precision :: taver,saver,zaver,vaver
	double precision, allocatable, save :: hydro(:,:)
	character*80 format

	if( icall == 0 ) then
	  allocate(hydro(nbox,nhydro))
	  call scalar_out_header(iubase+1,'water temperature [C]')
	  call scalar_out_header(iubase+2,'salinity [psu]')
	  call scalar_out_header(iubase+3,'current speed [m/s]')
	  call scalar_out_header(iubase+4,'water level [m]')

	  write(iuaver,'(a)')  '# averages of hydro data over all boxes'
          write(iuaver,'(4a)') '#      date_and_time'                   &
     &          // '   temperature'                                     &
     &          // '      salinity'                                     &
     &          // '   water_level'                                     &
     &          // '        volume'

	  icall = 1
	end if

	hydro(:,1) = temp(:)
	hydro(:,2) = salt(:)
	hydro(:,3) = cspeed(:)
	hydro(:,4) = znew(:)

	call scalar_out(iubase,nhydro,hydro)

	call box_average(temp,taver)
	call box_average(salt,saver)
	call box_average(znew,zaver)
	call box_average(vnew,vaver)

	format='(a,3f14.3,e14.6)'
        write(iuaver,format) pline,sngl(taver),sngl(saver)              &
     &                          ,sngl(zaver),sngl(vaver)

	end

!************************************************************

	subroutine meteo_out

	use boxes

	implicit none

	character*20 aline

	integer, save :: iubase = 104
	integer, save :: iuaver = 912
	integer, save :: icall = 0
	integer, save :: nmeteo = 5
	double precision :: taver,saver,haver,waver,caver
	double precision, allocatable, save :: lmeteo(:,:)
	character*80 format

	if( icall == 0 ) then
	  allocate(lmeteo(nbox,nmeteo))
	  call scalar_out_header(iubase+1,'solar radiation [W/m**2]')
	  call scalar_out_header(iubase+2,'air temperature [C]')
	  call scalar_out_header(iubase+3,'relative humidity [%]')
	  call scalar_out_header(iubase+4,'wind speed [m/s]')
	  call scalar_out_header(iubase+5,'cloud cover [0-1]')

	  write(iuaver,'(a)')  '# averages of meteo data over all boxes'
          write(iuaver,'(4a)') '#      date_and_time'                   &
     &          // '       srad'                                        &
     &          // '      atemp'                                        &
     &          // '       rhum'                                        &
     &          // '     wspeed'                                        &
     &          // '     ccover'

	  icall = 1
	end if

	lmeteo(:,1) = meteo(1,:)
	lmeteo(:,2) = meteo(2,:)
	lmeteo(:,3) = meteo(3,:)
	lmeteo(:,4) = meteo(4,:)
	lmeteo(:,5) = meteo(5,:)

	call scalar_out(iubase,nmeteo,lmeteo)

	call box_average(lmeteo(:,1),saver)
	call box_average(lmeteo(:,2),taver)
	call box_average(lmeteo(:,3),haver)
	call box_average(lmeteo(:,4),waver)
	call box_average(lmeteo(:,5),caver)

	format='(a,5f11.3)'
        write(iuaver,format) pline,sngl(saver),sngl(taver)              &
     &                  ,sngl(haver),sngl(waver),sngl(caver)

	end

!************************************************************

	subroutine fluxes_out

	use boxes

	implicit none

	character*20 aline

	integer, save :: iubase = 700
	integer, save :: icall = 0
	integer is,if,it,iu
	double precision :: taver,saver,zaver,vaver
	double precision, allocatable, save :: fluxes(:,:)
	character*80 format
	character*80 text

	if( icall == 0 ) then
	  allocate(fluxes(3,nsect))
	  do is=1,nsect
	    iu = iubase + is
            if = isect(1,is)      ! box from
            it = isect(2,is)      ! box to
            write(text,'(a,i5,a,i5)') 'fluxes between boxes [m**3/s]'   &
     &                  //' from box ',if,' to box ',it
	    call sect_out_header(iu,text)
	  end do

	  icall = 1
	end if

	call sect_out(iubase)

	end

!************************************************************
!************************************************************
!************************************************************

        subroutine sect_out_header(iu,what)

! writes header for scalar output

        use boxes

        implicit none

        integer iu
        character*(*) what
        integer nscal

        integer i,is
        character*80 formhead
        character*12, save :: descrp(3) =                               &
     &          (/'       total'                                        &
     &           ,'    positive'                                        &
     &           ,'    negative'/)

        formhead = ' '
        write(formhead,'(a,i2,a)') '(a,',nbox,'a)'

        write(iu,'(a,a)') '# ',trim(what)
        write(iu,formhead) '#      date_and_time'                       &
     &          ,(descrp(i),i=1,3)

        end

!************************************************************

        subroutine sect_out(iubase)

! outputs scalar

        use boxes

        implicit none

        integer iubase          !base unit to write to

        integer is,iu
        character*80 format

        write(format,'(a,i4,a)') '(a20,',3,'f12.4)'

        do is=1,nsect
          iu = iubase + is
          write(iu,format) pline,flux(:,is)
        end do

        end


!************************************************************
!************************************************************
!************************************************************

	subroutine box_average(val,valaver)

	use boxes

	implicit none

	double precision val(nbox)
	double precision valaver

	integer ib
	double precision aver,vol

	aver = 0.
	vol = 0.

	do ib=1,nbox
	  aver = aver + val(ib)*vnew(ib)
	  vol = vol + vnew(ib)
	end do

	valaver = aver / vol
	
	end

!************************************************************

