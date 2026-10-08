
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016,2018-2020  Georg Umgiesser
!    Copyright (C) 2016,2019  Erik Pascolo
!    Copyright (C) 2016  Leslie Aveytua
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

! routines for BFM module
!
! revision log :
!
! 22.02.2016	ggu&erp	new bfm routines created from newconz
! 22.03.2016	ggu	changed VERS_7_5_6
! 06.06.2016	ggu	initialization from file changed
! 09.06.2016 	laa	modified (Search LAA)
! 10.06.2016	ggu	changed VERS_7_5_13
! 17.06.2016	ggu	changed VERS_7_5_15
! 28.06.2016	ggu	initialize bfmv, new routine bfm_check()
! 09.09.2016	ggu	changed VERS_7_5_17
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 01.05.2019	erp	restructured for new BFM
! 17.10.2019	ggu	restructured and adapted
! 09.04.2020	ggu	new routine bfm_run()

!==================================================================
        module mod_bfm
!==================================================================

        implicit none

!#include "BFM_var_list.h"      

        integer,  save :: nst_bfm = 0
        integer,  save :: nkn_bfm = 0
        integer,  save :: nlv_bfm = 0
        integer,  save :: nst_bfm_d = 42 ! number of variables among the 114 variables of d to be saved in bfm_d

        real, allocatable, save :: bfmv(:,:,:)
        real, allocatable, save :: bfm_d(:,:,:)  ! here we save nst_bfm_d variables among the 114 of d
        real, allocatable, save :: ave_bfmv(:,:,:)
        real, allocatable, save :: ave_bfm_d(:,:,:)  ! here we save nst_bfm_d variables among the 114 of d
        real, allocatable, save :: ave_bfm_benthic(:,:,:)  ! here we save nst_bfm_d variables among the 114 of d
        real, allocatable, save :: wsink_bfm(:,:,:)  ! here we save nst_bfm_d variables among the 114 of d
        real, allocatable, save :: rload_bfm(:,:,:)  ! here we save nst_bfm_d variables among the 114 of d
!       Benthic variables
        real, allocatable, save :: bfm_benthic(:,:,:)

	integer, save :: ibfm = 0
	integer, save :: ibfm_state = 54	!number of state variables
        integer, save :: icall_bfm = 0
        integer, save :: ninfo_bfm = 0
	integer, save :: ibfm_state_benthic = 9!number of benthic state variables

        logical, private, save :: bdebug = .false.
        logical, save :: binfo_bfm = .true.

        real, save :: cref,rkpar,difmol

        real, save :: counter_ave = 0.

        real, save, allocatable :: bfminit(:)
        real, save, allocatable :: bfmbound(:)

	character*3, save :: what = 'bfm'
        logical,save :: mod_bfm_initialized=.False. 
        integer, save :: ieco_av=0

	integer, save :: nkn = 1		!number of boxes
	integer, save :: nlvdi = 1		!number of layers

!==================================================================
	contains
!==================================================================

        subroutine mod_bfm_init(nst,nkn,nlv)

	implicit none

        integer nst
        integer nkn
        integer nlv

        if( nst == nst_bfm .and. nkn == nkn_bfm .and. nlv == nlv_bfm ) return

        if( nst > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( nst == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nst,nkn,nlv: ',nst,nkn,nlv
            stop 'error stop mod_bfm_init: incompatible parameters'
          end if
        end if

        if( nkn_bfm > 0 ) then
          deallocate(bfmv)
          deallocate(bfm_d)
        end if

        nst_bfm = nst
        nkn_bfm = nkn
        nlv_bfm = nlv

        if( nkn == 0 ) return

	write(6,*) 'bfm mod init: ',nst,nkn,nlv
        allocate(bfmv(nlv,nkn,nst))
	bfmv = 0.
        allocate(bfm_d(nlv,nkn,nst_bfm_d))
	bfm_d = 0.
        allocate(wsink_bfm(nlv,nkn,nst))
        wsink_bfm(:,:,:)=0.0
        allocate(rload_bfm(nlv,nkn,nst))
        rload_bfm(:,:,:)=0.0
        
!       Benthic module
        allocate(bfm_benthic(1,nkn,ibfm_state_benthic))
        bfm_benthic(1,:,1)=3.0    ! PO4
        bfm_benthic(1,:,2)=30.0   ! NO3
        bfm_benthic(1,:,3)=200.0  ! NH4
        bfm_benthic(1,:,4)=8.0    ! SiO3
        bfm_benthic(1,:,5)=48000. ! DIC
        bfm_benthic(1,:,6)=.0     ! R6c
        bfm_benthic(1,:,7)=.0     ! R6n
        bfm_benthic(1,:,8)=.0     ! R6p
        bfm_benthic(1,:,9)=.0     ! R6s

        end subroutine mod_bfm_init

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine bfm_run
      
	implicit none
      
        call bfm_reactor

	end subroutine bfm_run

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine bfm_init

! initializes bfm computation

	implicit none

!-------------------------------------------------------------
! initialization of module
!-------------------------------------------------------------

        if( .not. mod_bfm_initialized ) then 

          call mod_bfm_init(ibfm_state,nkn,nlvdi)
	  call bfm_init_internal

          write(6,*) 'bfm initialized in bfm_init: ',ibfm_state,nkn,nlvdi
          mod_bfm_initialized = .true.

        end if

	end subroutine bfm_init

!*********************************************************************

      subroutine bfm_reactor
      
      implicit none
      
      call bfm_reactor_internal
      call bfm_benthic_internal

      end subroutine bfm_reactor

!==================================================================
        end module mod_bfm
!==================================================================

