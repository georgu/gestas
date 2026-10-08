
!============================================================
	module bfm_estas
!============================================================

	implicit none

	integer, parameter :: kb = 2

	real, save :: dz(kb-1) = (/1./)
	real, save :: ZZ(kb-1) = (/0.5/)

	real, save :: depth_of_box = 3.

!============================================================
	end module bfm_estas
!============================================================

!************************************************************

	subroutine bfm_estas_init

	use mem
	!use global_mem, ONLY:ZERO,RLEN
	use api_bfm
	use init_var_bfm_local
	use time
	use bfm_estas

	implicit none

	NO_BOXES_X=1
	NO_BOXES_Y=1
	NO_BOXES_Z=KB-1
	NO_BOXES_XY=NO_BOXES_X*NO_BOXES_Y
	NO_BOXES=NO_BOXES_XY*NO_BOXES_Z

!   -----ALLOCATE AND SET LOGICAL MASKS FOR:-----
!
!           --WHOLE WATER COLUMN---
!	allocate(SEAmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
!           ---BOTTOM---
!	allocate(BOTmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
!           ---SURFACE--
!	allocate(SRFmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
!
!	SEAmask = .TRUE.
!	BOTmask = .TRUE.
!	SRFmask = .TRUE.

!	allocate(ZEROS(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
!	ZEROS = ZERO

!  -----PELAGIC+BENTHIC SETUP-----

! something strange below ??
! NO_STATES seems wrong
! I would change it to:
!
!        NO_BOXES_Z_BEN  = 1
!        NO_BOXES_BEN = NO_BOXES_XY * NO_BOXES_Z_BEN
!        NO_STATES_BEN = NO_BOXES_BEN * NO_D2_BOX_STATES_BEN
!        NO_STATES_PEL = NO_BOXES * NO_D3_BOX_STATES 
!        NO_STATES = NO_STATES_PEL + NO_STATES_BEN

	NO_STATES   = NO_D3_BOX_STATES * NO_BOXES + NO_D2_BOX_STATES_BEN
	NO_BOXES_Z_BEN  = 1
	NO_BOXES_BEN = NO_BOXES_XY * NO_BOXES_Z_BEN
	NO_STATES_BEN = NO_BOXES_BEN * NO_D2_BOX_STATES_BEN

!  -----INITIALISE LAND (0) SEA (1) MASK. IN 1D ONLY SEA....-----
!
!	allocate(ocepoint(NO_BOXES))
!	ocepoint = 1

	allocate(BOTindices(NO_BOXES_XY)); BOTindices = KB-1
	allocate(SRFindices(NO_BOXES_XY)); SRFindices = 1

	call init_bfm
	call set_var_info_bfm

	call init_var_bfm(bio_setup)
	call ini_organic_quotas()

	! where to set depth of boxes?

	!do k = 1 , NO_BOXES_Z
	!  Depth(k) = dz(k)*h
	!end do

	! D3STATEB = D3STATE !PELAGIC
	! D2STATEB_BEN = D2STATE_BEN ! BENTHIC

	timestep = 3600		! timestep of BFM in [s] ??

	end

!************************************************************

	subroutine bfm_estas_run

	!use mem
	!use global_mem, ONLY:ZERO,RLEN
	!use api_bfm
	!use init_var_bfm_local
	use bfm_estas

	implicit none

	integer ib
	real ts,ss,sw,ws,cs

	!call env_forcing_pom_bfm_1d
	
	ib = 1	!here we have to loop over boxes
	!call bfm_get_forcing(ib,ts,ss,sw,ws,cs)    !get forcings from estas
	call bfm_estas_forcing(ts,ss,sw,ws,cs)

	call EcologyDynamics
	!call lf1d		!what is this ?? (in POM tree) do we need it ??

        !call pom_dia_bfm(TT)
        !call ResetFluxes

	end

!************************************************************

	subroutine bfm_estas_forcing(ts,ss,sw,ws,cs)

	use mem
	use bfm_estas

	implicit none

	real ts		!water temperature [C]
	real ss		!salinity [psu]
	real sw		!solar radiation [W/m2]
	real ws		!wind speed [m/s]
	real cs		!current speed [m/s]

	integer k
	real, parameter :: rhoa = 1.25
	real, parameter :: rho0 = 1024.0
	real, parameter :: Cd = 1.40E-3
	real, parameter :: GRAV = 9.81
	real taub,tauw,h

	do k=1,kb-1
	  h = depth(k)
	  ETW(k) = ts		!Celsius ??
	  ESW(k) = ss		!psu ??
	  ERHO(k) = rho0
	  EPR(k) = -ERHO(k)*GRAV*ZZ(k)*H*1.e-5
	  ESS(k) = 0.		!INORGANIC SUSPENDED MATTER ??
	end do

	!EIR(1) = (-ONE)*SWRAD*rcp	! SHORTWAVE RADIATION
	EIR(1) = -sw			! SW is in W/M2	! should be positive ??

	!tauw = sqrt(wusurf**2+wvsurf**2)*rho0
	tauw = ws*rho0
	EWIND = sqrt(tauw/(rhoa*Cd))	!why no index ??

	!taub=sqrt(WUBOT**2+WVBOT**2)*rho0
	taub = cs*rho0
	ETAUB = taub			!why no index ??

!        -----TRANSFER PHYSICAL VARIABLES INTO BFM-----
!
!         call pom_to_bfm
!
!        -----COMPUTE VERTICAL EXTINCTION COEFFICIENTS-----
!
         call CalcVerticalExtinction
!
!        -----COMPUTE THE IRRADIANCE VERTICAL PROFILE-----
!
         call CalcLightDistribution

	end

!************************************************************

	subroutine CalcLightDistribution

! this subroutine is in the POM tree and is not included in standalone
! we replicate it here to avoid including POM routines

	use mem
	use mem_par
	use constants

	implicit none

	integer box_no		! why is BoxNumberZ a global symbol ??
	real dep

        EIR(1) = EIR(1)*p_PAR/ E2W

	do BoxNumberZ=2,NO_BOXES_Z

	  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	  !     The calculations below are for the lower boxes only
	  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          box_no = BoxNumberZ
	  dep = Depth(box_no-1)

          EIR(box_no) = EIR(box_no-1)*exp( -xEPS(box_no-1) * dep )

	end do

	end

!************************************************************
!************************************************************
!************************************************************

	subroutine bfm_estas_test

	write(6,*) 'start running bfm_estas_init'
	call bfm_estas_init
	write(6,*) 'end running bfm_estas_init'

	end

!************************************************************

	program bfm_estas_main
	call bfm_estas_test
	end

!************************************************************

