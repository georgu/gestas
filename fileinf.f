
!***************************************************************
!
! fileinf - routine to check one and compare two general files
!
! file format:
!
!	time_record
!	time_record
!	...
!
! time_record:
!
!	nskip
!	...nskip lines to skip
!	nline nc1 nc2
!	line 1
!	line 2
!	...
!	line nline
!
! line:
!
!	col1 col2 ... col_n2
!
!	only columns from nc1 to nc2 are considered (nc1 <= nc2)
!
!***************************************************************

	module mod_fileinf

	logical, save :: bsilent = .false.
	logical, save :: bquiet = .false.
	logical, save :: bnodiff = .false.
	logical, save :: bverbose = .false.
	logical, save :: bnostop = .false.
	logical, save :: bstop = .true.
	double precision, save :: maxdiff = 0.

	double precision, save :: dstart,dend,dt
	character*20, save :: astart,aend

        INTERFACE
	  subroutine read_record(iu,n,values,ierr)
	  integer iu
	  integer n
	  double precision, allocatable :: values(:)
	  integer ierr
          END subroutine
        END INTERFACE

	end module mod_fileinf

!***************************************************************

	program fileinf

	use clo
	use mod_fileinf

	implicit none

	integer nc,ierr

	call fileinf_init

        nc = clo_number_of_files()

        if( nc == 0 ) then
          call clo_usage
        else if( nc == 1 ) then
          call read_file_file
        else if( nc == 2 ) then
          call compare_file_files(ierr)
        else
          write(6,*) 'nc = ',nc
          stop 'error stop fileinf: wrong number of files'
        end if

        if( ierr > 0 ) then
          if( ierr == 99 ) ierr = 100   !terrible hack - FIXME
          call exit(ierr)
        else
          call exit(99)
        end if

	end

!***************************************************************

	subroutine read_file_file

	use mod_fileinf
	use clo

	implicit none

	logical bwrite				!write single time records
	integer irec,n
	integer iu,ierr,ios
	character*80 name_one

	double precision, allocatable :: values(:)

	bwrite = bverbose
	irec = 0

!--------------------------------------------------------------
! open file
!--------------------------------------------------------------

        call clo_get_file(1,name_one)

	iu = 1
	open(iu,file=name_one,status='old',form='formatted',iostat=ios)

        if( ios /= 0 ) then
	  write(6,*) 'cannot open file: ',trim(name_one)
	  stop 'error opening file'
	end if

!--------------------------------------------------------------
! time loop
!--------------------------------------------------------------

	do
	  call read_record(iu,n,values,ierr)
	  if( ierr /= 0 ) exit
	  irec = irec + 1

	  if( bwrite ) write(6,*) irec
	end do

!--------------------------------------------------------------
! end of time loop
!--------------------------------------------------------------

	if( ierr > 0 ) stop 'error stop: read error'

	write(6,*) 'final message:'
	write(6,*) 'records read: ',irec

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!***************************************************************

	subroutine compare_file_files(ierr)

	use mod_fileinf
	use clo

	implicit none

	integer ierr

	logical bwrite				!write single time records
	double precision verr,verrmax
	integer ios,irec,iu1,iu2
	integer n,n1,n2
	character*80 name_one,name_two

	double precision, allocatable :: values1(:),values2(:)

	bwrite = bverbose
	irec = 0
	ierr = 0
	verrmax = 0.

!--------------------------------------------------------------
! open files
!--------------------------------------------------------------

        call clo_get_file(1,name_one)
        call clo_get_file(2,name_two)

	iu1 = 1
	open(iu1,file=name_one,status='old',form='formatted',iostat=ios)

        if( ios /= 0 ) then
	  write(6,*) 'cannot open file: ',trim(name_one)
	  stop 'error opening file'
	end if

	iu2 = 2
	open(iu2,file=name_two,status='old',form='formatted',iostat=ios)

        if( ios /= 0 ) then
	  write(6,*) 'cannot open file: ',trim(name_two)
	  stop 'error opening file'
	end if

!--------------------------------------------------------------
! loop over records
!--------------------------------------------------------------

	do
	  call read_record(iu1,n1,values1,ierr)
	  if( ierr /= 0 ) exit
	  call read_record(iu2,n2,values2,ierr)
	  if( ierr /= 0 ) exit
	  if( n1 /= n2 ) stop 'error stop: n1/= n2'
	  n = n1
	  irec = irec + 1

	  verr = maxval( abs( values1 - values2 ) )
	  verrmax = max(verrmax,verr)

	  if( bwrite ) write(6,*) irec,verr
	end do

!--------------------------------------------------------------
! end loop over records
!--------------------------------------------------------------

	if( ierr > 0 ) stop 'error stop: read error'

	write(6,*) 'final message:'
	write(6,*) 'records read: ',irec
	write(6,*) 'maximum difference:    ',verrmax

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*****************************************************

	subroutine read_record(iu,n,values,ierr)

	implicit none

	integer iu
	integer n
	double precision, allocatable :: values(:)
	integer ierr

	logical balloc
	integer i,nskip,nline,nc1,nc2,nc,ip
	double precision, allocatable :: val(:)

	read(iu,*,end=9) nskip

	do i=1,nskip
	  read(iu,*)
	end do

	read(iu,*,end=9) nline,nc1,nc2
	nc = nc2-nc1+1
	n = nline * nc
	allocate(val(nc2))

	balloc = .false.
	if( allocated(values) ) then
	  if( size(values) < n ) then
	    deallocate(values)
	    balloc = .true.
	  end if
	else
	  balloc = .true.
	end if
	if( balloc ) allocate(values(n))

	ip = 0
	do i=1,nline
	  read(iu,*) val(:)
	  values(ip+1:ip+nc) = val(nc1:nc2)
	  ip = ip + nc
	end do

	ierr = 0

	return
    9	continue
	ierr = -1
	return
	end

!*****************************************************

        subroutine fileinf_init

        use mod_fileinf
        use clo

        implicit none

        logical baux
        character*80 version

        version = '1.0'

        call clo_init('fileinf','file-file(s)',trim(version))

        call clo_add_info('checks box files')

        call clo_add_sep('general options:')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('nodiff',.false.,'do not show differences')
        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('nostop',.false.,'do not stop at error')
        call clo_add_option('maxdiff',0.,'maximum tolerated difference')

        call clo_parse_options

        call clo_get_option('silent',bsilent)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('nodiff',bnodiff)
        call clo_get_option('verbose',bverbose)
        call clo_get_option('nostop',baux)
        call clo_get_option('maxdiff',maxdiff)

        if( baux ) bstop = .false.
        if( bsilent ) bquiet = .true.
        if( bquiet ) bverbose = .false.

        end

!*****************************************************

