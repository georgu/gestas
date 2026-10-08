
! this is a shympi fake module that does nothing

!================================================================
	module shympi
!================================================================

	integer, save :: my_id = 0

        INTERFACE shympi_gather_and_sum
        MODULE PROCEDURE  shympi_gather_and_sum_i &
     &                   ,shympi_gather_and_sum_r 
        END INTERFACE

!================================================================
	contains
!================================================================

	function shympi_is_master()
	logical shympi_is_master
	shympi_is_master = .true.
	end

	function shympi_is_parallel()
	logical shympi_is_parallel
	shympi_is_parallel = .false.
	end

	function shympi_sum(i)
	integer shympi_sum
	integer i
	shympi_sum = i
	end

	subroutine shympi_gather_and_sum_i(i)
	integer i(:)
	end

	subroutine shympi_gather_and_sum_r(r)
	real r(:)
	end

	function id_node(node)		!fake id_node to be a function
	integer id_node
	integer node
	id_node = 0
	end

	subroutine shympi_barrier
	end

!================================================================
	end module shympi
!================================================================

