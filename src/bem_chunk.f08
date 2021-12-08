	! define a BEM chunk
	module chunk_
	!use hdf5
	!use h5interface_
	implicit none
	private
	public BEM_chunk
	!
	integer, parameter :: dp=kind(0.d0)
	!
	!
	!	
	type BEM_chunk 
		! A unit chunk or part of BEM layer or intity
		! metadata
		integer :: label ! label for BEM chunk
		integer :: bound ! boundary type
		integer :: num_bnode, num_belm, num_inode, num_elcon ! number of element connectivity
		!
		! members
		integer, allocatable, dimension(:) :: b_node ! list of node indices in boundary
		integer, allocatable, dimension(:) :: b_element ! list of element indices in boundary
		integer, allocatable, dimension(:) :: i_node ! list of internal node indices
		integer :: material ! material reference id
		!
		!
		contains
		!
		procedure, pass(this) :: coord_create => chunk_create_node_coordinates
		!
	end type BEM_chunk
	!
	!
	contains
	!
	!
	subroutine chunk_create_node_coordinates(this, coord, gnode)
		!
		! creates array of coordinates for this bem_chunk
		!
		class(BEM_chunk), intent(in) :: this
		real(dp), dimension(:,:), intent(in) :: gnode ! global node coordinates array
		integer :: cdim, nnode ! cartesian coordinate dimensions and num nodes
		integer :: indx, ii ! node index and loop bariables
		real(dp), intent(inout) :: coord(:,:) ! Nodal coordinates for chunk
		!
		cdim = size(gnode,1)
		nnode = size(this%b_node)
		!
		! allocation of nodal coordinate space
		!if (allocated(coord)) deallocate(coord)
		!allocate(coord(cdim, nnode))
		!
		do ii = 1, nnode
			!
			! current index
			!
			indx = this%b_node(ii)
			coord(:,ii) = gnode(:,indx)
			!
		end do
		! 
		!
	end subroutine chunk_create_node_coordinates
	!
	!
	end module chunk_
