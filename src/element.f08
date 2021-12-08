	! computational element parameters are included in this 
	! 
	module element_
    !use hdf5
	implicit none
	private
    public element
	integer, parameter :: dp=kind(0.d0)
	!
    !   
    type Element ! a single element definition
        !
        integer :: label ! element label
        integer :: eltype ! the element type recognition for computation
        integer :: num_node ! number of nodes in the element
        integer :: start_con ! start index of first node in connectivity list
        integer :: end_con ! end index of the last node in connectivity list
        integer, allocatable, dimension(:) :: nodeset ! list of global node index
        real(dp), allocatable, dimension(:,:) :: coord ! intrinsic coordinates of element nodes
        !
        contains
        !
        procedure, pass(this) :: assign_coord => elm_assign_coord_to_element
        !
    end type Element
    !
	!
	contains
	!
    !
    subroutine elm_assign_coord_to_element(this, node_coord)
        !
        ! assigns coordinate values to respective element
        class(Element), intent(inout) :: this ! this element
        real(dp), intent(in), dimension(:,:) :: node_coord ! current mesh to which the element belongs to
        integer :: num_node, dim, ii
        !
        ! allocating space to the coord parameter of element
        num_node = size(this%nodeset) ! number of nodes in the element
        dim = size(node_coord,1) ! number of coordinate spaces for each nodes
        if (.not. allocated(this%coord)) allocate(this%coord(dim, num_node))
        !
        !now assigning the coordinates
        do ii=1, num_node
            ! the coordinates values assigned from node indes in the nodeset of the element
            this%coord(:,ii) = node_coord(:,this%nodeset(ii))
            this%num_node = num_node
        end do
        !
    end subroutine elm_assign_coord_to_element
    !
    !
    end module element_