	! define and create mesh in the model
	module mesh_
	use element_
	use hdf5
	use h5interface_
	implicit none
	private
	public Mesh, Node
	integer, parameter :: dp=kind(0.d0)
	!
	! Defining different objects in mesh
	!
	type Node ! array of nodes
		!
		integer, allocatable, dimension(:) :: label ! list of element label
		real(dp), allocatable, dimension(:,:) :: coord ! list of coordinates
		!
	end type Node
	!
	!	
	type Intset ! a set of integer elements like nodeset elementset etc
		!
		character(:), allocatable :: name ! name of the nodeset
		integer, allocatable, dimension(:) :: set ! list of global node index
		!
	end type Intset
	!
	!
	type Mesh 
		!
		! metadata
		integer :: dim ! number of coordinate space
		integer :: num_node ! number of nodes
		integer :: num_element ! number of elements 
		integer :: num_elcon ! number of element connectivity
		integer :: num_nodeset ! number of nodesets
		integer :: num_elementset ! number of element sets 
		integer :: num_chunks ! number of chunks in the mesh
		!
		type(Node) :: node
		type(Element), allocatable, dimension(:) :: element
		integer, allocatable, dimension(:) :: elm_nodeConnectivity ! node list by elements used for G matrices
		!
		type(Intset), allocatable, dimension(:) :: nodeset ! set of node indices
		type(Intset), allocatable, dimension(:) :: elementset ! set of element indices
		!
		integer, dimension(4) :: chunk ! reference to element, b_node, i_node and bound
		!
		contains
		!
		procedure, pass(this) :: create_h5 => mesh_create_HDFInput
		!
	end type Mesh
	!
	!
	contains
	!
	subroutine mesh_create_HDFInput(this, mesh_id)
		! 
		class(Mesh), intent(inout) :: this ! mesh class
		integer(HID_T), intent(in) :: mesh_id ! the reference id number of the node array from hdf5 input
		integer(HID_T) :: node_id, elm_id, nset_id, elset_id, bem_id ! reference for node
		integer :: dim, num_node, num_element, num_nodeset, num_elementset, connectivity ! number of nodes and number of coordinates
		integer :: hdferr ! integer hdf5error reference
		!
		! reading attributes of the node array from HDF5 input (num_nodes and num_dim)
		! model id attribute
		call get_scalar_attribute_hdf5(attr_value=num_node, data_id=mesh_id, attr_name='num_node')
		call get_scalar_attribute_hdf5(attr_value=num_element, data_id=mesh_id, attr_name='num_element')
		call get_scalar_attribute_hdf5(attr_value=num_nodeset, data_id=mesh_id, attr_name='num_nodeset')
		call get_scalar_attribute_hdf5(attr_value=num_elementset, data_id=mesh_id, attr_name='num_elementset')
		call get_scalar_attribute_hdf5(attr_value=dim, data_id=mesh_id, attr_name='dim')
		call get_scalar_attribute_hdf5(attr_value=connectivity, data_id=mesh_id, attr_name='num_connectivity')
		!
		this%dim = dim
		this%num_node = num_node
		this%num_element = num_element
		this%num_elcon = connectivity
		this%num_nodeset = num_nodeset
		this%num_elementset = num_elementset
		!
		! open, read and close node group
		call h5gopen_f(mesh_id, 'node', node_id, hdferr)
		call mesh_create_nodes_HDFInput(this, node_id) ! calling local subroutine to assign nodes
		call h5gclose_f(node_id, hdferr)
		!
		! open, read and close element group
		call h5gopen_f(mesh_id, 'element', elm_id, hdferr)
		call mesh_create_elements_HDFInput(this, elm_id) ! calling local subroutine to assign elements
		call h5gclose_f(elm_id, hdferr)
		!
		! open, read and close element set group
		call h5gopen_f(mesh_id, 'nodeset', nset_id, hdferr)
		if (.not. allocated(this%nodeset)) allocate(this%nodeset(num_nodeset))
		call mesh_create_nodeset_HDFInput(this, nset_id) ! calling local subroutine to assign elements
		call h5gclose_f(nset_id, hdferr)
		!
		! open, read and close element set group
		call h5gopen_f(mesh_id, 'elementset', elset_id, hdferr)
		if (.not. allocated(this%elementset)) allocate(this%elementset(num_elementset))
		call mesh_create_elementset_HDFInput(this, elset_id) ! calling local subroutine to assign elements
		call h5gclose_f(elset_id, hdferr)
		!
		!
	end subroutine mesh_create_HDFInput
	!
	!
	subroutine mesh_create_nodes_HDFInput(this, node_id)
		!
		class(Mesh), intent(inout) :: this
		integer(HID_T), intent(in) :: node_id ! the reference id number of the node array from hdf5 input
		integer(HID_T) :: label_id, coord_id
		integer :: hdferr ! integer hdf5error reference
		integer(HSIZE_T), dimension(2) :: data_dims ! dimension of the data
		! allocation of node space
		if (.not. allocated(this%node%label)) allocate(this%node%label(this%num_node))
		if (.not. allocated(this%node%coord)) allocate(this%node%coord(this%dim,this%num_node))
		!
		! Assigning nodal coordinate values
		data_dims(1) = this%num_node
  		data_dims(2) = 1 ! one dimensional array of node labels
  		call h5dopen_f(node_id, 'label', label_id, hdferr) ! Opening the node label data
  		call h5dread_f(label_id, H5T_NATIVE_INTEGER, this%node%label, data_dims, hdferr)
  		call h5dclose_f(label_id, hdferr)
  		!
		data_dims(1) = this%dim
  		data_dims(2) = this%num_node
  		call h5dopen_f(node_id, 'coord', coord_id, hdferr)
		call h5dread_f(coord_id, H5T_NATIVE_DOUBLE, this%node%coord, data_dims, hdferr)
		call h5dclose_f(coord_id, hdferr)
		!this%node%coord = reshape(this%node%coord,(/this%dim, this%num_node/))
		!
		!
	end subroutine mesh_create_nodes_HDFInput
	!
	!
	subroutine mesh_create_elements_HDFInput(This, elm_id)
		!
		class(Mesh), intent(inout) :: this
		integer(HID_T), intent(in) :: elm_id ! the reference id number of the element array
		integer :: hdferr ! hdf5 error
		integer(HSIZE_T), dimension(2) :: data_dims ! dimension of the data
		integer(HID_T) :: index_id, elcon_id ! index of elements and connectivity
		! 
		integer, dimension(:,:), allocatable :: el_index
	
		!
		! allocation of elements
		if (.not. allocated(this%element)) allocate(this%element(this%num_element))
		!
		! Assigning element index array
		data_dims(1) = 4
  		data_dims(2) = this%num_element ! label, eltype, start node , end node
  		if (.not. allocated(el_index)) allocate(el_index(data_dims(1),data_dims(2)))
  		call h5dopen_f(elm_id, 'index', index_id, hdferr) ! Opening the node label data
  		call h5dread_f(index_id, H5T_NATIVE_INTEGER, el_index, data_dims, hdferr)
  		call h5dclose_f(index_id, hdferr)
  		el_index = reshape(el_index,(/4,this%num_element/))
  		write(*,*) 'element_index', el_index
  		!
		!
		! Assigning element connectivity array
		data_dims(1) = 1
  		data_dims(2) = this%num_elcon ! label, eltype, start node , end node
  		if (.not. allocated(this%elm_nodeConnectivity)) allocate(this%elm_nodeConnectivity(this%num_elcon))
  		call h5dopen_f(elm_id, 'connectivity', elcon_id, hdferr) ! Opening the node label data
  		call h5dread_f(elcon_id, H5T_NATIVE_INTEGER, this%elm_nodeConnectivity, data_dims, hdferr)
  		call h5dclose_f(elcon_id, hdferr)
  		write(*,*) 'element connectivity', this%elm_nodeConnectivity
		!
		! now we get element index and element connectivity in hand
		call createElement(this,el_index,this%elm_nodeConnectivity)
		!
		! deallocate temporary arrays
		if (allocated(el_index)) deallocate(el_index)
	
		!
	end subroutine mesh_create_elements_HDFInput
	!
	subroutine createElement(this,el_index, connectivity)
		!
		class(Mesh), intent(inout) :: this
		integer, dimension(:,:), intent(in) :: el_index ! contains array of label, eltype, start node, end node
		integer, dimension(:), intent(in) :: connectivity ! one dimensional array of nodes in the element
		integer :: ii
		!
		!
		do ii = 1, this%num_element
			this%element(ii)%label = el_index(1,ii) ! element label 
			this%element(ii)%eltype = el_index(2,ii) ! element type
			this%element(ii)%start_con = el_index(3,ii) ! start element connectivity index
			this%element(ii)%end_con = el_index(4,ii) ! end element connectivity index
			if (.not. allocated(this%element(ii)%nodeset)) allocate(this%element(ii)%nodeset(el_index(4,ii)-el_index(3,ii)+1))
			this%element(ii)%nodeset = connectivity(el_index(3,ii):el_index(4,ii))
			!
		end do 
	end subroutine createElement
	!
	!
	subroutine mesh_create_nodeset_HDFInput(this, group_id)
		!
		class(Mesh), intent(inout) :: this
		integer(HID_T), intent(in) :: group_id ! hdf5 id for node or element set
		integer(HID_T) :: set_id
		integer(HSIZE_T), dimension(2) :: data_dims ! dimension of the data
		character(:), allocatable :: data_name
		integer :: set_length, ii, hdferr
		character(:), allocatable :: set_name
		character(7) :: fmt
		!
		do ii = 1,this%num_nodeset
			!
			if (ii<10) then ! dset name 
				fmt = '(A4,I1)'
				data_name = 'nset0'
				write(data_name,fmt) 'nset', ii ! name of the model group
			else ! model10 to model99
				fmt = '(A4,I2)'
				data_name = 'nset00'
				write(data_name,fmt) 'nset', ii ! name of the model group
			end if
			!

			! now read from the dataset
			call h5dopen_f(group_id, data_name, set_id, hdferr) ! Opening the nodeset dataset
			call get_scalar_attribute_hdf5(attr_value=set_length, data_id=set_id, attr_name='len')
			if (.not. allocated(this%nodeset(ii)%set)) allocate(this%nodeset(ii)%set(set_length))
			!call get_scalar_attribute_hdf5(attr_value=set_name, data_id=set_id, attr_name='name')
			!this%nodeset(ii)%name = set_name
			data_dims(1) = set_length
  			data_dims(2) = 1
  			call h5dread_f(set_id, H5T_NATIVE_INTEGER, this%nodeset(ii)%set, data_dims, hdferr)
  			call h5dclose_f(set_id, hdferr)

		end do 

	end subroutine mesh_create_nodeset_HDFInput
	!
	!
	subroutine mesh_create_elementset_HDFInput(this, group_id)
		!
		class(Mesh), intent(inout) :: this
		integer(HID_T), intent(in) :: group_id ! hdf5 id for node or element set
		integer(HID_T) :: set_id
		integer(HSIZE_T), dimension(2) :: data_dims ! dimension of the data
		character(:), allocatable :: data_name
		integer :: set_length, ii, hdferr
		character(:), allocatable :: set_name
		character(7) :: fmt
		!
		do ii = 1,this%num_elementset
			!
			if (ii<10) then ! dset name 
				fmt = '(A5,I1)'
				data_name = 'elset0'
				write(data_name,fmt) 'elset', ii ! name of the model group
			else ! model10 to model99
				fmt = '(A5,I2)'
				data_name = 'elset00'
				write(data_name,fmt) 'elset', ii ! name of the model group
			end if
			!

			! now read from the dataset
			call h5dopen_f(group_id, data_name, set_id, hdferr) ! Opening the elementset dataset
			call get_scalar_attribute_hdf5(attr_value=set_length, data_id=set_id, attr_name='len')
			if (.not. allocated(this%elementset(ii)%set)) allocate(this%elementset(ii)%set(set_length))
			!call get_scalar_attribute_hdf5(attr_value=set_name, data_id=set_id, attr_name='name')
			!this%elementset(ii)%name = set_name
			data_dims(1) = set_length
  			data_dims(2) = 1
  			call h5dread_f(set_id, H5T_NATIVE_INTEGER, this%elementset(ii)%set, data_dims, hdferr)
  			call h5dclose_f(set_id, hdferr)

		end do 

	end subroutine mesh_create_elementset_HDFInput
	!
	
	end module mesh_