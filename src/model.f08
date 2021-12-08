	! A numerical model for computation 
	module model_
	use hdf5
	use h5interface_
	use mesh_
	use material_
	use chunk_
	!
	implicit none
	private
	integer, parameter :: dp=kind(0.d0)
	!
	type Model_step
		! gathers information for each step in the model
		integer :: stype ! code for computation type
		real(dp) :: frequency ! 
		!
	end type Model_step
	!
	! defining a Model class
	type, public :: Model
		!
		! metadata
		character(:), allocatable :: name ! name of the model
		integer :: id ! unique id given to model 
		integer :: dim ! model dimension (2 for 2D and so on)
		integer :: num_material ! number of materials in the model
		integer :: num_bchunk ! number of BEM superelement
		integer :: num_step ! number of steps
		!
		! members
		type(Mesh) :: mesh ! only one mesh is assumed to exist in a model
		type(Material), allocatable, dimension(:) :: material ! may contain multiple materials
		type(BEM_chunk), allocatable, dimension(:) :: bchunk
		type(Model_step), allocatable, dimension(:) :: step
		!
		contains
		!
		procedure, pass(this) :: create => model_create
		procedure, pass(this) :: create_h5 => model_create_hdfInput
		!
	end type Model
	!
	contains
	!
	subroutine model_create(this, id, dim, nMat, nStep, name)
		class(model), intent(inout) :: this
		integer, intent(in) :: id, dim ! model id and model dimension
		character(*), optional, intent(in) :: name ! name of the model
		integer, intent(in), optional :: nMat, nStep ! number of material and step
		character(7) :: string_format
		!
		this%id = id
		this%dim = dim
		!
		! formatting model id for default model name
		if (id<10) then
			string_format = '(A5,I1)'
		else ! id = 10 to 99
			string_format = '(A5,I2)'
		end if
		!
		! writing name for the model
		if (present(name)) then
			this%name = name
		else
			write(this%name, string_format) 'Model', id
		end if
		! 
		! assigning number of materials atrribute
		if (present(nMat)) then
			this%num_material = nMat
		else
			this%num_material = 1
		end if
		!
		! assigning number of step attributes
		if (present(nStep)) then
			this%num_step = nStep
		else
			this%num_step = 1
		end if
		!
	end subroutine model_create
	!
	!
	subroutine model_create_hdfInput(this, file_name)
		!
		class(Model), intent(inout) :: this
		character(*), intent(in) :: file_name
		character(:), allocatable :: model_name
		character(:), allocatable :: hdf_char
		character(:), allocatable :: mat_name ! name of the material index
		character(:), allocatable :: step_name
		integer(HID_T) :: model_id, mesh_id, material_group_id, material_id
		integer(HID_T) :: step_group_id, step_id
		integer :: hdf_int
		integer :: hdferr, ii
		character(7) :: fmt ! material name format
		!
		! initiate fortran interface for hdf5 file
		call h5open_f(hdferr)
		!
		! opening existing hdf5 file (input file)
		call h5fopen_f (file_name, H5F_ACC_RDONLY_F, model_id, hdferr)
		!
		! now reading attributes of each model
		!
		! model name attribute
		!if (allocated(hdf_char)) deallocate(hdf_char)
		!hdf_char = 'name'
		!call get_scalar_attribute_hdf5(attr_value=model_name, data_id=model_id, attr_name=hdf_char)
		!this%name = model_name
		!
		! model dimension attribute
		call get_scalar_attribute_hdf5(attr_value=hdf_int, data_id=model_id, attr_name='id')
		this%dim = hdf_int
		!
		! model dimension attribute
		call get_scalar_attribute_hdf5(attr_value=hdf_int, data_id=model_id, attr_name='dim')
		this%dim = hdf_int
		!
		! number of BEM chunks
		call get_scalar_attribute_hdf5(attr_value=hdf_int, data_id=model_id, attr_name='num_BEM_chunk')
		this%num_bchunk = hdf_int
		!
		! number of material
		call get_scalar_attribute_hdf5(attr_value=hdf_int, data_id=model_id, attr_name='num_material')
		this%num_material = hdf_int
		!
		!
		! number of steps
		call get_scalar_attribute_hdf5(attr_value=hdf_int, data_id=model_id, attr_name='num_step')
		this%num_step = hdf_int
		!
		!
		! opening mesh group within model
		call h5gopen_f(model_id, 'mesh', mesh_id, hdferr)
		!
		! now assigning mesh attributes to the model (node, elements, BEM chunks etc.)
		call this%mesh%create_h5(mesh_id)
		!
		! closing open mesh group
		call h5gclose_f(mesh_id, hdferr)
		!
		! 
		! creating bchunk groups
		if (.not. allocated(this%bchunk)) allocate(this%bchunk(this%num_bchunk)) ! allocate space for BEM chunks
		!
		call model_create_BEMchunk_HDFInput(this, model_id)
		!
		!
		! reading material in the model
		!
		if (.not. allocated(this%material)) allocate(this%material(this%num_material)) ! allocating materials
		!
		! opening material group
		call h5gopen_f(model_id, 'material', material_group_id, hdferr)
		!
		! opening each material folder
		do ii = 1,this%num_material
			!
			if (ii<10) then ! dset name 
				fmt = '(A3,I1)'
				mat_name = 'mat0'
				write(mat_name,fmt) 'mat', ii ! name of the material with index
			else ! model10 to model99
				fmt = '(A3,I2)'
				mat_name = 'mat00'
				write(mat_name,fmt) 'mat', ii ! name of the material with index
			end if
			!
			! opening material name mat1, mat2 etc
			call h5gopen_f(material_group_id, mat_name, material_id, hdferr)
			! now calling the material create module from the material itself
			!call this%material(ii)%create_elastic_h5(material_id)
			call create_elastic_material_hdf5(this%material(ii), material_id)
			call h5gclose_f(material_id, hdferr)
		!
		end do ! material loop
		!
		!
		call h5gclose_f(material_group_id, hdferr)
		!
		!
		if (.not. allocated(this%step)) allocate(this%step(this%num_step)) ! allocating materials
		! opening step group
		call h5gopen_f(model_id, 'step', step_group_id, hdferr)
		!
		! opening each material folder
		do ii = 1,this%num_step
			!
			if (ii<10) then ! dset name 
				fmt = '(A4,I1)'
				step_name = 'step0'
				write(step_name,fmt) 'step', ii ! name of the step with index
			else if (ii<100) then ! model10 to model99
				fmt = '(A4,I2)'
				step_name = 'step00'
				write(step_name,fmt) 'step', ii ! name of the step with index
			else ! model10 to model99
				fmt = '(A4,I3)'
				step_name = 'step000'
				write(step_name,fmt) 'step', ii ! name of the step with index
			end if
			!
			! opening step name mat1, mat2 etc
			call h5gopen_f(step_group_id, step_name, step_id, hdferr)
			! 
			write(*,*) 'Check for step names:', step_name
			!
			call create_step_hdf5(this%step(ii), step_id)
			call h5gclose_f(step_id, hdferr)
		!
		end do ! step loop
		!
		call h5gclose_f(step_group_id, hdferr)
		!
		call h5fclose_f(model_id, hdferr) !closing current model
		!
	end subroutine model_create_hdfInput
	!
	!
	subroutine model_create_BEMchunk_HDFInput(this, model_id)
		!
		class(Model), intent(inout) :: this
		integer(HID_T), intent(in) :: model_id ! hdf5 id for node or element set
		integer(HID_T) :: bem_group_id
		character(:), allocatable :: bem_name
		integer :: ii, hdferr, hdf_int
		character(7) :: fmt
		integer(HID_T) :: bem_id, set_id ! id of the bem group and then node or element sets inside
		integer(HSIZE_T), dimension(2) :: data_dims ! dimension of the data
		!
		!
		! opening bchunk group within model
		call h5gopen_f(model_id, 'bem_chunk', bem_group_id, hdferr)
		!
		do ii = 1,this%num_bchunk
			!
			if (ii<10) then ! dset name 
				fmt = '(A3,I1)'
				bem_name = 'bem0'
				write(bem_name,fmt) 'bem', ii ! name of the model group
			else ! model10 to model99
				fmt = '(A3,I2)'
				bem_name = 'bem00'
				write(bem_name,fmt) 'bem', ii ! name of the model group
			end if
			!
			!call this%bchunk(ii)%create_h5(bem_name, bem_group_id)

			! now for each bem group
		call h5gopen_f(bem_group_id, bem_name, bem_id, hdferr) ! Opening the bem group
		!
		call get_scalar_attribute_hdf5(attr_value=hdf_int, data_id=bem_id, attr_name='bound')
		this%bchunk(ii)%bound = hdf_int
		call get_scalar_attribute_hdf5(attr_value=hdf_int, data_id=bem_id, attr_name='num_bnode')
		this%bchunk(ii)%num_bnode = hdf_int
		call get_scalar_attribute_hdf5(attr_value=hdf_int, data_id=bem_id, attr_name='num_belm')
		this%bchunk(ii)%num_belm = hdf_int
		call get_scalar_attribute_hdf5(attr_value=hdf_int, data_id=bem_id, attr_name='num_inode')
		this%bchunk(ii)%num_inode = hdf_int
		call get_scalar_attribute_hdf5(attr_value=hdf_int, data_id=bem_id, attr_name='material')
		this%bchunk(ii)%material = hdf_int

		! reading boundary nodes in the chunk
		if (.not. allocated(this%bchunk(ii)%b_node)) allocate(this%bchunk(ii)%b_node(this%bchunk(ii)%num_bnode))
		data_dims(1) = this%bchunk(ii)%num_bnode
  		data_dims(2) = 1
		call h5dopen_f(bem_id, 'bnode', set_id, hdferr) ! Opening the boundary node dataset
		call h5dread_f(set_id, H5T_NATIVE_INTEGER, this%bchunk(ii)%b_node, data_dims, hdferr)
		call h5dclose_f(set_id, hdferr)
		!
		! reading boundary elements in the chunk
		if (.not. allocated(this%bchunk(ii)%b_element)) allocate(this%bchunk(ii)%b_element(this%bchunk(ii)%num_belm))
		data_dims(1) = this%bchunk(ii)%num_belm
  		data_dims(2) = 1
		call h5dopen_f(bem_id, 'belm', set_id, hdferr) ! Opening the boundary node dataset
		call h5dread_f(set_id, H5T_NATIVE_INTEGER, this%bchunk(ii)%b_element, data_dims, hdferr)
		call h5dclose_f(set_id, hdferr)
		!
		if (this%bchunk(ii)%num_inode > 0) then
			! reading internal nodes in the chunk
			if (.not. allocated(this%bchunk(ii)%i_node)) allocate(this%bchunk(ii)%i_node(this%bchunk(ii)%num_inode))
			data_dims(1) = this%bchunk(ii)%num_inode
  			data_dims(2) = 1
			call h5dopen_f(bem_id, 'inode', set_id, hdferr) ! Opening the boundary node dataset
			call h5dread_f(set_id, H5T_NATIVE_INTEGER, this%bchunk(ii)%i_node, data_dims, hdferr)
			call h5dclose_f(set_id, hdferr)
		end if
		!
		! closing each bem group
		call h5gclose_f(bem_id, hdferr)

			!
			!			
		end do 
		!
		! closing open groups
		call h5gclose_f(bem_group_id, hdferr)
		!
	end subroutine model_create_BEMchunk_HDFInput
	!
	subroutine create_elastic_material_hdf5(this, mat_id)
        !
        class(Material), intent(inout) :: this
        integer(HID_T), intent(in) :: mat_id
        integer(HID_T) :: prop_id ! material properties id
        integer(HSIZE_T), dimension(2) :: data_dims ! dimension of the data
        integer :: hdferr ! integer hdf5error reference
        !
        ! reading elastic material properties
        data_dims(1) = 1
        data_dims(2) = 1
        call h5dopen_f(mat_id, 'density', prop_id, hdferr)
        call h5dread_f(prop_id, H5T_NATIVE_DOUBLE, this%RHO, data_dims, hdferr)
        call h5dclose_f(prop_id, hdferr)
        write(*,*) 'density in subroutine:', this%RHO
        !
        call h5dopen_f(mat_id, 'E', prop_id, hdferr)
        call h5dread_f(prop_id, H5T_NATIVE_DOUBLE, this%E, data_dims, hdferr)
        call h5dclose_f(prop_id, hdferr)
        !
        call h5dopen_f(mat_id, 'nu', prop_id, hdferr)
        call h5dread_f(prop_id, H5T_NATIVE_DOUBLE, this%nu, data_dims, hdferr)
        call h5dclose_f(prop_id, hdferr)
        !
        ! calculating lami's parameters
        this%MU = this%E/(2.0_dp+2.0_dp*this%nu)
        this%LAMBDA = this%MU*2.0_dp*this%nu/(1.0_dp-2.0_dp*this%nu)
        !
        call this%complex_elastic()
        !
    end subroutine create_elastic_material_hdf5
    !
    !
    subroutine create_step_hdf5(this, step_id)
    	!
    	class(Model_step), intent(inout) :: this
    	integer(HID_T), intent(in) :: step_id
    	integer(HID_T) :: prop_id
    	integer(HSIZE_T), dimension(2) :: data_dims ! dimension of the data
        integer :: hdferr ! integer hdf5error reference
        !
        ! reading a step (frequency and step type)
        data_dims(1) = 1
        data_dims(2) = 1
        call h5dopen_f(step_id, 'step_type', prop_id, hdferr)
        call h5dread_f(prop_id, H5T_NATIVE_INTEGER, this%stype, data_dims, hdferr)
        call h5dclose_f(prop_id, hdferr)

        call h5dopen_f(step_id, 'frequency', prop_id, hdferr)
        call h5dread_f(prop_id, H5T_NATIVE_DOUBLE, this%frequency, data_dims, hdferr)
        call h5dclose_f(prop_id, hdferr)
        !
    end subroutine create_step_hdf5
    !
	!
	end module model_