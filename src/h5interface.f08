	! hdf5 interfaces
	module h5interface_
	use hdf5
	implicit none
	private
	public get_scalar_attribute_hdf5
	!
	!
	interface get_scalar_attribute_hdf5
		module procedure get_scalar_attribute_hdf5_integer 
		module procedure get_scalar_attribute_hdf5_char
	end interface get_scalar_attribute_hdf5
	!
	contains
	subroutine get_scalar_attribute_hdf5_integer(attr_value, data_id, attr_name)
		!
		integer, intent(out) :: attr_value ! the value of the attribute to be returned
		character(*), intent(in) :: attr_name ! name of the attribute
		integer(HID_T), intent(in) :: data_id ! File identifier
		integer :: hdferr ! error flag
		logical :: attr_exists ! .TRUE. if exists, .FALSE. otherwise
		integer(HID_T) :: attr_id ! Attribute identifier
		integer(HID_T) :: atype_id ! Attribute type identifier
		integer(HID_T) :: aspace_id ! Attribute space identifier
		integer, parameter :: arank = 1 ! rank of the attribute
		integer(HSIZE_T), parameter, dimension(1) :: adims = (/1/) ! dimension of attribute data
		! check if attribute exists
		call h5aexists_f(data_id, attr_name, attr_exists, hdferr)
		!
		if (attr_exists) then
			! open the 'num_models' attribute
			call h5aopen_f(data_id, attr_name, attr_id, hdferr)
			!
			! Create scalar data space for the attribute.
  			call h5screate_simple_f(arank, adims, aspace_id, hdferr)
  			!
			! Get dataspace and allocate memory for read buffer.
  			call h5aget_space_f(attr_id, aspace_id, hdferr)
  			!
  			! reading the attribute
  			call h5aread_f(attr_id, H5T_NATIVE_INTEGER, attr_value, adims, hdferr) 
  			!
  			! closing the attribute
  			call h5aclose_f(attr_id, hdferr)
		else
			! through error
			write(*,*) 'Attribute ''' // attr_name // ''' does not exist.'
		end if
		!
	end subroutine get_scalar_attribute_hdf5_integer
	!
	!
	subroutine get_scalar_attribute_hdf5_char(attr_value, data_id, attr_name)
		!
		character(:), allocatable, intent(out) :: attr_value ! the value of the attribute to be returned
		character(*), intent(in) :: attr_name ! name of the attribute
		integer(HID_T), intent(in) :: data_id ! File identifier
		integer :: hdferr ! error flag
		logical :: attr_exists ! .TRUE. if exists, .FALSE. otherwise
		integer(HID_T) :: attr_id ! Attribute identifier
		integer(HID_T) :: atype_id ! Attribute type identifier
		integer(HID_T) :: aspace_id ! Attribute space identifier
		integer, parameter :: arank = 1 ! rank of the attribute
		integer(HSIZE_T), parameter, dimension(1) :: adims = (/1/) ! dimension of attribute data
		! check if attribute exists
		call h5aexists_f(data_id, attr_name, attr_exists, hdferr)
		!
		if (attr_exists) then
			! open the 'num_models' attribute
			call h5aopen_f(data_id, attr_name, attr_id, hdferr)
			!
			! Create scalar data space for the attribute.
  			call h5screate_simple_f(arank, adims, aspace_id, hdferr)
  			!
			! Get dataspace and allocate memory for read buffer.
  			call h5aget_space_f(attr_id, aspace_id, hdferr)
  			!
  			! reading the attribute
  			call h5aread_f(attr_id, H5T_NATIVE_CHARACTER, attr_value, adims, hdferr) 
  			write(*,*) 'The character attribute is:', attr_value
  			!
  			! closing the attribute
  			call h5aclose_f(attr_id, hdferr)
		else
			! through error
			write(*,*) 'Attribute ''' // attr_name // ''' does not exist.'
		end if
		!
	end subroutine get_scalar_attribute_hdf5_char

	end module h5interface_