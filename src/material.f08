	! material used in the model
	module material_
    use hdf5
	implicit none
	private
    public Material
	integer, parameter :: dp=kind(0.d0)
	!
	type Material ! currently elastic material only
		! 
		integer :: label ! maetrial id for ease of use
		character(10) :: name ! name of the material
		!
		! general properties
        real(dp) :: RHO ! density
        !
        ! Mechanical Properties
        !
        ! Elasticity
        real(dp):: E  ! Young's Modulus
        real(dp) :: NU ! poisson's ratio
        !
        ! Damping
        real(dp):: ZETA  ! Damping factor (which damping is it??)
        !	                
        ! Lami's constants
        real(dp)::MU
        real(dp)::LAMBDA  
        !
        ! complex parameters
        complex(dp) :: c1_c, c2_c, mu_c, lambda_c! complex wave velocities and lami's constant mu
        !
        contains
        !
        procedure, pass(this) :: create_elastic => create_elastic_material
        procedure, pass(this) :: create_elastic_h5 => create_elastic_material_hdf5
        procedure, pass(this) :: complex_elastic => elastic_complexVelocity
        !
    end type Material
    !
    contains
    !
    subroutine create_elastic_material(this, Rho, E, nu, zeta)
        !
        class(Material), intent(inout) :: this ! current material
        real(dp), intent(in) :: Rho, E, nu! density, young's modulus, poisso's ratio
        real(dp), optional, intent(in) :: zeta ! damping
        !
        this%RHO = Rho
        this%E = E
        this%nu = nu
        !
        if (present(zeta)) then
            !
            this%zeta = zeta ! damping 
        else 
            !
            this%zeta = 0.0_dp ! no damping
        end if
        !
        !
        ! calculating lami's parameters
        this%MU = this%E/(2.0_dp+2.0_dp*this%nu)
        this%LAMBDA = this%MU*2.0_dp*this%nu/(1.0_dp-2.0_dp*nu)
        !
    end subroutine create_elastic_material
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
    end subroutine create_elastic_material_hdf5
    !
    !
    subroutine elastic_complexVelocity(this)
        ! Computes complex wave velocity and shear modulus
        !
        Class(Material), intent(inout) :: this
        !
        this%mu_c = this%MU +(0.0_dp,2.0_dp)*this%ZETA
        this%lambda_c = this%mu_c*2.0_dp*this%NU/(1.0_dp-2.0_dp*this%NU)
        this%c2_c = sqrt(this%mu_c/this%RHO)
        this%c1_c = sqrt((this%mu+2.0_dp*this%lambda_c)/this%RHO)
        !
    end subroutine elastic_complexVelocity
    !
    end module material_


