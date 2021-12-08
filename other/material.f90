	module material_
	! defines material and 
	! computes material properties
	!
	implicit none
	private
	integer, parameter:: dp=kind(0.d0)
	!
	!
	! -----------------------
	! definition of material
	! -----------------------
	type, public :: material
        !
        ! metadata
        integer :: modelid ! id of containing model
        integer :: id ! id number of the material
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
        contains
        !
        procedure, pass(this) :: read_data => read_material_data
        !
	!
	end type material 
	!
	!
	! ------------------
	! module procedures
	! ------------------
	!
	contains
	!
	!
	!
	subroutine read_material_data(this,cmid,mtrlid)
	! read material data from database
	! /////////currently defined locally //////
	!
		class(material), intent(INOUT) :: this
		integer, intent(IN) :: cmid ! model id
		integer, intent(IN) :: mtrlid ! material id
		! records data from respective material index
		!////temporary variables ////
		!
		this%RHO  = 2000.
		this%E = 2.E+10
		this%NU = 0.33333
		this%ZETA = 0.00
		!
	!
	end subroutine read_material_data
	!
	!
	end module material_

