	! computational element parameters are included in this 
	! 
	module element_calc
    use element_sub
	implicit none
	private
    public tri_rect_func, tri_jac
	integer, parameter :: dp=kind(0.d0)
	!
    interface shape_fun
        module procedure lagrangian_1D_standard
        module procedure lagrangian_2D_standard
    end interface
	!
	contains
	!
    subroutine tri_rect_func(Ni, xsi, eta)
        !
        ! Computes shape functions for sub triangle to rectangle 
        ! The boundary element method with programming Eq. 6.51
        !
        real(dp), dimension(:), intent(out) :: Ni ! Array of Shape functions
        real(dp) :: xsi, eta ! intrinsic coordinate
        !
        !
        Ni(1) = 0.25*(1+xsi)*(1-eta)
        Ni(2) = 0.25*(1+xsi)*(1+eta)
        Ni(3) = 0.5*(1-xsi)
        !
        !
    end subroutine tri_rect_func
    !
    !
    !
    subroutine tri_rect_deriv(dNi, xsi, eta)
        !
        ! Computes shape functions for sub triangle to rectangle 
        ! The boundary element method with programming Eq. 6.51
        !
        real(dp), dimension(:,:), allocatable, intent(out) :: dNi ! Array of SF derivatives
        real(dp) :: xsi, eta ! intrinsic coordinate
        !
        if (.not. allocated(dNi)) allocate(dNi(3,2))
        ! taking derivatives of the functions
        ! wrt xsi
        dNi(1,1) = 0.25*(1-eta)
        dNi(2,1) = 0.25*(1+eta)
        dNi(3,1) = -0.5
        !
        ! wrt eta
        dNi(1,2) = -0.25*(1+xsi)
        dNi(2,2) = 0.25*(1+xsi)
        dNi(3,2) = 0.0
        !
        !
    end subroutine tri_rect_deriv
    !
    subroutine tri_jac(jac, xsi, eta, coord)
        !
        ! computes jacobian and unit normal
        !
        real(dp), intent(out) :: jac ! jacobian of transformation
        real(dp), intent(in) :: xsi, eta ! intrinsic coordinates
        real(dp), dimension(:,:), intent(in) :: coord ! nodal coordinates of element
        !
        real(dp), allocatable, dimension(:,:) :: dNi ! SF derivative
        real(dp) :: J11, J12, J21, J22 ! jacobian matrix entities
        !
        integer:: ii

        ! calling the SF derivatives
        call tri_rect_deriv(dNi, xsi, eta)
        !
        J11 = sum(dNi(:,1)*coord(1,:))
        J12 = sum(dNi(:,1)*coord(2,:))
        J21 = sum(dNi(:,2)*coord(1,:))
        J22 = sum(dNi(:,2)*coord(2,:))
        !
        jac = J11*J22-J12*J21
        !
        if (allocated(dNi)) deallocate(dNi)
        !
    end subroutine tri_jac
    !
    ! shape functions
    subroutine lagrangian_1D_standard(Ni, xsi, num_node)
        !
        ! computes one dimentional shape function
        ! standard means normal nominclature [-1, 1] or [-1, 0, 1]
        !
        real(dp), dimension(:), intent(out) :: Ni ! Array of Shape functions
        real(dp) :: xsi ! intrinsic coordinate
        integer, intent(in) :: num_node ! number of nodes
        !
        select case(num_node)
        case(2) ! linear element
            Ni(1) = 0.5*(1-xsi)
            Ni(2) = 0.5*(1+xsi)
            !
        case(3) ! quadratic element
            Ni(1) = 0.5*xsi*(1-xsi)
            Ni(2) = 0.5*xsi*(1+xsi)
            Ni(3) = 1.0 - xsi*xsi
        case default
            !call Error_message('Number of nodes mismatch')
            !
        end select
        !
    end subroutine lagrangian_1D_standard
    !
    !
    subroutine lagrangian_2D_standard(Ni, xsi, eta, num_node)
        !
        ! computes one dimentional shape function
        ! standard means normal nominclature [(-1, 1] or [-1, 0, 1]
        ! if number of nodes are 8, then it is serendipity element
        !
        real(dp), dimension(:), intent(out) :: Ni ! Array of Shape functions
        real(dp) :: xsi, eta ! intrinsic coordinate
        real(dp) :: mxs, pxs, met, pet, qxs, qet ! internal variables
        integer, intent(in) :: num_node ! number of nodes
        !
        select case(num_node)
        ! rectangular elements (4, 8, 9)
        case(4) ! rectangular bilinear elements
            mxs=1.0-xsi
            pxs=1.0+xsi
            met=1.0-eta
            pet=1.0+eta
            !
            Ni(1)= 0.25*mxs*met
            Ni(2)= 0.25*pxs*met
            Ni(3)= 0.25*pxs*pet
            Ni(4)= 0.25*mxs*pet
            !
        case(9) ! rectangular biquadratic lagrangian elements
            mxs=1.0-xsi
            pxs=1.0+xsi
            met=1.0-eta
            pet=1.0+eta
            qxs = 1.0-xsi*xsi
            qet = 1.0-eta*eta
            !
            Ni(1) = 0.25*xsi*mxs*eta*met
            Ni(2) = 0.25*xsi*pxs*eta*met
            Ni(3) = 0.25*xsi*pxs*eta*pet
            Ni(4) = 0.25*xsi*mxs*eta*pet
            Ni(5) = 0.5*qxs*eta*met
            Ni(6) = 0.5*qet*xsi*pxs
            Ni(7) = 0.5*qxs*eta*pet
            Ni(8) = 0.5*qet*xsi*mxs
            Ni(9) = qxs*qet
            !
        case(8) ! rectangular biquadratic serendipity elements (no middle node)
            mxs=1.0-xsi
            pxs=1.0+xsi
            met=1.0-eta
            pet=1.0+eta
            qxs = 1.0-xsi*xsi
            qet = 1.0-eta*eta
            !
            Ni(5)= 0.5*qxs*met
            Ni(6)= 0.5*qet*pxs
            Ni(7)= 0.5*qxs*pet
            Ni(8)= 0.5*qet*mxs
            !
            Ni(1)= 0.25*mxs*met-0.5*Ni(8)-0.5*Ni(5)
            Ni(2)= 0.25*pxs*met-0.5*Ni(5)-0.5*Ni(6)
            Ni(3)= 0.25*pxs*pet-0.5*Ni(6)-0.5*Ni(7)
            Ni(4)= 0.25*mxs*pet-0.5*Ni(7)-0.5*Ni(8)
            !
        case(3) ! triangular bilinear element
            Ni(1) = 1.0-xsi-eta
            Ni(2) = xsi
            Ni(3) = eta
            !
        case(6) ! triangular biquadratic element
            Ni(4) = 4.0*xsi*(1.0-xsi-eta)
            Ni(5) = 4.0*xsi*eta
            Ni(6) = 4*eta*(1-xsi-eta)
            !
            Ni(1) = 1.0-xsi-eta - 0.5*Ni(6)-0.5*Ni(4)
            Ni(2) = xsi - 0.5*Ni(4)-0.5*Ni(5)
            Ni(3) = eta - 0.5*Ni(5)-0.5*Ni(6)
            !
        case default
            !call Error_message('Number of nodes mismatch')
            !
        end select
    end subroutine lagrangian_2D_standard
    !
    !
    end module element_calc