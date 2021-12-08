	! discrete forms like shape function jacobian etc
	module discrete_
	!
	use element_sub
	implicit none
	!
	private
	public Jac_normal
	!
	integer, private, parameter :: dp=kind(0.d0)
	!
	!
	contains
	!
	!
	subroutine Jac_normal(Jac, Un, xi, eta, coord)
		!
		! Calculates Jacobian and unit normal using the derivatives of shape function
		!
		real(dp), intent(out) :: Jac ! jacobian of the transformation
		real(dp), intent(out) :: Un(3) ! unit normal vector
		real(dp), intent(in) :: xi, eta ! intrinsic coordinates
		real(dp), intent(in) :: coord(:,:) ! cartesian nodal coordinate of the element
		!
		real(dp) :: dNi(2,4) ! derivative of the shape funciton
		real(dp) :: dNx(2,3) ! dot product of coordinates and gradient shape functionos
		!
		
		dNi = dSF(xi, eta)
		!
		dNx = matmul(dNi, coord)
		write(*,*) 'dNx:'
		write(*,*) dNx
		!
		Un = vector_ex(dNx(1,:), dNx(2,:))
		write(*,*) 'V3:'
		write(*,*) Un
		!
		call vector_norm(Un, Jac)
		!
	end subroutine Jac_normal
	!
	!
	function SF(xi, eta)
		! returns the shape function for given element type
		! initially only one type of shape function is given
		!
		real(dp), intent(in) :: xi, eta ! intrinsic coordinates
		real(dp) :: SF(4)
		!
		! shape function for four node quadrilateral ([-1,-1], [1,-1], [1,1], [-1,1])
		SF(1) = 0.25*(xi -1.0)*(eta -1.0)
		SF(2) = -0.25*(xi +1.0)*(eta -1.0)
		SF(3) = 0.25*(xi +1.0)*(eta +1.0)
		SF(4) = -0.25*(xi -1.0)*(eta +1.0)

	end function SF
	!
	!
	function dSF(xi, eta)
		! returns the derivatives of the shape function for given element type
		! initially only for one type of element
		!
		real(dp), intent(in) :: xi, eta ! intrinsic coordinates
		real(dp) :: dSF(2,4)
		!
		! shape function for four node quadrilateral ([-1,-1], [1,-1], [1,1], [-1,1])
		!
		dSF(1,1) = 0.25*(eta -1.0)
		dSF(1,2) = -0.25*(eta -1.0)
		dSF(1,3) = 0.25*(eta +1.0)
		dSF(1,4) = -0.25*(eta +1.0)
		!
		dSF(2,1) = 0.25*(xi -1.0)
		dSF(2,2) = -0.25*(xi +1.0)
		dSF(2,3) = 0.25*(xi +1.0)
		dSF(2,4) = -0.25*(xi -1.0)
		!
		write(*,*) 'The derivatives of SF are:'
		write(*,*) dSF
		!
	end function dSF
	!
	!
	end module discrete_
