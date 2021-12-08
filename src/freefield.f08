	! free field motion 
	module freefield_
	!
	use element_
	use material_
	use element_sub
	!
	implicit none
	private
	public Uff3p, Tff3p
	!
	integer, parameter :: dp=kind(0.d0)
	!
	contains
	!
	subroutine Tff3p(Tff, elm, mtrl, Ap, theta, theta_h, omega)
		!
		! calculates 3D free field motion (displacement and traction)
		!
		Class(Element), pointer, intent(in) :: elm ! element
		Class(Material), pointer, intent(in) :: mtrl
		complex(dp), dimension(:), intent(out) :: Tff
		!
		! wave parameters
		real, intent(in) :: Ap ! Amplitude of the wave
		real, intent(in) :: theta, theta_h ! respective angles of the wave propagation
		complex(dp), intent(in) :: omega ! Complex frequency
		!
		! material parameters
		complex(dp) :: c1, c2, mu, lambda! wave velocities and lami constants
		real(dp), allocatable, dimension(:) :: Un !Unit normal 
		real(dp) :: Jac ! jacobian of the transformation
		!
		!
		real(dp), allocatable, dimension(:,:) :: elcoord
		integer :: eltype
		real(dp) :: inci(2,4) ! coordinates of transformed element
		integer :: ii 
		!
		! extracting element type and element coordinates
		! 
		elcoord = elm%coord ! element coordinates
		eltype = elm%eltype ! element type index
		!
		! material parameters
		c1 = mtrl%c1_c
		c2 = mtrl%c2_c
		mu = mtrl%mu_c
		lambda = mtrl%lambda_c
		!
		! for rectangular element
		inci = reshape((/ -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0 /), shape=(/ 2,4 /))
		!
		!
		do ii = 1, size(elcoord,2) ! for all nodes in the element
			!
			call Normal_Jac(Un, Jac, inci(1,ii), inci(2,ii), 2, 4, elm%coord)
			!
			Tff(3*(ii-1)+1:3*ii) = Tff3ps(Ap, theta, theta_h, omega, c1, c2, lambda, mu, elcoord(:,ii), Un)
			!write(*,*) 'elsize', size(elcoord,2), elcoord(:,ii)
		end do 
		!
	end subroutine Tff3p
	!
	!
	function Uff3p(Ap, theta, theta_h, omega, Cp, Cs, coord)
		!
		! Ap : amplitude of the incident wave
		! theta : angle made by the wave to half plane (0: horizontal right to left, 90: vertical upwards)
		! theta_h : angle made by the wave to XZ plane in clockwise direction (0: X to -X direction, 90: Y to -Y direction)
		! omega: angular frequency of the wave
		! Cp, Cs : primary and secondary wave velocities of the medium
		!
		!
		complex(dp) :: Uff3p(3) ! Free field motion due to P wave
		real, intent(in) :: Ap ! Amplitude of the wave
		real, intent(in) :: theta, theta_h ! respective angles of the wave propagation
		complex(dp), intent(in) :: omega ! angular frequency
		complex(dp), intent(in) :: Cp, Cs ! complex wave velocity of the medium
		real(dp), intent(in) :: coord(3) ! coordinates of the point
		!
		! temporary variables
		complex(dp) :: k, kp, ks, dk, v, v1, temp1, temp2, RPP, RPS, x, z
		real(dp) :: T(3,3) ! coordinate transformation matrix
		!
		!
		kp = omega/Cp
		ks = omega/Cs
		k = kp*cos(theta)

		v = (0.0, 1.0)*kp*sin(theta)
		v1 = sqrt(k*k - ks*ks)
		temp1 = (2.0*k*k-ks*ks)**2 
		temp2 = 4.0*k*k*v*v1
		dk = temp1 - temp2
		RPP = (-temp1 - temp2)/dk
		RPS = -(0.0, 4.0) *k*v*temp1/dk
		x = coord(1)
		z = coord(3)
		!
		! free field displacement
		Uff3p(1) = Ap*((k/kp)*exp(-v*z - (0.0, 1.0)*k*x) &
				  +RPP*(k/kp)*exp(v*z - (0.0, 1.0)*k*x) &
				  +RPS*((0.0, 1.0)*v1/kp)*exp(v1*z - (0.0, 1.0)*k*x))

		Uff3p(2) = 0.0_dp

		Uff3p(3) = Ap*((-(0.0,1.0)*v/kp)*exp(-v*z - (0.0, 1.0)*k*x) &
				  +RPP*((0.0, 1.0)*v/kp)*exp(v*z - (0.0, 1.0)*k*x) &
				  +RPS*(-k/kp)*exp(v1*z - (0.0, 1.0)*k*x))
		!
		! Transformation of axes
		T = reshape((/ cos(theta_h), -sin(theta_h), &
			0.0, sin(theta_h), cos(theta_h), 0.0, &
			0.0, 0.0, 1.0 /), shape=(/ 3, 3/))
		!
		Uff3p = MATMUL(T, Uff3p)
		!
	end function Uff3p
	!
	function Tff3ps(Ap, theta, theta_h, omega, Cp, Cs, lambda, mu, coord, u_norm)
		!
		! Ap : amplitude of the incident wave
		! theta : angle made by the wave to half plane (0: horizontal right to left, 90: vertical upwards)
		! theta_h : angle made by the wave to XZ plane in clockwise direction (0: X to -X direction, 90: Y to -Y direction)
		! omega: angular frequency of the wave
		! Cp, Cs : primary and secondary wave velocities of the medium
		!
		!
		complex(dp) :: Tff3ps(3) ! Free field motion due to P wave
		real, intent(in) :: Ap ! Amplitude of the wave
		real, intent(in) :: theta, theta_h ! respective angles of the wave propagation
		complex(dp), intent(in) :: omega ! angular frequency
		complex(dp), intent(in) :: Cp, Cs ! complex wave velocity of the medium
		complex(dp), intent(in) :: lambda, mu ! lami constants
		real(dp), intent(in) :: coord(3) ! coordinates of the point
		real(dp), intent(in) :: u_norm(3) ! unit normal of the element at the point
		!
		! temporary variables
		complex(dp) :: k, kp, ks, dk, v, v1, temp1, temp2, RPP, RPS, x, z, Kai
		complex(dp) :: Ux_x, Ux_z, Uz_x, Uz_z
		real(dp) :: T(3,3)
		!
		kp = omega/Cp
		ks = omega/Cs
		k = kp*cos(theta)

		v = (0.0, 1.0)*kp*sin(theta)
		v1 = sqrt(k*k - ks*ks)
		temp1 = (2.0*k*k-ks*ks)**2 
		temp2 = 4.0*k*k*v*v1
		dk = temp1 - temp2
		RPP = (-temp1 - temp2)/dk
		RPS = -(0.0, 4.0) *k*v*temp1/dk
		x = coord(1)
		z = coord(3)
		!
		Ux_x = Ap*(0.0, -1.0)*k*((k/kp)*exp(-v*z - (0.0, 1.0)*k*x) &
				  +RPP*(k/kp)*exp(v*z - (0.0, 1.0)*k*x) &
				  +RPS*((0.0, 1.0)*v1/kp)*exp(v1*z - (0.0, 1.0)*k*x))


		Uz_x = Ap*(0.0, -1.0)*k*((-(0.0,1.0)*v/kp)*exp(-v*z - (0.0, 1.0)*k*x) &
				  +RPP*((0.0, 1.0)*v/kp)*exp(v*z - (0.0, 1.0)*k*x) &
				  +RPS*(-k/kp)*exp(v1*z - (0.0, 1.0)*k*x))

		!temp
		Ux_z = Ap*(-v*(k/kp)*exp(-v*z - (0.0, 1.0)*k*x) &
				  +RPP*v*(k/kp)*exp(v*z - (0.0, 1.0)*k*x) &
				  +RPS*v1*((0.0, 1.0)*v1/kp)*exp(v1*z - (0.0, 1.0)*k*x))


		Uz_z = Ap*(((0.0,1.0)*v*v/kp)*exp(-v*z - (0.0, 1.0)*k*x) &
				  +RPP*v*((0.0, 1.0)*v/kp)*exp(v*z - (0.0, 1.0)*k*x) &
				  +RPS*v1*(-k/kp)*exp(v1*z - (0.0, 1.0)*k*x))
		!
		Kai = lambda+2*mu
		!
		! as Uy = 0 and constant expressions in terms of y all U with y or wrt y are zero
		!
		Tff3ps(1) = (Kai*Ux_x +  lambda*Uz_z) * u_norm(1) &
				+ mu*(Ux_z + Uz_x) * u_norm(3)

		Tff3ps(2) = lambda*(Ux_x + Uz_z)*u_norm(2)

		Tff3ps(3) = mu*(Ux_z + Uz_x) * u_norm(1) &
					+ (Kai*Ux_x + lambda*Uz_z) * u_norm(3)
		!
		! Transformation of axes
		T = reshape((/ cos(theta_h), -sin(theta_h), &
			0.0, sin(theta_h), cos(theta_h), 0.0, &
			0.0, 0.0, 1.0 /), shape=(/ 3, 3/))
		!
		Tff3ps = MATMUL(T, Tff3ps)
		!
	end function Tff3ps
	!
	end module freefield_

