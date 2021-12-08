	! fundamental solutions
	module funda_el
	!
	implicit none
	!
	!
	private
	public funEl_freq
	!
	integer, parameter :: dp=kind(0.d0)
	integer, parameter :: pi=22.0_dp/7.0_dp
	!
	contains
	!
	subroutine funEl_freq(Uf, Tdf, Ts, r, dxr, Unorm, c1, c2, G, nu, omega)
		!
		! Computes required fundamental solutions for frequency domain elastodynamics
		!
		complex(dp), intent(out), dimension(:,:) :: Uf, Tdf ! fundamental Disp. and dynamic Traction
		real(dp), intent(out) :: Ts(:,:) ! Static Traction
		!
		real(dp), intent(in) :: r ! radius
		real(dp), intent(in), dimension(:) :: dxr, Unorm ! radius derivatives and unit normal
		complex(dp), intent(in) :: c1, c2, G ! wave velocities and shear modulus
		real(dp), intent(in) :: nu ! Poisson's ratio
		complex(dp), intent(in) :: omega ! Complex frequency
		!
		!
		integer :: Cdim ! cartesian dimension
		complex(dp) :: psi, kapa, dpsi, dkapa ! psi kapa functions
		!
		! finding out the dimensions
		Cdim = size(dxr)
		!
		! computing psi and kapa functions
		if (Cdim == 3) then
			call psi_kapa_3D(psi, kapa, dpsi, dkapa, r, omega, c1, c2) 
		else
			write(*,*) 'Cartesian dimension mismatch!'
		end if
		!
		! displacement fundamental solution
		Uf = Uf_eldyn(r, dxr, G, psi, kapa, Cdim)
		!
		! traction fundamental solution (Total)
		Tdf = Tf_eldyn(r, dxr, Unorm, c1, c2, psi, kapa, dpsi, dkapa, Cdim)
		!
		! Static traction fundamental solution
		Ts = Ts_elstat(r, dxr, Unorm, nu, Cdim)
		!
		! Using only dynamic part of total fundamental traction
		!Tdf = Tdf - Ts
		!
		!write(*,*) 'Motions: ', Uf(1,1), Tdf(1,1), Ts(1,1)
	end subroutine funEl_freq
	!
	!
	function Uf_eldyn(r, dxr, G, psi, kapa, Cdim)
		!
		! 2D or 3D fundamental solution for displacement
		!
		real(dp), intent(in) :: r ! 
		real(dp), intent(in) :: dxr(:) ! direction derivatives to radius dx/dr
		complex(dp), intent(in) :: G ! shear modulus (Material parameter)
		complex(dp), intent(in) :: psi, kapa ! psi and kapa functions
		integer, intent(in) :: Cdim ! cartesian dimension
		integer :: a ! alpha variable in solution
		integer :: l, k ! displacement in k direction, load in l direction
		complex(dp) :: Uf_eldyn(Cdim,Cdim) ! Fundamental displacement
		!
		complex(dp) :: mul ! temporary variable
		!
		! determinint alpha 'a' variable
		if (Cdim==2) then 
			a = 2
		else if (Cdim==3) then
			a = 4
		else 
			write(*,*) 'Dimension mismatch for the solution'
		end if
		!
		! Computation of displacement fundamental solution
		mul = a*pi*G
		! 
		do l=1, Cdim
			do k=1, Cdim
				! Dominguez equation (2.199)
				Uf_eldyn(l, k) = (psi*KRON(l,k)-kapa*dxr(l)*dxr(k))/mul
			end do
		end do
		!
	end function Uf_eldyn
	!
	!
	function Tf_eldyn(r, dxr, Unorm, c1, c2, psi, kapa, dpsi, dkapa, Cdim)
		!
		! 2D or 3D fundamental solution for displacement
		!
		real(dp), intent(in) :: r ! 
		real(dp), intent(in) :: dxr(:) ! direction derivatives to radius dx/dr
		real(dp), intent(in) :: Unorm(:) ! Unit normal to the element
		complex(dp), intent(in) :: c1, c2 ! wave velocities
		complex(dp), intent(in) :: psi, kapa, dpsi, dkapa ! psi and kapa functions
		integer, intent(in) :: Cdim ! cartesian dimension
		integer :: alpha ! alpha variable in solution
		integer :: l, k ! displacement in k direction, load in l direction
		complex(dp) :: Tf_eldyn(Cdim,Cdim) ! Fundamental displacement
		!
		complex(dp) :: A, B, C, mul ! temp functions 
		real(dp) :: drn ! radius derivative to unit normal

		!
		! determinint alpha 'a' variable
		if (Cdim==2) then 
			alpha = 2
		else if (Cdim==3) then
			alpha = 4
		else 
			write(*,*) 'Dimension mismatch for the solution'
		end if
		!
		! Computation of dr/dn 
		drn = DOT_PRODUCT(dxr, Unorm)
		!
		! computing temporary variables
		!
		A = dpsi-kapa/r
		B = 4.0_dp*kapa/r - 2.0_dp*dkapa
		C = ((c1*c1)/(c2*c2)-2.0_dp)*(dpsi-dkapa-2.0_dp*kapa/r) - 2.0_dp*psi/r
		mul = alpha*pi

		! Computation of Traction fundamental solution
		! 
		do l=1, Cdim
			do k=1, Cdim
				! Dominguez equation (2.199)
				Tf_eldyn(l, k) = A*(drn*kron(l,k)+dxr(k)*Unorm(l)) &
								 + B*(dxr(k)*dxr(l)*drn) &
								 + C*(dxr(l)*Unorm(k))
			end do
		end do
		!
	end function Tf_eldyn
	!
	!
	function Ts_elstat(r, dxr, Unorm, nu, Cdim)
		! computes static traction for elastostatics
		real(dp), intent(in) :: r ! 
		real(dp), intent(in) :: dxr(:) ! direction derivatives to radius dx/dr
		real(dp), intent(in) :: Unorm(:) ! Unit normal to the element
		real(dp), intent(in) :: nu ! Poisson's ratio
		integer, intent(in) :: Cdim ! cartesian dimension
		real(dp):: Ts_elstat(Cdim, Cdim) ! Static Traction
		!
		integer :: l, k ! displacement in k direction, load in l direction
		complex(dp) :: drn ! radius derivative to unit normal
		real(dp) :: A1, A2, A3 ! temporary variables
		!
		!
		! determinint alpha 'a' variable
		if (Cdim==2) then 
			A1 = 2.0
		else if (Cdim==3) then
			A1 = 3.0
		else 
			write(*,*) 'Dimension mismatch for the solution'
		end if
		!
		! Computation of dr/dn 
		drn = DOT_PRODUCT(dxr, Unorm)
		!
		! computing temporary variables
		!
		A2 = 1.0/(4.0*pi*(1.0-nu)*r*r) ! eqv to C2/r^2 in book
		A3 = 1.0 - 2.0*nu 



		! Computation of Traction fundamental solution
		! 
		do l=1, Cdim
			do k=1, Cdim
				! Dominguez equation (2.199)
				Ts_elstat(l, k) = A2*(A3*kron(l,k)+A1*dxr(l)*dxr(k)*drn &
								  - A3*(1.0-kron(l,k)) &
								  *(Unorm(k)*dxr(l)-Unorm(l)*dxr(k)))
			end do
		end do
		!
	end function Ts_elstat
	!
	!
	! computation of psi and kappa functions
	!
	subroutine psi_kapa_3D(psi, kapa, dpsi, dkapa, r, omega, c1, c2) 
		!
		! computes psi and kapa functions and their derivatives
		!
		complex(dp), intent(out) :: psi, kapa, dpsi, dkapa ! output functions
		real(dp) :: r ! radius
		complex(dp) :: omega ! angular frequency (Omega)
		complex(dp) :: c1, c2 ! wave velocities
		!
		complex(dp) :: i, k1, k2, e1, e2, b1, b2, vr
		!
		i = (0.0_dp, 1.0_dp)
		k1 = i*omega/c1
		k2 = i*omega/c2
		!
		e1 = EXP(-k1*r)/r
		e2 = EXP(-k2*r)/r
		!
		b1 = 1.0_dp/(k1*k1*r*r)
		b2 = 1.0_dp/(k2*k2*r*r)
		!
		vr = c2*c2/(c1*c1)
		!
		! Dominguez equation (2.200)
		psi = e2 + (b2+b2*k2*r)*e2 - vr*(1.0_dp/(k1*k1)+k1*r*b1)*e1
		!
		kapa = (3.0_dp*b2 + 3.0_dp*b2*k2*r + 1.0_dp)*e2 &
				-vr*(3.0_dp*b1 + 3.0_dp*b1*k1*r + 1.0_dp)*e1
		!
		! Dominguez equation (3.80)
		dpsi = (-2.0_dp/r-k2-3.0_dp*b2*k2-3.0_dp*b2/r)*e2 &
				+ vr*(1.0_dp/r + 3.0_dp*b1*k1 + 3.0_dp*b1/r)*e1
		!
		! Dominguez equation (3.81)
		dkapa = (-4.0_dp/r - k2 - 9.0_dp*b2*k2 - 9.0_dp*b2/r)*e2 & 
				+ vr*(4.0_dp/r + k1 + 9.0_dp*b1*k1 + 9.0_dp*b1/r)*e1
		!
		!
	end subroutine psi_kapa_3D
	!
	!
	function kron(i,j)
		! Kronecker delta function
		integer, intent(in) :: i, j
		integer :: kron 
		!
		kron = 0
		if (i==j) kron = 1
	end function kron
	!
	end module funda_el



