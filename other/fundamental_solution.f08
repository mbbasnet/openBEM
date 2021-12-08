	!
    module fundamental_solution
	!
	! fundamental solution 
	! for 1. 2D plane wave propagation 
	! in frequency domain
	!
    use special_functions
    implicit none
	private
	public psi_kappa, reciprocal_displacement, reciprocal_traction, reciprocal_traction_static
	integer, parameter :: dp = kind(0.d0)
	real(kind=dp), parameter :: pi = 22._dp/7._dp ! Pi
	real(kind=dp), parameter :: eug= pi/(2._dp*exp(1._dp)) ! euler gamma number
	! --------------------------------------
	! Modified bessel function of second kind
	! for small arguments
	! ---------------------------------------
	!
	contains
	!
    subroutine cik01_small_arg(z, cbk0m,cbk1m,cbk0,cbk1)
	!
	! z : complex argument
	!  k0m, k1m : modified bessel function 2nd kind (small argument)	
	!  k0, km : modified bessel function 2nd kind
	! 
	! Discussion:
	! 
	!  This module used the subroutine for modified bessel functions
	!  and calculates the modified bessel function of second kind for 
	!  small arguments
	!
	complex(kind=dp), intent(IN) :: z, cbk0, cbk1
	complex(kind=dp), intent(OUT) :: cbk0m, cbk1m
	!
	!
		cbk0m = cbk0 - cdlog(z/(2.0_dp,0.0_dp)) + eug
		cbk1m = cbk1 - (1.0_dp,0.0_dp)/z - z/(2.0_dp,0.0_dp)*(cdlog(z/(1.0_dp,0.0_dp)) &
		             + eug -(1.0_dp, 0.0_dp)/(2.0_dp, 0.0_dp))
	return
	end subroutine cik01_small_arg
	!
	!
	!
	! ---------------------------------------
	! computation of psi and kappa parameters
	! ---------------------------------------
	!
 	subroutine psi_kappa(r, omega, c1, c2, psi, kap, kapm, dkap, dpsi)
	!
	! computes psi and kappa functions
	! required for 2D plane wave propagation problems
	!
	! real r : radius (distance between integratio pt. and collocation pt.)
	! complex omega : complex angular frequency
	! complex c1 : P-wave velocity
	! complex c2 : S-wave velocity
	! kap : KAPPA function
	! psi : PSI function
	! kapm : KAPPA function for small arguments
	! dkap : derivative for kappa for small arg
	! dpsi : derivative of psi for small arg
	!
	real(kind=dp), intent(IN) :: r
	complex(kind=dp), intent(IN) :: omega, c1, c2
	complex(kind=dp), intent(OUT) :: psi,kap,kapm,dkap,dpsi
	!
	complex(kind=dp) :: z1 ! argument
	complex(kind=dp) :: k0z1 ! modified bessel function k0(z1)
	complex(kind=dp) :: k1z1 ! modified bessel function k1(z1)
	complex(kind=dp) :: k2z1 ! modified bessel function k2(z1)
	complex(kind=dp) :: k0mz1 ! k0(z1) for small argument
	complex(kind=dp) :: k1mz1 ! k1(z1) for small argument
	complex(kind=dp) :: k2mz1 ! k2(z1) for small argument
	!
	complex(kind=dp) :: z2 ! argument
	complex(kind=dp) :: k0z2 ! modified bessel function k0(z2)
	complex(kind=dp) :: k1z2 ! modified bessel function k1(z2)
	complex(kind=dp) :: k2z2 ! modified bessel function k2(z2)
	complex(kind=dp) :: k0mz2 ! k0(z2) for small argument
	complex(kind=dp) :: k1mz2 ! k1(z2) for small argument
	complex(kind=dp) :: k2mz2 ! k2(z2) for small argument
	!
	complex(kind=dp) :: AA,BB ! additional terms for dpsi and dkap
	!
	complex(kind=dp) :: cbi0,cdi0,cbi1,cdi1,cdk0,cdk1 ! unused (bessel first kind)
	!
	!
	z1 = (0.0_dp, 1.0_dp)*omega*r/c1
	z2 = (0.0_dp, 1.0_dp)*omega*r/c2
	!
	! call modified bessel function for normal argument (z1)
	call cik01(z1,cbi0,cdi0,cbi1,cdi1,k0z1,cdk0,k1z1,cdk1)
	! cbi0,cdi0,cbi1,cdi1,cdk0 and cdk1 are not used here 
	!
	! call modified bessel function for small argument (z1)
	call cik01_small_arg(z1,k0mz1,k1mz1,k0z1,k1z1)
	!
	!
	! call modified bessel function for normal argument (z2)
	call cik01(z2,cbi0,cdi0,cbi1,cdi1,k0z2,cdk0,k1z2,cdk1)
	! cbi0,cdi0,cbi1,cdi1,cdk0 and cdk1 are not used here 
	!
	! call modified bessel function for small argument (z2)
	call cik01_small_arg(z2,k0mz2,k1mz2,k0z2,k1z2)
	!	
	! calculate modified bessel funciton
	! second kind second order
	!
	k2z1 = k0z1+(2.0_dp,1.0_dp)*k1z1/z1
	k2mz1 = k0mz1+(2.0_dp,1.0_dp)*k1mz1/z1
	k2z2 = k0z2+(2.0_dp,1.0_dp)*k1z2/z1
	k2mz2 = k0mz2+(2.0_dp,1.0_dp)*k1mz2/z1
	!
	!
	! -------------------------------
	! calculate parameter PSI & KAPPA
	! --------------------------------
	!
	! for normal arguments
	!
	! Dominguez equation 2.209 [page 135]
	psi = k0z2 + (k1z2-c2/c1*k1z1)/z2
	kap = k2z2 - (c2/c1)**2*k2z1 
	!
	!
	! for small arguments
	!
	!psim = k0mz2 + (k1mz2-c2/c1*k1mz1)/z2
	kapm = k0mz2 +(2.0_dp,0.0_dp)*k1mz2/z2 &
	    -(c2/c1)**2*(k0mz2 + (2.0_dp,0.0_dp)*k1mz1/z1) 
	!
	! derivatives of PSI and KAPPA for small agruments
	! with addition of significant parts AA and BB
	! Dominguez equation 2.308
	AA = z2**2*(cdlog(z2/2.0_dp) + eug - 0.5_dp)
	BB = (c2/c1)**2*z1**2/(2.0_dp*r)*(cdlog(z1/2.0_dp) + eug -0.5_dp)
	!
	!
	! Dominguez equation 2.266 and 2.267 [page 152] 
	! with addition of significant part for larger arguments
	!
	dpsi = (0.0_dp,1.0_dp)*omega*k1mz2*&
	       ((2.0_dp,0.0_dp)*c2/(omega*r)**2-(1.0_dp,0.0_dp))&
	       + k1mz1*((2.0_dp,0.0_dp)*c2**2/(c1*(0.0_dp,1.0_dp)*omega))&
	       - k0mz2/r + k0mz1*(c2/c1)**2/r-AA
	!
	dkap = -k1z2*(0.0_dp,1.0_dp)*omega/c2-(2.0_dp,0.0_dp)*k2mz2 &
	        + (c2/c1)**2*((0.0_dp,1.0_dp)*omega*k1mz1/c1 &
	        + (2.0_dp,0.0_dp)*k2z1/r)-AA+BB 
	!
	!
	end subroutine psi_kappa
	!
	!
	
    subroutine reciprocal_displacement(U, r_i, mu, psi, kappa)
    ! calculation of reciprocal traction and displacement
	! calculates reciprocal displacement
	complex(dp), intent(out), dimension(:,:) :: U ! reciprocal displacement
	real(dp), intent(in), dimension(:) :: r_i ! radius derivatives
	complex(dp), intent(in) :: mu, psi, kappa ! shear modulus
	
	U(1,1) = (psi - kappa * r_i(1)*r_i(1))/(6.2821852 * mu)
	U(1,2) = (- kappa * r_i(1)*r_i(2))/(6.2821852 * mu)
	U(2,1) = U(1,2)
	U(2,2) = (psi - kappa * r_i(2)*r_i(2))/(6.2821852 * mu)
	!
	end subroutine reciprocal_displacement
	

    subroutine reciprocal_traction_static(P_st, v_norm, r, r_i, r_n, nu)
	!
	complex(dp), intent(out), dimension(:,:) :: P_st ! reciprocal Tr. static
	real(dp), intent(in), dimension(:) :: r_i, v_norm ! r der and unit norm
	real(dp), intent(in) :: r, r_n ! radius and derivitave to normal
	complex(dp), intent(in) :: nu ! poisson's ratio
	!
	! temporary variables
	complex(dp) :: S1, S2
	S2 = -1./(4.*pi*(1.-nu)*r)
    S1 = 1.-2.*nu
        
    P_st(1,1) =  (S2 * r_n * (S1 + 2. * r_i(1) *r_i(1)))
    P_st(1,2) = S2*(r_n*(2.*r_i(1)*r_i(2)) &
						+ S1*(v_norm(1)*r_i(2)-v_norm(2)*r_i(1)))
	P_st(2,1) = S2*(r_n*(2.*r_i(2)*r_i(1)) &
						+ S1*(v_norm(2)*r_i(1)-v_norm(1)*r_i(2)))
	P_st(2,2) =  (S2 * r_n * (S1 + 2. * r_i(2) *r_i(2)))
	
	end subroutine reciprocal_traction_static
	
	
    subroutine reciprocal_traction(P, v_norm, r, r_i, r_n, nu, c1, c2, kap, dkap, dpsi)
	!
	complex(dp), intent(out), dimension(:,:) :: P! reciprocal Traction
	real(dp), intent(in), dimension(:) :: r_i, v_norm ! r der and unit norm
	real(dp), intent(in) :: r, r_n ! radius and derivitave to normal
	complex(dp), intent(in) :: c1, c2 ! poisson's ratio & velocities
	complex(dp), intent(in) :: nu ! poisson's ratio
	complex(dp), intent(in) :: kap, dkap, dpsi
	!
	! Temporary variables
	real(dp) :: T0
	complex(dp) :: T1, T2, T3, T4, sq_velRatio
	!
	T0 = 1./(2.*pi)
    T1 = dpsi - kap/r 
	T2 = (2./r)*kap
 	T3 = 2*dkap*r_n
    T4 = dpsi - dkap - 2.*kap/(2.*r)
    sq_velRatio = (c1/c2)*(c1/c2) ! ration P_vel to SV_vel
                
    P(1,1) = T0*(T1*(r_n+r_i(1)*v_norm(1)) &
                -T2*(v_norm(1)*r_i(1)-2*r_i(1)*r_i(1)*r_n) &
                -T3*r_i(1)*r_i(1)+((sq_velRatio)-2) &
                *T4*r_i(1)*v_norm(1))
        
    P(1,2) = T0*(T1*(r_i(2)*v_norm(1)) &
                -T2*(v_norm(2)*r_i(1)-2*r_i(1)*r_i(2)*r_n) &
                -T3*r_i(1)*r_i(2)+((sq_velRatio)-2) &
                *T4*r_i(1)*v_norm(2))
        
    P(2,1) = T0*(T1*(r_i(1)*v_norm(2)) &
                -T2*(v_norm(1)*r_i(2)-2*r_i(2)*r_i(1)*r_n) &
                -T3*r_i(2)*r_i(1)+((sq_velRatio)-2) &
                *T4*r_i(2)*v_norm(1))
        
    P(2,2) = T0 *(T1*(r_n+r_i(2)*v_norm(2)) &
                -T2*(v_norm(2)*r_i(2)-2*r_i(2)*r_i(2)*r_n) &
                -T3*r_i(2)*r_i(2)+((sq_velRatio)-2) &
                *T4*r_i(2)*v_norm(2))
	
	end subroutine reciprocal_traction
	
	
	end module fundamental_solution

