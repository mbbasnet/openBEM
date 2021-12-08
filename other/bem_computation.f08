	! computation of BEM 
	! currently steady state dynamic steps only
	! 
    module bem_computation
    use element_sub
    use fundamental_solution
    implicit none
	private 
	public ssd_elem
	!
	integer, parameter :: dp = kind(0.d0)
	!
	contains
	!
	! lets make subroutine for each gauss point 
	! for inplane wave propagation problem
	
    subroutine ssd_elem(Ni, v_norm, jac, coord)
	real(dp), intent(out) :: v_norm(2)
	real(dp), intent(out) :: Jac
	real(dp) :: xsi, eta
	integer :: ldim, nodes
	real(dp), intent(in) :: coord(2,3)
	real(dp),allocatable, dimension(:) :: ipt ! coordinates of integration point
	
	real(dp) :: Ni(3)
	complex(dp) :: c1, c2, mu, nu, omega
	integer :: node_dim
	real(dp) :: node(2)
	node = (/ 100.,15. /)

	node_dim = 2
	c1 = (1.0, 2.0)
	c2 = c1
	mu = (7.0d+7, 0.)
	nu = (0.3, 0.)
	omega = (10.,0.)
	
	xsi = 0.4567
	eta = 0.0
	ldim = 1
	nodes = 3

	! Getting unit normal and jacobian
	call Normal_Jac(v_norm,Jac,xsi,eta,ldim,nodes,coord)
	
	! getting shape function at the point
	call serendip_func_1D(Ni, xsi, eta, ldim, nodes)
	
	! calculating coordinate at gauss point
	if(.not. allocated(ipt)) allocate(ipt(2))
	call integration_pt(ipt,Ni,coord)
	!
	call ssd_elem_node(Ni, Jac, v_norm, ipt,node_dim, node, c1, c2, mu, nu, omega)
    write(*,*) "Coordinates of integration point: ", ipt
	if(allocated(ipt)) deallocate(ipt)
	!
	end subroutine ssd_elem
	
    subroutine ssd_elem_node(Ni, Jac, v_norm, ipt,node_dim, node, c1, c2, mu, nu, omega)
	! now lets calculate radius parameters
	! this is dependent upon node external to the element
	! inputs (integration point and coordinate of source point)
	real(dp), intent(in), dimension(:) :: ipt, node, v_norm
	real(dp), intent(in), dimension(:) :: Ni ! shape function
	integer, intent(in) :: Jac, node_dim ! Jacobian & coordinate space
	!
	! material parameters
	complex(dp) :: c1, c2, mu, nu, omega
	
	! position vectors
	real(dp) :: r, r_n! radius (positional vector)
	real(dp), allocatable, dimension(:) :: r_i

	! parameters used in fundamental solutions
	complex(dp) :: psi, kap, kapm, dkap, dpsi
	!
	! calculate position vectors
	call position_vectors(r, r_i, r_n, ipt, node_dim, node, v_norm)
	!
	! calculating psi and kappa variables
	call psi_kappa(r, omega, c1, c2, psi, kap, kapm, dkap, dpsi)
	! now calculating reciprocal traction (static and dynamic)
	!call reciprocal_traction_static(P_st, v_norm, r, r_i, r_n, nu)
	
	
	end subroutine ssd_elem_node
	
	
	end module bem_computation
