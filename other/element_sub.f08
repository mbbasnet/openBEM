	! calculates element shape functions, jacobian, unit normal etc
	! referenced from:
	! Beer, G., Smith, I., & Duenser, C. (2008). 
	! The boundary element method with programming: 
	! for engineers and scientists. 
	! Springer Science & Business Media.
	!
	module element_sub
	implicit none
	private
	public serendip_func_1D, Normal_jac
	public integration_pt, position_vectors
	integer, parameter :: dp=kind(0.d0)
	!
	!
	contains
	
	
	subroutine vector_norm(V, V_len)
	! ------------------
	! Normalise a vector
	! ------------------
	! 
	real(kind=dp), intent(inout) :: V(:) ! vector to be normalized
	real(kind=dp), intent(out) :: V_len ! length of vector
	!
	V_len = sqrt(sum(V*V)) ! magnitude of vector
	!
	if(abs(V_len)<1.0e-10) return
	V = V/V_len ! normalize 
	return
	end subroutine vector_norm
	
	
	!
	real(dp) function vector_ex(V1, V2)
	! ----------------------------
	! cross product of two vectors
	! ----------------------------
	!
	real(dp), intent(in) :: V1(3), V2(3)
	real(dp) :: V_cp(3)
	!
	V_cp(1) = V1(2) * V2(3) - V2(2) * V1(3)
	V_cp(2) = V1(3) * V2(1) - V2(3) * V1(1)
	V_cp(3) = V1(1) * V2(2) - V2(1) * V1(2)
	!
	return
	end function vector_ex
	
	
	subroutine serendip_func_1D(Ni, xsi, eta, ldim, nodes, inci)
	! ----------------------------------------------------------
	! computes serindipity shape functions Ni(xsi, eta)
	! for 1D and 2D (linear/parabolic) finite/Boundary elements
	! ----------------------------------------------------------
	!
	real(dp), intent(out) :: Ni(:) ! shape funcitons
	real(dp), intent(in) :: xsi, eta ! intrinsic coordinates
	integer, intent(in) :: ldim, nodes ! element dimension & num nodes
	!
	real, intent(in), optional :: inci(:) 
	! intrinsic coordinates, default is [-1,1,0] and so on
	! not included in current version. Assumes default
	
	! computation for different case
	select case(ldim)
	case(1) ! one-dimensional element
		if (nodes == 2) then ! linear elements
			Ni(1) = 0.5*(1.0-xsi)
			Ni(2) = 0.5*(1.0 +xsi)
		else if (nodes == 3) then ! quadratic elements
			Ni(1) = -xsi * 0.5*(1.0-xsi)
			Ni(2) = xsi * 0.5*(1.0 +xsi)
			Ni(3) = 1.0 - xsi * xsi ! quadratic elements
		end if
	case default ! error message
		! call Error_message (Only 1D elements are supported)
	end select
	end subroutine serendip_func_1D
	
	
		
	subroutine serendip_deriv_1D(DNi, xsi, eta, ldim, nodes, inci)
	! ----------------------------------------------------------
	! computes serindipity shape functions derivatives DNi(xsi, eta)
	! for 1D and 2D (linear/parabolic) finite/Boundary elements
	! ----------------------------------------------------------
	!
	real(dp), intent(out) :: DNi(:) ! shape funcitons
	real(dp), intent(in) :: xsi, eta ! intrinsic coordinates
	integer, intent(in) :: ldim, nodes ! element dimension & num nodes
	!
	real, intent(in), optional :: inci(:) 
	! intrinsic coordinates, default is [-1,1,0] and so on
	! not included in current version. Assumes default
	
	! computation for different case
	select case(ldim)
	case(1) ! one-dimensional element
		if (nodes==2) then ! linear elements
			DNi(1) = -0.5
			DNi(2) = 0.5
		else if (nodes ==3) then
			DNi(1) = xsi - 0.5
			DNi(2) = xsi + 0.5
			DNi(3) = -2.0 * xsi
		end if
	case default ! error message
	! call Error_message(Only 1D elements are supported)
	end select
	end subroutine serendip_deriv_1D
		
	
	
	SUBROUTINE Normal_Jac(v3,Jac,xsi,eta,ldim,nodes,coords,inci)
	!-------------------------------------------------------
	! Computes normal vector and Jacobian
	!-------------------------------------------------------
	REAL(kind=dp),INTENT(OUT) :: v3(:)
	! Vector normal to point
	REAL(kind=dp),INTENT(OUT) :: Jac
	! Jacobian
	REAL(kind=dp), INTENT(IN) :: xsi,eta
	! intrinsic coords of point
	INTEGER,INTENT(IN):: ldim
	! element dimension
	INTEGER,INTENT(IN):: nodes
	! number of nodes
	INTEGER,INTENT(IN), optional :: inci(:)
	! element incidences
	REAL(KIND=DP), INTENT(IN) :: coords(:,:)! node coordinates
	REAL(KIND=DP),ALLOCATABLE :: DNi(:)
	! Derivatives of Ni
	REAL(KIND=DP),ALLOCATABLE :: v1(:),v2(:)! Vectors in xsi,eta dir
	INTEGER :: Cdim ! Cartesian dimension
	integer :: i
	Cdim= ldim+1
	!Allocate temporary arrays
	ALLOCATE (DNi(nodes),V1(Cdim),V2(Cdim))
	!Compute derivatives of shape function
	Call Serendip_deriv_1D(DNi,xsi,eta,ldim,nodes)
	!
	! Compute vectors in xsi (eta) direction(s)
	DO I=1,Cdim
		V1(I)= DOT_PRODUCT(DNi(:),COORDS(I,:))
	END DO
	!Compute normal vector
	IF(ldim == 1) THEN
		v3(1)= V1(2)
		v3(2)= -v1(1)
	ELSE
		V3 = vector_ex(v1,v2)
	END IF
	!Normalise
	CAll Vector_norm(V3,Jac)
	DEALLOCATE (DNi,V1,V2)
	RETURN
	END SUBROUTINE Normal_Jac
	!
	! 
	subroutine integration_pt(ipt,Ni,coords)
	! calculate integration points
	! to the global coordinates
	!
	REAL(KIND=DP),intent(out) :: ipt(2)
	REAL(kind=dp),INTENT(in) :: Ni(:) ! Array with shape function
	!
	REAL(KIND=DP), INTENT(IN) :: coords(:,:)! node coordinates
	!
	! compute coordinates at integration point
	ipt(1) = DOT_PRODUCT(Ni,COORDS(1,:))
	ipt(2) = DOT_PRODUCT(Ni,COORDS(2,:))
	!
	end subroutine integration_pt
	
	subroutine position_vectors(r, r_i, r_n, ipt, node_dim, node, v_norm)
	! calculates position vectors (radius parameters)
	! r, dr/dxi, dr/dnorm
	! node is nodal coordinates, v_norm is unit normal vector
	real(dp), intent(out) :: r, r_n
	real(dp), intent(out) :: r_i(:)
	integer, intent(in) :: node_dim
	real(dp), intent(in), dimension(:) :: ipt, node, v_norm
	
	integer :: ii
	!
	do ii = 1 ,node_dim
		r_i(ii) = ipt(ii) - node(ii) ! position vector (not derivatives now)
	end do
	r = sqrt (sum(r_i*r_i)) ! radius ! r_i is ri until now
	r_i = r_i/r ! derivatives 
	r_n = sum(r_i * v_norm) ! derivative to unit normal
	!
	end subroutine position_vectors
	end module element_sub
	
	
	
	
