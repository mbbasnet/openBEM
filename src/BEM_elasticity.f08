	module bem_el
	! computes BEM for each element and available nodesets
	use element_sub
	use element_calc
	use element_
	use mesh_
	use material_
	use gaussquad_
	use funda_el
	use discrete_
	!
	implicit none
	private
	public bem_elasticity, rigid_body_motion_c, surface_area! bemMat
	!
	integer, parameter :: dp=kind(0.d0)
	!
	!
	contains
	!
	!
	subroutine bem_elasticity(G, Hd, Hs, elm, bnode, ncoord, mtrl, omega, ndof)
		!
		! computes G, H, Hs and M matrices for 3D elastodynamic problems in freq Domain
		! G - G matrix of BEM system
		! Hd - Dynamic part of H matrix in the system
		! Hs - Static part of H matrix in the system
		! M - traction to force conversion matrix for the element
		!
		class(Element), pointer, intent(in) :: elm
		integer, intent(in) :: bnode(:) ! list of boundary nodes
		real(dp), dimension(:,:), intent(in) :: ncoord ! nodal coordinates array
		Class(Material), pointer, intent(in) :: mtrl
		complex(dp), dimension(:,:), intent(out) :: G, Hd
		real(dp), dimension(:,:), intent(out) :: Hs
		integer, intent(in) :: ndof 
		!real(dp), dimension(:,:), optional, intent(inout) :: M
		!
		complex(dp), intent(in) :: omega ! Complex frequency
		!
		! material parameters
		complex(dp) :: c1, c2, mu ! wave velocities and shear modulus
		real(dp) :: nu ! Poisson's ratio
		!
		!
		real(dp), allocatable, dimension(:,:) :: elcoord
		integer :: eltype
		integer, allocatable, dimension(:) :: snode_array
		integer :: ii, jj, n1, n2
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
		nu = mtrl%nu
		! 
		! preparing list of singular nodes by the array number
		if (.not. allocated(snode_array)) allocate(snode_array(size(elm%nodeset)))
		snode_array = 0 
		!
		! now calling integration scheme
		! initially for all nodes. The singular parts will be replaced later
		! 
		call gauss_integration_3Del(G, Hd, Hs, elcoord, eltype, ncoord, 10, &
											c1, c2, mu, nu, omega, ndof)
		! 
		! the values of singular nodes is incorrect
		! thus to change the values of singular nodes lets use singular integration
		!
		!
		! listing out the singular node index from the node array
!		do ii = 1, size(elm%nodeset) ! for each nodes in the element
!			do jj = 1, size(bnode) ! for each node in the node array
				!
				! filter out the singular node
				! there will be same number of singular nodes as the element size
!				if (bnode(jj) == elm%nodeset(ii)) snode_array(ii) = jj
!			end do ! jj loop
!		end do ! ii loop
!		write(*,*) 'Singular node array', snode_array
		!
		!
!		do ii = 1, size(elm%nodeset) ! for ii th singular node in the element
!			n2 = snode_array(ii)*ndof ! respective matrix indices
!			n1 = n2-ndof+1
			!
!			call gauss_integration_3D_partition(G(n1:n2,:),elcoord, eltype, &
!									ncoord(:,ii), ii, 10, c1, c2, mu, nu, omega, ndof)
			!
!		end do ! ii loop
		!
		!
!		if (allocated(snode_array)) deallocate(snode_array)
		!
	end subroutine bem_elasticity
	!
	!
	subroutine gauss_integration_3D_partition(G, elcoord, eltype, ncoord, snode, ngauss, &
										c1, c2, mu, nu, omega, ndof)
		!
		! gauss integration for all required parameters in 3D elasticity
		!
		integer, intent(in) :: ngauss ! number of gauss points
		real(dp), intent(in) :: elcoord(:,:) ! element coordinates
		integer, intent(in) :: eltype ! the element type code
		real(dp), intent(in) :: ncoord(:) ! coordinates of single node
		integer, intent(in) :: snode ! node index in element with singularity
		!
		! matrices
		complex(dp), dimension(:,:), intent(out) :: G! G and dynamic part of H matrices
		!
		complex(dp), intent(in) :: c1, c2, mu ! wave velocities and shear modulus
		real(dp), intent(in) :: nu ! Poisson's ratio
		complex(dp), intent(in) :: omega ! Complex frequency

		integer :: gx, nn, en, e1, e2, ii ! loop variables
		integer :: Cdim, eldim ! cartesian dimensions, element level dimension
		!
		! Element functions
		real(dp), allocatable, dimension(:,:) :: gp! gauss points
		real(dp), allocatable, dimension(:) :: gw ! gauss weights
		real(dp) :: Jac, Jact ! Jacobian
		real(dp), allocatable, dimension(:) :: Ni, Un, ipt ! SF, Unit normal and integration point
		!
		! radius parameters
		real(dp) :: r ! radius 
		real(dp), allocatable :: dxr(:) ! radius derivatives
		!
		! fundamental solutions
		complex(dp), allocatable, dimension(:,:) :: Uf, Tdf 
		real(dp), allocatable, dimension(:,:) :: Ts

		! Temporary matrices
		complex(dp), dimension(:,:), allocatable :: G_temp, Hd_temp ! G and dynamic part of H matrices
		real(dp), dimension(:,:), allocatable :: Hs_temp ! Static part of H matrix

		real(dp), allocatable :: melcoord(:,:) ! modified element coordinates
		real(dp), allocatable :: lcoord(:) ! local coordinates of first transformation array of(xsi, eta)
		integer :: nnode_el ! number of nodes in the element
		integer :: ndof ! number of degrees of freedom
		!
		!
		! making subdivision of the element
		! currently only for four nodes element
		
		! partitioning the element always as scheme 1 of figure 6.14 in BEM for programming
		! rearranging the elements (for 4 node quad element)
		nnode_el = size(elcoord,2)
		if (.not. allocated(melcoord)) allocate(melcoord, mold=elcoord) ! modified element coordinates
		do ii = 1, nnode_el ! for all nodes in the element
			!
			nn = snode+ii-1
			if (nn .gt. nnode_el) then
				nn = nn - nnode_el
			end if
			melcoord(:,ii) = elcoord(:,nn)
		end do ! ii The nodes are rearranged now
		!
		! calling gauss quadrature scheme for rectangle
		Cdim = size(ncoord)
		eldim = Cdim-1
		!
		! shape function unit normal and integration point in the element
		if (.not. allocated(Ni)) allocate(Ni(nnode_el))
		if (.not. allocated(Un)) allocate(Un(Cdim))
		if (.not. allocated(ipt)) allocate(ipt(Cdim))
		!
		! radius derivatives
		if (.not. allocated(dxr)) allocate(dxr(ndof))
		!
		! fundamental solutions (allocated for all nodes) and also setting them to zero
		if (.not.allocated(Uf)) allocate(Uf(ndof, ndof)) 
		Uf = (0.0_dp, 0.0_dp)
		if (.not. allocated(Tdf)) allocate(Tdf(ndof, ndof))
		!
		! BEM matrix allocataion
		if (.not. allocated(G_temp)) allocate(G_temp(ndof, nnode_el*ndof))
		!
		G_temp = (0.0_dp, 0.0_dp)
		!
		! allocating gauss points
		if (allocated(gp)) deallocate(gp)
		if (allocated(gw)) deallocate(gw)
		if (.not. allocated(gp)) allocate(gp(eldim, ngauss**eldim))
		if (.not. allocated(gw)) allocate(gw(ngauss**eldim))
		!
		call gaussquad_rect(gp, gw, ngauss)
		!
		!
		! now working with partitioned triangles
		triang: do ii = 1, 2
		!
		gauss_coord: do gx = 1, size(gw) ! gauss coord in x direction
			! for two triangles
			!
			call tri_transform(jact, lcoord, gp(1, gx), gp(2,gx), ii)
			!
			call Serendip_func(Ni, lcoord(1), lcoord(2), 2, 4)
			call Normal_Jac(Un, Jac, lcoord(1), lcoord(2), 2, 4, melcoord)
			call integration_pt(ipt, Ni, melcoord)
			!
			!
			call radius_param(r, dxr, ncoord, ipt)
			!
			!
			! now calling fundamental solutions for 3D elasticity
			call funEl_freq(Uf, Tdf, Ts, r, dxr, Un, c1, c2, mu, nu, omega)
			!
			
			! gauss integration of the function within element
			withinElement: do en = 1, nnode_el
				!
				e1= (en-1)*ndof+1
				e2 = en*ndof
				!
				G_temp(:,e1:e2) = G_temp(:,e1:e2) + Ni(en)*Uf*jact*Jac*gw(gx)
				!
			end do withinElement
			!
		end do gauss_coord 
		end do triang ! partitioned triangles
		!
		! The obtained G and H matrices are arranged as per the modified element
		! Thus changing the arrangement of G and H matrices
		!
		do ii = 1, size(elcoord,2) ! for all nodes in the element
			!
			nn = snode+ii-1
			if (nn .gt. nnode_el) nn = nn - nnode_el
			G(:,nn:nn+ndof-1) = G_temp(:,ii:ii+ndof-1)
			!
		end do ! ii The nodes are rearranged now
		!

		!
		if (allocated(Ni)) deallocate(Ni)
		if (allocated(ipt)) deallocate(ipt)
		if (allocated(Un)) deallocate(Un)
		!
		if (allocated(dxr)) deallocate(dxr)
		! 
		if (allocated(Uf)) deallocate(Uf)
		if (allocated(Tdf)) deallocate(Tdf)
		if (allocated(Ts)) deallocate(Ts)
		!
		if (allocated(gp)) deallocate(gp)
		if (allocated(gw)) deallocate(gw)
		!
		if (allocated(G_temp)) deallocate(G_temp)
		!
		if (allocated(lcoord)) deallocate(lcoord)
		if (allocated(melcoord)) deallocate(melcoord)
		!
	end subroutine gauss_integration_3D_partition
	!
	!
	subroutine gauss_integration_3Del(G, Hd, Hs, elcoord, eltype, ncoord, ngauss, &
										c1, c2, mu, nu, omega, ndof)
		!
		! gauss integration for all required parameters in 3D elasticity
		!
		integer, intent(in) :: ngauss ! number of gauss points
		real(dp), intent(in) :: elcoord(:,:) ! element coordinates
		integer, intent(in) :: eltype ! the element type code
		real(dp), intent(in) :: ncoord(:,:) ! coordinates of all nodes
		integer, intent(in) :: ndof ! number of degrees of freedom
		!
		! matrices
		complex(dp), dimension(:,:), intent(out) :: G, Hd ! G and dynamic part of H matrices
		real(dp), dimension(:,:), intent(out) :: Hs ! Static part of H matrix
		!real(dp), optional, dimension(:,:) :: M 
		!
		complex(dp), intent(in) :: c1, c2, mu ! wave velocities and shear modulus
		real(dp), intent(in) :: nu ! Poisson's ratio
		complex(dp), intent(in) :: omega ! Complex frequency
		!
		integer :: gx, nn, en, n1, n2, e1, e2 ! loop variables
		integer :: Cdim, nnode, nel_node, eldim ! cartesian dimensions, number of nodes, number of nodes in element
		!
		! Element functions
		real(dp), allocatable, dimension(:,:) :: gp! gauss points
		real(dp), allocatable, dimension(:) :: gw ! gauss weights
		real(dp) :: Jac ! Jacobian
		real(dp), allocatable, dimension(:) :: Ni, Un, ipt ! SF, Unit normal and integration point
		!
		! radius parameters
		real(dp) :: r ! radius 
		real(dp), allocatable :: dxr(:) ! radius derivatives
		!
		! fundamental solutions
		complex(dp), allocatable, dimension(:,:) :: Uf, Tdf 
		real(dp), allocatable, dimension(:,:) :: Ts
		! 
		real(dp) :: H1 ! foe check of H
		!
		H1 = 0.0_dp
		! allocation of allocatable variables
		Cdim = size(elcoord, 1)
		eldim = Cdim - 1 ! one degree down
		nnode = size(ncoord, 2)
		nel_node = size(elcoord, 2) 
		!
		! shape function unit normal and integration point in the element
		if (.not. allocated(Ni)) allocate(Ni(nel_node))
		if (.not. allocated(Un)) allocate(Un(Cdim))
		if (.not. allocated(ipt)) allocate(ipt(Cdim))
		!
		! radius derivatives
		if (.not. allocated(dxr)) allocate(dxr(Cdim))
		!
		! fundamental solutions (allocated for all nodes) and also setting them to zero
		if (.not.allocated(Uf)) allocate(Uf(ndof*nnode, ndof))
		Uf = (0.0_dp, 0.0_dp)
		if (.not. allocated(Tdf)) allocate(Tdf(ndof*nnode, ndof))
		Tdf = (0.0_dp, 0.0_dp)
		if (.not. allocated(Ts)) allocate(Ts(ndof*nnode, ndof))
		Ts = 0.0_dp
		!
		!
		! allocating gauss points
		if (.not. allocated(gp)) allocate(gp(eldim, ngauss**eldim))
		if (.not. allocated(gw)) allocate(gw(ngauss**eldim))
		!
		! calling gauss quadrature scheme for rectangle
		call gaussquad_rect(gp, gw, ngauss)
		!
		! for each gauss points
		gauss_coord: do gx = 1, size(gw) ! gauss coord in x direction
			!
			! compute element parameters SF, Jacobian, Unit Normal, integration point coord
			! 
			call Serendip_func(Ni, gp(1, gx), gp(2,gx), 2, 4)
			call Normal_Jac(Un, Jac, gp(1,gx), gp(2,gx), 2, 4, elcoord)
			call integration_pt(ipt, Ni, elcoord)
			!write(*,*) 'We are calculating shape functions'
			!write(*,*) Ni
			!
			! now for each node 
			colloc : do nn = 1, size(ncoord,2) ! for all nodes in the array
				n2 = nn*ndof ! end node array number
				n1 = n2-ndof+1 ! start node array number
				! calculate radius parameters
				! Now prepare the functions to be integrated
				! has to be taken from fundamental solutions
				call radius_param(r, dxr, ncoord(:,nn), ipt)
				!write(*,*) 'radius is: ', r, dxr
				!
				! now calling fundamental solutions for 3D elasticity
				call funEl_freq(Uf(n1:n2,:), Tdf(n1:n2,:), Ts(n1:n2,:), &
									r, dxr, Un, c1, c2, mu, nu, omega)
				
				!
				! 
				! Performing the integration to form G and H matrices
				!el_node : do en = 1, size(elcoord,2) ! for each node in the element
			end do colloc
			!
			! At this point the fundamental solution for all nodes and each element is calculated
			! Now, it needs to be integrated with the shape function and jacobian
			within_element : do en = 1, nel_node
				! for all nodes in the element
				e2 = en*ndof
				e1 = e2-ndof+1
				! gauss integration of the function
				G(:, e1:e2) = G(:, e1:e2) + Ni(en)*Uf*Jac*gw(gx)
				Hd(:, e1:e2) = Hd(:, e1:e2) + Ni(en)*Tdf*Jac*gw(gx)
				Hs(:, e1:e2) = Hs(:, e1:e2) + Ni(en)*Ts*Jac*gw(gx)
				
				H1 = H1 +Jac*gw(gx)
	
			end do within_element
				!
		end do gauss_coord
		!
		if (allocated(Ni)) deallocate(Ni)
		if (allocated(ipt)) deallocate(ipt)
		if (allocated(Un)) deallocate(Un)
		!
		if (allocated(dxr)) deallocate(dxr)
		! 
		if (allocated(Uf)) deallocate(Uf)
		if (allocated(Tdf)) deallocate(Tdf)
		if (allocated(Ts)) deallocate(Ts)
		!
		if (allocated(gp)) deallocate(gp)
		if (allocated(gw)) deallocate(gw)
		!
	end subroutine gauss_integration_3Del
	!
	!
	subroutine radius_param(r, dxr, p0, p1)
		!
		real(dp), intent(out) :: r ! radius
		real(dp), intent(out) :: dxr(:) ! radius derivatives
		!
		real(dp), intent(in), dimension(:) :: p0, p1 ! point coordinates
		!
		! calculating the differences
		dxr = p1 - p0
		r = sqrt(sum(dxr*dxr))
		dxr = dxr/r ! dx/dr
		!
	end subroutine radius_param
	!
	!
	subroutine gaussquad_rect(gx, gw, ngauss)
		! gives coordinates of gauss point and weights for 2D rectangle [(-1,-1),(1,1)]
		integer, intent(in) :: ngauss
		real(dp), dimension(2, ngauss*ngauss), intent(out) :: gx ! rectangular gauss coordinates
		real(dp), dimension(ngauss*ngauss), intent(out) :: gw ! gauss weight
		!
		real(dp), allocatable, dimension(:) :: xsi, wt ! gauss coords and weights in 1D
		integer :: ii, jj, count ! loop variables
		!
		!
		if (allocated(xsi)) deallocate(xsi)
		if (allocated(wt)) deallocate(wt)
		allocate(xsi(ngauss))
		allocate(wt(ngauss))
		!
		! now getting gauss point and gauss weights
		call cpquad(n=ngauss,alpha=1.0_dp,name="Legendre",w=wt,x=xsi)
		!
		count = 1 ! initializing the loop
		axis1 : do ii = 1, ngauss
		axis2 : do jj = 1, ngauss
			! filling two dimensional gauss points
			gx(1, count) = xsi(ii)
			gx(2, count) = xsi(jj)
			gw(count)= wt(ii)*wt(jj)
			count = count + 1
			
		end do axis2
		end do axis1
		!
		if (allocated(xsi)) deallocate(xsi)
		if (allocated(wt)) deallocate(wt)
		!
	end subroutine gaussquad_rect
	!
	!
	subroutine rigid_body_motion_c(c, Hs, ndof, dt)
		!
		! combines Hd and Hs matrix, depending upon the type of domain
		real(dp), allocatable, dimension(:,:), intent(out) :: c ! dynamic part of H matrix
		real(dp), dimension(:,:), intent(in) :: Hs ! static part of H matrix
		integer :: ndof ! number of degrees of freedom
		integer, intent(in) :: dt 	! 0= finite, 1 = semiinfiinte, 2 = infinite
									! 3= constant elements
		integer :: ii, jj, kk ! loop variables
		!
		if (.not. allocated(c)) allocate(c, mold=Hs)
		c = 0.0_dp
		!
		!
		select case(dt)
		case (0:1)
			do ii = 1, size(c,1) ! a square matrix
				c(ii, ii) = 0.5*dt 
			end do
			!
			do ii = 1, size(Hs, 1), ndof
				!
				do jj = 1, ndof
					do kk = 1, size(Hs, 2), ndof
						c(ii+jj-1,ii+jj-1) = c(ii+jj-1,ii+jj-1) &
						                        -Hs(ii+jj-1, kk+jj-1) 
					end do

				end do
				!
			end do
		case(3)
			do ii = 1, size(c,1) ! a square matrix
				c(ii, ii) = 0.5
			end do
		end select
		!
	end subroutine rigid_body_motion_c
	!
	subroutine tri_transform(jac, lcoord, xsit, etat, tr)
		!
		! gives jacobian of transformation from rectangle to triangle
		! also gives coordinates of rectangle
		!
		real(dp), intent(out) :: jac ! jacobian of the transformation
		real(dp), allocatable, intent(out) :: lcoord(:) ! intrinsic coordinates for the rectangle
		real(dp), intent(in) :: xsit, etat ! intrinsic coordinates for the triangle
		integer, intent(in) :: tr ! index of the triangle
		real(dp) :: elcoord(2,3) ! the coordinates of the triangle element
		real(dp), allocatable :: Ni(:), dNi(:,:), Un(:)	
		!
		! assigning intrinsic coordinates of the nodes in triangle
		select case (tr) ! index of the triangle (1 or 2)
		case(1) ! lower (first triangle)
			! nodes 2, 3, 1
			elcoord = reshape((/1.0, -1.0, 1.0, 1.0, -1.0, -1.0/), shape=(/2,3/))
		 	
		case(2)
			! nodes 3, 4, 1
			elcoord = reshape((/1.0, 1.0, -1.0, 1.0, -1.0, -1.0/), shape=(/2,3/))
		case default
			! error message
		end select
		!
		if (.not. allocated(Ni)) allocate(Ni(3))
		if (.not. allocated(dNi)) allocate(dNi(3,2))
		!
		if (.not. allocated(lcoord)) allocate(lcoord(2))
		!
		! computing Shape function and jacobian
		call tri_rect_func(Ni, xsit, etat) ! shape function 
		call tri_jac(jac, xsit, etat, elcoord)
		!
		lcoord(1) = sum(Ni(:)*elcoord(1,:)) ! local coordinates
		lcoord(2) = sum(Ni(:)*elcoord(2,:)) ! local coordinates
		!
		!
		if (allocated(Ni)) deallocate(Ni)
		if (allocated(dNi)) deallocate(dNi)
		!
	end subroutine tri_transform
	!
	!
	! example calculating surface area of the element
	subroutine surface_area(SF, ngauss)
		!
		! Element functions
		integer, intent(in) :: ngauss
		real(dp), allocatable, dimension(:,:) :: gp! gauss points
		real(dp), allocatable, dimension(:) :: gw ! gauss weights
		real(dp) :: Jac ! Jacobian
		real(dp), allocatable, dimension(:) :: Ni, Un, ipt ! SF, Unit normal and integration point
		!
		real(dp) :: elcoord(3,4)! coordinates of the elements
		real(dp), intent(out) :: SF
		integer :: gx
		!
		elcoord = reshape((/ -1.0, -1.0, 0.0, 1.0, -1.0, 0.0, 1.0, 1.0, 0.0, -1.0, 1.0,  0.0/),shape=(/ 3, 4 /))
		!
		
		if (.not. allocated(gp)) allocate(gp(2, ngauss*ngauss))
	
		if (.not. allocated(gw)) allocate(gw(ngauss*ngauss))
		
		! calling gauss quadrature scheme for rectangle
		call gaussquad_rect(gp, gw, ngauss)
		!
		if (.not. allocated(Un)) allocate(Un(2))
		gauss_coord: do gx = 1, size(gw) ! gauss coord in x direction
			!
			! compute element parameters SF, Jacobian, Unit Normal, integration point coord
			! 
			call Normal_Jac(Un, Jac, gp(1,gx), gp(2,gx), 2, 4, elcoord)

			!call Jac_normal(Jac, Un, gp(1,gx), gp(2,gx), elcoord)

			SF= SF + jac*gw(gx)
			write(*,*) 'jac:', jac
		end do gauss_coord
		write(*,*) 'Surface:', SF
		deallocate(gp, gw)
		!
	end subroutine surface_area
	!
	!
	end module bem_el