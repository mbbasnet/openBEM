	! 
	module comp_model
	!
	use chunk_
	use mesh_
	use element_
	use model_
	use material_
	use bem_el
	use mat_ops
	use freefield_
	!
	implicit none
	private
	public Step, bemMat, initialize_step
	!
	integer, parameter :: dp = kind(0.d0)
	real, parameter :: pi = 3.1415926535897932_dp
	!
	type Step
		!
		integer :: step_num
		integer :: ndof
		complex(dp) :: omega
		!
		type(bemMat) :: gmat ! global matrix array
		type(bemMat), allocatable, dimension(:) :: cmat ! chunk matrix arrays
		type(bemMat) :: elmat ! element level matrix
		!
		contains
		! 
		procedure, pass(this) :: compute_bem => step_bem_computation
		procedure, pass(this) :: chunk_compute_bem => step_chunk_bem_computation
		!
	end type Step
	!
	contains
	!
	subroutine step_bem_computation(this, mdl)
		!
		class(Step), target, intent(inout) :: this ! this particular step
		class(Model), pointer, intent(in) :: mdl ! Discrete model over which computation is done
		!
		! 
		integer :: ii ! loop variable
		!
		!
		! allocate global matrices
		call this%gmat%mat_alloc(this%ndof, mdl%mesh%num_node, mdl%mesh%num_elcon, 1)
		!
		! allocating space for chunk matrices
		if (.not. allocated(this%cmat)) allocate(this%cmat(mdl%num_bchunk))
		!
		! start bem computation for the model
		do ii= 1, mdl%num_bchunk
			! compute bem for each chunk
			call this%chunk_compute_bem(mdl, ii)
			!
		end do 
		!write(*,*) 'Global G matrix', shape(this%gmat%Hd)
		!write(*,*) sum(this%gmat%Gd(13:,:))

		!write(*,*) 'Global H matrix', shape(this%gmat%Hd)
		!write(*,*) sum(this%gmat%Hd)
		!
	end subroutine step_bem_computation
	!
	!
	subroutine step_chunk_bem_computation(this, mdl, cn)
		!
		! computes bem for each chunk
		!
		class(Step), target, intent(inout) :: this ! this particular step
		class(Model), target, intent(in) :: mdl ! mesh within the model
		integer, intent(in) :: cn ! this particular chunk number
		!
		class(Material), pointer :: mtrl! material within the model
		class(Element), pointer :: elm ! points to the element
		!
		real(dp), allocatable :: coord(:,:) ! Nodal coordinates
		integer :: ii, ee, loc , nodeNum, indx! loop variables
		logical :: occr ! occurance of the number in the set
		!
		! element level matrices in 3d elasticity
		class(bemMat), pointer :: elmat ! element level matrix
		integer, allocatable :: nmap(:), gmap(:) ! element mapping to the global matrix
		!
		! allocation of bem matrices required for the chunk
		call this%cmat(cn)%mat_alloc(this%ndof, mdl%bchunk(cn)%num_bnode, mdl%bchunk(cn)%num_elcon, 2)
		!
		! preparing nodal coordinates for the nodes within the chunk
		coord = chunk_coord(mdl%bchunk(cn), mdl%mesh%node%coord)
		!
		! now computing BEM within chunk
		! 
		! Getting material for the chunk
		mtrl => mdl%material(mdl%bchunk(cn)%material) ! material from chunk index
		!
		!
		! Calculating free field displacement 
		do ii = 1, mdl%bchunk(cn)%num_bnode
			nodeNum = mdl%bchunk(cn)%b_node(ii)
			indx = this%ndof*nodeNum
			! Uff is directly calculated on the global matrix
			this%gmat%Uff(indx-this%ndof+1:indx) = Uff3p(Ap=1.0, theta=pi/2.0, theta_h=0.0, omega=this%omega, &
													 Cp=mtrl%c1_c, Cs=mtrl%c2_c, coord=mdl%mesh%node%coord(:, nodeNum))
		end do
		!
		!
		do ii = 1, mdl%bchunk(cn)%num_belm ! for elements inside the 
			!
			! BEM computation for each element
			!
			write(*,*) 'Element No.:', ii
			!
			if (associated(elm)) nullify(elm)
			elm => mdl%mesh%element(mdl%bchunk(cn)%b_element(ii))
			!
			!
			if (associated(elmat)) nullify(elmat)
			elmat =>this%elmat
			!
			! allocation of element level G, Hd and Hs matrices
			call elmat%mat_alloc(this%ndof, mdl%bchunk(cn)%num_bnode, elm%num_node, 3)
			!
			! now calling BEM elasticity prodecure for each element
			call bem_elasticity(elmat%Gd, elmat%Hd, elmat%Hs, elm, mdl%bchunk(cn)%b_node, coord, mtrl, this%omega, this%ndof)
			!
			!
			! computation of free field motion
			! currently assume only incident P wave
			call Tff3p(elmat%Tff, elm, mtrl, 1.0, pi/2.0, 0.0, this%omega)
			!write(*,*) 'Tff: ', ii
			!write(*,*) elmat%Tff
			!
			! we got element level G, Hd and Hs matrices
			! now update it to chunk level matrices
			!
			! node mapping for updating columns of H matrix in chunk 
			if (allocated(nmap)) deallocate(nmap)
			allocate(nmap(elm%num_node))
			!
			loc = 0
			occr = .false.
			do ee = 1, size(nmap)
				do while (occr .eqv. .false.) 
					loc = loc+1
					if (mdl%bchunk(cn)%b_node(loc) ==  elm%nodeset(ee)) occr = .true.
				end do
				nmap(ee) = loc ! locatin of the node in b_nodes of the chunk
				occr = .false.
				loc = 0
			end do ! ee
			!
			! element mapping to update element G matrix to Global G matrix :: element connectivity start to end
			!
			if (allocated(gmap)) deallocate(gmap)
			allocate(gmap(elm%num_node)) ! the size is same as number of nodes in the element
			!
			gmap = (/ (ee, ee=elm%start_con, elm%end_con) /) ! start of the elm connectivity to end of elm connectivity
			!
			!
			call matrix_mapping(this%gmat%Gd, elmat%Gd, this%ndof, gmap, mdl%bchunk(cn)%b_node)
			call matrix_mapping(this%cmat(cn)%Hd, elmat%Hd, this%ndof, nmap)
			call matrix_mapping(this%cmat(cn)%Hs, elmat%Hs, this%ndof, nmap)
			! 
			! Updating free field motion (Traction)
			call vector_mapping(this%gmat%Tff, elmat%Tff, this%ndof, gmap)
			!
		end do ! ii for each elements in the chunk
		!
		! computing c matrix for the chunk
		!if (this%step_num == 1) then ! only once in the whole model
			! using rigid body motion or smooth boudary condition
			call rigid_body_motion_c(this%cmat(cn)%c, this%cmat(cn)%Hs, this%ndof, mdl%bchunk(cn)%bound)
			!write(*,*) 'C matrix is'
			!write(*,*) this%cmat(cn)%c
		!end if
		!
		! combining dynamic and static part for complete H matrix for the chunk
		this%cmat(cn)%Hd = this%cmat(cn)%c + this%cmat(cn)%Hs + this%cmat(cn)%Hd
		!
		! now updating chunk H matrix to the global H matrix
		!
		call matrix_mapping(this%gmat%Hd, this%cmat(cn)%Hd, this%ndof, mdl%bchunk(cn)%b_node, mdl%bchunk(cn)%b_node)
		!call matrix_mapping(this%gmat%Gd, this%cmat(cn)%Gd, this%ndof, mdl%bchunk(cn)%b_node)
		!

	end subroutine step_chunk_bem_computation
	!
	!
	function chunk_coord(chunk, gnode)
		!
		! creates array of coordinates for this bem_chunk
		!
		class(BEM_chunk), intent(in) :: chunk
		real(dp), dimension(:,:), intent(in) :: gnode ! global node coordinates array
		integer :: cdim, nnode ! cartesian coordinate dimensions and num nodes
		integer :: indx, ii ! node index and loop bariables
		real(dp), allocatable :: chunk_coord(:,:) ! Nodal coordinates for chunk
		!
		cdim = size(gnode,1)
		nnode = size(chunk%b_node)
		!
		! allocation of nodal coordinate space
		if (allocated(chunk_coord)) deallocate(chunk_coord)
		allocate(chunk_coord(cdim, nnode))
		!
		!
		do ii = 1, nnode
			!
			! current index
			!
			indx = chunk%b_node(ii)
			chunk_coord(:,ii) = gnode(:,indx)
			!
		end do
		! 
		!
	end function chunk_coord
	!
	!
	subroutine initialize_step(this, mdl)
		!
		! takes info from the step class in the model
		!
		class(Step), allocatable, intent(inout) :: this(:)
		class(Model), intent(in) :: mdl
		integer :: ii
		!
		if (.not. allocated(this)) allocate(this(mdl%num_step))
		!
		do ii = 1, mdl%num_step
			if (mdl%step(ii)%stype == 1) this(ii)%ndof = 1
			if (mdl%step(ii)%stype == 2) this(ii)%ndof = 2
			if (mdl%step(ii)%stype == 3) this(ii)%ndof = 3
			!
			this(ii)%omega = mdl%step(ii)%frequency *(2.0_dp,0.0_dp)*pi
			this(ii)%step_num = ii
			!
		end do
		!
	end subroutine initialize_step
	!
	!
	end module comp_model