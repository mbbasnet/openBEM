	! BEM matrix operations
	module mat_ops
	!
	implicit none 
	private
	public bemMat
	public matrix_mapping, vector_mapping
	!
	integer, parameter :: dp = kind(0.d0)
	!
	!
	interface matrix_mapping
		module procedure matrix_mapping_real
		module procedure matrix_mapping_complex
	end interface matrix_mapping
	!
	!
	interface vector_mapping
		module procedure vector_mapping_real
		module procedure vector_mapping_complex
	end interface vector_mapping
	!
	!
	type bemMat
		! List of unallocated matrices, that may have to be used for bem computation
		complex(dp), dimension(:,:), allocatable :: Gd, Hd
		real(dp), dimension(:,:), allocatable :: Gs, Hs, M
		real(dp), dimension(:,:), allocatable :: c ! constant part of H matrix
		complex(dp), dimension(:), allocatable :: U, t ! Total displacement and traction
		complex(dp), dimension(:), allocatable :: Uff, Tff ! free field motions
		!
		contains
		!
		procedure, pass(this) :: mat_alloc => bem_matrix_allocate
		procedure, pass(this) :: mat_dealloc => bem_matrix_allocate
		!
	end type bemMat
	!
	!
	contains
	!
	subroutine bem_matrix_allocate(this, ndof, nnode, nel_con, idn)
		!
		! allocates the required matrix for the matrix group chunkwise
		!
		class(bemMat), intent(inout) :: this
		integer, intent(in) :: ndof ! Number of DOF
		integer, intent(in) :: nnode ! number of nodes
		integer, intent(in) :: nel_con ! number of nodes in element connectivity
		integer, intent(in) :: idn ! allocation identifier
		!
		!
		! allocate the required matrices (now default allocation)
		select case(idn)
		case(1) ! global matrices
		if (allocated(this%Gd)) deallocate(this%Gd)
		if (allocated(this%Hd)) deallocate(this%Hd)
		!
		allocate(this%Gd(ndof*nnode, ndof*nel_con))
		this%Gd = (0.0_dp, 0.0_dp)
		allocate(this%Hd(ndof*nnode, ndof*nnode))
		this%Hd = (0.0_dp, 0.0_dp)
		!
		! T and U vectors
		if (allocated(this%U)) deallocate(this%U)
		if (allocated(this%t)) deallocate(this%t)
		!
		allocate(this%U(ndof*nnode))
		allocate(this%t(ndof*nel_con))
		!
		this%U = (0.0_dp, 0.0_dp)
		this%t = (0.0_dp, 0.0_dp)
		!
		! free field motion
		if (allocated(this%Uff)) deallocate(this%Uff)
		if (allocated(this%Tff)) deallocate(this%Tff)
		!
		allocate(this%Uff(ndof*nnode))
		allocate(this%Tff(ndof*nel_con))
		!
		this%Uff = (0.0_dp, 0.0_dp)
		this%Tff = (0.0_dp, 0.0_dp)
		!
		!
		case(2) ! local chunk matrics for elastodynamics
		! currently only H matrices required in part level for rigid body motion
		if (allocated(this%Gd)) deallocate(this%Gd)
		if (allocated(this%Hd)) deallocate(this%Hd)
		if (allocated(this%Hs)) deallocate(this%Hs)
		!
		allocate(this%Gd(ndof*nnode, ndof*nel_con))
		allocate(this%Hd(ndof*nnode, ndof*nnode))
		this%Hd = (0.0_dp, 0.0_dp)
		allocate(this%Hs(ndof*nnode, ndof*nnode))
		this%Hs = 0.0_dp
		!
		!
		case(3) ! element level matrix allocation 
		if (allocated(this%Gd)) deallocate(this%Gd)
		if (allocated(this%Hd)) deallocate(this%Hd)
		if (allocated(this%Hs)) deallocate(this%Hs)
		!
		allocate(this%Gd(ndof*nnode, ndof*nel_con))
		this%Gd = (0.0_dp, 0.0_dp)
		allocate(this%Hd(ndof*nnode, ndof*nel_con))
		this%Hd = (0.0_dp, 0.0_dp)
		allocate(this%Hs(ndof*nnode, ndof*nel_con))
		this%Hs = 0.0_dp
		!
		!
		! free field motion
		if (allocated(this%Uff)) deallocate(this%Uff)
		if (allocated(this%Tff)) deallocate(this%Tff)
		!
		allocate(this%Uff(ndof*nel_con))
		allocate(this%Tff(ndof*nel_con))
		!
		this%Uff = (0.0_dp, 0.0_dp)
		this%Tff = (0.0_dp, 0.0_dp)
		!
		case default ! error message
			!
			write(*,*) 'Identifier mismatch for matrix allocation'
		end select
		!
	end subroutine bem_matrix_allocate
	!
	subroutine bem_matrix_deallocate(this)
		!
		! deallocate the matrices within the group of matrices
		!
		class(bemMat), intent(inout) :: this
		!
		if (allocated(this%Gd)) deallocate(this%Gd)
		if (allocated(this%Hd)) deallocate(this%Hd)
		if (allocated(this%Hs)) deallocate(this%Hs)
		if (allocated(this%Gs)) deallocate(this%Gs)
		if (allocated(this%Uff)) deallocate(this%Uff)
		if (allocated(this%Tff)) deallocate(this%Tff)
		!
		! T and U vectors
		if (allocated(this%U)) deallocate(this%U)
		if (allocated(this%t)) deallocate(this%t)
	end subroutine bem_matrix_deallocate
	!
	!
	subroutine matrix_mapping_real(gmat, lmat, ndof, cmap, rmap)
		!
		! updates G type matrix (node rows and element columns)
		!
		real(dp), dimension(:,:), intent(inout) :: gmat ! global matrix
		real(dp), dimension(:,:), intent(in) :: lmat ! local matrix
		integer, intent(in), dimension(:) :: cmap !  column mapping mapping
		integer, optional, intent(in), dimension(:) :: rmap ! rows mapping(optional)
		! if node mapping is not present it is understood that the whole nodes are passed to the array
		integer, intent(in) :: ndof ! number of degrees of freedom
		integer :: col, row ! temporary indices for rows and columns
		integer :: ii, jj ! loop indices
		!
		do ii= 1, size(cmap)
			col = (cmap(ii)-1)*ndof ! the column number of last entry in the global matrix
			if (.not. present(rmap)) then
				
				!write(*,*) 'column numbers', ii , col, (ii-1)*ndof+1, col+1, ii*ndof, col+ndof
				gmat(:,col+1:col+ndof) = gmat(:,col+1:col+ndof)+ lmat(:, (ii-1)*ndof+1: ii*ndof)
			else
				do jj = 1, size(rmap)
					! updating the global matrix node wise as well
					row = (rmap(jj)-1)*ndof ! the row number for last entry in the global matrix
					gmat(row+1:row+ndof,col+1:col+ndof) = gmat(row+1:row+ndof,col+1:col+ndof) + lmat((jj-1)*ndof+1:jj*ndof, (ii-1)*ndof+1: ii*ndof)
				end do ! jj
			end if 
		end do ! ii
		!
	end subroutine matrix_mapping_real
	!
	!
	subroutine matrix_mapping_complex(gmat, lmat, ndof, cmap, rmap)
		!
		! updates G type matrix (node rows and element columns)
		!
		complex(dp), dimension(:,:), intent(inout) :: gmat ! global matrix
		complex(dp), dimension(:,:), intent(in) :: lmat ! local matrix
		integer, intent(in), dimension(:) :: cmap !  column mapping mapping
		integer, optional, intent(in), dimension(:) :: rmap ! rows mapping(optional)
		! if node mapping is not present it is understood that the whole nodes are passed to the array
		integer, intent(in) :: ndof ! number of degrees of freedom
		integer :: col, row ! temporary indices for rows and columns
		integer :: ii, jj ! loop indices
		!
		do ii= 1, size(cmap)
			col = (cmap(ii)-1)*ndof ! the column number of last entry in the global matrix
			if (.not. present(rmap)) then
				
				!write(*,*) 'column numbers', ii , col, (ii-1)*ndof+1, col+1, ii*ndof, col+ndof
				gmat(:,col+1:col+ndof) = gmat(:,col+1:col+ndof)+ lmat(:, (ii-1)*ndof+1: ii*ndof)
			else
				do jj = 1, size(rmap)
					! updating the global matrix node wise as well
					row = (rmap(jj)-1)*ndof ! the row number for last entry in the global matrix
					gmat(row+1:row+ndof,col+1:col+ndof) = gmat(row+1:row+ndof,col+1:col+ndof) + lmat((jj-1)*ndof+1:jj*ndof, (ii-1)*ndof+1: ii*ndof)
				end do ! jj
			end if 
		end do ! ii
		!
	end subroutine matrix_mapping_complex
	!
	!
	subroutine vector_mapping_real(gmat, lmat, ndof, rmap)
		!
		! updates G type matrix (node rows and element columns)
		!
		real(dp), dimension(:), intent(inout) :: gmat ! global matrix
		real(dp), dimension(:), intent(in) :: lmat ! local matrix
		integer, intent(in), dimension(:) :: rmap !  row mapping mapping
		! if node mapping is not present it is understood that the whole nodes are passed to the array
		integer, intent(in) :: ndof ! number of degrees of freedom
		integer :: row ! temporary indices for rows and rowumns
		integer :: ii, jj ! loop indices
		!
		do ii= 1, size(rmap)
			row = (rmap(ii)-1)*ndof ! the rowumn number of last entry in the global matrix
			gmat(row+1:row+ndof) = gmat(row+1:row+ndof)+ lmat((ii-1)*ndof+1: ii*ndof)
		end do ! ii
		!
	end subroutine vector_mapping_real
	!
	!
	subroutine vector_mapping_complex(gmat, lmat, ndof, rmap)
		!
		! updates G type matrix (node rows and element rowumns)
		!
		complex(dp), dimension(:), intent(inout) :: gmat ! global matrix
		complex(dp), dimension(:), intent(in) :: lmat ! local matrix
		integer, intent(in), dimension(:) :: rmap !  rowumn mapping mapping
		! if node mapping is not present it is understood that the whole nodes are passed to the array
		integer, intent(in) :: ndof ! number of degrees of freedom
		integer :: row ! temporary indices for rows and rowumns
		integer :: ii ! loop indices
		!
		do ii= 1, size(rmap)
			row = (rmap(ii)-1)*ndof ! the rowumn number of last entry in the global matrix
			gmat(row+1:row+ndof) = gmat(row+1:row+ndof)+ lmat((ii-1)*ndof+1: ii*ndof)
		end do ! ii
		!
	end subroutine vector_mapping_complex
	!
	end module mat_ops