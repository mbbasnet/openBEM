	! solver modules
	module solver_
	!
	implicit none
	private
	public bem_solver
	!
	integer, parameter :: dp = kind(0.d0)
	!
	!
	!external cgesv
	contains
	!
	subroutine bem_solver(U, G, H, t, Uff, Tff)
		!
		! solves system for U in HU - Gt = Phi
		! H and G are BEM matrices
		! U : Total displacement matrix
		! t : Total traction matrix
		! Phi : Body traction matrix of appropriate unit
		!
		complex(dp), intent(out), dimension(:) :: U
		complex(dp), intent(in), dimension(:,:) :: G, H
		complex(dp), intent(in), dimension(:) :: t, Uff, Tff
		!
		complex(dp), allocatable, dimension(:) :: Phi, tsc
		integer, allocatable, dimension(:,:) :: ipiv
		!
		integer :: n, nhrs
		integer :: lda, ldb
		integer :: info
		!
		!
		n = size(H,1)
		nhrs = 1
		lda = n
		ldb = n
		!info = 0
		!
		if (.not. allocated(ipiv)) allocate(ipiv(n,n))
		if (.not. allocated(Phi)) allocate(Phi, mold=Uff)
		if (.not. allocated(Phi)) allocate(tsc, mold=Tff)
		tsc = t - Tff
		Phi = MATMUL(H, Uff) + MATMUL(G, tsc)  
		!
		!write(*,*) 'Phi:'
		!write(*,*) Phi
		!
		! LU factorization of matrix H
		!
		! solving the matrix for A*x = B for x
		!
		!write(*,*) 'H matrix:', H
		call zgetrf( N, N, H, LDA, IPIV, INFO )
		!
		write(*,*) 'info:', info
		! solving the system of eqation
		!
		call zgetrs('N', n, nhrs, H, lda, ipiv, Phi, ldb, info )
		!
		!write(*,*) 'H matrix:', H
		write(*,*) 'info:', info
		U = Phi
		!	
		!write(*,*) 'Utotal:'
		!write(*,*) U
		!
		if (allocated(ipiv)) deallocate(ipiv)
		if (allocated(Phi)) deallocate(Phi)
		!
	end subroutine bem_solver
	!
	end module solver_
		
