	! main executable file 
	program main
	use model_
	use material_
	use gaussquad_
	use comp_model
	use solver_
	use bem_el
	!
	implicit none
	!
	integer, parameter :: dp= kind(0.d0)
	integer, parameter :: n = 2
	integer :: ii, jj, kk ! do loop variables
	real(dp) :: w(n), x(n)
	!
	type(Model), target :: mdl ! model space
	class(Step), allocatable :: cstep(:) ! computational step
	!
	! Temporary variables
	complex(dp) :: A(2,2), B(2), C(2), D(2)
	real(dp) :: SA
	!real(dp), allocatable :: tempCoord(:)
	!
	!nodes
	call mdl%create_h5('./inp/firstTry.h5')
	call initialize_step(cstep, mdl)

	! ------------------------------------------------
	! writing on the terminal screen
	!
	write(*,*) '*Node', shape(mdl%mesh%node%coord)
	do ii= 1,mdl%mesh%num_node
		write(*,*) mdl%mesh%node%label(ii),  mdl%mesh%node%coord(:,ii)
	end do

	! elements
	write(*,*) '*elements'
	do ii= 1,mdl%mesh%num_element
			write(*,*)  mdl%mesh%element(ii)%label, mdl%mesh%element(ii)%nodeset
	end do

	! node set
	do ii= 1,mdl%mesh%num_nodeset
		write(*,*) '*nset'
		write(*,*)  mdl%mesh%nodeset(ii)%set
	end do

	! element set
	do ii= 1,mdl%mesh%num_elementset
		write(*,*) '*elset'
		write(*,*)  mdl%mesh%elementset(ii)%set
	end do

	!material
	do ii = 1,mdl%num_material
		write(*,*) '*material'
		write(*,*) '*density=', mdl%material(ii)%rho
		write(*,*) '*E=', mdl%material(ii)%E
		write(*,*) '*nu=', mdl%material(ii)%nu
	end do

	! BEM chunk
	write(*,*) 'The number of BEM chunk are', mdl%num_bchunk
	do ii = 1, mdl%num_bchunk
		write(*,*) '*BEM, bound=',  mdl%bchunk(ii)%bound, 'material=', mdl%bchunk(ii)%material
		write(*,*) '*bnode'
		write(*,*)  mdl%bchunk(ii)%b_node
		write(*,*) '*belm'
		write(*,*)  mdl%bchunk(ii)%b_element
		if (mdl%bchunk(ii)%num_inode>0) then
			write(*,*) '*inode'
			write(*,*)  mdl%bchunk(ii)%i_node
		end if
	end do
	!
	!
	write(*,*) 'The number of steps are', mdl%num_step
	do ii = 1, mdl%num_step
		write(*,*) 'Step type:', mdl%step(ii)%stype, 'frequency', mdl%step(ii)%frequency
		write(*,*) '-------------------------------------------------------------------'
	end do
	! ------------------------------------------------------------------------
	!
	! Temporarily changing y and z coordinates in the model
	!if (.not. allocated(tempCoord)) allocate(tempCoord(mdl%mesh%num_node))
	!tempCoord = mdl%mesh%node%coord(2,:)
	!mdl%mesh%node%coord(2,:) = mdl%mesh%node%coord(3,:)
	!mdl%mesh%node%coord(3,:) = tempCoord
	!if (allocated(tempCoord)) deallocate(tempCoord)
	!
	!
	! assigning coordinates to each elements
	do ii = 1, mdl%mesh%num_element
		!write(*,*) ii, 'The element coordinates are'
		call mdl%mesh%element(ii)%assign_coord(mdl%mesh%node%coord)
		!write(*,*) mdl%mesh%element(ii)%coord
		!call bem(mdl%mesh%element(ii), mdl%mesh%node, mdl%material(1), (1.0_dp, 0.0_dp))!, G, Hd, Hs)
	end do
	!
	!
	! Using computational step for the moment
	!

	call cstep(1)%compute_bem(mdl)
	!
	!
	! solving the system
	call bem_solver(cstep(1)%gmat%U, cstep(1)%gmat%Gd, cstep(1)%gmat%Hd, cstep(1)%gmat%t, cstep(1)%gmat%Uff, cstep(1)%gmat%Tff)
	!
	!
	write(*,*) 'Free field traction:'
	do ii = 1, mdl%mesh%num_node
				write(*,*) mdl%mesh%node%coord(1,ii), mdl%mesh%node%coord(2,ii), mdl%mesh%node%coord(3,ii), &
				abs(cstep(1)%gmat%Tff(ii*3-2)), abs(cstep(1)%gmat%Tff(ii*3-1)), cstep(1)%gmat%Tff(ii*3)
	
	end do
	!
	!
	write(*,*) 'Solved displacement vector:'
	do ii = 1, mdl%mesh%num_node
				write(*,*) mdl%mesh%node%coord(1,ii), mdl%mesh%node%coord(2,ii), mdl%mesh%node%coord(3,ii), &
						abs(cstep(1)%gmat%U(ii*3-2)), abs(cstep(1)%gmat%U(ii*3-1)), REAL(cstep(1)%gmat%U(ii*3)), aimag(cstep(1)%gmat%U(ii*3))
	end do
	!
	!
	! check the solver
!	A = reshape(&
!    (/ (1.0, 0.0), (3.0,0.0), (2.0,0.0), (4.0,0.0) /), shape=(/ 2, 2 /))
!    B = (/ (2.0,0.0), (3.0,0.0) /)
!    C = (/ (0.0,0.0), (0.0,0.0) /)
    !
!    call bem_solver(U=D, G=A, H=A, t=C, Phi=B)
!    write(*,*) 'The solved vector is'
!    write(*,*) D
	!
	!call surface_area(SA, 10)
	!write(*,*) 'Surface Area:', SA
	!
	!
	end program main