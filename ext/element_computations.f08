	module element_sub
	! contains the element discretization functions
	! from the Book : Boundary Element Method with Programming
	!
	implicit none
	!private
	public Serendip_func, Normal_Jac, integration_pt, vector_norm, vector_ex
	!
	integer, private, parameter :: dp=kind(0.d0)
	!e
	!
	contains
	!
	!
	SUBROUTINE VECTOR_NORM(V,VLEN)
	!----------------------------------------
	! Normalise vector
	!----------------------------------------
	REAL(kind=dp), INTENT(INOUT) :: V(:)
	! Vector to be normalised
	REAL(kind=dp), INTENT(OUT) :: VLEN ! Length of vector
	VLEN= SQRT( SUM(V*V))
	IF(ABS(VLEN)<1.E-10) RETURN
	V= V/VLEN
	RETURN
	END SUBROUTINE VECTOR_NORM
	!
	!
	!
	FUNCTION VECTOR_EX(V1,V2)
	!----------------------------------------
	! Returns vector x-product v1xv2
	! where v1 and v2 are dimension 3
	!----------------------------------------
	REAL(kind=dp), INTENT(IN) :: V1(3),V2(3)
	!
	REAL(kind=dp) :: VECTOR_EX(3)
	!
	VECTOR_EX(1)=V1(2)*V2(3)-V2(2)*V1(3)
	VECTOR_EX(2)=V1(3)*V2(1)-V1(1)*V2(3)
	VECTOR_EX(3)=V1(1)*V2(2)-V1(2)*V2(1)
	RETURN
	END FUNCTION VECTOR_EX
	!
	!
	!calculates array with shape functions 
	SUBROUTINE Serendip_func(Ni,xsi,eta,ldim,nodes,inci)
	!---------------------------------
	! Computes Serendipity shape functions Ni(xsi,eta)
	! for one and two-dimensional (linear/parabolic) finite
	! boundary elements
	!---------------------------------
	REAL(kind=dp),INTENT(OUT) :: Ni(:) ! Array with shape function
	REAl(kind=dp),INTENT(IN) :: xsi,eta! intrinsic coordinates
	INTEGER,INTENT(IN):: ldim ! element dimension
	INTEGER,INTENT(IN):: nodes ! number of nodes
	INTEGER,optional, INTENT(IN):: inci(:)! element incidences
	REAL(kind=dp) :: mxs,pxs,met,pet
	! temporary variables
	SELECT CASE (ldim)
	CASE(1)! one-dimensional element
		Ni(1)= 0.5*(1.0 - xsi)
		Ni(2)= 0.5*(1.0 + xsi)
		IF(nodes == 2) RETURN! linear element finished
		Ni(3)= 1.0 - xsi*xsi
		Ni(1)= Ni(1) - 0.5*Ni(3)
		Ni(2)= Ni(2) - 0.5*Ni(3)
	CASE(2)! two-dimensional element
		mxs=1.0-xsi
		pxs=1.0+xsi
		met=1.0-eta
		pet=1.0+eta
		Ni(1)= 0.25*mxs*met
		Ni(2)= 0.25*pxs*met
		Ni(3)= 0.25*pxs*pet
		Ni(4)= 0.25*mxs*pet
		IF(nodes == 4) RETURN! linear element finished
!		IF(Inci(5) > 0) THEN !zero node = node missing
!			Ni(5)= 0.5*(1.0 -xsi*xsi)*met
!			Ni(1)= Ni(1) - 0.5*Ni(5)
!			Ni(2)= Ni(2) - 0.5*Ni(5)
!		END IF
!		IF(Inci(6) > 0) THEN
!			Ni(6)= 0.5*(1.0 -eta*eta)*pxs
!			Ni(2)= Ni(2) - 0.5*Ni(6)
!			Ni(3)= Ni(3) - 0.5*Ni(6)
!		END IF
!		IF(Inci(7) > 0) THEN
!			Ni(7)= 0.5*(1.0 -xsi*xsi)*pet
!			Ni(3)= Ni(3) - 0.5*Ni(7)
!			Ni(4)= Ni(4)- 0.5*Ni(7)
!		END IF
!		IF(Inci(8) > 0) THEN
!			Ni(8)= 0.5*(1.0 -eta*eta)*mxs
!			Ni(4)= Ni(4) - 0.5*Ni(8) 
!			Ni(1)= Ni(1) - 0.5*Ni(8)
!		END IF
	CASE DEFAULT ! error message
		!CALL Error_message('Element dimension not 1 or 2')
	END SELECT
	RETURN
	END SUBROUTINE Serendip_func
	!
	!
	!
	SUBROUTINE Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
	!---------------------------------
	! Computes Derivatives of Serendipity shape functions
	! for one and two-dimensional (linear/parabolic)
	! finite boundary elements
	!---------------------------------
	REAL(kind=dp),INTENT(OUT) :: DNi(:,:) ! Derivatives of Ni
	REAL(kind=dp), INTENT(IN) :: xsi,eta ! intrinsic coordinates
	INTEGER,INTENT(IN):: ldim
	! element dimension
	INTEGER,INTENT(IN):: nodes
	! number of nodes
	INTEGER,optional, INTENT(IN):: inci(:) ! element incidences
	REAL(kind=dp):: mxs,pxs,met,pet
	! temporary variables
	SELECT CASE (ldim)
	CASE(1) ! one-dimensional element
		DNi(1,1)= -0.5
		DNi(2,1)= 0.5
		IF(nodes == 2)RETURN ! linear element finished
		DNi(3,1)= -2.0*xsi
		DNi(1,1)= DNi(1,1) - 0.5*DNi(3,1)
		DNi(2,1)= DNi(2,1) - 0.5*DNi(3,1)
	CASE(2) ! two-dimensional element
		mxs= 1.0-xsi
		pxs= 1.0+xsi
		met= 1.0-eta
		pet= 1.0+eta
		DNi(1,1)= -0.25*met
		DNi(1,2)= -0.25*mxs
		DNi(2,1)= 0.25*met
		DNi(2,2)= -0.25*pxs
		DNi(3,1)= 0.25*pet
		DNi(3,2)= 0.25*pxs
		DNi(4,1)= -0.25*pet
		DNi(4,2)= 0.25*mxs
		IF(nodes == 4) RETURN ! linear element finshed
!		IF(Inci(5) > 0) THEN ! zero node = node missing
!			DNi(5,1)= -xsi*met
!			DNi(5,2)= -0.5*(1.0 -xsi*xsi)
!			DNi(1,1)= DNi(1,1) - 0.5*DNi(5,1)
!			DNi(1,2)= DNi(1,2) - 0.5*DNi(5,2)
!			DNi(2,1)= DNi(2,1) - 0.5*DNi(5,1)
!			DNi(2,2)= DNi(2,2) - 0.5*DNi(5,2)
!		END IF
!		IF(Inci(6) > 0) THEN
!			DNi(6,1)= 0.5*(1.0 -eta*eta)
!			DNi(6,2)= -eta*pxs
!			DNi(2,1)= DNi(2,1) - 0.5*DNi(6,1)
!			DNi(2,2)= DNi(2,2) - 0.5*DNi(6,2)
!			DNi(3,1)= DNi(3,1) - 0.5*DNi(6,1)
!			DNi(3,2)= DNi(3,2) - 0.5*DNi(6,2)
!		END IF
!		IF(Inci(7) > 0) THEN
!			DNi(7,1)= -xsi*pet
!			DNi(7,2)= 0.5*(1.0 -xsi*xsi)
!			DNi(3,1)= DNi(3,1) - 0.5*DNi(7,1)
!			DNi(3,2)= DNi(3,2) - 0.5*DNi(7,2)
!			DNi(4,1)= DNi(4,1) - 0.5*DNi(7,1)
!			DNi(4,2)= DNi(4,2) - 0.5*DNi(7,2)
!		END IF
!		IF(Inci(8) > 0) THEN
!			DNi(8,1)= -0.5*(1.0-eta*eta)
!			DNi(8,2)= -eta*mxs
!			DNi(4,1)= DNi(4,1) - 0.5*DNi(8,1)
!			DNi(4,2)= DNi(4,2) - 0.5*DNi(8,2)
!			DNi(1,1)= DNi(1,1) - 0.5*DNi(8,1)
!			DNi(1,2)= DNi(1,2) - 0.5*DNi(8,2)
!		END IF
	CASE DEFAULT ! error message
		! CALL Error_message('Element dimension not 1 or 2' )
	END SELECT
	RETURN
	END SUBROUTINE Serendip_deriv
	!
	!
	SUBROUTINE Normal_Jac(v3,Jac,xsi,eta,ldim,nodes,coords, inci)
	!-------------------------------------------------------
	! Computes normal vector and Jacobian
	!-------------------------------------------------------
	REAL(kind=dp),allocatable, INTENT(OUT) :: v3(:)
	! Vector normal to point
	REAL(kind=dp),INTENT(OUT) :: Jac
	! Jacobian
	REAL(kind=dp), INTENT(IN) :: xsi,eta
	! intrinsic coords of point
	INTEGER,INTENT(IN):: ldim
	! element dimension
	INTEGER,INTENT(IN):: nodes
	! number of nodes
	INTEGER,optional, INTENT(IN):: inci(:)
	! element incidences
	REAL(KIND=DP), INTENT(IN) :: coords(:,:)! node coordinates
	REAL(KIND=DP),ALLOCATABLE :: DNi(:,:)
	! Derivatives of Ni
	REAL(KIND=DP),ALLOCATABLE :: v1(:),v2(:)! Vectors in xsi,eta dir
	INTEGER :: Cdim ! Cartesian dimension
	integer :: i
	Cdim= ldim+1
	!Allocate temporary arrays
	ALLOCATE (DNi(nodes,Cdim),V1(Cdim),V2(Cdim))
	!Compute derivatives of shape function
	Call Serendip_deriv(DNi,xsi,eta,ldim,nodes)
	!
	! Compute vectors in xsi (eta) direction(s)
	DO I=1,Cdim
		V1(I)= DOT_PRODUCT(DNi(:,1),COORDS(I,:))
		IF(ldim == 2) THEN
			V2(I)= DOT_PRODUCT(DNi(:,2),COORDS(I,:))
		END IF
	END DO
	!Compute normal vector
	IF(ldim == 1) THEN
		v3(1)= V1(2)
		v3(2)= -v1(1)
	ELSE
		V3= Vector_ex(v1,v2)
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
	REAL(KIND=DP),intent(out) :: ipt(:)
	REAL(kind=dp),INTENT(in) :: Ni(:) ! Array with shape function
	!
	REAL(KIND=DP), INTENT(IN) :: coords(:,:)! node coordinates
	!
	! compute coordinates at integration point
	ipt(1) = DOT_PRODUCT(Ni,COORDS(1,:))
	ipt(2) = DOT_PRODUCT(Ni,COORDS(2,:))

	if (size(Coords,1) == 3) then
	ipt(3) = DOT_PRODUCT(Ni, coords(3,:))
	end if
	!
	end subroutine integration_pt
	!
	end module element_sub
	
