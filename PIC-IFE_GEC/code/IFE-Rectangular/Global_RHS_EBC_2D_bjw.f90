SUBROUTINE Global_RHS_EBC_2D_bjw (A_stiff_xt, p_basic, t_basic_int, e_basic, node_type, U_full, bnd_elem_index, rhs_EBC, delta)

! WSY REVISE FOR SIDG 2022/7/22
! Essential boundary conditions contribution to IFE RHS vector
! F_EBC = Int[Ud_j*Grad(Epsi_i).Grad(Epsi_j) dOmega]

USE IFE_MAIN_PARAM
USE IFE_Boundary
USE Gauss_Data
USE IFE_INTERFACE, ONLY: Userdef_EBC_Value_bjw

IMPLICIT NONE
TYPE(SPARSE), DIMENSION(:), POINTER	::	A_stiff_xt
REAL(8), DIMENSION(:,:), POINTER	  ::	p_basic 
INTEGER, DIMENSION(:), POINTER		  ::	bnd_elem_index
INTEGER, DIMENSION(:,:), POINTER	  ::	node_type, t_basic_int, e_basic
REAL(8), DIMENSION(:), POINTER		  ::	U_full, rhs_EBC

INTEGER		num_bound_elem, n_nodes_in_elem, i, j	
INTEGER   num_of_nodes, m_unknowns, delta
			
REAL(8), DIMENSION(:,:), POINTER			::	EBC_Value
REAL(8), DIMENSION(:), ALLOCATABLE      ::	EBC_Value_xt

WRITE(*,*) '================='
num_of_nodes = SIZE(p_basic,2)

num_bound_elem	= SIZE(bnd_elem_index)
n_nodes_in_elem = 4 ! tangular element

!ALLOCATE(EBC_Value(n_nodes_in_elem,SIZE(t_basic_int,2)))
!ALLOCATE(EBC_Value_xt(num_of_nodes))
!EBC_Value_xt=Zero
!EBC_Value=Zero
!CALL Userdef_EBC_Value_bjw(delta, p_basic, t_basic_int, node_type, EBC_Value, EBC_Value_xt)    ! For IFE TEST
!
!DO i=1,num_of_nodes
!    
!    U_full(i) = EBC_Value_xt(i)   
!
!END DO

!LY REVISE, 2021-11-30
!USE ZWZ's METHOD
DO i = 1, SIZE(e_basic,2)   !$ Traverse all the edges
    IF (bc_index(e_basic(4,i))==1) THEN     !$ Dirichlet Boundary
        U_full(e_basic(1,i)) = bc_value(e_basic(4,i))   !$ Dirichlet boundary edge first node
        !print*,bc_value(e_basic(4,i))
        IF (node_type(1,e_basic(1,i)) /= 1) THEN        !$ IF the first node is not Dirichlet node
            PRINT*,'node_type(1,e_basic(1,i)) /= 1', e_basic(1,i)
            STOP
        ENDIF
    ENDIF
    IF (bc_index(e_basic(5,i))==1) THEN     !$ Dirichlet Boundary
        U_full(e_basic(2,i)) = bc_value(e_basic(5,i))   !$ Dirichlet boundary edge second node
        IF (node_type(1,e_basic(2,i)) /= 1) THEN        !$ IF the first node is not Dirichlet node
            PRINT*,'node_type(1,e_basic(2,i)) /= 1', e_basic(2,i)
            STOP
        ENDIF
    ENDIF
ENDDO

U_full(68) =0.0   !四边诺伊曼，内部一点已知电势

m_unknowns = MAXVAL(node_type)
ALLOCATE(rhs_EBC(m_unknowns))
DO i =1, SIZE(rhs_EBC)
	rhs_EBC(i)	= Zero
END DO

DO i=1,num_of_nodes 
    IF (node_type(2, i)>0)THEN
        DO j=A_stiff_xt(node_type(2, i))%SROW, A_stiff_xt(node_type(2, i)+1)%SROW-1
            !WRITE(*,*) i ,j
            rhs_EBC(node_type(2, i))=rhs_EBC(node_type(2, i))+U_full(A_stiff_xt(j)%JCOL)*A_stiff_xt(j)%K
        END DO
    END IF
END DO




END


