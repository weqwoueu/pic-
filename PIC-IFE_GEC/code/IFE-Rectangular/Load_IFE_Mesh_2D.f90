SUBROUTINE Load_IFE_Mesh_2D(	IFE_Mesh_Filename,		    &
							t_basic_int, p_basic, e_basic,  &
                            n_elements,information_1,information_2,information_3,information_3_D,&
                            element_index, Node_Index, edge_index, edge_index_D, &
                            HP, HT, HE, E_Dirichlet, P_average, P_flag)

! Purpose:		Load IFE mesh from IFE mesh file(SIDG)
! Last Update:	7/22/2022 8:34 PM
USE IFE_Boundary
IMPLICIT NONE

CHARACTER(*)							IFE_Mesh_Filename
INTEGER, DIMENSION(:,:), POINTER	::	t_basic_int, e_basic, HT, HE, E_Dirichlet
REAL(8), DIMENSION(:,:), POINTER	::	p_basic, HP, P_average, P_flag
!INTEGER, DIMENSION(:), POINTER		::	Node_Index
INTEGER, DIMENSION(:,:),POINTER   ::  information_1
REAL(8), DIMENSION(:,:), POINTER  ::  information_2, information_3, information_3_D
INTEGER, DIMENSION(:), POINTER		::	element_index
INTEGER, DIMENSION(:), POINTER		::	Node_Index, edge_index, edge_index_D


INTEGER :: n_elements, n_nodes, n_int_elements, i, j, N_Objects, N_Boundary, n_bc_edge, &
          n_elements_A, n_nodes_DG, n_elements_DG, n_edge_inf3, n_edge_inf3_D, n_edge, &
          n_edge_D, n_nodes_P_average, n_nodes_P_flag

! Mesh and Interface elements data file
OPEN(1,ACTION='READ',FILE='ife.msh')
READ(1,*)
READ(1,*) N_Objects
READ(1,*)
READ(1,*) N_Boundary

READ(1,*)
READ(1, *) n_nodes
ALLOCATE(p_basic(2,n_nodes))
DO	i=1,n_nodes
	READ(1,*) p_basic(1,i), p_basic(2,i)!, element_index(i)
END DO

READ(1,*)
READ(1, *) n_elements

ALLOCATE(t_basic_int(4,n_elements))	! The fifth entry is used to store hardwiring information

t_basic_int= -1

DO	i=1,n_elements
	READ(1,*) t_basic_int(1,i),t_basic_int(2,i),t_basic_int(3,i),t_basic_int(4,i) !t_c(4,n) the index of of every grid
END DO


READ(1,*)
READ(1,*) n_elements_A                  !LY REVISE,2021-11-27
ALLOCATE(element_index(n_elements_A))   !LY REVISE,2021-11-27
DO i=1, n_elements_A
  READ(1,*) element_index(i)
ENDDO

!------------------------------------LY REVISE,2021-11-27-------------------------------
READ(1,*)
READ(1,*) n_nodes_DG
ALLOCATE(HP(2, n_nodes_DG))
DO	i=1, n_nodes_DG
	READ(1,*) HP(1,i), HP(2,i)
END DO

READ(1,*)
READ(1,*) n_elements_DG
ALLOCATE(HT(5, n_elements_DG))
DO	i=1,n_elements_DG
	READ(1,*) HT(1,i),HT(2,i),HT(3,i),HT(4,i),HT(5,i)
END DO

READ(1,*)
READ(1,*) n_edge
ALLOCATE(HE(8, n_edge))
DO	i=1,n_edge
	READ(1,*) HE(1,i),HE(2,i),HE(3,i),HE(4,i),HE(5,i),HE(6,i),HE(7,i),HE(8,i)
END DO 

READ(1,*)
READ(1,*) n_edge_D
ALLOCATE(E_Dirichlet(7, n_edge_D))
DO	i=1,n_edge_D
	READ(1,*) E_Dirichlet(1,i),E_Dirichlet(2,i), &
                E_Dirichlet(3,i),E_Dirichlet(4,i), E_Dirichlet(5, i), E_Dirichlet(6, i), E_Dirichlet(7, i)
END DO

READ(1,*)
READ(1,*) n_edge_inf3
ALLOCATE(information_3(7, n_edge_inf3))
DO	i=1,n_edge_inf3
	READ(1,*) information_3(1,i),information_3(2,i),information_3(3,i),information_3(4,i), &
                information_3(5,i),information_3(6,i),information_3(7,i)
END DO


READ(1,*)
READ(1,*) n_edge_inf3_D
ALLOCATE(information_3_D(7, n_edge_inf3_D))
DO	i=1,n_edge_inf3_D
	READ(1,*) information_3_D(1,i),information_3_D(2,i),information_3_D(3,i), &
                information_3_D(4,i), information_3_D(5,i),information_3_D(6,i), &
                information_3_D(7,i)
END DO

READ(1,*)
ALLOCATE(edge_index(n_edge))
DO i=1, n_edge
    READ(1,*) edge_index(i)
ENDDO

READ(1,*)
ALLOCATE(edge_index_D(n_edge_D))
DO i=1, n_edge_D
    READ(1,*) edge_index_D(i)
ENDDO

READ(1,*)
READ(1,*) n_nodes_P_average
ALLOCATE(P_average(6, n_nodes_P_average))
DO	i=1, n_nodes_P_average
	READ(1,*) P_average(1,i), P_average(2,i), P_average(3,i), P_average(4,i), &
            P_average(5,i), P_average(6,i)
END DO

READ(1,*)
READ(1,*) n_nodes_P_flag
ALLOCATE(P_flag(5, n_nodes_P_flag))
DO	i=1, n_nodes_P_flag
	READ(1,*) P_flag(1,i), P_flag(2,i), P_flag(3,i), P_flag(4,i), P_flag(5,i)
END DO
!---------------------------------------------------------------------------------------

IF (N_Objects /= 0  ) THEN
   READ(1,*)
   READ(1, *) n_int_elements
   IF (n_int_elements /= 0) THEN
   ALLOCATE(information_1(18,n_int_elements), information_2(8,n_int_elements))
   READ(1,*)
   DO i=1, n_int_elements
    READ(1,*) information_1(:,i)
    READ(1,*) information_2(:,i)
   ENDDO
   ELSEIF (n_int_elements == 0) THEN
   ALLOCATE(information_1(18,n_elements), information_2(8,n_elements))
 
	information_1 = 0
	information_2 = 0

ENDIF


   !READ(1,*)
   !ALLOCATE(Node_Index(n_nodes))
   ! DO	i=1,n_nodes
	  ! READ(1,*) Node_Index(i)
   ! END DO
   READ(1,*)
   ALLOCATE(Node_Index(SIZE(HP,2)))     !LY REVISE, 2021-11-27
    DO	i=1,SIZE(HP,2)
	   READ(1,*) Node_Index(i)
    END DO

    
ENDIF

IF (N_Boundary /= 0 ) THEN
   READ(1,*)
   ALLOCATE(bc_index(N_Boundary), bc_value(N_Boundary), bc_point_1(2,N_Boundary), bc_point_2(2,N_Boundary))
   DO	i=1,N_Boundary
    READ(1, *)bc_index(i), bc_value(i)
    READ(1, *)bc_point_1(:,i)
	  READ(1, *)bc_point_2(:,i)
   ENDDO

   READ(1, *)
   READ(1, *) n_bc_edge
   ALLOCATE(e_basic(5,n_bc_edge))
   DO	i=1,n_bc_edge
     READ(1, *) e_basic(:,i)
   ENDDO
ENDIF


CLOSE(1)


END 