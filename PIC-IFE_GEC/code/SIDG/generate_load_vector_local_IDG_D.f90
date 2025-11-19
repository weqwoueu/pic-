SUBROUTINE generate_load_vector_local_IDG_D(Global_Beta, node_type, num_of_unknowns, &
                                            information_1, information_2, information_3, HP, element_index, edge_index, &
                                            HT, HE, &
                                            Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
                                            test_basis_type, delta, con_penalty, vector)
  
! -------------------LY add vector of IDG edge for IDG solver 2021-9-28-------------------

USE IFE_MAIN_PARAM
USE Gauss_Data
!=========LY modification, 2022-7-25=========
Use IFE_Data, Only: E_DirichletValue
!=========LY modification, 2022-7-25=========
USE IFE_INTERFACE, ONLY: generate_Gauss_local_1D_Linear, gauss_integration_local_vector_DG_FE_D, &
                         gauss_integration_local_vector_IDG_D

IMPLICIT NONE


REAL(8), DIMENSION(:), POINTER    :: Global_Beta
INTEGER, DIMENSION(:,:), POINTER  :: node_type
INTEGER, INTENT(IN)               :: num_of_unknowns
INTEGER, DIMENSION(:,:), POINTER  :: information_1
REAL(8), DIMENSION(:,:), POINTER  :: information_2, information_3
REAL(8), DIMENSION(:,:), POINTER  :: HP
INTEGER, DIMENSION(:,:), POINTER  :: HT, HE
INTEGER, DIMENSION(:), POINTER    :: element_index, edge_index
REAL(8), DIMENSION(:)             :: Gauss_coefficient_reference_1D
REAL(8), DIMENSION(:)             :: Gauss_point_reference_1D
INTEGER                           :: test_basis_type
INTEGER, INTENT(IN)               :: delta
REAL(8)                           :: con_penalty
REAL(8), DIMENSION(:), POINTER    :: vector

INTEGER                           :: out_index, in_index
INTEGER                           :: i, j, k, m, n
REAL(8), DIMENSION(:), POINTER    :: Gauss_coefficient_local
REAL(8), DIMENSION(:,:), POINTER  :: Gauss_point_local

INTEGER                           :: information_vector_1(18)
REAL(8)                           :: information_vector_2(8)
REAL(8)                           :: vertices_linear(2,2)
REAL(8)                           :: vertices_Rectangular(2,4)

INTEGER,DIMENSION(4,2)            :: this_flag1

!-----------------------use the sparse matrix to store the Global stiff-----------------------
INTEGER                           ::  NZ
INTEGER                           ::  n_nodes_in_elem, m_unknowns, loc1, loc2
REAL(8)                           ::  int_value
REAL(8)                           ::  El_Vector(4)
!DATA                                  El_Vector /0.0, 0.0, 0.0, 0.0/

m_unknowns = MAXVAL(node_type)
ALLOCATE(vector(m_unknowns))
vector = 0.0
this_flag1 = 0
n_nodes_in_elem = 4

out_index = -1
in_index = -2

DO n = 1, SIZE(HE, 2)

  !LY Add for Boundary Condition, 2021-11-18
  IF (HE(6,n) == 1) THEN    !This boundary edge is Dirichlet boundary edge.
    El_Vector = 0
    vertices_linear(:, 1) = HP(1:2, HE(1,n))
    vertices_linear(:, 2) = HP(1:2, HE(2,n))
    vertices_Rectangular(:, :) = HP(:, HT(1:4, HE(5, n)))
    this_flag1(2:3,1) = HE(3:4,n)
    this_flag1(2:3,2) = HE(3:4,n)
    this_flag1(4,:) = 1
    
    CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, vertices_linear(:,1), &
                                        vertices_linear(:,2), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                        Gauss_Point_Local_1D_Linear_Eight)
    
    IF (element_index(HE(5,n)) == out_index) THEN
    
      DO k = 1, 4
      
        CALL gauss_integration_local_vector_DG_FE_D(delta, Global_Beta(1), Gauss_Coefficient_Local_1D_Linear_Eight, &
                              Gauss_Point_Local_1D_Linear_Eight, &
              this_flag1, vertices_Rectangular, test_basis_type, k, con_penalty, E_DirichletValue(n), int_value, HE(5,n))
        
        El_Vector(k) = El_Vector(k) + int_value
        
      END DO
    
    ELSE IF (element_index(HE(5,n)) == in_index) THEN
    
     DO k = 1, 4
      
        CALL gauss_integration_local_vector_DG_FE_D(delta, Global_Beta(2), Gauss_Coefficient_Local_1D_Linear_Eight, &
                              Gauss_Point_Local_1D_Linear_Eight, &
              this_flag1, vertices_Rectangular, test_basis_type, k, con_penalty, E_DirichletValue(n), int_value, HE(5,n))
        
        El_Vector(k) = El_Vector(k) + int_value
        
      END DO 
    
    ELSE IF (element_index(HE(5,n)) > 0) THEN
    
      information_vector_1 = information_1(:, element_index(HE(5, n)))
      information_vector_2 = information_2(:, element_index(HE(5, n)))
      
      IF (edge_index(n) < 0) THEN
      
        DO k = 1, 4
        
          CALL gauss_integration_local_vector_IDG_D(delta, Global_Beta, information_vector_1, information_vector_2, &
                                        Gauss_Coefficient_Local_1D_Linear_Eight, Gauss_Point_Local_1D_Linear_Eight, &
                                              this_flag1, vertices_Rectangular, test_basis_type, k, edge_index(n), &
                                                con_penalty, E_DirichletValue(n), int_value, HE(5,n))
          
          El_Vector(k) = El_Vector(k) + int_value
          
        END DO
      
      ELSE IF (edge_index(n) > 0) THEN
      
        CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices_linear(:,1),&
                                            information_3(4:5, edge_index(n)), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                            Gauss_Point_Local_1D_Linear_Eight)
        
        loc1 = information_3(6, edge_index(n))
        
        DO k = 1, 4
        
          CALL gauss_integration_local_vector_IDG_D(delta, Global_Beta, information_vector_1, information_vector_2, &
                                      Gauss_Coefficient_Local_1D_Linear_Eight, Gauss_Point_Local_1D_Linear_Eight, &
                                                    this_flag1, vertices_Rectangular, test_basis_type, k, loc1, &
                                                    con_penalty, E_DirichletValue(n), int_value, HE(5,n))
          
          El_Vector(k) = El_Vector(k) + int_value
          
        END DO
        
        CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
                                            information_3(4:5, edge_index(n)), &
                                            vertices_linear(:,2), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                            Gauss_Point_Local_1D_Linear_Eight)
        
        loc2 = information_3(7, edge_index(n))
        
        DO k = 1, 4
        
          CALL gauss_integration_local_vector_IDG_D(delta, Global_Beta, information_vector_1, information_vector_2, &
                                        Gauss_Coefficient_Local_1D_Linear_Eight, Gauss_Point_Local_1D_Linear_Eight, &
                                                    this_flag1, vertices_Rectangular, test_basis_type, k, loc2, &
                                                    con_penalty, E_DirichletValue(n), int_value, HE(5,n))
          
          El_Vector(k) = El_Vector(k) + int_value
          
        END DO
        
      ELSE
        WRITE(6,*) 'error edge_index in generate_load_vector_local_IDG_D'
        STOP
      END IF
      
    ELSE
      WRITE(6,*) 'error element_index in generate_load_vector_local_IDG_D'
      STOP
    END IF
    
      !=======================STORE VECTOR================================LY: not definite
      DO k = 1, n_nodes_in_elem
  
        IF (node_type(2, HT(k, HE(5, n)))>0) THEN
    
          vector(node_type(2, HT(k, HE(5, n)))) = vector(node_type(2, HT(k, HE(5, n)))) + El_Vector(k)
          
        END IF
        
      END DO
      !=======================STORE VECTOR================================
  END IF
  
  
END DO

END SUBROUTINE