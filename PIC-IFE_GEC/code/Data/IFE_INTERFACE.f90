MODULE IFE_INTERFACE

IMPLICIT NONE

! Inject
Interface
    SUBROUTINE InjectBeams_2D(it, delta, PB)
        Use ModuleParticleBundle
        INTEGER		it, delta
        Type(ParticleBundle), intent(inout) :: PB(:)
    End SUBROUTINE
    Subroutine generate_Elementmap(Element_Map, repeat_refinement)
    Integer, Allocatable, Intent(InOut) ::  Element_Map(:,:)
    Integer, Intent(In)                 ::  repeat_refinement
End Subroutine  
End Interface


! Drivers
! =======
INTERFACE
Subroutine ParticlePositioning(part_x, part_y, n_element_old, traverse_flag, EdgeArray)
    !========= Input and Output ================
    Integer, Intent(InOut) :: n_element_old
    !================================

    !========= Input ================
    Real(8), Intent(In) :: part_x, part_y
    !================================

    !========= Output ================
    Integer, Intent(Out) :: traverse_flag
    Integer, Allocatable, Intent(Out) :: EdgeArray(:)
    !================================
End Subroutine ParticlePositioning

    Subroutine Improve_HT(HT, HP,HE, repeat_refinement)
        REAL(8), DIMENSION(:,:), POINTER    ::	HP
        INTEGER, DIMENSION(:,:), POINTER    ::	HT
        INTEGER, DIMENSION(:,:), POINTER    ::  HE
        INTEGER                             ::  repeat_refinement
    End Subroutine
    
	SUBROUTINE IFE_Start_2D(    PIC_phi_bkgd, PIC_Te_bkgd, PIC_Rho_bkgd, nbkgd,		&
							    Phi, nnx, nny, delta, time	)

		INTEGER					nbkgd
		REAL(8)					PIC_phi_bkgd(nbkgd), PIC_Te_bkgd(nbkgd), PIC_Rho_bkgd(nbkgd)
		INTEGER					nnx, nny, delta
        REAL(8)					Phi(0:nnx+1,0:nny+1)  
        REAL(8),INTENT(IN) :: time !$ ab.ZWZ for AC objects
	END SUBROUTINE
	
	SUBROUTINE Input_2D(xmin, xmax, ymin, ymax, nnx, nny, N_Objects, N_Boundary, objects, bc_point_1, bc_point_2, bc_index, bc_value)
		USE Object_Data_2D
		REAL(8)		                                    ::  xmin, xmax, ymin, ymax
		INTEGER		                                    ::  nnx, nny, N_Objects, N_Boundary
		TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
		REAL(8), DIMENSION(:,:), POINTER	            ::	bc_point_1, bc_point_2
		INTEGER, DIMENSION(:), POINTER	                ::	bc_index
		REAL(8), DIMENSION(:), POINTER		            ::	bc_value
	END SUBROUTINE
	
	SUBROUTINE Setup_IFE_Mesh_2D(delta, xmin, xmax, ymin, ymax, nnx, nny, N_Objects, N_Boundary, objects)
		USE Object_Data_2D
		INTEGER		delta 
		REAL(8)		xmin, xmax, ymin, ymax
		INTEGER		nnx, nny, N_Objects, N_Boundary
		TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
	END SUBROUTINE
	
END INTERFACE

!------------------ DG PART------
INTERFACE

SUBROUTINE generate_Adaptive_DGT_IFE(DGP, DGT, objects, DGT_IFE_partition, area_flag)

	USE Object_Data_2D
	USE IFE_MAIN_PARAM
		
	REAL(8), DIMENSION(:,:), POINTER        ::  DGP
	INTEGER, DIMENSION(:,:), POINTER        ::  DGT, DGT_IFE_partition
	TYPE(ObjectType), DIMENSION(:), INTENT(IN) :: objects
	INTEGER :: area_flag
		
	END SUBROUTINE
  
  SUBROUTINE Adaptive_interface_function(x, y, objects, temp)

    USE Object_Data_2D
        
    TYPE(ObjectType), INTENT(IN) :: objects
    REAL(8), INTENT(IN) :: x, y
    REAL(8)             :: temp
    
  END SUBROUTINE
  
  SUBROUTINE generate_HP_HT_2D(P, T, DGT_IFE_partition, HP, HT, P_average)
  
    USE IFE_MAIN_PARAM
        
    REAL(8), DIMENSION(:,:), POINTER 	    ::	P, HP, HP_temp, P_average
    INTEGER, DIMENSION(:,:), POINTER	    ::	T, HT
    INTEGER, DIMENSION(:,:), POINTER      ::  DGT, DGT_IFE_partition
    
  END SUBROUTINE
  
  SUBROUTINE generate_HE_2D(HP, HT, ele_col_number, HE)
  
    REAL(8), DIMENSION(:,:), POINTER 	    ::	HP
    INTEGER, DIMENSION(:,:), POINTER	    ::	HT, HE
    INTEGER                               ::  ele_col_number
    
  END SUBROUTINE
  
  SUBROUTINE generate_DGP_DGT_2D(dimensions, n_nodes, P, T, DGP, DGT)
        
    REAL(8), DIMENSION(2,2), INTENT(IN) ::	dimensions
    INTEGER, DIMENSION(2), INTENT(IN)		::	n_nodes
    REAL(8), DIMENSION(:,:), POINTER 		::	P, DGP
		INTEGER, DIMENSION(:,:), POINTER		::	T, DGT
        
  END SUBROUTINE
  
  SUBROUTINE generate_DGE_2D(DGP, DGT, DGE, dimensions)
        
    REAL(8), DIMENSION(:,:), POINTER 	    ::	DGP
    INTEGER, DIMENSION(:,:), POINTER	    ::	DGT, DGE
    REAL(8), DIMENSION(2,2), INTENT(IN)		::	dimensions
        
  END SUBROUTINE
  
  SUBROUTINE generate_Adaptive_DG_P_T_E(HP, HT, HE, DGT_IFE, AHP, AHT, AHE, HT_flag, &
                                        CellMesh_Old, CellMesh_New, N_Tier_Old)
    
    !==========LY Add for Multigrid Store, 2022-1-17==========
    USE Cell_Data_2D
    INTEGER :: N_Tier_Old
    TYPE(CellDataType), DIMENSION(:), POINTER :: CellMesh_Old, CellMesh_New
    !==========LY Add for Multigrid Store, 2022-1-17==========
    REAL(8), DIMENSION(:,:), POINTER :: HP, AHP
    INTEGER, DIMENSION(:,:), POINTER :: HT, HE, DGT_IFE, AHT, AHE, HT_flag
  
  END SUBROUTINE
  
  SUBROUTINE generate_Adaptive_DG_P_T_E_repeat(DGP, DGT, DGE, DGT_IFE, repeat_index, ADGP, ADGT, ADGE, HT_flag, &
                                                CellMesh_Old, CellMesh_New, N_Tier_Old)
    
    !==========LY Add for Multigrid Store, 2022-1-18==========
    USE Cell_Data_2D
    INTEGER :: N_Tier_Old
    TYPE(CellDataType), DIMENSION(:), POINTER :: CellMesh_Old, CellMesh_New
    !==========LY Add for Multigrid Store, 2022-1-18==========
    REAL(8), DIMENSION(:,:), POINTER :: DGP, ADGP
    INTEGER, DIMENSION(:,:), POINTER :: DGT, DGE, DGT_IFE, ADGT, ADGE, HT_flag
    INTEGER                          :: repeat_index
    
  END SUBROUTINE
  
  SUBROUTINE generate_Adaptive_DG_P_average(HP, HP_average)
  
    !USE qsort_c_module_2D
    !USE ShellSort_2D
    
    REAL(8), DIMENSION(:,:), INTENT(IN)   :: HP
    REAL(8), DIMENSION(:,:), POINTER      :: HP_average
    
  END SUBROUTINE
  
  SUBROUTINE generate_DG_node_flag(P, T, P_average, P_flag)
    REAL(8), DIMENSION(:,:), POINTER :: P, P_flag, P_average
    INTEGER, DIMENSION(:,:), POINTER :: T
  END SUBROUTINE
  
  SUBROUTINE generate_E_Dirichlet_2D(HP, HT, E_Dirichlet)
  
    USE IFE_Boundary
    !=========LY modification for dealing with Floating Dirichlet boundary condition, 2022-7-25=========
    Use IFE_Data, Only: E_DirichletValue
    !=========LY modification for dealing with Floating Dirichlet boundary condition, 2022-7-25=========
    REAL(8), DIMENSION(:,:), POINTER 	    ::	HP
    INTEGER, DIMENSION(:,:), POINTER	    ::	HT, E_Dirichlet

  END SUBROUTINE
  
  SUBROUTINE generate_DG_edge_flag(DGE, dimensions, DGP, DGT, DG_flag)
        
    REAL(8), DIMENSION(:,:), POINTER 	    ::	DGP
    REAL(8), DIMENSION(2,2), INTENT(IN)		::	dimensions
    INTEGER, DIMENSION(:,:), POINTER	    ::	DGT, DGE, DG_flag
        
  END SUBROUTINE
  
  SUBROUTINE Update_Edge_Objects_Intersection_2D(DGE, DGP, DGT, information_2, &
                                                    node_index, element_index, objects, information_3, edge_index)
        
    USE Object_Data_2D
    USE IFE_MAIN_PARAM !LY REVISE, 2022-1-6, Float number error
    
    INTEGER, DIMENSION(:), POINTER          ::  edge_index, node_index, element_index
    REAL(8), DIMENSION(:,:), POINTER        ::	information_2, information_3, information_3_temp
    REAL(8), DIMENSION(:,:), POINTER 	      ::	DGP
    INTEGER, DIMENSION(:,:), POINTER	      ::	DGT, DGE
    INTEGER                                 ::  number_of_edges, i, j, N_Objects, xc, yc, r
    TYPE(ObjectType), DIMENSION(:), POINTER ::	objects
    
  END SUBROUTINE
  
  SUBROUTINE generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, start_nodes, &
                                    end_nodes, Gauss_coefficient_local_1D_Linear, Gauss_point_local_1D_Linear)
  
    REAL(8), DIMENSION(:), INTENT(IN)		                ::	Gauss_coefficient_reference_1D
    REAL(8), DIMENSION(:), INTENT(IN)                   ::	Gauss_point_reference_1D
    REAL(8), DIMENSION(:), INTENT(IN)                   ::	start_nodes, end_nodes	
    REAL(8), DIMENSION(:)                               ::	Gauss_coefficient_local_1D_Linear
    REAL(8), DIMENSION(:, :)                            ::	Gauss_point_local_1D_Linear
    
  END SUBROUTINE
  
  SUBROUTINE gauss_integration_local_stiffness_DG_FE_D(delta, function_coefficient,Gauss_coefficient_local_1D_linear, &
                                                        Gauss_point_local_1D_linear, this_flag1,vertices,trial_basis_type, &
                                                        trial_basis_index, test_basis_type, test_basis_index, &
                                                        con_penalty, int_value, element_Gauss)
    !=========LY modification, 2022-7-25=========
    Use IMPIC_Data_2D
    !=========LY modification, 2022-7-25=========
        REAL(8)                                           function_coefficient
        INTEGER                                           trial_basis_type, test_basis_type
        REAL(8)                                           int_value

        REAL(8),DIMENSION(:)		                ::    Gauss_coefficient_local_1D_linear
        REAL(8),DIMENSION(:,:)			            ::    Gauss_point_local_1D_linear
        REAL(8),DIMENSION(2,4)                      ::    vertices
        INTEGER,DIMENSION(4,2)                      ::    this_flag1
        REAL(8)                                     ::    con_penalty

        INTEGER                                     ::    trial_basis_index
        INTEGER                                     ::    test_basis_index
        INTEGER                                     ::	  delta
        Integer :: element_Gauss

  END SUBROUTINE
  
  SUBROUTINE gauss_integration_local_stiffness_IDG_D(delta, Global_Beta, information1_vector, information2_vector, &
                                                Gauss_coefficient_local_1D_linear,Gauss_point_local_1D_linear, &
	                                            this_flag1,vertices,trial_basis_type, trial_basis_index, &
											    test_basis_type, test_basis_index, location, con_penalty, int_value, element_Gauss)
        !=========LY modification, 2022-7-25=========
        Use IMPIC_Data_2D
        !=========LY modification, 2022-7-25=========
        USE Gauss_Data
        USE IFE_MAIN_PARAM !LY REVISE, 2022-1-6, Float number error
        
        INTEGER                     :: delta
        REAL(8), DIMENSION(:), POINTER	        	::	  Global_Beta
        REAL(8), DIMENSION(:)       :: Gauss_coefficient_local_1D_linear
        REAL(8), DIMENSION(:, :)    :: Gauss_point_local_1D_linear
        INTEGER, DIMENSION(4, 2)    :: this_flag1
        REAL(8), DIMENSION(2, 4)    :: vertices
        INTEGER                     :: trial_basis_type, trial_basis_index, test_basis_index, test_basis_type, location
        REAL(8)                     :: con_penalty, int_value
        INTEGER                     :: information1_vector(18)
        REAL(8)                     :: information2_vector(8)
        Integer :: element_Gauss


  END SUBROUTINE
  
  SUBROUTINE generate_stiffness_matrix_local_IDG_D( Global_Beta, node_type,num_of_unknowns, &
                                                information_1, information_2, information_3_D, HP, &
                                                element_index, edge_index_D,&
                                                HT, E_Dirichlet, HE, dimensions, &
                                                Gauss_coefficient_reference_1D, Gauss_point_reference_1D,&
                                                trial_basis_type, &
                                                test_basis_type, matrix, matrix_xt, delta, con_penalty)

        USE IFE_MAIN_PARAM
        
        INTEGER, DIMENSION(:,:), POINTER        ::	node_type
        REAL(8), DIMENSION(:), POINTER          ::  Global_Beta
        INTEGER                                 ::  num_of_unknowns
        INTEGER,DIMENSION(:,:),POINTER          ::  information_1
        REAL(8),DIMENSION(:,:),POINTER          ::  information_2, information_3_D
        REAL(8), DIMENSION(:,:), POINTER        ::  HP
        INTEGER,DIMENSION(:), POINTER           ::  element_index, edge_index_D
        INTEGER, DIMENSION(:,:), POINTER        ::  HT, E_Dirichlet, HE
        REAL(8)                                 ::  dimensions(2,2)
        REAL(8),DIMENSION(:)                    ::  Gauss_coefficient_reference_1D
        REAL(8),DIMENSION(:)                    ::  Gauss_point_reference_1D
        INTEGER                                 ::  trial_basis_type, test_basis_type
        TYPE(SPARSE), DIMENSION(:), POINTER     ::	matrix
        TYPE(SPARSE), DIMENSION(:), ALLOCATABLE	::	matrix_xt
        INTEGER, INTENT(IN)                     ::	delta
        REAL(8)                                 ::  con_penalty
    
  END SUBROUTINE
  
  SUBROUTINE gauss_integration_local_stiffness_DG_FE_O(delta, function_coefficient,Gauss_coefficient_local_1D_linear, &
                                                        Gauss_point_local_1D_linear, &
	                                                    this_flag1, vertices, brother_vertices, trial_basis_type, trial_basis_index, &
											            test_basis_type, test_basis_index, con_penalty, int_value, element_Gauss)
        !=========LY modification, 2022-5-20=========
        Use IMPIC_Data_2D
        !=========LY modification, 2022-5-20=========
        REAL(8)                                           function_coefficient
        INTEGER                                           trial_basis_type, test_basis_type
        REAL(8)                                           int_value

        REAL(8),DIMENSION(:)		                ::    Gauss_coefficient_local_1D_linear
        REAL(8),DIMENSION(:,:)			            ::    Gauss_point_local_1D_linear
        REAL(8),DIMENSION(2,4)                      ::    vertices, brother_vertices
        INTEGER,DIMENSION(4,2)                      ::    this_flag1
        REAL(8)                                     ::    con_penalty

        INTEGER                                     ::    trial_basis_index
        INTEGER                                     ::    test_basis_index
        INTEGER                                     ::	  delta
        Integer :: element_Gauss

  END SUBROUTINE
  
  SUBROUTINE gauss_integration_local_stiffness_IDG_O(delta, Global_Beta, information1_vector, information2_vector, &
                                                Gauss_coefficient_local_1D_linear,Gauss_point_local_1D_linear, &
	                                            this_flag1,vertices, trial_basis_type, trial_basis_index, &
											    test_basis_type, test_basis_index, location, con_penalty, int_value, element_Gauss)
        !=========LY modification, 2022-7-25=========
        Use IMPIC_Data_2D
        !=========LY modification, 2022-7-25=========
        USE Gauss_Data
        USE IFE_MAIN_PARAM !LY REVISE, 2022-1-6, Float number error
  
        INTEGER                     :: delta
        REAL(8), DIMENSION(:), POINTER	        	::	  Global_Beta
        REAL(8), DIMENSION(:)       :: Gauss_coefficient_local_1D_linear
        REAL(8), DIMENSION(:, :)    :: Gauss_point_local_1D_linear
        INTEGER, DIMENSION(4, 2)    :: this_flag1
        REAL(8), DIMENSION(2, 4)    :: vertices
        INTEGER                     :: trial_basis_type, trial_basis_index, test_basis_index, test_basis_type, location
        REAL(8)                     :: con_penalty, int_value
        INTEGER                     :: information1_vector(18)
        REAL(8)                     :: information2_vector(8)
        Integer :: element_Gauss


  END SUBROUTINE
  
  SUBROUTINE gauss_integration_local_stiffness_Mixed_O(delta, Global_Beta, information1_vector, information2_vector, &
                                                Gauss_coefficient_local_1D_linear,Gauss_point_local_1D_linear, &
	                                            this_flag1, vertices1, vertices2, trial_basis_type, trial_basis_index, &
											    test_basis_type, test_basis_index, location, this_basic, next_basic, &
                                                con_penalty, int_value, element_Gauss)
        !=========LY modification, 2022-4-25=========
        Use IMPIC_Data_2D
        !=========LY modification, 2022-4-25=========
        USE Gauss_Data
        USE IFE_MAIN_PARAM !LY REVISE, 2022-1-6, Float number error
  
        INTEGER                     :: delta
        REAL(8), DIMENSION(:), POINTER	        	::	  Global_Beta
        REAL(8), DIMENSION(:)       :: Gauss_coefficient_local_1D_linear
        REAL(8), DIMENSION(:, :)    :: Gauss_point_local_1D_linear
        INTEGER, DIMENSION(4, 2)    :: this_flag1
        REAL(8), DIMENSION(2, 4)    :: vertices1, vertices2
        INTEGER                     :: trial_basis_type, trial_basis_index, test_basis_index, &
                                        test_basis_type, location,  this_basic, next_basic
        REAL(8)                     :: con_penalty, int_value 
        INTEGER                     :: information1_vector(18)
        REAL(8)                     :: information2_vector(8)
        Integer :: element_Gauss


  END SUBROUTINE
  
  SUBROUTINE gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information1_vector, information2_vector, &
                                                information1_vector_next, information2_vector_next, &
                                                Gauss_coefficient_local_1D_linear,Gauss_point_local_1D_linear, &
	                                            this_flag1,vertices1, vertices2, trial_basis_type, trial_basis_index, &
											    test_basis_type, test_basis_index, location, con_penalty, int_value,element_Gauss_trial, element_Gauss_test)
        !=========LY modification, 2022-7-25=========
        Use IMPIC_Data_2D
        !=========LY modification, 2022-7-25=========
        USE Gauss_Data
        USE IFE_MAIN_PARAM !LY REVISE, 2022-1-6, Float number error
  
        INTEGER                     :: delta
        REAL(8), DIMENSION(:), POINTER	        	::	  Global_Beta
        REAL(8), DIMENSION(:)       :: Gauss_coefficient_local_1D_linear
        REAL(8), DIMENSION(:, :)    :: Gauss_point_local_1D_linear
        INTEGER, DIMENSION(4, 2)    :: this_flag1
        REAL(8), DIMENSION(2, 4)    :: vertices1, vertices2
        INTEGER                     :: trial_basis_type, trial_basis_index, test_basis_index, test_basis_type, &
                                        location,  this_basic, next_basic
        REAL(8)                     :: con_penalty, int_value 
        INTEGER                     :: information1_vector(18), information1_vector_next(18)
        REAL(8)                     :: information2_vector(8), information2_vector_next(8)
        Integer :: element_Gauss_trial, element_Gauss_test


  END SUBROUTINE
  
  SUBROUTINE generate_stiffness_matrix_local_IDG_O( Global_Beta, node_type,num_of_unknowns, &
                                                    information_1, information_2, information_3, HP, element_index, edge_index,&
                                                    HT, HE, dimensions, &
                                                    Gauss_coefficient_reference_1D, Gauss_point_reference_1D,&
                                                    trial_basis_type, &
                                                    test_basis_type, matrix, matrix_xt, delta, con_penalty)

        USE IFE_MAIN_PARAM
        
        INTEGER, DIMENSION(:,:), POINTER        ::	node_type
        REAL(8), DIMENSION(:), POINTER          ::  Global_Beta
        INTEGER                                 ::  num_of_unknowns
        INTEGER,DIMENSION(:,:),POINTER          ::  information_1
        REAL(8),DIMENSION(:,:),POINTER          ::  information_2, information_3
        REAL(8), DIMENSION(:,:), POINTER        ::  HP
        INTEGER,DIMENSION(:), POINTER           ::  element_index, edge_index
        INTEGER, DIMENSION(:,:), POINTER        ::  HT, HE
        REAL(8)                                 ::  dimensions(2,2)
        REAL(8),DIMENSION(:)                    ::  Gauss_coefficient_reference_1D
        REAL(8),DIMENSION(:)                    ::  Gauss_point_reference_1D
        INTEGER                                 ::  trial_basis_type, test_basis_type
        TYPE(SPARSE), DIMENSION(:), POINTER     ::	matrix
        TYPE(SPARSE), DIMENSION(:), ALLOCATABLE	::	matrix_xt
        INTEGER, INTENT(IN)                     ::	delta
        REAL(8)                                 ::  con_penalty
    
  END SUBROUTINE
  
  SUBROUTINE gauss_integration_local_vector_DG_FE_D(delta, function_coefficient, Gauss_coefficient_local_1D_linear, &
                                                        Gauss_point_local_1D_linear, this_flag1, vertices, test_basis_type, & 
                                                        test_basis_index, con_penalty, boundary_value, int_value, element_Gauss)
        USE Gauss_Data
        Use IMPIC_Data_2D
        
        INTEGER                 :: delta
        REAL(8)                 :: function_coefficient
        REAL(8), DIMENSION(:)   :: Gauss_coefficient_local_1D_linear
        REAL(8), DIMENSION(:,:) :: Gauss_point_local_1D_linear
        INTEGER, DIMENSION(4,2) :: this_flag1
        REAL(8), DIMENSION(2,4) :: vertices
        INTEGER                 :: test_basis_type, test_basis_index
        Real(8) :: boundary_value
        Integer :: element_Gauss
        REAL(8)                 :: con_penalty, int_value
     
  END SUBROUTINE
  
  SUBROUTINE gauss_integration_local_vector_IDG_D(delta, Global_Beta, information1_vector, information2_vector, &
                                                    Gauss_coefficient_local_1D_linear, Gauss_point_local_1D_linear, &
                                                    this_flag1, vertices, test_basis_type, test_basis_index, location, & 
                                                    con_penalty, boundary_value, int_value, element_Gauss)
    
      USE Gauss_Data
      Use IMPIC_Data_2D
      USE IFE_MAIN_PARAM !LY REVISE, 2022-1-6, Float number error
      
      INTEGER                        :: delta
      REAL(8), DIMENSION(:), POINTER :: Global_Beta
      INTEGER                        :: information1_vector(18)
      REAL(8)                        :: information2_vector(8)
      REAL(8), DIMENSION(:)          :: Gauss_coefficient_local_1D_linear
      REAL(8), DIMENSION(:, :)       :: Gauss_point_local_1D_linear
      INTEGER, DIMENSION(4, 2)       :: this_flag1
      REAL(8), DIMENSION(2, 4)       :: vertices
      INTEGER                        :: test_basis_type, test_basis_index, location
      Real(8) :: boundary_value
      Integer :: element_Gauss
      REAL(8)                        :: con_penalty, int_value
      
  END SUBROUTINE
  
  SUBROUTINE generate_load_vector_local_IDG_D(Global_Beta, node_type, num_of_unknowns, &
                                            information_1, information_2, information_3, HP, element_index, edge_index, &
                                            HT, HE, &
                                            Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
                                            test_basis_type, delta, con_penalty, vector) 
        USE IFE_MAIN_PARAM
        USE Gauss_Data
        !=========LY modification, 2022-7-25=========
        Use IFE_Data, Only: E_DirichletValue
        !=========LY modification, 2022-7-25=========
  
        REAL(8), DIMENSION(:), POINTER      :: Global_Beta
        INTEGER, DIMENSION(:,:), POINTER    :: node_type
        INTEGER, INTENT(IN)                 :: num_of_unknowns
        INTEGER, DIMENSION(:,:), POINTER    :: information_1
        REAL(8), DIMENSION(:,:), POINTER    :: information_2, information_3
        REAL(8), DIMENSION(:,:), POINTER    :: HP
        INTEGER, DIMENSION(:,:), POINTER    :: HT, HE
        INTEGER, DIMENSION(:), POINTER      :: element_index, edge_index
        REAL(8), DIMENSION(:)               :: Gauss_coefficient_reference_1D
        REAL(8), DIMENSION(:)               :: Gauss_point_reference_1D
        INTEGER                             :: test_basis_type
        INTEGER, INTENT(IN)                 :: delta
        REAL(8)                             :: con_penalty
        REAL(8), DIMENSION(:), POINTER      :: vector
        
  END SUBROUTINE
     
  SUBROUTINE generate_Tier_initial_2D(P, T, CellMesh)
      USE Cell_Data_2D
      
      REAL(8), DIMENSION(:,:), POINTER :: P
      INTEGER, DIMENSION(:,:), POINTER :: T
      TYPE(CellDataType), DIMENSION(:), POINTER :: CellMesh
  END SUBROUTINE
  
  !=========LY modification for Periodic Boundary Condition,2022-7-25=========
  Subroutine Periodic_pair_boundary_nodes(nnx, nny, bc_index, p_basic, bc_point_1, bc_point_2, PairNodes)
  
    Integer :: nnx, nny
    Integer, Dimension(:), Pointer :: bc_index
    Real(8), Dimension(:,:), Pointer :: p_basic
    Real(8), Dimension(:,:), Pointer :: bc_point_1, bc_point_2
    Integer, Dimension(:,:), Pointer :: PairNodes
  
  End Subroutine
  !=========LY modification for Periodic Boundary Condition,2022-7-25=========
  !================================================================================================== 
END INTERFACE

! ------- PPR part----------
Interface
    SUBROUTINE BRINV(A,N,L,IS,JS)
	    !DIMENSION A(N,N),IS(N),JS(N)
	    !DOUBLE PRECISION A,T,D
        IMPLICIT NONE  
        REAL(8), DIMENSION(N,N) :: A
        INTEGER                    N, L, I, J, K
        INTEGER, DIMENSION(N)   :: IS, JS
        REAL(8)                    T, D  
    END SUBROUTINE
    
    SUBROUTINE Bubble_Sort_new(A,N) 
        IMPLICIT NONE
        Integer :: N
        Integer :: A(N)
        Integer :: i, j
        Integer :: insertVal, insertIndex
    END SUBROUTINE
  
    SUBROUTINE Linear_FE_Basis_Coeff_triangle_2D(vert, coef)
    USE IFE_MAIN_PARAM
    REAL(8), INTENT(IN)		                ::	vert(2,3)
    REAL(8), INTENT(OUT)	                ::	coef(3,3)
    END SUBROUTINE
    
    SUBROUTINE   Generate_Gauss_local_ppr( Gauss_coefficient_reference, Gauss_point_reference, &
                                      left_lower_point, h_partition, Gauss_coefficient_local, Gauss_point_local)
	    REAL(8),DIMENSION(:)		                ::    Gauss_coefficient_reference
	    REAL(8),DIMENSION(:,:)		                ::    Gauss_point_reference
	    REAL(8)                           			::    left_lower_point(2)
        REAL(8)	                                          h_partition(2)
        REAL(8),DIMENSION(:)		                ::    Gauss_coefficient_local
        REAL(8),DIMENSION(:,:)				        ::    Gauss_point_local
    END SUBROUTINE 

End Interface
!--------------------------


! Mesh Generation
! =======
INTERFACE
	SUBROUTINE Cubic_Partition_2D(dimensions, n_nodes, P, T)

		REAL(8), DIMENSION(2,2), INTENT(IN)		::	dimensions
		INTEGER, DIMENSION(2), INTENT(IN)		::	n_nodes
		REAL(8), DIMENSION(:,:), POINTER 		::	P
		INTEGER, DIMENSION(:,:), POINTER		::	T
		
	END SUBROUTINE

	SUBROUTINE Boundary_Condition_info_2D  ( n_grid, t_basic, p_basic, e_basic)

         INTEGER, DIMENSION(:,:), POINTER			::	t_basic
         INTEGER, DIMENSION(:,:), POINTER	        ::	e_basic
         REAL(8), DIMENSION(:,:), POINTER			::	p_basic
		 INTEGER                                    ::  n_grid
	END SUBROUTINE

	SUBROUTINE	modify_interface_elements(h, t_c, p_basic, element_index, information_1, information_2,	&
										modified_element_index, modified_information_1, modified_information_2)
		REAL(8), INTENT(IN)								::	h(2)
		INTEGER, DIMENSION(:,:), INTENT(IN)				::	t_c
		REAL(8), DIMENSION(:,:), INTENT(IN)				::	p_basic
		INTEGER, DIMENSION(:), INTENT(IN)				::	element_index
		INTEGER, DIMENSION(:,:), INTENT(IN)				::  information_1
		REAL(8), DIMENSION(:,:), INTENT(IN)				::  information_2
		INTEGER, DIMENSION(:), POINTER					::	modified_element_index
		INTEGER, DIMENSION(:,:), POINTER				::  modified_information_1
		REAL(8), DIMENSION(:,:), POINTER				::  modified_information_2
	END SUBROUTINE

	SUBROUTINE	Mesh_Objects_Intersection_info_2D_BL(t_basic, p_basic, 	&
													objects, Int_El_Frac,	&
													element_index, information_1, information_2, &
													node_index_el)

		USE Object_Data_2D
		IMPLICIT NONE

		REAL(8),				 INTENT(IN)						::	Int_El_Frac
		INTEGER, DIMENSION(:,:), INTENT(IN)						::	t_basic
		REAL(8), DIMENSION(:,:), INTENT(IN)						::	p_basic
		TYPE(ObjectType), DIMENSION(:), INTENT(IN)				::	objects
		INTEGER, DIMENSION(:,:), POINTER						::	p_int_f_tmp
		INTEGER, DIMENSION(:), POINTER							::	element_index
		INTEGER, DIMENSION(:,:), POINTER			            ::  information_1
		REAL(8), DIMENSION(:,:), POINTER			            ::  information_2
		INTEGER, DIMENSION(:), INTENT(IN)						::	node_index_el
	END SUBROUTINE

	SUBROUTINE El_Object_Intersection_info_2D_BL(vert, object, P_intrs, el_type, &
	                                                el_region, information_1_el, &
	                                                information_2_el, node_index_el)
			
		USE Object_Data_2D
		REAL(8), DIMENSION(2,4), INTENT(IN)			::	vert
		TYPE(ObjectType), INTENT(IN)				::	object
		REAL(8), DIMENSION(4,6)						::	P_intrs
		INTEGER											el_region, el_type
		INTEGER,DIMENSION(:),POINTER                ::  information_1_el
		REAL(8),DIMENSION(:),POINTER                ::  information_2_el
		INTEGER, DIMENSION(:), INTENT(IN)			::	node_index_el

	END SUBROUTINE

	SUBROUTINE Load_IFE_Mesh_2D(	IFE_Mesh_Filename,&
							t_c, p_basic, e_basic,  &
                            n_elements,information_1,information_2,information_3,information_3_D,&
                            element_index, Node_Index, edge_index, edge_index_D, &
                            HP, HT, HE, E_Dirichlet, P_average, P_flag)

		CHARACTER(*)							IFE_Mesh_Filename
		INTEGER, DIMENSION(:,:), POINTER	::	t_c, e_basic, HE, HT, E_Dirichlet
		REAL(8), DIMENSION(:,:), POINTER	::	p_basic, HP, P_average, P_flag
		INTEGER                             ::  n_elements
  		INTEGER, DIMENSION(:,:), POINTER    ::  information_1
        REAL(8), DIMENSION(:,:), POINTER    ::  information_2 ,information_3, information_3_D
		INTEGER, DIMENSION(:), POINTER		::	element_index
		INTEGER, DIMENSION(:), POINTER		::	Node_Index, edge_index, edge_index_D

    END SUBROUTINE  

	SUBROUTINE Partition_Cube_2D(	p_basic, e_basic, dimensions, 	&
								node_type, num_of_unknowns)

		REAL(8), DIMENSION(:,:), POINTER	::	p_basic
        INTEGER, DIMENSION(:,:), POINTER		::	e_basic
		REAL(8), DIMENSION(2,2), INTENT(IN)		::	dimensions
		INTEGER, DIMENSION(:,:), POINTER	::	node_type
		INTEGER, INTENT(OUT)					::	num_of_unknowns

	END SUBROUTINE

	SUBROUTINE Get_Boundary_Elements_2D( n_elements ,t_basic, node_type, node_loc, bnd_elem_index)

		INTEGER, DIMENSION(:,:), POINTER		::	t_basic, node_type
		INTEGER, DIMENSION(:), POINTER			::	node_loc
		INTEGER, DIMENSION(:), POINTER			::	bnd_elem_index
		INTEGER                                 ::  n_elements 

	END	SUBROUTINE

END INTERFACE

! IFE Assembler
! =======
INTERFACE
	SUBROUTINE Generate_Gauss_reference_triangle( Gauss_point_number, &
	                                                Gauss_coefficient_reference_triangle, &
	                                                Gauss_point_reference_triangle)

		IMPLICIT NONE
        INTEGER                                              Gauss_point_number
		REAL(8),DIMENSION(:)							::   Gauss_coefficient_reference_triangle
		REAL(8),DIMENSION(:,:)							::   Gauss_point_reference_triangle

	END SUBROUTINE

	SUBROUTINE Generate_Gauss_reference( Gauss_point_number, Gauss_coefficient_reference, Gauss_point_reference)
		
        INTEGER                                              Gauss_point_number
		REAL(8),DIMENSION(:)		                   ::    Gauss_coefficient_reference
		REAL(8),DIMENSION(:,:)			               ::    Gauss_point_reference

	END SUBROUTINE

	SUBROUTINE Generate_stiffness_matrix_local_IFE( minus_coefficient_function_name,Global_Beta, node_type, num_of_unknowns,&
                                                information_1, information_2, p_basic, element_index, &
									            t_c, nnx, nny, dimensions, &
												Gauss_coefficient_reference, &
												Gauss_point_reference, Gauss_coefficient_reference_triangle, &
												Gauss_point_reference_triangle, &
                                                trial_basis_type, &
												trial_derivative_degree_x,trial_derivative_degree_y, &
												test_basis_type,test_derivative_degree_x,test_derivative_degree_y, matrix, matrix_xt, delta)

		USE IFE_MAIN_PARAM
		INTEGER,DIMENSION(:,:),POINTER              ::    information_1
		REAL(8),DIMENSION(:,:),POINTER              ::    information_2
		INTEGER,DIMENSION(:), POINTER	        	::	  element_index
		REAL(8),DIMENSION(:)						::    Gauss_coefficient_reference
		REAL(8),DIMENSION(:,:)						::    Gauss_point_reference
		REAL(8),DIMENSION(:)						::    Gauss_coefficient_reference_triangle
		REAL(8),DIMENSION(:,:)		                ::    Gauss_point_reference_triangle
		REAL(8),DIMENSION(:,:),POINTER              ::    r
		INTEGER, DIMENSION(:,:), POINTER			::    t_c
		REAL(8)                                           minus_coefficient_function_name
		REAL(8)	                                          dimensions(2,2)
		INTEGER                                           nnx,nny
		INTEGER                                           trial_basis_type, trial_derivative_degree_x,trial_derivative_degree_y, &
														  test_basis_type,test_derivative_degree_x,test_derivative_degree_y
		REAL(8), DIMENSION(:,:), POINTER	        ::	  p_basic
		TYPE(SPARSE), DIMENSION(:), POINTER	        ::	  matrix
        TYPE(SPARSE), DIMENSION(:), ALLOCATABLE	    ::	matrix_xt
		INTEGER, DIMENSION(:,:), POINTER	        ::	  node_type
		INTEGER, INTENT(IN)					        ::	  num_of_unknowns, delta
		REAL(8), DIMENSION(:), POINTER	        	::	  Global_Beta

	END SUBROUTINE

	SUBROUTINE Classify_Nodes_2D(	p_basic, node_type,	node_loc, node_permute )

		REAL(8), DIMENSION(:,:), POINTER	::	p_basic
		INTEGER, DIMENSION(:), POINTER		::	node_loc
		INTEGER, DIMENSION(:,:), POINTER	::	node_type, node_permute

	END SUBROUTINE

	SUBROUTINE Global_RHS_EBC_2D (	information_1,information_2, element_index, Coeff_FUN, U_full, p_basic, t_basic_int, e_basic,	&
							    el_type1, p_int_x, p_int_y, beta, node_type,	&
							    bnd_elem_index, Stiff_HW, vector, delta)

		EXTERNAL								Coeff_FUN
		INTEGER, DIMENSION(:), POINTER		::  element_index
		INTEGER, DIMENSION(:), POINTER		::	bnd_elem_index, el_type1
		INTEGER, DIMENSION(:,:), POINTER	    ::	t_basic_int, node_type, e_basic
		REAL(8), DIMENSION(:), POINTER		::	U_full, vector, beta
		REAL(8), DIMENSION(:,:), POINTER 	::	p_basic, p_int_x, p_int_y
		REAL(8)									Stiff_HW(2,4,4)
		INTEGER                                 delta
		INTEGER, DIMENSION(:,:),POINTER             ::  information_1
		REAL(8), DIMENSION(:,:), POINTER            ::  information_2
	
    END SUBROUTINE
                                
    !$=============== ab.ZWZ 2021/7/9 ====================\\
    SUBROUTINE Global_RHS_EBC_2D_bjw(A_stiff_xt, p_basic, t_basic_int, e_basic, node_type, U_full, bnd_elem_index, rhs_EBC, &
																		delta)
  
        USE IFE_MAIN_PARAM
        IMPLICIT NONE
        TYPE(SPARSE), DIMENSION(:), POINTER	::	A_stiff_xt
        REAL(8), DIMENSION(:,:), POINTER	::	p_basic 
        INTEGER, DIMENSION(:), POINTER		::	bnd_elem_index
        INTEGER, DIMENSION(:,:), POINTER	::	node_type, t_basic_int, e_basic
        REAL(8), DIMENSION(:), POINTER		::	U_full, rhs_EBC
				!=========LY modification, 2022-4-26=========
				Integer, Intent(In)                 ::  delta
				!=========LY modification, 2022-4-26=========

        INTEGER		num_bound_elem, n_nodes_in_elem, i, j	
			
        REAL(8), DIMENSION(:,:), POINTER			::	EBC_Value
        REAL(8), DIMENSION(:), ALLOCATABLE      ::	EBC_Value_xt
	
	END SUBROUTINE
    !$=============== ab.ZWZ 2021/7/9 ====================//            
    !$=============== ab.ZWZ 2021/7/9 ====================//                     

	SUBROUTINE Global_RHS_NBC_2D_dj (element_index, el_type1, information_1, information_2, p_basic, t_basic_int, e_basic,	&
					  	        p_int_x, p_int_y,  beta,	&
						        node_type, bnd_elem_index, vector, delta)
	
		REAL(8), DIMENSION(:,:), POINTER	::	p_basic, p_int_x, p_int_y
		INTEGER, DIMENSION(:,:), POINTER	::	t_basic_int, node_type, e_basic
		INTEGER, DIMENSION(:), POINTER		::	bnd_elem_index
		REAL(8), DIMENSION(:), POINTER		::	vector, beta
		INTEGER  								delta
		INTEGER, DIMENSION(:,:), POINTER             ::  information_1
		REAL(8), DIMENSION(:,:), POINTER             ::  information_2
		INTEGER, DIMENSION(:), POINTER		         ::  el_type1
		INTEGER, DIMENSION(:), POINTER		         ::  element_index

	END SUBROUTINE

	SUBROUTINE   Generate_Gauss_local( Gauss_coefficient_reference, Gauss_point_reference, &
                                   left_lower_point, h_partition, Gauss_coefficient_local, Gauss_point_local)

		REAL(8),DIMENSION(:)		                ::    Gauss_coefficient_reference
		REAL(8),DIMENSION(:,:)		                ::    Gauss_point_reference
		REAL(8)                           			::    left_lower_point(2)
		REAL(8)	                                          h_partition(2)
		REAL(8),DIMENSION(:)		                ::    Gauss_coefficient_local
		REAL(8),DIMENSION(:,:)		                ::    Gauss_point_local
	END SUBROUTINE

	SUBROUTINE Inter_Tetra_Partition_2D(tri_sign, pointer_reference_to_local, el_flag, vert, Intrs_pts, t_int_pt, p_int_pt)

		REAL(8), DIMENSION(:,:), INTENT(IN)	::	vert, Intrs_pts
		INTEGER, DIMENSION(:), INTENT(IN)	::	tri_sign(2)
		REAL(8), DIMENSION(:,:), POINTER 	::	p_int_pt
		INTEGER, DIMENSION(:,:), POINTER	::	t_int_pt
		INTEGER									el_flag
		INTEGER                                 pointer_reference_to_local(4)

	END SUBROUTINE

	SUBROUTINE Check_Sub_Nodes_2D(Subvert1, Subvert2, BCvert1, BCvert2, IN1, IN2 )

		REAL(8), INTENT(IN)		::	Subvert1(2), Subvert2(2), BCvert1(2), BCvert2(2)
		INTEGER, INTENT(OUT)	::	IN1, IN2

	END SUBROUTINE

	SUBROUTINE Global_RHS_PDE_2D(	RHS_FUN,RHS_FUN1, U_full, R_full, p_basic, t_c1,		&
							    h_partition, node_type1, num_of_unknowns, RHS_HW, Interpol_HW, vector, delta,element_index,information_1,information_2)

		EXTERNAL								RHS_FUN,RHS_FUN1
		REAL(8), DIMENSION(:,:), POINTER	::	p_basic
		REAL(8), DIMENSION(:),	INTENT(IN)	::	U_full, R_full
		REAL(8), DIMENSION(:),	POINTER		::	beta
		INTEGER, DIMENSION(:,:), POINTER	::	t_c1, node_type1
		INTEGER, INTENT(IN)					::	num_of_unknowns, delta
		REAL(8), DIMENSION(:), POINTER		::	vector
		REAL(8)									RHS_HW(2,4,4), Interpol_HW(4,4)
		REAL(8)	                                          h_partition(2)
		INTEGER, DIMENSION(:), POINTER				::	  element_index
		INTEGER,DIMENSION(:,:),POINTER              ::    information_1
		REAL(8),DIMENSION(:,:),POINTER              ::    information_2

	END SUBROUTINE

	SUBROUTINE	Global_Mass( RHS_FUN, U_full, R_full, p_basic, t_basic_int, &
    					 h_partition, node_type, num_of_unknowns, Mass_HW, Interpol_HW, matrix, element_index,information_1,information_2)
		
		USE IFE_MAIN_PARAM
		IMPLICIT NONE
		EXTERNAL								  RHS_FUN
		REAL(8), DIMENSION(:),	INTENT(IN)	::	  U_full, R_full
		INTEGER, DIMENSION(:,:), POINTER	::	  t_basic_int, node_type
		INTEGER, INTENT(IN)					::	  num_of_unknowns
		TYPE(SPARSE), DIMENSION(:), POINTER	::	  matrix
		REAL(8)								::	  Mass_HW(4,4,4), Interpol_HW(4,4)
		REAL(8)	                                  h_partition(2)
		INTEGER,DIMENSION(:,:),POINTER      ::    information_1
		REAL(8),DIMENSION(:,:),POINTER      ::    information_2
		REAL(8), DIMENSION(:,:), POINTER	::	  p_basic
		INTEGER, DIMENSION(:), POINTER		::	  element_index

	END SUBROUTINE

	SUBROUTINE Gauss_Nodes_Elem_2D(vert, gnodes)

		REAL(8), INTENT(IN)		::	vert(2,3)
		REAL(8), INTENT(OUT)	::	gnodes(2,3)

	END SUBROUTINE

	SUBROUTINE Global_HardWire_2D(	Coeff_FUN, p_basic, t_basic_int,element_index	,				&
							    G_Stiff_HW, G_RHS_HW, G_Mass_HW, G_Interpol_HW, delta)

		EXTERNAL								Coeff_FUN
		REAL(8), DIMENSION(:,:), POINTER	::	p_basic
		INTEGER, DIMENSION(:,:), POINTER	::	t_basic_int
		REAL(8)									G_Stiff_HW(2,4,4), G_RHS_HW(2,4,4), &
												G_Mass_HW(2,4,4), G_Interpol_HW(4,4)
		INTEGER                             ::	delta
		INTEGER, DIMENSION(:), POINTER		::	element_index

	END SUBROUTINE

	SUBROUTINE IFE_RHS_2D(	RHS_FUN, U_val, R_val, vert, Intrs_pts, &
					    el_region, n_nodes_in_elem_1, vector, delta, information_vector_1, information_vector_2)

		EXTERNAL					RHS_FUN
		REAL(8), INTENT(IN)		::	vert(2,4)
		INTEGER, INTENT(IN)		::	el_region(2), n_nodes_in_elem_1, delta
		REAL(8), INTENT(IN)		::	U_val(n_nodes_in_elem_1), R_val(n_nodes_in_elem_1)	
		REAL(8), INTENT(OUT)	::	vector(4)
		INTEGER                     information_vector_1(18)
		REAL(8)                     information_vector_2(6), Intrs_pts(2,2)

	END SUBROUTINE

	SUBROUTINE FE_RHS_HW_2D( vert, gwght, gnodes, n_nodes_in_elem_1, hardwire, delta)

		REAL(8), INTENT(IN)						::	vert(2,4), gnodes(2,4),gwght(4)
		INTEGER, INTENT(IN)						::	n_nodes_in_elem_1, delta
		REAL(8), INTENT(OUT)					::	hardwire(4,4)

	END SUBROUTINE

	SUBROUTINE IFE_Interpol_2D(vert, u_val, x_int, y_int, &
						   Intrs_pts, el_type, el_beta, sub_region_ind, dind, f)

		REAL(8), INTENT(IN)					::	vert(2,3), el_beta(2)
		REAL(8), DIMENSION(:,:), INTENT(IN)	::	Intrs_pts
		REAL(8), DIMENSION(:), INTENT(IN)	::	u_val, x_int, y_int
		INTEGER, INTENT(IN)					::	dind(2), el_type, sub_region_ind
		REAL(8), DIMENSION(:), INTENT(OUT)	::	f

	END SUBROUTINE

	SUBROUTINE	FE_Mass(RHS_FUN, U_val, R_val, el_region, &
					n_nodes_in_elem_1, n_nodes_in_elem_2, mass_hw, interpol_hw, matrix)

		EXTERNAL					RHS_FUN
		INTEGER, INTENT(IN)		::	el_region,	&
									n_nodes_in_elem_1, n_nodes_in_elem_2
		REAL(8), INTENT(IN)		::	U_val(n_nodes_in_elem_1), R_val(n_nodes_in_elem_1)
		REAL(8), INTENT(OUT)	::	matrix(4,4)
		REAL(8), INTENT(IN)		::	mass_hw(4,4,4), interpol_hw(4,4)
	END SUBROUTINE

	SUBROUTINE	IFE_Mass(	RHS_FUN, U_val, R_val, vert, Intrs_pts,					&
							el_region, n_nodes_in_elem_1, n_nodes_in_elem_2, matrix, information_vector_1, information_vector_2)

		EXTERNAL				RHS_FUN
		REAL(8), INTENT(IN)		::	vert(2,4) !, el_beta(2)
		INTEGER, INTENT(IN)		::	n_nodes_in_elem_1, n_nodes_in_elem_2,	&
				                    el_region(2)
		REAL(8), INTENT(IN)		::	Intrs_pts(2, 2)
		REAL(8), INTENT(IN)		::	U_val(n_nodes_in_elem_1), R_val(n_nodes_in_elem_1)
		REAL(8), INTENT(OUT)	::	matrix(4,4)
		INTEGER                     information_vector_1(18)
		REAL(8)                     information_vector_2(6)

	END SUBROUTINE


END INTERFACE

! Linear solver
! ========
INTERFACE
	SUBROUTINE Diagonal_Preconditioner(K_stiff, K_diag)

		USE IFE_MAIN_PARAM
		TYPE(SPARSE), DIMENSION(:), INTENT(IN)	::	K_stiff
		REAL(8), DIMENSION(:), POINTER			::	K_diag

	END SUBROUTINE

	SUBROUTINE Sparse_Structure(t_basic_int, HE, node_type, K_VROW, SMatrix, NZ)

		USE IFE_MAIN_PARAM
		INTEGER, DIMENSION(:,:), INTENT(IN)		::	t_basic_int, node_type
        INTEGER, DIMENSION(:,:), POINTER        ::	HE
		INTEGER, DIMENSION(:,:), POINTER		::	K_VROW
		TYPE(SPARSE), DIMENSION(:), POINTER		::	SMatrix
		INTEGER, INTENT(OUT)					::	NZ

    END SUBROUTINE
    
    SUBROUTINE Sparse_Structure_xt(t_basic_int, HE, node_type, K_VROW, SMatrix, NZ)

		USE IFE_MAIN_PARAM
		INTEGER, DIMENSION(:,:), INTENT(IN)		::	t_basic_int, node_type
        INTEGER, DIMENSION(:,:), POINTER        ::	HE
		INTEGER, DIMENSION(:,:), POINTER		::	K_VROW
		TYPE(SPARSE), DIMENSION(:), ALLOCATABLE		::	SMatrix
		INTEGER, INTENT(OUT)					::	NZ

	END SUBROUTINE
    
	SUBROUTINE MY_PCG(ANEL, A, AI, AJ, KNEL, K, KI, KJ, N, R, X, RELTOL, MAXITER, ITER, IERR, MATVEC, MSOLVE)

		INTEGER								ANEL, KNEL, N 
		REAL(8), DIMENSION(:)			::	A, K, R, X
		INTEGER, DIMENSION(:)			::	AI, AJ, KI, KJ
		REAL(8)								RELTOL
		INTEGER								MAXITER, ITER, IERR
		EXTERNAL							MATVEC, MSOLVE

    END SUBROUTINE
    
    !$=============== ab.ZWZ 2021/7/9 ====================\\
    SUBROUTINE MY_JPCG_Solver(matrix, R, X, M, NZ, N, ITER)
        USE IFE_MAIN_PARAM   
        TYPE(SPARSE), DIMENSION(NZ), TARGET		::	matrix
        REAL(8), DIMENSION(N)					::	X
        REAL(8), DIMENSION(N)					::	R
        REAL(8), DIMENSION(N), TARGET			::	M
        INTEGER									::	NZ, N
        INTEGER, INTENT(OUT)					::	ITER
    END SUBROUTINE
    !$=============== ab.ZWZ 2021/7/9 ====================//

END INTERFACE

! Integration
! ======
INTERFACE
	
	SUBROUTINE Gauss_integration_local_stiffness(delta, coefficient_function_name,Gauss_coefficient_local,Gauss_point_local, &
	                                         left_lower_point,h_partition,trial_basis_type, trial_basis_index, &
											 trial_derivative_degree_x,trial_derivative_degree_y, test_basis_type, &
											 test_basis_index,test_derivative_degree_x,test_derivative_degree_y, r, element_Gauss)

		REAL(8)                                           coefficient_function_name
		INTEGER                                           trial_basis_type, trial_derivative_degree_x,trial_derivative_degree_y, &
														  test_basis_type,test_derivative_degree_x,test_derivative_degree_y
		REAL(8)								        	  r
		REAL(8)	                                          h_partition(2)
		REAL(8),DIMENSION(:)						::    Gauss_coefficient_local
		REAL(8),DIMENSION(:,:)		                ::    Gauss_point_local
		REAL(8)                           			::    left_lower_point(2)
		INTEGER                                           trial_basis_index
		INTEGER                                           test_basis_index
		INTEGER                                     ::	  delta
    Integer :: element_Gauss

	END SUBROUTINE

	SUBROUTINE  Gauss_integration_local_stiffness_IFE(delta, Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
			                                      vertices,information_vector_1,information_vector_2,trial_basis_type, &
												  trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y,test_basis_type, &
												  test_basis_index,test_derivative_degree_x,test_derivative_degree_y,r, element_Gauss)
		
		REAL(8),DIMENSION(:)					    ::    Gauss_coefficient_reference_triangle
		REAL(8),DIMENSION(:,:)		                ::    Gauss_point_reference_triangle
		INTEGER                                           information_vector_1(18)
		REAL(8)                                           information_vector_2(8)
		REAL(8)                                           vertices(2,4)
		INTEGER                                           trial_basis_type, trial_derivative_degree_x,trial_derivative_degree_y, &
														  test_basis_type,test_derivative_degree_x,test_derivative_degree_y
		INTEGER                                           trial_basis_index
		INTEGER                                           test_basis_index
		REAL(8)								        	  r
		INTEGER                                     ::	  delta
    Integer :: element_Gauss

	END SUBROUTINE

	SUBROUTINE  Retangular_local_basis( x, y, left_lower_point, h_partition, basis_type, &
	                                    basis_index, derivative_degree_x, derivative_degree_y, r)

		REAL(8)                                               x, y
		REAL(8)	                                              h_partition(2) 
		INTEGER                                               basis_type, basis_index, derivative_degree_x, derivative_degree_y
		REAL(8)                           	        	      left_lower_point(2)  
		REAL(8)                                               r

	END SUBROUTINE

	SUBROUTINE  Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                      vertices_triangle,Gauss_coefficient_local_triangle,Gauss_point_local_triangle) 

		REAL(8),DIMENSION(:),INTENT(IN)                ::    Gauss_coefficient_reference_triangle
		REAL(8),DIMENSION(:,:),INTENT(IN)              ::    Gauss_point_reference_triangle
		REAL(8),DIMENSION(:)			               ::    Gauss_coefficient_local_triangle
		REAL(8),DIMENSION(:,:)						   ::    Gauss_point_local_triangle
		REAL(8)                                           vertices_triangle(2,3)

	END SUBROUTINE

	SUBROUTINE  Retangular_local_basis_IFE(x,y,vertices,information_vector_1,information_vector_2,piece_flag,basis_type, &
		                               basis_index,derivative_degree_x,derivative_degree_y,r,element_Gauss) 

		REAL(8)                                           x, y
		INTEGER                                           information_vector_1(18)
		REAL(8)                                           information_vector_2(8)
		REAL(8)                                           vertices(2,4)
		INTEGER                                           basis_type, basis_index, derivative_degree_x, derivative_degree_y  
		REAL(8)								        	  r
		INTEGER                                           piece_flag
    Integer :: element_Gauss
	
	END SUBROUTINE

	SUBROUTINE FE_Stiff_Eval_2D( el_beta, n_nodes_in_elem_1, n_nodes_in_elem_2, stiff_hw, matrix )

		INTEGER, INTENT(IN)						::	n_nodes_in_elem_1, n_nodes_in_elem_2
		REAL(8), INTENT(IN)						::	el_beta, stiff_hw(4,4)
		REAL(8), INTENT(OUT)					::	matrix(4,4)

	END SUBROUTINE

	SUBROUTINE FE_Stiff_HW_2D( Coeff_FUN, vert, n_nodes_in_elem_1, n_nodes_in_elem_2, hardwire, delta )

		EXTERNAL									Coeff_FUN
		REAL(8), INTENT(IN)						::	vert(2,4)
		INTEGER, INTENT(IN)						::	n_nodes_in_elem_1, n_nodes_in_elem_2
		REAL(8), INTENT(OUT)					::	hardwire(4,4)
		INTEGER                                 ::  delta               
                 
	END SUBROUTINE

	SUBROUTINE IFE_Stiff_2D(	tri_sign, pointer_reference_to_local, Coeff_FUN, vert, Intrs_pts, el_type, el_beta, &
					  	    n_nodes_in_elem_1, n_nodes_in_elem_2, matrix, delta )

		EXTERNAL									Coeff_FUN
		REAL(8), INTENT(IN)						::	vert(2,4), el_beta(2)
		INTEGER, INTENT(IN)						::	n_nodes_in_elem_1, n_nodes_in_elem_2,	&
													el_type, delta,tri_sign(2)
		REAL(8), INTENT(IN)						::	Intrs_pts(2, el_type)
		REAL(8), INTENT(OUT)					::	matrix(4,4)
		INTEGER                                 ::  pointer_reference_to_local(4)

	END SUBROUTINE

	SUBROUTINE Linear_FE_Basis_Eval_2D(coef, x, y, basis_ind, d_ind_x, d_ind_y, phi)

		REAL(8), INTENT(IN)					::	coef(3,3)
		REAL(8), DIMENSION(:), INTENT(IN)	::	x, y
		INTEGER, INTENT(IN)					::	basis_ind, d_ind_x, d_ind_y
		REAL(8), DIMENSION(:), INTENT(OUT)	::	phi

	END SUBROUTINE

	SUBROUTINE Gauss_Nodes_2D(p_basic, t_basic, g_x, g_y)

		REAL(8), DIMENSION(:,:), INTENT(IN)		::	p_basic
		INTEGER, DIMENSION(:,:), INTENT(IN)		::	t_basic
		REAL(8), DIMENSION(:,:), POINTER		::	g_x, g_y

	END SUBROUTINE

	SUBROUTINE Linear_IFE_Basis_Eval_2D(coef, x, y, region_ind,				&
								    basis_ind, d_ind_x, d_ind_y, phi)

		REAL(8), INTENT(IN)					::	coef(3,6)
		REAL(8), DIMENSION(:), INTENT(IN)	::	x, y
		INTEGER, INTENT(IN)					::	basis_ind, d_ind_x, d_ind_y, &
												region_ind
		REAL(8), DIMENSION(:), INTENT(OUT)	::	phi

	END SUBROUTINE

	SUBROUTINE Linear_IFE_Basis_Coeff_2D(vert, intrs, confg_ind, ElBeta, coef)

		REAL(8), INTENT(IN)					::	vert(2,3), intrs(2,2), ElBeta(2)										
		INTEGER									confg_ind
		REAL(8), INTENT(OUT)				::	coef(3,6)

	END SUBROUTINE

	SUBROUTINE Gauss_Nodes_Edge_2D(vert1, vert2, gnodes)

		REAL(8), INTENT(IN)		::	vert1(2), vert2(2)
		REAL(8), INTENT(OUT)	::	gnodes(2,3)

	END SUBROUTINE

	SUBROUTINE Linear_FE_Basis_Coeff_2D(vert, coef)

		REAL(8), INTENT(IN)		::	vert(2,4)
		REAL(8), INTENT(OUT)	::	coef(4,4)

	END SUBROUTINE

END INTERFACE

! Test
! =====
INTERFACE
	
	SUBROUTINE	Userdef_EBC_Value(p_basic, t_basic_int, node_type, EBC_Value)
		REAL(8), DIMENSION(:,:), INTENT(IN)			::	p_basic
		INTEGER, DIMENSION(:,:), INTENT(IN)			::	t_basic_int, node_type
		REAL(8), DIMENSION(:,:), INTENT(OUT)		::	EBC_Value
	END SUBROUTINE

	
END INTERFACE

! Surface Jump
! ==========
INTERFACE

	SUBROUTINE	generate_Gauss_reference_1D(Gauss_point_number,	Gauss_coefficient_reference_1D,	Gauss_point_reference_1D)

		INTEGER, INTENT(IN)						::	Gauss_point_number
		REAL(8), DIMENSION(:)					::	Gauss_coefficient_reference_1D, Gauss_point_reference_1D

	END SUBROUTINE

	SUBROUTINE	generate_load_vector_local_IFE_on_interface(interface_coefficient_function_name, information_1, information_2,	&
															p_basic, element_index, t_c, nnx, nny,	&
															Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
															test_basis_type, r, node_type, node_index)
		
		EXTERNAL										interface_coefficient_function_name
		INTEGER, DIMENSION(:,:), INTENT(IN)	            ::  information_1
		REAL(8), DIMENSION(:,:), INTENT(IN)	            ::  information_2
		REAL(8), DIMENSION(:,:), INTENT(IN)		        ::	p_basic
		INTEGER, DIMENSION(:), INTENT(IN)	        	::	element_index
		INTEGER, DIMENSION(:,:), INTENT(IN)				::	t_c
		INTEGER, INTENT(IN)                             ::  nnx,nny
		REAL(8), DIMENSION(:), INTENT(IN)	            ::	Gauss_coefficient_reference_1D
		REAL(8), DIMENSION(:), INTENT(IN)	            ::	Gauss_point_reference_1D
		INTEGER, INTENT(IN)								::	test_basis_type
		REAL(8), DIMENSION(:), POINTER		            ::	r
		INTEGER, DIMENSION(:,:), INTENT(IN)				::	node_type
		INTEGER, DIMENSION(:), INTENT(IN)				::	node_index

	END SUBROUTINE

	SUBROUTINE	generate_load_vector_local_IFE_nonhomogeneous(flux_jump_function_name, Global_Beta, information_1, information_2,	&
															p_basic, element_index, t_c, nnx, nny, &
															Gauss_coefficient_reference_1D, Gauss_point_reference_1D,	&
															Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle, &
															test_basis_type, test_derivative_degree_x, test_derivative_degree_y, &
															nonhomogeneous_trial_basis_type, nonhomogeneous_trial_derivative_degree_x, &
															nonhomogeneous_trial_derivative_degree_y, r, node_type)

		EXTERNAL										flux_jump_function_name
		REAL(8), DIMENSION(:), INTENT(IN)        	::	Global_Beta
		INTEGER, DIMENSION(:,:), INTENT(IN)         ::  information_1
		REAL(8), DIMENSION(:,:), INTENT(IN)         ::  information_2
		REAL(8), DIMENSION(:,:), INTENT(IN)	        ::	p_basic
		INTEGER, DIMENSION(:), INTENT(IN)        	::	element_index
		INTEGER, DIMENSION(:,:), INTENT(IN)			::	t_c
		INTEGER, INTENT(IN)                         ::  nnx,nny
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_1D
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_point_reference_1D
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_triangle
		REAL(8), DIMENSION(:,:), INTENT(IN)         ::	Gauss_point_reference_triangle
		INTEGER, INTENT(IN)                         ::  nonhomogeneous_trial_basis_type, nonhomogeneous_trial_derivative_degree_x, &
		                                                nonhomogeneous_trial_derivative_degree_y
		INTEGER, INTENT(IN)							::	test_basis_type, test_derivative_degree_x, test_derivative_degree_y
		REAL(8), DIMENSION(:), POINTER	            ::	r
		INTEGER, DIMENSION(:,:), INTENT(IN)			::	node_type

	END SUBROUTINE

	SUBROUTINE	Gauss_integration_local_load_IFE_on_interface(interface_coefficient_function_name, Gauss_coefficient_reference_1D, &
	                                                            Gauss_point_reference_1D, vertices, information_vector_1, &
                                                information_vector_2, test_basis_type, test_basis_index, r, element_Gauss)
		EXTERNAL										interface_coefficient_function_name
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_1D
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_point_reference_1D
		REAL(8), INTENT(IN)							::	vertices(2,4)
		INTEGER, INTENT(IN)							::	information_vector_1(18)
		REAL(8), INTENT(IN)							::	information_vector_2(8)
		INTEGER, INTENT(IN)							::	test_basis_type
		INTEGER, INTENT(IN)							::	test_basis_index
		REAL(8)										::	r
    Integer :: element_Gauss

	END SUBROUTINE

	SUBROUTINE	Gauss_integration_local_load_IFE_on_interface_Rectangular(interface_coefficient_function_name, &
	                                                                        Gauss_coefficient_reference_1D, &
	                                                                        Gauss_point_reference_1D, vertices, test_basis_type, &
	                                                                        test_basis_index, r, count_x, count_y, x_start, y_start)
		EXTERNAL										interface_coefficient_function_name
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_1D
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_point_reference_1D
		REAL(8), INTENT(IN)							::	vertices(2,4)
		INTEGER, INTENT(IN)							::	test_basis_type
		INTEGER, INTENT(IN)							::	test_basis_index
		REAL(8)										::	r, x_start, y_start
		INTEGER, INTENT(IN)							::	count_x, count_y

	END SUBROUTINE

	SUBROUTINE	generate_Gauss_local_1D(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, lower_bound, &
	                                    upper_bound, Gauss_coefficient_local_1D, Gauss_point_local_1D)
		REAL(8), DIMENSION(:), INTENT(IN)	        ::	Gauss_coefficient_reference_1D
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_point_reference_1D
		REAL(8), INTENT(IN)							::	lower_bound, upper_bound	
		REAL(8), DIMENSION(:)						::	Gauss_coefficient_local_1D
		REAL(8), DIMENSION(:)						::	Gauss_point_local_1D
	END SUBROUTINE

	SUBROUTINE	Gauss_integration_of_flux_jump_on_interface(flux_jump_function_name, Gauss_coefficient_reference_1D, &
	                                                        Gauss_point_reference_1D, vertices, information_vector_1, &
	                                                        information_vector_2, r)

		EXTERNAL										flux_jump_function_name
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_1D
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_point_reference_1D
		REAL(8), INTENT(IN)							::	vertices(2,4)
		INTEGER, INTENT(IN)							::	information_vector_1(18)
		REAL(8), INTENT(IN)							::	information_vector_2(6)
		REAL(8)										::	r

	END SUBROUTINE

	SUBROUTINE	Gauss_integration_local_load_IFE_nonhomogeneous(Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle,	&
																	vertices, information_vector_1, information_vector_2, test_basis_type,	&
																	test_basis_index, test_derivative_degree_x, test_derivative_degree_y, &
																	nonhomogeneous_trial_basis_type,	&
																	nonhomogeneous_trial_derivative_degree_x, &
																	nonhomogeneous_trial_derivative_degree_y, r, element_Gauss)
		REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_triangle
		REAL(8), DIMENSION(:,:), INTENT(IN)         ::	Gauss_point_reference_triangle
		REAL(8), INTENT(IN)							::	vertices(2,4)
		INTEGER, INTENT(IN)							::	information_vector_1(18)
		REAL(8), INTENT(IN)							::	information_vector_2(6)
		INTEGER, INTENT(IN)							::	test_basis_type, test_derivative_degree_x, test_derivative_degree_y
		INTEGER, INTENT(IN)							::	test_basis_index
		INTEGER, INTENT(IN)                         ::  nonhomogeneous_trial_basis_type, nonhomogeneous_trial_derivative_degree_x, &
		                                                nonhomogeneous_trial_derivative_degree_y
		REAL(8)										::	r
    Integer :: element_Gauss
	END SUBROUTINE

	SUBROUTINE	Matrix_Inverse(A,N)
	
		REAL(8), DIMENSION(N,N), INTENT(INOUT)		:: A
		INTEGER, INTENT(IN)							:: N
	
	END SUBROUTINE

	SUBROUTINE	Left_Divide(y, h, g)

		REAL(8), DIMENSION(:,:), INTENT(INOUT)		::	y
		REAL(8), DIMENSION(:,:), INTENT(IN)		::	h, g

	END SUBROUTINE

	SUBROUTINE periodic_boundary_conditions(nnx,nny,bc_index,p_basic,node_type,bc_point_1,bc_point_2,ALQ,A_stiff,rhs)
	    USE IFE_MAIN_PARAM

		INTEGER nnx,nny
		INTEGER, DIMENSION(:),POINTER	        ::	bc_index
		REAL(8), DIMENSION(:,:),POINTER	        ::	p_basic
		INTEGER, DIMENSION(:,:), POINTER			::	node_type
		TYPE(SPARSE), DIMENSION(:),POINTER		::	A_stiff
		REAL(8), DIMENSION(:,:),POINTER	        ::	bc_point_1, bc_point_2
		REAL(8), DIMENSION(:), POINTER	        ::  rhs
		REAL(8),DIMENSION(:,:),	POINTER		::	ALQ
     
     END SUBROUTINE

	SUBROUTINE periodic_boundary_corresponding_nodes(nnx,nny,bc_index,p_basic,node_type,ALQ,bc_point_1,bc_point_2,unknow_nodes)
			
		INTEGER nnx,nny
		INTEGER, DIMENSION(:),POINTER	        ::	bc_index
		REAL(8), DIMENSION(:,:),POINTER	        ::	p_basic
		REAL(8)	,DIMENSION(:,:),POINTER 	    ::	unknow_nodes
		REAL(8), DIMENSION(:,:),POINTER	        ::	bc_point_1, bc_point_2
		REAL(8),DIMENSION(:,:),	POINTER		::	ALQ
		INTEGER, DIMENSION(:,:), POINTER        ::  node_type
			
	ENDSUBROUTINE

	SUBROUTINE solve_periodic_boundary_conditions(i,j,AL,Q,num_unknows,A_stiff,f )
		USE IFE_MAIN_PARAM

		REAL(8) AL,Q
		INTEGER i,j
		INTEGER		num_unknows
		TYPE(SPARSE), DIMENSION(:),POINTER		::	A_stiff
		REAL(8), DIMENSION(:), POINTER	        ::  f

	ENDSUBROUTINE

END INTERFACE


! PIC Boundary
! ==========
INTERFACE
    !$ ===================== ab.ZWZ 2021/7/9 =========================\\
    SUBROUTINE AdjustBoundary_2D(time, dt, delta, N_Objects, objects)
        USE Object_Data_2D
        REAL(8)    :: time, dt
        INTEGER    :: delta
        TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
        INTEGER    :: N_Objects
    END SUBROUTINE
    !$ ===================== ab.ZWZ 2021/7/9 =========================//

    SUBROUTINE dielectric_2D(i_part, ispe, inters_point, boundary, nnp, newloss, i_Boundary)

        USE Boundary_Data_2D

        INTEGER        i_part, ispe
        REAL(8), DIMENSION(2)            ::  inters_point
        TYPE(BoundaryType)                 ::  boundary
        INTEGER        nnp
        INTEGER, DIMENSION(:)            ::  newloss
        INTEGER                            ::  i_Boundary

    END SUBROUTINE

END INTERFACE


! Huygens_Wave
! =======
INTERFACE

	SUBROUTINE Interp1(A,B,nnx,c,inte)
	INTEGER			i, j,nnx
	REAL(8)			A(1:nnx),B(1:nnx)
	REAL(8)			amax,amin,da,d,di,db,ai,c,inte,delt
	END SUBROUTINE

	SUBROUTINE Get_1st_Point_2D(nodes, xp, yp,	yc)
	INTEGER		np, count ,i
	REAL(8)		k2 , yc
	REAL(8)		xyp(2,2)
	REAL(8), DIMENSION(:,:), POINTER		::	nodes
	REAL(8), DIMENSION(:), POINTER			::	xp, yp
	END SUBROUTINE
	
	SUBROUTINE Get_includeflag_2D(nor1, x1, y1, nor2, x2, y2,include_flag)
		REAL(8), INTENT(IN)		::		nor1(2), nor2(2)
		REAL(8), INTENT(IN)		::		x1, y1, x2, y2
		INTEGER, INTENT(OUT)	::		include_flag
		REAL(8)							t1,t2
	END SUBROUTINE
	
	SUBROUTINE Get_Wall_Points_2D(vxs, vys, num_angle,			&
							start_angle, final_angle,kp)
		INTEGER		num_angle, count ,i,ii
		REAL(8)		k1, k3 , start_angle, final_angle,kp,vys_i,s_angle,f_angle,t1,t2,t3,t4,rsangle
		REAL(8)		xt(3), yt(3)
		REAL(8), DIMENSION(:), POINTER		::	vxs, vys
		REAL(8), DIMENSION(:), POINTER		::	vxs_tmp, vys_tmp
	END SUBROUTINE

	SUBROUTINE Get_includearea_2D(kp,deltat,left_theta, right_theta , &
							num_angle,start_angle,final_angle, included_vx,included_vy,walljp)
	INTEGER		num_angle,i,kk,nang
	REAL(8)		kp,deltat,left_theta, right_theta,walljp
	REAL(8)		includearea,half_alpha,start_angle,final_angle,dangle,sec
	REAL(8), DIMENSION(:), POINTER			::	all_angle,included_vx,included_vy
	END SUBROUTINE

	SUBROUTINE Get_includepoint_2D(theta, kp, mj, vxp, vyp)
		REAL(8), INTENT(IN)		::		theta, kp, mj
		REAL(8)							sec,A,Y,diffY,vxp,vyp
	END SUBROUTINE

	SUBROUTINE Get_normalpoint_2D(theta, kp, vnx, vny)
	REAL(8), INTENT(IN)		::		theta, kp
	REAL(8)							sec,A,Y,diffY,vnx,vny
	END SUBROUTINE

	SUBROUTINE Get_Wave_Direction_2D(n_elements, p_x_tmp, 	&
							p_y_tmp, xp, yp , xc, yc)
	REAL(8)									rx1, ry1, rx2, ry2, rz, xc, yc
	REAL(8), DIMENSION(:,:), POINTER	::	p_x_tmp, p_y_tmp
	REAL(8), DIMENSION(:), POINTER		::	xp, yp
	END SUBROUTINE

	SUBROUTINE Curve_tracing_2D(n_int_elements, xp,  yp, p_int_x_tmp2, p_int_y_tmp2)
	INTEGER		n_int_elements, i
	REAL(8)			xp(3),yp(3)
	REAL(8), DIMENSION(:,:), POINTER		::	p_int_x_tmp2,p_int_y_tmp2
	REAL(8), DIMENSION(:,:), POINTER		::	p_int_x_tmp3,p_int_y_tmp3
	END SUBROUTINE

	SUBROUTINE Huygens_Wave_At_Twopoints_2D(xpi, ypi, xp2, yp2, vxp, vyp, nxp, num_angle,xp_right,yp_right,deltat)
		INTEGER		nxp, num_angle ,ii ,xpint,ypint
		REAL(8)		xpi,ypi,xp_right,yp_right,deltat
		REAL(8)		xpreal,ypreal,xptmin,xptmax,yptmin,yptmax,xi,yi
		REAL(8)		xyp(2),yxp(2)
		REAL(8), DIMENSION(:), POINTER		::	xptmp,yptmp
		REAL(8), DIMENSION(:), POINTER		::	xp2, yp2
		REAL(8), DIMENSION(:), POINTER		::	vxp,vyp	
	END SUBROUTINE

	SUBROUTINE Huygens_Wave_2D(it, objects, alphax, E_par,nc_par,Update)		!,alphax,E_par
		USE Object_Data_2D 
        INTEGER		it
        REAL(8), DIMENSION(:,:,:), POINTER	::	alphax,E_par		!,nc_par
        INTEGER(4), DIMENSION(:,:,:), POINTER			::	nc_par
        TYPE(ObjectType), DIMENSION(:), INTENT(INOUT)		::	objects
        INTEGER		nnxp,iwall,n_objects,include_flag,Update

	END SUBROUTINE

	SUBROUTINE Huygens_Wave_Rotate_2D(normalp,beita, vxp,vyp, num_angle)
		INTEGER		num_angle,ii
		REAL(8)		beita,Comega,Somega,Cbeita,Sbeita,r	
		REAL(8)		normalp(2)
		REAL(8)		d1(2),d2(2),r_t(2)
		REAL(8), DIMENSION(:), POINTER	::	vxp,vyp
	END SUBROUTINE

	SUBROUTINE Update_IFE_Start_2D(it,Phi, nnx, nny, delta, objects)
		USE Object_Data_2D
		REAL(8)		xmin, xmax, ymin, ymax 
		INTEGER		it,nnx, nny, N_Objects, N_Boundary,nnxp
		INTEGER		ni, nj, nk, nindex, delta, ix, jy
		INTEGER		nodes(2),count
		INTEGER		 num_of_nodes, n_int_elements, i, j, FID,num_of_unknowns
		REAL(8)		dimensions(2,2)
		REAL(8), DIMENSION(:), POINTER				::	rhs_EBC, rhs_BGS, U_fe_full, U_fe 
		REAL(8), DIMENSION(:), POINTER				::	rhs_NBC
		EXTERNAL										FUN_2D_ONE
		INTEGER, DIMENSION(:,:), POINTER			::	t_basic
		INTEGER, DIMENSION(:), POINTER				::	element_index , node_index,bnd_elem_index,t_iel_tmp 
		REAL(8), DIMENSION(:,:), POINTER			::	p_int_x_tmp, p_int_y_tmp
		INTEGER, DIMENSION(:,:), POINTER			::	p_int_f_tmp, p_int_f
		TYPE(ObjectType), DIMENSION(:), INTENT(INOUT)		::	objects
		REAL(8)			Phi(0:nnx+1,0:nny+1)  
	END SUBROUTINE

	SUBROUTINE Setup_IFE_Wall_Mesh_2D(t_c, p_basic, information_2, objects)
		USE Object_Data_2D
		INTEGER			N_Boundary
		INTEGER		 num_of_nodes, n_int_elements, i, j,num_of_waindex,ii
		REAL(8)		dimensions(2,2),xo,yo,xc,yc,rx1,ry1,rx2,ry2,rz
		INTEGER		count,nnx,nny
		INTEGER		nodes(2)
		INTEGER, DIMENSION(:,:), POINTER			::	t_c
		REAL(8), DIMENSION(:,:), POINTER			::	p_basic
		REAL(8), DIMENSION(:), POINTER			::	xp_tmp,yp_tmp
		REAL(8), DIMENSION(:,:), POINTER		::	p_int_x_tmp2,p_int_y_tmp2
		REAL(8), DIMENSION(:,:), POINTER		::	p_int_x_tmp1,p_int_y_tmp1
		REAL(8), DIMENSION(:,:), POINTER		::	information_2
		REAL(8), DIMENSION(:,:,:), POINTER		::	wallindex_tmp
		TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
	END SUBROUTINE

	SUBROUTINE Update_IFE_Wall_Mesh_2D(it, t_c, p_basic, &
							information_2, objects,node_index)
		USE Object_Data_2D
		INTEGER			N_Objects
		CHARACTER*30	fname, filename
		INTEGER		it, num_of_nodes, n_int_elements, i, j,num_of_waindex, num_dnode
		REAL(8)		xyp(2,2),xo,yo,xc,yc,rx1,ry1,rx2,ry2,rz
		INTEGER		np,count,nnx,nny,n_xo
		REAL(8)		k2,rp
		INTEGER, DIMENSION(:,:), POINTER			::	t_c
		REAL(8), DIMENSION(:,:), POINTER			::	p_basic
		REAL(8), DIMENSION(:), POINTER			::	xp_tmp,yp_tmp,num_p
		REAL(8), DIMENSION(:,:), POINTER		::	p_int_x_tmp2,p_int_y_tmp2
		REAL(8), DIMENSION(:,:), POINTER		::	p_int_x_tmp1,p_int_y_tmp1
		REAL(8), DIMENSION(:,:), POINTER			::	information_2,nodes		!p_int_x_tmp, p_int_y_tmp
		REAL(8), DIMENSION(:,:,:), POINTER		::	wallindex_tmp
		TYPE(ObjectType), DIMENSION(:), INTENT(INOUT)		::	objects
		INTEGER, DIMENSION(:), POINTER				::    node_index
	END SUBROUTINE

	SUBROUTINE Update_Mesh_Objects_Intersection_2D  (t_basic, p_basic, objects,		&
												element_index, information_1, information_2, node_index)
	USE IFE_MAIN_PARAM
	USE Object_Data_2D
		INTEGER, DIMENSION(:,:), POINTER			::	t_basic
		REAL(8), DIMENSION(:,:), POINTER			::	p_basic
		TYPE(ObjectType), DIMENSION(:), INTENT(IN)		::	objects
		INTEGER, DIMENSION(:), POINTER				::	element_index
		INTEGER, DIMENSION(:,:),POINTER							::  information_1
		INTEGER, DIMENSION(:,:),POINTER							::  information_1_temp
		INTEGER, DIMENSION(:),POINTER							::  information_1_el
		REAL(8), DIMENSION(:,:),POINTER							::  information_2
		REAL(8), DIMENSION(:,:),POINTER							::  information_2_temp
		REAL(8), DIMENSION(:),POINTER							::  information_2_el
		INTEGER, DIMENSION(:), INTENT(IN)						::	node_index
		INTEGER		n_elements, n_nodes, ee, ie, n_int_elements,n_objects
		INTEGER		el_type, el_region, i, j		
		REAL(8)		vert(2,4)
		REAL(8), DIMENSION(:,:), POINTER			::	ObjLowBound, ObjUpBound
		INTEGER, DIMENSION(4)						::	node_index_el
		LOGICAL											OutOfBound
		REAL(8)											P_intrs(4,6)
	END SUBROUTINE

	SUBROUTINE	Delete_Used_Elements_2D (n_elements, i, p_x_1, p_y_1)
	REAL(8), DIMENSION(:,:), POINTER		::	p_x_1,p_y_1
	REAL(8), DIMENSION(:,:), POINTER		::	p_x_2,p_y_2
	END SUBROUTINE

	SUBROUTINE Delete_Improper_Points_2D (	n_p, i, p_x_1, p_y_1)
	INTEGER		n_p, i, count
	REAL(8), DIMENSION(:), POINTER		::	p_x_1,p_y_1
	REAL(8), DIMENSION(:), POINTER		::	p_x_2,p_y_2
	END SUBROUTINE

	SUBROUTINE Normalization_2D(Eion,alphax, nc_ion, objects,kx)
	USE Object_Data_2D
!	USE Erosion_PARAM
	REAL(8), DIMENSION(:,:,:), POINTER				::	Eion,alphax
	INTEGER(4), DIMENSION(:,:,:), POINTER			::	nc_ion
	REAL(8), DIMENSION(:,:), POINTER				::	kx
	TYPE(ObjectType), DIMENSION(:), INTENT(IN)		::	objects
	INTEGER		iwall,n_objects,n_sect,i
	REAL(8)		s,wallN
	END SUBROUTINE

	SUBROUTINE AdjustBoundary_Sputter_2D(dt, delta,it)
	USE PIC_MAIN_PARAM_2D
	USE Domain_2D
	REAL(8)		dt
	INTEGER		delta
	END SUBROUTINE

	SUBROUTINE Update_Locate_Node_2D	(	p_node, object, node_region,Phi_p)
	USE IFE_MAIN_PARAM
	USE Object_Data_2D
	REAL(8), DIMENSION(2), INTENT(IN)	::	p_node
	TYPE(ObjectType), INTENT(IN)		::	object
	INTEGER									i,node_region, coord1, coord2 
	INTEGER									n_section,n_nodes,count,count2
	REAL(8), DIMENSION(2)				::	X1, X2, Xb, Xt, Xc, Xtb, Xpb, Xpt, Xpc
	REAL(8)									y_node,xi,Phi_p
	REAL(8), DIMENSION(2,2)				::	xyp
	REAL(8)		Norm_2 
	END SUBROUTINE

	SUBROUTINE El_Curve_Intersection_info_2D_BL(vert, object, P_intrs, el_type, el_region,     &
	                            information_1_el, information_2_el, node_index_el)
	USE IFE_MAIN_PARAM
	USE Object_Data_2D
	IMPLICIT NONE
	REAL(8), DIMENSION(2,4), INTENT(IN)		::	vert
	TYPE(ObjectType), INTENT(IN)			::	object
	REAL(8), DIMENSION(4,6)					::	P_intrs
	INTEGER										el_region
	INTEGER, DIMENSION(:), POINTER	        ::  information_1_el
	REAL(8), DIMENSION(:), POINTER	        ::  information_2_el
	INTEGER, DIMENSION(:), INTENT(IN)		::	node_index_el
	INTEGER										el_type
	END SUBROUTINE

	SUBROUTINE Output_To_Tecplot_Wall_2D(	it, objects)
	USE Object_Data_2D
	INTEGER						::	it	
	TYPE(ObjectType), DIMENSION(:), INTENT(IN)		::	objects
	INTEGER			i, j, N_Objects
	END SUBROUTINE
	
	SUBROUTINE Setup_IFE_Wall_Mesh_2D_NEW(objects)
		USE Object_Data_2D
		INTEGER			N_Boundary
		INTEGER		 num_of_nodes, n_int_elements, i, j,num_of_waindex,ii
		REAL(8)		dimensions(2,2),xo,yo,xc,yc,rx1,ry1,rx2,ry2,rz
		INTEGER		count,nnx,nny
		INTEGER		nodes(2)
		TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
    END SUBROUTINE

    
    
    
    SUBROUTINE	GTOP(i_part,parapre,para)

        USE PIC_MAIN_PARAM_2D
        USE Domain_2D
        !USE Particle_2D
        USE Field_2D
        IMPLICIT NONE
        REAL(8)		            ::	para(0:nx+1,0:ny+1)
        INTEGER					::	i_part,i,j,delta
        REAL(8)					::	parapre, xp, yp, dx, dy
        REAL(8)		            xcellmdx, ycellmdy, R1, R2, den
        
    END SUBROUTINE
    
    
    
END INTERFACE

! IFE error --- ab.ZWZ 2021/7/9
! ==========
INTERFACE
    SUBROUTINE Compute_IFE_Norm_error(delta, Gauss_coefficient_reference, Gauss_point_reference, Gauss_point_reference_triangle, &
                                     Gauss_coefficient_reference_triangle, derivative_degree_x, derivative_degree_y, error_norm)  
  
        INTEGER, INTENT(IN) :: delta      
        REAL(8),DIMENSION(:),POINTER                ::  Gauss_coefficient_reference
        REAL(8),DIMENSION(:,:),POINTER              ::  Gauss_point_reference
        REAL(8),DIMENSION(:,:),POINTER              ::  Gauss_point_reference_triangle
        REAL(8),DIMENSION(:),POINTER                ::  Gauss_coefficient_reference_triangle
        INTEGER     :: derivative_degree_x, derivative_degree_y
        REAL(8)     :: error_norm
    END SUBROUTINE
                                     
    SUBROUTINE  Gauss_integration_local_error_IFE(delta, Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
			                                          vertices,information_vector_1,information_vector_2,elem_index, &
												      test_basis_type,test_derivative_degree_x,test_derivative_degree_y,U_fe,error, element_Gauss)

        INTEGER                                     ::	  delta
        REAL(8),DIMENSION(:)				        ::    Gauss_coefficient_reference_triangle
        REAL(8),DIMENSION(:,:)		                ::    Gauss_point_reference_triangle
        REAL(8)                                           vertices(2,4)
        INTEGER                                           information_vector_1(18)
        REAL(8)                                           information_vector_2(8)
        INTEGER                                           elem_index
        INTEGER                                           test_basis_type,test_derivative_degree_x,test_derivative_degree_y
        REAL(8),DIMENSION(:),ALLOCATABLE            ::    U_fe 
        REAL(8)								        	  error
        Integer :: element_Gauss
                                            
    END SUBROUTINE      
                                     
END INTERFACE

! Test --- ab.ZWZ 2021/7/9
! =====
INTERFACE

    SUBROUTINE	Userdef_EBC_Value_bjw(delta, p_basic, t_basic_int, node_type, EBC_Value, EBC_Value_xt)
        INTEGER, INTENT(IN) :: delta    !$ ab.ZWZ 2021/7/7
		REAL(8), DIMENSION(:,:), INTENT(IN)			::	p_basic
		INTEGER, DIMENSION(:,:), INTENT(IN)			::	t_basic_int, node_type
		REAL(8), DIMENSION(:,:), INTENT(OUT)		::	EBC_Value
        REAL(8), DIMENSION(:), INTENT(OUT)		    ::	EBC_Value_xt
    END SUBROUTINE
    
    SUBROUTINE CheckTrail(N_Objects,objects)
        USE Object_Data_2D
        TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
        INTEGER    :: N_Objects
    END SUBROUTINE
END INTERFACE
END MODULE