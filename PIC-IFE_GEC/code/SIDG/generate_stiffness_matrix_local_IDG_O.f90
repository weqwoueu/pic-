SUBROUTINE generate_stiffness_matrix_local_IDG_O( Global_Beta, node_type,num_of_unknowns, &
                                                information_1, information_2, information_3, HP, element_index, edge_index,&
									            HT, HE, dimensions, &
                                                Gauss_coefficient_reference_1D, Gauss_point_reference_1D,&
                                                trial_basis_type, test_basis_type, matrix, matrix_xt, delta, con_penalty)
   
USE IFE_MAIN_PARAM
USE Gauss_Data
USE IFE_INTERFACE, ONLY: generate_Gauss_local_1D_Linear, Sparse_Structure, Sparse_Structure_xt, &
                            gauss_integration_local_stiffness_DG_FE_O, gauss_integration_local_stiffness_IDG_O, &
                            gauss_integration_local_stiffness_Mixed_O, gauss_integration_local_stiffness_IDG_both_O
    
IMPLICIT NONE

INTEGER,DIMENSION(:,:),POINTER        ::    information_1
REAL(8),DIMENSION(:,:),POINTER        ::    information_2, information_3
INTEGER,DIMENSION(:), POINTER         ::    element_index, edge_index
REAL(8),DIMENSION(:)                  ::    Gauss_coefficient_reference_1D
REAL(8),DIMENSION(:)                  ::    Gauss_point_reference_1D
INTEGER, DIMENSION(:,:), POINTER      ::    HT, HE
REAL(8)                               ::    dimensions(2,2)
INTEGER                               ::    out_index, in_index
INTEGER                               ::    trial_basis_type, test_basis_type
REAL(8)                               ::    con_penalty
    
INTEGER                               ::    i, j ,k, m, n , ii!, NZ
REAL(8),DIMENSION(:),POINTER          ::    Gauss_coefficient_local
REAL(8),DIMENSION(:,:),POINTER        ::    Gauss_point_local
REAL(8), DIMENSION(:,:), POINTER      ::    HP 
INTEGER                               ::    information_vector_1(18), information_vector_1_next(18)
REAL(8)                               ::    information_vector_2(8), information_vector_2_next(8)
REAL(8)                               ::    vertices_linear(2,2)
REAL(8), DIMENSION(2, 4)              ::    vertices_Rectangular, brother_vertices_Rectangular
REAL(8), DIMENSION(:), POINTER        ::	  Global_Beta
INTEGER,DIMENSION(4,2)                ::    this_flag1, this_flag2

!-----------------------use the sparse matrix to store the Global stiff-----------------------
INTEGER, DIMENSION(:,:), POINTER	      ::	VROW
TYPE(SPARSE), DIMENSION(:), POINTER	    ::	matrix
TYPE(SPARSE), DIMENSION(:), ALLOCATABLE	::	matrix_xt
INTEGER, DIMENSION(:,:), POINTER	      ::	node_type
INTEGER                                 ::  NZ
INTEGER, INTENT(IN)					            ::	num_of_unknowns, delta
INTEGER                                 ::  n_nodes_in_elem, loc1, loc2
REAL(8)	                                ::  El_Stiff(4,4), El_Stiff_brother(4,4), &
                                            El_Stiff_b_b(4,4), El_Stiff_b_t(4,4)
REAL(8)                                 ::  int_value, beta_fun_coe
INTEGER                                 ::  this_element, next_element, this_basic, next_basic

out_index = -1
in_index = -2
n_nodes_in_elem = 4
int_value = 0
this_flag1 = 0

CALL Sparse_Structure(HT, HE, node_type, VROW, matrix, NZ)
CALL Sparse_Structure_xt(HT, HE, node_type, VROW, matrix_xt, NZ)

do n=1,SIZE(HE,2)
    
    if (HE(6, n) /= 0 .And. HE(8, n) == 1) then
        
        El_Stiff_brother = 0
        El_Stiff = 0
        El_Stiff_b_b = 0
        El_Stiff_b_t = 0
        this_element = HE(5, n)
        next_element = HE(6, n)
        ! write(*,*) 'element_index(this_element) = ', element_index(this_element), &
        ! 'element_index(next_element) = ', element_index(next_element)
        vertices_linear(:,1) = HP(:, HE(1, n))
        vertices_linear(:,2) = HP(:, HE(2, n)) 
        vertices_Rectangular = HP(:,HT(1:4,HE(5,n)))
        brother_vertices_Rectangular = HP(:, HT(1:4, HE(6, n)))
      CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices_linear(:,1),&
                                    vertices_linear(:,2), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
        
        if (element_index(this_element) == out_index .AND. element_index(next_element) == out_index) then
            
            do j = 1, n_nodes_in_elem
                do k = 1 ,n_nodes_in_elem
                    
                    this_flag1(2:3,1) = HE(3:4, n)
                    this_flag1(2:3,2) = HE(3:4, n)
                    if (HE(5, n) > HE(6, n)) then
                        this_flag1(4, 1) = 1
                        this_flag1(4, 2) = 1
                    else
                        this_flag1(4, 1) = -1
                        this_flag1(4, 2) = -1
                    end if
                    
                CALL gauss_integration_local_stiffness_DG_FE_O(delta,Global_Beta(1),Gauss_Coefficient_Local_1D_Linear_Eight,&
                          Gauss_Point_Local_1D_Linear_Eight, &
                          this_flag1, vertices_Rectangular, vertices_Rectangular, trial_basis_type, j, &
											    test_basis_type, k, con_penalty, int_value, HE(5,n))
                    
                    El_Stiff(j,k) = El_Stiff(j,k) + int_value
                    
                enddo
            enddo
            
            
            do j = 1, n_nodes_in_elem
                do k = 1 ,n_nodes_in_elem
                    
                    this_flag1(2:3,1) = HE(3:4, n)
                    this_flag1(2:3,2) = HE(3:4, n)
                    if (HE(5, n) > HE(6, n)) then
                        this_flag1(4, 1) = 1
                        this_flag1(4, 2) = -1
                    else
                        this_flag1(4, 1) = -1
                        this_flag1(4, 2) = 1
                    end if
                    
                CALL gauss_integration_local_stiffness_DG_FE_O(delta,Global_Beta(1),Gauss_Coefficient_Local_1D_Linear_Eight,&
                                                Gauss_Point_Local_1D_Linear_Eight, &
	                                            this_flag1, vertices_Rectangular, &
                                                brother_vertices_Rectangular, trial_basis_type, j, &
											    test_basis_type, k, con_penalty, int_value, HE(5,n))
                    
                    El_Stiff_brother(j,k) = El_Stiff_brother(j,k) + int_value
                    
                enddo
            enddo

            if  (HT(5, next_element) == 0) then
                
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = -1
                        else
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = 1
                        end if
                    
                CALL gauss_integration_local_stiffness_DG_FE_O(delta, Global_Beta(1),Gauss_Coefficient_Local_1D_Linear_Eight,&
                                                  Gauss_Point_Local_1D_Linear_Eight, &
                            this_flag1, brother_vertices_Rectangular, brother_vertices_Rectangular, trial_basis_type, j, &
											        test_basis_type, k, con_penalty, int_value, HE(6,n))
                    
                        El_Stiff_b_b(j,k) = El_Stiff_b_b(j,k) + int_value
                    
                    enddo
                enddo
                
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = 1
                        else
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = -1
                        end if
                    
                CALL gauss_integration_local_stiffness_DG_FE_O(delta,Global_Beta(1),Gauss_Coefficient_Local_1D_Linear_Eight,&
                                                    Gauss_Point_Local_1D_Linear_Eight, &
	                                                this_flag1, brother_vertices_Rectangular, &
                                                    vertices_Rectangular, trial_basis_type, j, &
											        test_basis_type, k, con_penalty, int_value, HE(6,n))
                    
                        El_Stiff_b_t(j,k) = El_Stiff_b_t(j,k) + int_value
                    
                    enddo
                enddo
                
            end if
        elseif (element_index(this_element) == in_index .AND. element_index(next_element) == in_index) then
            
            do j = 1, n_nodes_in_elem
                do k = 1 ,n_nodes_in_elem
                    
                    this_flag1(2:3,1) = HE(3:4, n)
                    this_flag1(2:3,2) = HE(3:4, n)
                    if (HE(5, n) > HE(6, n)) then
                        this_flag1(4, 1) = 1
                        this_flag1(4, 2) = 1
                    else
                        this_flag1(4, 1) = -1
                        this_flag1(4, 2) = -1
                    end if
                    
                    CALL gauss_integration_local_stiffness_DG_FE_O(delta, Global_Beta(2),&
                                                Gauss_Coefficient_Local_1D_Linear_Eight, &
                                                Gauss_Point_Local_1D_Linear_Eight, &
                                this_flag1, vertices_Rectangular, vertices_Rectangular, trial_basis_type, j, &
											    test_basis_type, k, con_penalty, int_value, HE(5,n))
                    
                    El_Stiff(j,k) = El_Stiff(j,k) + int_value
                    
                enddo
            enddo
            
            
            do j = 1, n_nodes_in_elem
                do k = 1 ,n_nodes_in_elem
                    
                    this_flag1(2:3,1) = HE(3:4, n)
                    this_flag1(2:3,2) = HE(3:4, n)
                    if (HE(5, n) > HE(6, n)) then
                        this_flag1(4, 1) = 1
                        this_flag1(4, 2) = -1
                    else
                        this_flag1(4, 1) = -1
                        this_flag1(4, 2) = 1
                    end if
                    
                    CALL gauss_integration_local_stiffness_DG_FE_O(delta, Global_Beta(2), &
                                                Gauss_Coefficient_Local_1D_Linear_Eight, &
                                                Gauss_Point_Local_1D_Linear_Eight, &
	                                            this_flag1, vertices_Rectangular, brother_vertices_Rectangular, &
                                                trial_basis_type, j, test_basis_type, k, con_penalty, int_value, HE(5,n))
                    
                    El_Stiff_brother(j,k) = El_Stiff_brother(j,k) + int_value
                    
                enddo
            enddo
            
            if (HT(5, next_element) == 0) then
                
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = -1
                        else
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = 1
                        end if
                    
              CALL gauss_integration_local_stiffness_DG_FE_O(delta, Global_Beta(2),Gauss_Coefficient_Local_1D_Linear_Eight,&
                                                    Gauss_Point_Local_1D_Linear_Eight, &
                          this_flag1, brother_vertices_Rectangular, brother_vertices_Rectangular, trial_basis_type, j, &
											        test_basis_type, k, con_penalty, int_value, HE(6,n))
                    
                        El_Stiff_b_b(j,k) = El_Stiff_b_b(j,k) + int_value
                    
                    enddo
                enddo
                
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = 1
                        else
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = -1
                        end if
                    
                        CALL gauss_integration_local_stiffness_DG_FE_O(delta, Global_Beta(2), &
                                                    Gauss_Coefficient_Local_1D_Linear_Eight, &
                                                    Gauss_Point_Local_1D_Linear_Eight, &
	                                                this_flag1, brother_vertices_Rectangular, vertices_Rectangular, &
                                                    trial_basis_type, j, test_basis_type, k, con_penalty, int_value, HE(6,n))
                    
                        El_Stiff_b_t(j,k) = El_Stiff_b_t(j,k) + int_value
                    
                    enddo
                enddo
                
                
            end if
            
        elseif (element_index(this_element) > 0 .AND. element_index(next_element) < 0) then 
            
            information_vector_1 = information_1(:, element_index(this_element))
            information_vector_2 = information_2(:, element_index(this_element))
            
            if (edge_index(n) < 0) then
            
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = 1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = -1
                        end if
                    
                        CALL gauss_integration_local_stiffness_IDG_O(delta, Global_Beta,information_vector_1 ,&
                                                    information_vector_2, &
                                                    Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
	                                                this_flag1, vertices_Rectangular, trial_basis_type, j, &
											        test_basis_type, k, edge_index(n), con_penalty, int_value, HE(5,n))
                    
                        El_Stiff(j,k) = El_Stiff(j,k) + int_value
                    
                    enddo
                enddo
            
            
                this_basic = 1; !IFE
                next_basic = 0; ! FE
            
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = -1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = 1
                        end if
                    
                        CALL gauss_integration_local_stiffness_Mixed_O(delta, Global_Beta, &
                                                information_vector_1, information_vector_2, &
                                                Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                                this_flag1, vertices_Rectangular, brother_vertices_Rectangular, &
                                                trial_basis_type, j, &
                                                test_basis_type, k, edge_index(n), this_basic, next_basic, con_penalty, int_value, HE(5,n))
                    
                        El_Stiff_brother(j,k) = El_Stiff_brother(j,k) + int_value
                    
                    enddo
                enddo
            
                if (HT(5, next_element) == 0) then
                
                    if ( element_index(next_element) == in_index ) then
                        beta_fun_coe = Global_Beta(2)
                    elseif ( element_index(next_element) == out_index ) then
                        beta_fun_coe = Global_Beta(1)
                    end if
                
                    do j = 1, n_nodes_in_elem
                        do k = 1 ,n_nodes_in_elem
                    
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = -1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = 1
                            end if
                    
                            CALL gauss_integration_local_stiffness_DG_FE_O(delta, beta_fun_coe, &
                                        Gauss_Coefficient_Local_1D_Linear_Eight,&
                                        Gauss_Point_Local_1D_Linear_Eight, this_flag1, brother_vertices_Rectangular, &
                                        brother_vertices_Rectangular, trial_basis_type, j, &
										test_basis_type, k, con_penalty, int_value, HE(6,n))
                    
                            El_Stiff_b_b(j,k) = El_Stiff_b_b(j,k) + int_value
                    
                        enddo
                    enddo
                
                    this_basic = 0; ! FE
                    next_basic = 1; ! IFE
                
                    do j = 1, n_nodes_in_elem
                        do k = 1 ,n_nodes_in_elem
                    
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = 1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = -1
                            end if
                    
                            CALL gauss_integration_local_stiffness_Mixed_O(delta, Global_Beta, &
                                    information_vector_1, information_vector_2, &
                                    Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                    this_flag1, brother_vertices_Rectangular, vertices_Rectangular, &
                                    trial_basis_type, j, &
                                    test_basis_type, k, edge_index(n), this_basic, next_basic, con_penalty, int_value, HE(5,n))
                    
                            El_Stiff_b_t(j,k) = El_Stiff_b_t(j,k) + int_value
                    
                        enddo
                    enddo
    
                end if
            else if (edge_index(n) > 0) Then
                

                
                CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, &
                                    Gauss_point_reference_1D, vertices_linear(:,1), &
                                    information_3(4:5, edge_index(n)), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
                
                loc1 = information_3(6, edge_index(n))
                loc1 = out_index
                
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = 1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = -1
                        end if
                    
                        CALL gauss_integration_local_stiffness_IDG_O(delta, Global_Beta,information_vector_1 ,&
                                    information_vector_2, &
                                    Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                    this_flag1, vertices_Rectangular, trial_basis_type, j, &
                                    test_basis_type, k, loc1, con_penalty, int_value, HE(5,n))
                    
                        El_Stiff(j,k) = El_Stiff(j,k) + int_value
                    
                    enddo
                enddo
                
                this_basic = 1; !IFE
                next_basic = 0; ! FE
            
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = -1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = 1
                        end if
                    
                        CALL gauss_integration_local_stiffness_Mixed_O(delta, Global_Beta, &
                                    information_vector_1, information_vector_2, &
                                    Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                    this_flag1, vertices_Rectangular, brother_vertices_Rectangular, &
                                    trial_basis_type, j, &
                                    test_basis_type, k, loc1, this_basic, next_basic, con_penalty, int_value, HE(5,n))
                    
                        El_Stiff_brother(j,k) = El_Stiff_brother(j,k) + int_value
                    
                    enddo
                enddo
                
                CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
                                    information_3(4:5, edge_index(n)),&
                                    vertices_linear(:,2), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
                
                loc2 = information_3(7, edge_index(n))
                loc2 = out_index
                
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = 1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = -1
                        end if
                    
                        CALL gauss_integration_local_stiffness_IDG_O(delta, Global_Beta,information_vector_1 ,&
                                    information_vector_2, &
                                    Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                    this_flag1, vertices_Rectangular, trial_basis_type, j, &
                                    test_basis_type, k, loc2, con_penalty, int_value, HE(5,n))
                    
                        El_Stiff(j,k) = El_Stiff(j,k) + int_value
                    
                    enddo
                enddo
                
                this_basic = 1; !IFE
                next_basic = 0; ! FE
            
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = -1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = 1
                        end if
                    
                        CALL gauss_integration_local_stiffness_Mixed_O(delta, Global_Beta, &
                                information_vector_1, information_vector_2, &
                                Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                this_flag1, vertices_Rectangular, brother_vertices_Rectangular, &
                                trial_basis_type, j, &
                                test_basis_type, k, loc2, this_basic, next_basic, con_penalty, int_value, HE(5,n))
                    
                        El_Stiff_brother(j,k) = El_Stiff_brother(j,k) + int_value
                    
                    enddo
                enddo
                
                If (HT(5, next_element) == 0) then
                    
                    CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, &
                                        Gauss_point_reference_1D, vertices_linear(:,1), &
                                        information_3(4:5, edge_index(n)), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                        Gauss_Point_Local_1D_Linear_Eight)
                
                    loc1 = information_3(6, edge_index(n))
                    loc1 = out_index

                    if ( element_index(next_element) == in_index ) then
                        beta_fun_coe = Global_Beta(2)
                    elseif ( element_index(next_element) == out_index ) then
                        beta_fun_coe = Global_Beta(1)
                    end if
                    
                    do j = 1, n_nodes_in_elem
                        do k = 1 ,n_nodes_in_elem
                    
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = -1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = 1
                            end if
                    
                            CALL gauss_integration_local_stiffness_DG_FE_O(delta, beta_fun_coe, &
                                            Gauss_Coefficient_Local_1D_Linear_Eight,&
                                            Gauss_Point_Local_1D_Linear_Eight, this_flag1, brother_vertices_Rectangular, &
                                            brother_vertices_Rectangular, trial_basis_type, j, &
                                            test_basis_type, k, con_penalty, int_value, HE(6,n))
                    
                            El_Stiff_b_b(j,k) = El_Stiff_b_b(j,k) + int_value
                    
                        enddo
                    enddo
                
                    this_basic = 0; ! FE
                    next_basic = 1; ! IFE
                
                    do j = 1, n_nodes_in_elem
                        do k = 1 ,n_nodes_in_elem
                    
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = 1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = -1
                            end if
                    
                            CALL gauss_integration_local_stiffness_Mixed_O(delta, Global_Beta, &
                                    information_vector_1, information_vector_2, &
                                    Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                    this_flag1, brother_vertices_Rectangular, vertices_Rectangular, &
                                    trial_basis_type, j, &
                                    test_basis_type, k, loc1, this_basic, next_basic, con_penalty, int_value, HE(5,n))
                    
                            El_Stiff_b_t(j,k) = El_Stiff_b_t(j,k) + int_value
                    
                        enddo
                    enddo
                    
                    CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
                                    information_3(4:5, edge_index(n)),&
                                    vertices_linear(:,2), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
                
                    loc2 = information_3(7, edge_index(n))
                    loc2 = out_index
                    
                    if ( element_index(next_element) == in_index ) then
                        beta_fun_coe = Global_Beta(2)
                    elseif ( element_index(next_element) == out_index ) then
                        beta_fun_coe = Global_Beta(1)
                    end if
                    
                    do j = 1, n_nodes_in_elem
                        do k = 1 ,n_nodes_in_elem
                    
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = -1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = 1
                            end if
                    
                            CALL gauss_integration_local_stiffness_DG_FE_O(delta, beta_fun_coe, &
                                    Gauss_Coefficient_Local_1D_Linear_Eight,&
                                    Gauss_Point_Local_1D_Linear_Eight, this_flag1, brother_vertices_Rectangular, &
                                    brother_vertices_Rectangular, trial_basis_type, j, &
                                    test_basis_type, k, con_penalty, int_value, HE(6,n))
                    
                            El_Stiff_b_b(j,k) = El_Stiff_b_b(j,k) + int_value
                    
                        enddo
                    enddo
                
                    this_basic = 0; ! FE
                    next_basic = 1; ! IFE
                
                    do j = 1, n_nodes_in_elem
                        do k = 1 ,n_nodes_in_elem
                    
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = 1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = -1
                            end if
                    
                            CALL gauss_integration_local_stiffness_Mixed_O(delta, Global_Beta, &
                                    information_vector_1, information_vector_2, &
                                    Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                    this_flag1, brother_vertices_Rectangular, vertices_Rectangular, &
                                    trial_basis_type, j, &
                                    test_basis_type, k, loc2, this_basic, next_basic, con_penalty, int_value, HE(5,n))
                    
                            El_Stiff_b_t(j,k) = El_Stiff_b_t(j,k) + int_value
                    
                        enddo
                    enddo
                    
                End If
                
                !write(*,*) 'not exist this setting pls add'
                !stop
                
            End if


        elseif (element_index(this_element) < 0 .AND. element_index(next_element) > 0) then 
            
            information_vector_1 = information_1(:, element_index(next_element))
            information_vector_2 = information_2(:, element_index(next_element))
            
            if ( element_index(this_element) == in_index ) then
                beta_fun_coe = Global_Beta(2)
            elseif ( element_index(this_element) == out_index ) then
                beta_fun_coe = Global_Beta(1)
            end if
            
            do j = 1, n_nodes_in_elem
                do k = 1 ,n_nodes_in_elem
                    
                    this_flag1(2:3,1) = HE(3:4, n)
                    this_flag1(2:3,2) = HE(3:4, n)
                    if (HE(5, n) > HE(6, n)) then
                        this_flag1(4, 1) = 1
                        this_flag1(4, 2) = 1
                    else
                        this_flag1(4, 1) = -1
                        this_flag1(4, 2) = -1
                    end if
                    
            CALL gauss_integration_local_stiffness_DG_FE_O(delta, beta_fun_coe,Gauss_Coefficient_Local_1D_Linear_Eight, &
                                                Gauss_Point_Local_1D_Linear_Eight, &
                                                this_flag1, vertices_Rectangular, vertices_Rectangular, trial_basis_type, j, &
											    test_basis_type, k, con_penalty, int_value, HE(5,n))
                    
                    El_Stiff(j,k) = El_Stiff(j,k) + int_value
                    
                enddo
            enddo
            
            
            this_basic = 0; !IFE
            next_basic = 1; ! FE
            
            do j = 1, n_nodes_in_elem
                do k = 1 ,n_nodes_in_elem
                    
                    this_flag1(2:3,1) = HE(3:4, n)
                    this_flag1(2:3,2) = HE(3:4, n)
                    if (HE(5, n) > HE(6, n)) then
                        this_flag1(4, 1) = 1
                        this_flag1(4, 2) = -1
                    else
                        this_flag1(4, 1) = -1
                        this_flag1(4, 2) = 1
                    end if
                    
                    CALL gauss_integration_local_stiffness_Mixed_O(delta, Global_Beta, information_vector_1, &
                                                information_vector_2, &
                                                Gauss_Coefficient_Local_1D_Linear_Eight, &
                                                Gauss_Point_Local_1D_Linear_Eight, &
                                          this_flag1, vertices_Rectangular, brother_vertices_Rectangular, &
                                                trial_basis_type, j, &
											    test_basis_type, k, edge_index(n), this_basic, next_basic, &
                                                con_penalty, int_value, HE(6,n))
                    
                    El_Stiff_brother(j,k) = El_Stiff_brother(j,k) + int_value
                    
                enddo
            enddo
            
            if (HT(5, next_element) == 0) then
                
                this_basic = 1; !IFE
                next_basic = 0; ! FE
            
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = 1
                        else
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = -1
                        end if
                    
                        CALL gauss_integration_local_stiffness_Mixed_O(delta, Global_Beta, information_vector_1, &
                                                    information_vector_2, &
                                                    Gauss_Coefficient_Local_1D_Linear_Eight, &
                                                    Gauss_Point_Local_1D_Linear_Eight, &
	                                                this_flag1, brother_vertices_Rectangular, vertices_Rectangular, &
                                                    trial_basis_type, j, &
											        test_basis_type, k, edge_index(n), this_basic, next_basic, &
                                                    con_penalty, int_value, HE(6,n))
                    
                        El_Stiff_b_t(j,k) = El_Stiff_b_t(j,k) + int_value
                    
                    enddo
                enddo
            
                do j = 1, n_nodes_in_elem
                    do k = 1 ,n_nodes_in_elem
                    
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = -1
                        else
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = 1
                        end if
                    
                        CALL gauss_integration_local_stiffness_IDG_O(delta, Global_Beta,information_vector_1 ,&
                                                    information_vector_2, &
                                                    Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
	                                                this_flag1, brother_vertices_Rectangular, trial_basis_type, j, &
											        test_basis_type, k, edge_index(n), con_penalty, int_value, HE(6,n))
                    
                        El_Stiff_b_b(j,k) = El_Stiff_b_b(j,k) + int_value
                    
                    enddo
                enddo
                
                !write(*,*) 'interface element should be a DG element'
                !stop
            end if
            
            
        elseif (element_index(this_element) > 0 .AND. element_index(next_element) > 0) then 
            
            information_vector_1 = information_1(:, element_index(this_element))
            information_vector_2 = information_2(:, element_index(this_element))
            information_vector_1_next = information_1(:, element_index(next_element))
            information_vector_2_next = information_2(:, element_index(next_element))
            
            if (edge_index(n) > 0) then
                
                CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, &
                                    Gauss_point_reference_1D, vertices_linear(:,1), &
                                    information_3(4:5, edge_index(n)), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
                
                loc1 = information_3(6, edge_index(n))
                
                do j = 1, 4
                    do k = 1, 4
                        
                        
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = 1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = -1
                        end if
                        CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1, &
                                            information_vector_2, information_vector_1, information_vector_2, &
                                            Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                            this_flag1,vertices_Rectangular, vertices_Rectangular, trial_basis_type, j, &
                                            test_basis_type, k, loc1, con_penalty, int_value, HE(5,n), HE(5,n))
                        
                        El_Stiff(j,k) = El_Stiff(j,k) + int_value   ! unsure
                    enddo
                enddo
                
                do j = 1, 4
                    do k = 1, 4
                        
                        
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = -1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = 1
                        end if
                        CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1, &
                                            information_vector_2, information_vector_1_next, information_vector_2_next, &
                                            Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                  this_flag1,vertices_Rectangular, brother_vertices_Rectangular, trial_basis_type, j, &
                                            test_basis_type, k, loc1, con_penalty, int_value, HE(5,n), HE(6,n))
                        
                        El_Stiff_brother(j,k) = El_Stiff_brother(j,k) + int_value   ! unsure
                    enddo
                enddo
                
                
                
                CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
                                    information_3(4:5, edge_index(n)),&
                                    vertices_linear(:,2), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
                
                loc2 = information_3(7, edge_index(n))
                
                do j = 1, 4
                    do k = 1, 4
                        
                        
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = 1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = -1
                        end if
                        CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1, &
                                                information_vector_2, information_vector_1, information_vector_2, &
                                          Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                        this_flag1,vertices_Rectangular, vertices_Rectangular, trial_basis_type, j, &
											    test_basis_type, k, loc2, con_penalty, int_value, HE(5,n), HE(5,n))
                        
                        El_Stiff(j,k) = El_Stiff(j,k) + int_value   ! unsure
                    enddo
                enddo
                
                do j = 1, 4
                    do k = 1, 4
                        
                        
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = -1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = 1
                        end if
                        CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1, &
                                            information_vector_2, information_vector_1_next, information_vector_2_next, &
                                            Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                  this_flag1,vertices_Rectangular, brother_vertices_Rectangular, trial_basis_type, j, &
                                            test_basis_type, k, loc2, con_penalty, int_value, HE(5,n), HE(6,n))
                        
                        El_Stiff_brother(j,k) = El_Stiff_brother(j,k) + int_value   ! unsure
                    enddo
                enddo
                
                if (HT(5, next_element) == 0) then
                    
                    CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, &
                                    Gauss_point_reference_1D, vertices_linear(:,1), &
                                    information_3(4:5, edge_index(n)), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
                
                    loc1 = information_3(6, edge_index(n))
                
                    do j = 1, 4
                        do k = 1, 4
                        
                        
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = -1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = 1
                            end if
                            CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1_next, &
                                                information_vector_2_next, information_vector_1_next, information_vector_2_next, &
                                                Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                                this_flag1,brother_vertices_Rectangular, brother_vertices_Rectangular, &
                                                trial_basis_type, j, &
                                                test_basis_type, k, loc1, con_penalty, int_value, HE(6,n), HE(6,n))
                        
                            El_Stiff_b_b(j,k) = El_Stiff_b_b(j,k) + int_value   ! unsure
                        enddo
                    enddo
                
                    do j = 1, 4
                        do k = 1, 4
                        
                        
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = 1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = -1
                            end if
                            CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1_next, &
                                                information_vector_2_next, information_vector_1, information_vector_2, &
                                                Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                                this_flag1,brother_vertices_Rectangular, vertices_Rectangular,&
                                                trial_basis_type, j, &
                                                test_basis_type, k, loc1, con_penalty, int_value, HE(6,n), HE(5,n))
                        
                            El_Stiff_b_t(j,k) = El_Stiff_b_t(j,k) + int_value   ! unsure
                        enddo
                    enddo
                    
                    CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
                                    information_3(4:5, edge_index(n)),&
                                    vertices_linear(:,2), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
                
                    loc2 = information_3(7, edge_index(n))
                        
                    
                    do j = 1, 4
                        do k = 1, 4
                        
                        
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = -1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = 1
                            end if
                            CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1_next, &
                                                information_vector_2_next, information_vector_1_next, information_vector_2_next, &
                                                Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                                this_flag1,brother_vertices_Rectangular, brother_vertices_Rectangular, &
                                                trial_basis_type, j, &
                                                test_basis_type, k, loc2, con_penalty, int_value, HE(6,n), HE(6,n))
                        
                            El_Stiff_b_b(j,k) = El_Stiff_b_b(j,k) + int_value   ! unsure
                        enddo
                    enddo
                
                    do j = 1, 4
                        do k = 1, 4
                        
                        
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = 1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = -1
                            end if
                            CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1_next, &
                                                information_vector_2_next, information_vector_1, information_vector_2, &
                                                Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                                this_flag1,brother_vertices_Rectangular, vertices_Rectangular, &
                                                trial_basis_type, j, &
                                                test_basis_type, k, loc2, con_penalty, int_value, HE(6,n), HE(5,n))
                        
                            El_Stiff_b_t(j,k) = El_Stiff_b_t(j,k) + int_value   ! unsure
                        enddo
                    enddo
                    
                    write(*,*) 'wsy revise'
                    !stop
                End IF
                
            elseif (edge_index(n) < 0) then
                
                CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
                                    vertices_linear(:,1), &
                                    vertices_linear(:,2), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
                
                do j = 1, 4
                    do k = 1, 4
                        
                        
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = 1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = -1
                        end if
                        CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1, &
                                                information_vector_2, information_vector_1, information_vector_2, &
                                                Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                            this_flag1,vertices_Rectangular, vertices_Rectangular, trial_basis_type, j, &
											    test_basis_type, k, edge_index(n), con_penalty, int_value, HE(5,n), HE(5,n))
                        
                        El_Stiff(j,k) = El_Stiff(j,k) + int_value   ! unsure
                    enddo
                enddo
                
                do j = 1, 4
                    do k = 1, 4
                        
                        
                        this_flag1(2:3,1) = HE(3:4, n)
                        this_flag1(2:3,2) = HE(3:4, n)
                        if (HE(5, n) > HE(6, n)) then
                            this_flag1(4, 1) = 1
                            this_flag1(4, 2) = -1
                        else
                            this_flag1(4, 1) = -1
                            this_flag1(4, 2) = 1
                        end if
                        CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1, &
                                            information_vector_2, information_vector_1_next, information_vector_2_next, &
                                            Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                  this_flag1,vertices_Rectangular, brother_vertices_Rectangular, trial_basis_type, j, &
                                            test_basis_type, k, edge_index(n), con_penalty, int_value, HE(5,n), HE(6,n))
                        
                        El_Stiff_brother(j,k) = El_Stiff_brother(j,k) + int_value   ! unsure
                    enddo
                enddo
                
                if (HT(5, next_element) == 0) then
                    
                    do j = 1, 4
                        do k = 1, 4
                        
                        
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = -1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = 1
                            end if
                            CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1_next, &
                                                    information_vector_2_next, information_vector_1_next,&
                                                    information_vector_2_next, &
                                                    Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
	                                                this_flag1,brother_vertices_Rectangular, brother_vertices_Rectangular, &
                                                    trial_basis_type, j, &
											        test_basis_type, k, edge_index(n), con_penalty, int_value, HE(6,n), HE(6,n))
                        
                            El_Stiff_b_b(j,k) = El_Stiff_b_b(j,k) + int_value   ! unsure
                        enddo
                    enddo
                    
                    do j = 1, 4
                        do k = 1, 4
                        
                        
                            this_flag1(2:3,1) = HE(3:4, n)
                            this_flag1(2:3,2) = HE(3:4, n)
                            if (HE(5, n) > HE(6, n)) then
                                this_flag1(4, 1) = -1
                                this_flag1(4, 2) = 1
                            else
                                this_flag1(4, 1) = 1
                                this_flag1(4, 2) = -1
                            end if
                            CALL gauss_integration_local_stiffness_IDG_both_O(delta, Global_Beta, information_vector_1_next, &
                                                information_vector_2_next, information_vector_1, information_vector_2, &
                                                Gauss_Coefficient_Local_1D_Linear_Eight,Gauss_Point_Local_1D_Linear_Eight, &
                                                this_flag1,brother_vertices_Rectangular, vertices_Rectangular, &
                                                trial_basis_type, j, &
                                                test_basis_type, k, edge_index(n), con_penalty, int_value, HE(6,n), HE(5,n))
                        
                            El_Stiff_b_t(j,k) = El_Stiff_b_t(j,k) + int_value   ! unsure
                        enddo
                    enddo
                    
                End If
                
            end if
            
        end if
            
        
        !=======================STORE MATRIX================================
        i = HE(5,n)

    	DO k=1,n_nodes_in_elem  
		    IF ( node_type(2,HT(k,i)) > 0 ) THEN	! ˵���õ����ڲ��ĵ� iΪ��ǰ������Ԫ���
		    ! This is an Unknown node           
    		    DO j=1,n_nodes_in_elem
    			    IF ( node_type(2,HT(j,i)) > 0 ) THEN	
				    ! This is an Unknown node
                        IF (size(matrix%K)==1)  Then
                            m=matrix( node_type(2,HT(k,i)) )%SROW
					        IF( matrix(m)%JCOL==node_type(2,HT(j,i)) )THEN
							        matrix(m)%K = matrix(m)%K+El_Stiff(k,j)
							        EXIT
					        END IF
					    ElSE
					        DO m=matrix( node_type(2,HT(k,i)) )%SROW,  &
							        matrix( node_type(2,HT(k,i))+1 )%SROW-1
						        IF( matrix(m)%JCOL==node_type(2,HT(j,i)) )THEN
							        matrix(m)%K = matrix(m)%K+El_Stiff(k,j)
							        EXIT
						        END IF
					        END DO
                        ENDIF
                    ELSEIF( node_type(2,HT(j,i)) < 0 ) THEN
                        ! This is an known node��1
                        IF (size(matrix_xt%K)==1)  Then
                            m=matrix_xt( node_type(2,HT(k,i)) )%SROW
                            IF( matrix_xt(m)%JCOL==node_type(2,HT(j,i)) )THEN
                                matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff(k,j)
							    EXIT
                            END IF
                        ElSE
                            DO m=matrix_xt( node_type(2,HT(k,i)) )%SROW, &
                                matrix_xt( node_type(2,HT(k,i))+1 )%SROW-1
                                IF( matrix_xt(m)%JCOL==HT(j,i) )THEN
                                    matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff(k,j)
                                    EXIT
                                END IF
                            END DO
                        ENDIF
				    END IF
			    END DO
		    END IF
        END DO
        
        
        !=======================STORE MATRIX(The element of brother edge)================================
        i = HE(6,n)
        ii = HE(5,n)
        
    	DO k=1,n_nodes_in_elem  
		    IF ( node_type(2,HT(k,ii)) > 0 ) THEN	! ˵���õ����ڲ��ĵ� iΪ��ǰ������Ԫ���
		    ! This is an Unknown node           
    		    DO j=1,n_nodes_in_elem
    			    IF ( node_type(2,HT(j,i)) > 0 ) THEN	
				    ! This is an Unknown node
                        IF (size(matrix%K)==1)  Then
                            m=matrix( node_type(2,HT(k,ii)) )%SROW
					        IF( matrix(m)%JCOL==node_type(2,HT(j,i)) )THEN
							        matrix(m)%K = matrix(m)%K+El_Stiff_brother(k,j)
							        EXIT
					        END IF
					    ElSE
					        DO m=matrix( node_type(2,HT(k,ii)) )%SROW,  &
							        matrix( node_type(2,HT(k,ii))+1 )%SROW-1
						        IF( matrix(m)%JCOL==node_type(2,HT(j,i)) )THEN
							        matrix(m)%K = matrix(m)%K+El_Stiff_brother(k,j)
							        EXIT
						        END IF
					        END DO
                        ENDIF
                    ELSEIF( node_type(2,HT(j,i)) < 0 ) THEN
                        ! This is an known node��1
                        IF (size(matrix_xt%K)==1)  Then
                            m=matrix_xt( node_type(2,HT(k,ii)) )%SROW
                            IF( matrix_xt(m)%JCOL==node_type(2,HT(j,i)) )THEN
                                matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff_brother(k,j)
							    EXIT
                            END IF
                        ElSE
                            DO m=matrix_xt( node_type(2,HT(k,ii)) )%SROW, &
                                matrix_xt( node_type(2,HT(k,ii))+1 )%SROW-1
                                IF( matrix_xt(m)%JCOL==HT(j,i) )THEN
                                    matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff_brother(k,j)
                                    EXIT
                                END IF
                            END DO
                        ENDIF
				    END IF
			    END DO
		    END IF
        END DO
        
    
        i = HE(6,n)
        ii = HE(6,n)
        
    	DO k=1,n_nodes_in_elem  
		    IF ( node_type(2,HT(k,ii)) > 0 ) THEN	! ˵���õ����ڲ��ĵ� iΪ��ǰ������Ԫ���
		    ! This is an Unknown node           
    		    DO j=1,n_nodes_in_elem
    			    IF ( node_type(2,HT(j,i)) > 0 ) THEN	
				    ! This is an Unknown node
                        IF (size(matrix%K)==1)  Then
                            m=matrix( node_type(2,HT(k,ii)) )%SROW
					        IF( matrix(m)%JCOL==node_type(2,HT(j,i)) )THEN
							        matrix(m)%K = matrix(m)%K+El_Stiff_b_b(k,j)
							        EXIT
					        END IF
					    ElSE
					        DO m=matrix( node_type(2,HT(k,ii)) )%SROW,  &
							        matrix( node_type(2,HT(k,ii))+1 )%SROW-1
						        IF( matrix(m)%JCOL==node_type(2,HT(j,i)) )THEN
							        matrix(m)%K = matrix(m)%K+El_Stiff_b_b(k,j)
							        EXIT
						        END IF
					        END DO
                        ENDIF
                    ELSEIF( node_type(2,HT(j,i)) < 0 ) THEN
                        ! This is an known node��1
                        IF (size(matrix_xt%K)==1)  Then
                            m=matrix_xt( node_type(2,HT(k,ii)) )%SROW
                            IF( matrix_xt(m)%JCOL==node_type(2,HT(j,i)) )THEN
                                matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff_b_b(k,j)
							    EXIT
                            END IF
                        ElSE
                            DO m=matrix_xt( node_type(2,HT(k,ii)) )%SROW, &
                                matrix_xt( node_type(2,HT(k,ii))+1 )%SROW-1
                                IF( matrix_xt(m)%JCOL==HT(j,i) )THEN
                                    matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff_b_b(k,j)
                                    EXIT
                                END IF
                            END DO
                        ENDIF
				    END IF
			    END DO
		    END IF
        END DO
    
        i = HE(5,n)
        ii = HE(6,n)
        
    	DO k=1,n_nodes_in_elem  
		    IF ( node_type(2,HT(k,ii)) > 0 ) THEN	! ˵���õ����ڲ��ĵ� iΪ��ǰ������Ԫ���
		    ! This is an Unknown node           
    		    DO j=1,n_nodes_in_elem
    			    IF ( node_type(2,HT(j,i)) > 0 ) THEN	
				    ! This is an Unknown node
                        IF (size(matrix%K)==1)  Then
                            m=matrix( node_type(2,HT(k,ii)) )%SROW
					        IF( matrix(m)%JCOL==node_type(2,HT(j,i)) )THEN
							        matrix(m)%K = matrix(m)%K+El_Stiff_b_t(k,j)
							        EXIT
					        END IF
					    ElSE
					        DO m=matrix( node_type(2,HT(k,ii)) )%SROW,  &
							        matrix( node_type(2,HT(k,ii))+1 )%SROW-1
						        IF( matrix(m)%JCOL==node_type(2,HT(j,i)) )THEN
							        matrix(m)%K = matrix(m)%K+El_Stiff_b_t(k,j)
							        EXIT
						        END IF
					        END DO
                        ENDIF
                    ELSEIF( node_type(2,HT(j,i)) < 0 ) THEN
                        ! This is an known node��1
                        IF (size(matrix_xt%K)==1)  Then
                            m=matrix_xt( node_type(2,HT(k,ii)) )%SROW
                            IF( matrix_xt(m)%JCOL==node_type(2,HT(j,i)) )THEN
                                matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff_b_t(k,j)
							    EXIT
                            END IF
                        ElSE
                            DO m=matrix_xt( node_type(2,HT(k,ii)) )%SROW, &
                                matrix_xt( node_type(2,HT(k,ii))+1 )%SROW-1
                                IF( matrix_xt(m)%JCOL==HT(j,i) )THEN
                                    matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff_b_t(k,j)
                                    EXIT
                                END IF
                            END DO
                        ENDIF
				    END IF
			    END DO
		    END IF
        END DO
    
    end if
    
enddo


!=======================STORE MATRIX================================



END