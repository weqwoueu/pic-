!Main
!SUBROUTINE IFE_Start_2D(	IFE_Mesh_Filename,									&
!						    PIC_BC_Index,										&
!						    eps_objects, nobjects,								&
!						    PIC_phi_bkgd, PIC_Te_bkgd, PIC_Rho_bkgd, nbkgd,		&
!						    i_grid_flag,										&
!						    Phi, nnx, nny, delta		)

SUBROUTINE IFE_Start_2D(	PIC_phi_bkgd, PIC_Te_bkgd, PIC_Rho_bkgd, nbkgd,		&
						    Phi, nnx, nny, delta, time		)

! Purpose:		Initialize the IFE solver. Define the mesh information, the
!				interface elements, the stiffenss matrix and the boundary conditions.
!				Modified to work with the F77 PIC
! Last Updated:	10/11/2003 05:15 PM

USE IMPIC_Data_2D
USE IFE_MAIN_PARAM
USE IFE_RHS_PARAM
USE IFE_Data
USE IFE_Boundary
USE Object_Data_2D
USE Gauss_Data
USE TimeControl
USE IFE_INTERFACE, ONLY: Load_IFE_Mesh_2D, Partition_Cube_2D, Generate_Gauss_reference_triangle, Generate_Gauss_reference, &
						 Generate_stiffness_matrix_local_IFE, Sparse_Structure, Classify_Nodes_2D, Get_Boundary_Elements_2D, &
						 Diagonal_Preconditioner, Global_RHS_EBC_2D, Global_RHS_NBC_2D_dj, generate_Gauss_reference_1D, &
             Sparse_Structure_xt, generate_stiffness_matrix_local_IDG_D, generate_stiffness_matrix_local_IDG_O
USE Domain_2D, ONLY: VertX                          !$ ab.ZWZ for check node_type field
USE Constant_Variable_2D, ONLY: Phi_ref, time_ref !$ ab.ZWZ

IMPLICIT NONE

INTEGER					nbkgd				!, nobjects

REAL(8)					PIC_phi_bkgd(nbkgd), PIC_Te_bkgd(nbkgd), PIC_Rho_bkgd(nbkgd)

! The following declarations are true in case of using ghost cells at 0 and nnx+1
INTEGER					nnx, nny

!REAL(8)					Phi(0:nnx+1,0:nny+1)
REAL(8) :: Phi(Field_Size, 1)     !WSY REVISE, 2022-7-22

INTEGER                 delta	 !  0: 2-D Cartesian coordinates, 1: cylindrical coordinate, axisymmetric situation
REAL(8),INTENT(IN) :: time !$ ab.ZWZ for AC objects

INTEGER                                              Gauss_point_number, Gauss_point_number_1D
REAL(8),DIMENSION(:),POINTER                   ::    Gauss_coefficient_reference
REAL(8),DIMENSION(:,:),POINTER                 ::    Gauss_point_reference
INTEGER                                              Gauss_point_number_triangle
REAL(8),DIMENSION(:),POINTER                   ::    Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:),POINTER                 ::    Gauss_point_reference_triangle
REAL(8),DIMENSION(:,:),POINTER                 ::    A1, A2, A
REAL(8),DIMENSION(:),POINTER                   ::    b
INTEGER,DIMENSION(:,:),POINTER                 ::    pointer_reference_to_local,BETA_SIGN


REAL(8)                                              function_beta_minus, function_beta_plus

REAL(8)                                              function_f_minus,function_f_plus


!INTEGER		ni, nj, num_of_nodes, num_of_unknowns, &
!			i, j, k, nindex
INTEGER		ni, nj, i, j, k, nindex
INTEGER, SAVE                                   ::   num_of_nodes, num_of_unknowns  !$ mb.ZWZ for calling IFE_Start multiple times
REAL(8), SAVE	                                ::   dimensions(2,2) !, CPUtime

INTEGER, DIMENSION(:), POINTER				::  el_type!, bnd_elem_index

REAL(8), DIMENSION(:), POINTER				::	rhs_BGS, U_fe_full, U_fe !, eps_objects
REAL(8), DIMENSION(:), POINTER				::	rhs_NBC

EXTERNAL										FUN_2D_ONE
INTEGER											nodes(3), n_elements
INTEGER,Save ::									N_Objects, object_region
TYPE(ObjectType), DIMENSION(:), POINTER, SAVE		::	objects


INTEGER, DIMENSION(:,:), POINTER	::	VROW			
TYPE(SPARSE), DIMENSION(:), POINTER,SAVE	::	matrix, matrix1, matrix2, matrix_D, matrix_O ! wsy add for SIDG
TYPE(SPARSE), DIMENSION(:), POINTER	::	matrix3, matrix4  !!! bjw add for impic 2019-6-3
TYPE(SPARSE), DIMENSION(:), ALLOCATABLE	,SAVE	::	matrix_xt, matrix5, matrix6, matrix7, matrix8, matrix_Dt, matrix_Ot ! wsy add for SIDG
INTEGER                                 NZ, e, n_nodes_in_elem, m_unknowns, be

INTEGER, SAVE :: IFE_Init_flag=0  !$ ab.ZWZ for adding IFE_start into loop
INTEGER, SAVE :: Recall_flag=0  !$ ab.ZWZ for adding IFE_start into loop
Logical, Save :: StiffUpdateFlag = .True. !> ab.ZWZ for updating dirhichlet value without updating the stiffness matrix


    !OPEN(1, ACTION = 'WRITE', FILE = './OUTPUT/Phi.dat')
    !WRITE(1,*) 'TITLE = "Field Plot"'
    !WRITE(1,*) 'VARIABLES = "x" "y" "Phi" '
    !WRITE(1,4000) nnx, nny
	   ! DO j=1,nny
		  !  DO i=1,nnx
			 !   WRITE(1,5000) REAL(i-1), REAL(j-1), Phi(i,j)
		  !  END DO
	   ! END DO
    !4000 FORMAT (' ZONE I = ',I6,', J= ',I6)
    !5000 FORMAT (E15.6,' ',E15.6,' ',E15.6)
    !CLOSE(1)

IF(IFE_Init_flag == 0 )THEN !$ ab.ZWZ temporily

    WRITE(6,*)
    WRITE(6,*) 'Initialize IFE Solver'
    WRITE(6,*) '===================== '

    !IFE_Mesh_Filename = 'ife.msh'
    ! Get IFE Main Parameters from ife.inp file
    OPEN(1, ACTION='READ', FILE = './INPUT/ife.inp')
    READ(1,*)	! Skip IFE mesh resave option
    READ(1,*)	! Skip Int_El_Frac
    READ(1,*)	NSolver
    READ(1,*)	SLSolver
    READ(1,*)	Blocking
    READ(1,*)	PCG_MaxIT
    READ(1,*)	PCG_Tol
    READ(1,*)	BGS_MaxIT
    READ(1,*)	BGS_Tol
    READ(1,*)	Newton_MaxIT
    READ(1,*)	Newton_Tol
    CLOSE(1)

    ! Get IFE objects Parameters from object.inp file
    OPEN(1, ACTION='READ', FILE='./INPUT/object.inp')
    READ(1,*) N_Objects,vacuum

    IF( N_Objects /= 0) THEN
      WRITE(6,*)
      WRITE(6,*)'Number of Object=',N_Objects,' vacuum=',vacuum
      ALLOCATE(objects(N_Objects))
      DO i = 1, N_Objects
        objects(i)%Shape		=0
        objects(i)%Axis			=0
        objects(i)%Dimensions	=0
        objects(i)%Locations	=0
        objects(i)%Regions		=0
        objects(i)%Phi			=0
        objects(i)%Eps			=0
        READ(1,*) objects(i)%Shape
        READ(1,*) objects(i)%Axis
        READ(1,*) objects(i)%Dimensions(:)
        READ(1,*) objects(i)%Locations(1,:)
        READ(1,*) objects(i)%Locations(2,:)
        !=========LY modification, 2022-6-13=========
        READ(1,*) objects(i)%Locations(3,:)
        READ(1,*) objects(i)%Locations(4,:)
        !=========LY modification, 2022-6-13=========
        READ(1,*) objects(i)%Regions
        READ(1,*) objects(i)%Direction
        READ(1,*) objects(i)%Phi
        READ(1,*) objects(i)%Eps
        READ(1,*) objects(i)%Erosion
        IF (objects(i)%Erosion > 0) THEN
          READ(1,*) objects(i)%Wall(1)%Shape
          READ(1,*) objects(i)%Wall(1)%Channelwall
          READ(1,*) objects(i)%Wall(1)%Limits(1,1:2)
          READ(1,*) objects(i)%Wall(1)%Limits(2,1:2)
          READ(1,*) objects(i)%Wall(1)%stepx
        ELSE
          READ(1,*)
          READ(1,*)
          READ(1,*)
          READ(1,*)
          READ(1,*)
        ENDIF
      END DO
    END IF
    CLOSE(1)
    
    !-----------------------LY REVISE, 2021-11-27--------------------------
    OPEN(1,ACTION='READ',FILE='./INPUT/IDG_inf.inp')
	    READ(1, *) con_penalty
    CLOSE(1)
    !----------------------------------------------------------------------

    
    ALLOCATE(Global_Beta(N_Objects+1))	! this ONE accounts for the vacuum epsilon
    !==========================test===================================
    !=========PIC=========
    Global_Beta(1)	= one
    !=========PIC=========
    !=========IFE=========
    !Global_Beta(1)	= 10
    !=========IFE=========
    !==========================test===================================
    IF( N_Objects /= 0) THEN
	    DO i = 1, N_Objects
		    object_region = -objects(i)%Regions(1)
		    Global_Beta(object_region) = objects(i)%Eps
	    END DO
    ENDIF

    IF (.NOT.ALLOCATED(phi_bkgd)) ALLOCATE(phi_bkgd(SIZE(PIC_phi_bkgd)))
    IF (.NOT.ALLOCATED(Te_bkgd)) ALLOCATE(Te_bkgd(SIZE(PIC_Te_bkgd)))
    IF (.NOT.ALLOCATED(Rho_bkgd)) ALLOCATE(Rho_bkgd(SIZE(PIC_Rho_bkgd)))
    !ALLOCATE(phi_bkgd(SIZE(PIC_phi_bkgd)), Te_bkgd(SIZE(PIC_Te_bkgd)), Rho_bkgd(SIZE(PIC_Rho_bkgd)))
    phi_bkgd = PIC_phi_bkgd
    Te_bkgd  = PIC_Te_bkgd
    Rho_bkgd = PIC_Rho_bkgd
    !
    WRITE(6,*) '111111111111111'
    !========================NEW ADD=============================================
    OPEN(1, ACTION='READ', FILE='./INPUT/mesh.inp')
    READ(1,*) , 
    READ(1,*) , 
    READ(1,*) nnx, nny
    CLOSE(1)
    !============================================================================


    WRITE(6,*) 'Mesh Loading ....'

    IF(irestart==0)THEN
    CALL Load_IFE_Mesh_2D(	'ife.msh',		    &
						    t_c, p_basic, e_basic,  &
                n_elements, information_1, information_2, information_3, information_3_D,&
                element_index, node_index, edge_index, edge_index_D, HP, HT, HE, E_Dirichlet, P_average, P_flag)
    ELSEIF(irestart==1)THEN
    CALL Load_IFE_Mesh_2D(	'update_ife.msh',		    &
                t_c, p_basic, e_basic,  &
                n_elements, information_1, information_2, information_3, information_3_D,&
                element_index, node_index, edge_index, edge_index_D, HP, HT, HE, E_Dirichlet, P_average, P_flag)

    ENDIF

 
dimensions(1,1) = MINVAL(HP(1,:))
dimensions(2,1) = MINVAL(HP(2,:))
dimensions(1,2) = MAXVAL(HP(1,:))
dimensions(2,2) = MAXVAL(HP(2,:))


num_of_nodes	= SIZE(HP,2)
n_elements    = SIZE(HT,2)

!IF(IFE_Init_flag == 0 )THEN !$ ab.ZWZ temporily

    CALL Partition_Cube_2D(HP, e_basic, dimensions, &
				  	       node_type, num_of_unknowns)

    CALL Get_Boundary_Elements_2D(n_elements, HT, node_type, node_index, bnd_elem_index)

    IF (N_Objects  /=0) THEN
      CALL Classify_Nodes_2D( HP, node_type, node_index, node_permute )
    ENDIF

    Gauss_point_number = 9

    Gauss_point_number_triangle = 9

    Gauss_point_number_1D = 8
    
    CALL Generate_Gauss_reference( Gauss_point_number, Gauss_Coefficient_Reference_Nine, Gauss_Point_Reference_Nine)

    CALL Generate_Gauss_reference_triangle( Gauss_point_number_triangle, Gauss_Coefficient_Reference_Triangle_Nine, &
                                            Gauss_Point_Reference_Triangle_Nine)

    CALL generate_Gauss_reference_1D(Gauss_point_number_1D, Gauss_coefficient_reference_1D_Eight, Gauss_point_reference_1D_Eight) !LY REVISE, 2021-11-27, Line Integration

    
    num_of_unknowns = MAXVAL(node_type(2,:))

    WRITE(6,*) 'No of Mesh Elements      = ',SIZE(HT,2)
    WRITE(6,*) 'No of Mesh Nodes         = ',SIZE(HP,2)
    WRITE(6,*) 'No of Boundary Elements  = ',SIZE(bnd_elem_index)

    WRITE(6,*) 'No of Unknown Nodes      = ',num_of_unknowns
    WRITE(6,*) 'No of Interface Elements = ',SIZE(information_1,2)
    WRITE(6,*)


    WRITE(6,*) 'Mesh Generation Done.'

    IFE_Init_flag = 1   !$ ab.ZWZ temporily
    
    
    
ELSE
    Recall_flag = 1
ENDIF   !$ ab.ZWZ temporily


!=========================================NEW=====================================================
ALLOCATE(U_fe_full(num_of_nodes), U_fe(num_of_unknowns))    !$ A U = R
U_fe_full = Zero

! wsy revise for SIDG
DO ni=1,SIZE(HP, 2)
  U_fe_full(ni) = Phi(ni,1)
END DO


!! Exract U_fe (unknown nodes) from U_fe_full (all nodes)
DO k=1,num_of_nodes
	IF (node_type(2, k) > 0) THEN	! Unknown node
		U_fe(node_type(2,k)) = U_fe_full(k)
	END IF
END DO

WRITE(6,*) "Max U_fe      = ",MAXVAL(U_fe)
WRITE(6,*) "Min U_fe      = ",MINVAL(U_fe)
WRITE(6,*) "Max U_fe_full = ",MAXVAL(U_fe_full)
WRITE(6,*) "Min U_fe_full = ",MINVAL(U_fe_full)
WRITE(6,*)

!======================================================================================================

function_beta_minus = Global_Beta(1)
!function_beta_plus = Global_Beta(1)
If (StiffUpdateFlag) Then !ab.ZWZ

    WRITE(6,*) 'Constructing Sitffness Matrix ....'

    !-------------------------------------initiate the sparse---------------------------------------------------
    !n_nodes_in_elem = 4
    !IF(IFE_init_flag == 0 ) Then
        CALL Sparse_Structure( HT, HE, node_type, VROW, matrix, NZ )
        CALL Sparse_Structure_xt( HT, HE, node_type, VROW, matrix_xt, NZ ) !$ ab.ZWZ 2021/7/9
    
        CALL Sparse_Structure( HT, HE, node_type, VROW, matrix1, NZ ) 
        CALL Sparse_Structure_xt( HT, HE, node_type, VROW, matrix5, NZ ) 
    
        CALL Sparse_Structure(  HT, HE,  node_type, VROW, matrix2, NZ ) 
        CALL Sparse_Structure_xt(  HT, HE,  node_type, VROW, matrix6, NZ ) 
    
        IF (Bfiled_index) Then
            CALL Sparse_Structure(  HT, HE,  node_type, VROW, matrix3, NZ ) 
            CALL Sparse_Structure_xt(  HT, HE,  node_type, VROW, matrix7, NZ ) 
        
            CALL Sparse_Structure(  HT, HE,  node_type, VROW, matrix4, NZ )
            CALL Sparse_Structure_xt(  HT, HE,  node_type, VROW, matrix8, NZ )
        End If
    
    !    IFE_init_flag = 1
    !End If

    !> -- ab.ZWZ to zero out all the sparse matrix value
    DO i=1,SIZE(matrix,1)
	    matrix(i)%K = 0.
	    matrix1(i)%K = 0.
        matrix2(i)%K = 0.
        IF (Bfiled_index) Then
	        matrix3(i)%K = 0.
	        matrix4(i)%K = 0.
        End If
    ENDDO

    DO i=1,SIZE(matrix_xt,1)
        matrix_xt(i)%K = 0.
        matrix5(i)%K = 0.
        matrix6(i)%K = 0.
        IF (Bfiled_index) Then
            matrix7(i)%K = 0.
            matrix8(i)%K = 0.
        End If
    ENDDO
    !> -- ab.ZWZ to zero out all the sparse matrix

    !-------------------------------------initiate the sparse---------------------------------------------------


    CALL Generate_stiffness_matrix_local_IFE( function_beta_minus, Global_Beta, node_type,num_of_unknowns, &
                                              information_1, information_2, HP, element_index, &
								   	                          HT, nnx, nny, dimensions, &
										                          Gauss_coefficient_reference_Nine, &
		   					                              Gauss_point_reference_Nine, Gauss_coefficient_reference_triangle_Nine, &
                                              Gauss_point_reference_triangle_Nine, &										  								      
									                            1, 1, 0, 1, 1, 0, matrix1, matrix5, delta)    !!! xx


    CALL Generate_stiffness_matrix_local_IFE( function_beta_minus, Global_Beta, node_type,num_of_unknowns, &
                                              information_1, information_2, HP, element_index, &
								   	                          HT, nnx, nny, dimensions, &
										                          Gauss_coefficient_reference_Nine, &
		   					                              Gauss_point_reference_Nine, Gauss_coefficient_reference_triangle_Nine, &
                                              Gauss_point_reference_triangle_Nine, &
									                            1, 0, 1, 1, 0, 1, matrix2, matrix6, delta)   !!! yy

    !=========================================WSY ADD, 2022-7-22=====================================
    CALL generate_stiffness_matrix_local_IDG_D( Global_Beta, node_type,num_of_unknowns, &
                                                information_1, information_2, information_3_D, HP, element_index, edge_index_D,&
                                                HT, E_Dirichlet, HE, dimensions, &
                                                Gauss_coefficient_reference_1D_Eight, Gauss_point_reference_1D_Eight,&
                                                1, 1, matrix_D, matrix_Dt, delta, con_penalty)          !LY modification.

    CALL generate_stiffness_matrix_local_IDG_O( Global_Beta, node_type,num_of_unknowns, &
                                                information_1, information_2, information_3, HP, element_index, edge_index,&
									                              HT, HE, dimensions, &
                                                Gauss_coefficient_reference_1D_Eight, Gauss_point_reference_1D_Eight,&
                                                1, 1, matrix_O, matrix_Ot, delta, con_penalty)          !LY modification.
    !===================================================================================================


    !!! ************************ bjw add for impic 2019-6-3 **********************************************
    IF (Bfiled_index .AND. IMPIC_index) THEN   !!!! 有磁场
        CALL Generate_stiffness_matrix_local_IFE( function_beta_minus, Global_Beta, node_type,num_of_unknowns, &
                                                  information_1, information_2, HP, element_index, &
								   	                              HT, nnx, nny, dimensions, &
										                              Gauss_coefficient_reference_Nine, &
		   					                                  Gauss_point_reference_Nine, Gauss_coefficient_reference_triangle_Nine, &
                                                  Gauss_point_reference_triangle_Nine, &										  								      
									                                1, 1, 0, 1, 0, 1, matrix3, matrix7, delta)  !!! xy

        CALL Generate_stiffness_matrix_local_IFE( function_beta_minus, Global_Beta, node_type,num_of_unknowns, &
                                                  information_1, information_2, HP, element_index, &
								   	                              HT, nnx, nny, dimensions, &
										                              Gauss_coefficient_reference_Nine, &
		   					                                  Gauss_point_reference_Nine, Gauss_coefficient_reference_triangle_Nine, &
                                                  Gauss_point_reference_triangle_Nine, &										  								      
									                                1, 0, 1, 1, 1, 0, matrix4, matrix8, delta)  !!! yx

         !!!constructing the Global IFE stiff
        DO i=1,SIZE(matrix,1)
	        matrix(i)%K = matrix1(i)%K + matrix2(i)%K + matrix3(i)%K + matrix4(i)%K + matrix_D(i)%K + matrix_O(i)%K
        ENDDO
    
        DO i=1,SIZE(matrix_xt,1)
            matrix_xt(i)%K = matrix5(i)%K + matrix6(i)%K + matrix7(i)%K + matrix8(i)%K + matrix_Dt(i)%K + matrix_Ot(i)%K
        ENDDO
    
        DEALLOCATE(matrix1, matrix2, matrix3, matrix4,matrix7,matrix8,matrix_D,matrix_Dt,matrix_O,matrix_Ot)
        !$ ============================ ab.ZWZ ============================== //
    
    ELSE  !!!! 无磁场
        !!constructing the Global IFE stiff
        ! the Improve SIDG method do not need to handle Dirichlet edge
        DO i=1,SIZE(matrix,1)
            matrix(i)%K = matrix1(i)%K + matrix2(i)%K + matrix_O(i)%K  ! + matrix_D(i)%K
        ENDDO
    
        DO i=1,SIZE(matrix_xt,1)
            matrix_xt(i)%K = matrix5(i)%K + matrix6(i)%K + matrix_Ot(i)%K ! + matrix_Dt(i)%K
        ENDDO
    
        DEALLOCATE(matrix1, matrix2, matrix5, matrix6, matrix_D,matrix_Dt,matrix_O,matrix_Ot)
    
    END IF

    !!! ************************ bjw add for impic 2019-6-3 **********************************************

    !!!constructing the Global IFE stiff
    !DO i=1,SIZE(matrix,1)
    !	matrix(i)%K = matrix1(i)%K + matrix2(i)%K
    !ENDDO

    !!! ************************ bjw add for impic 2019-6-3 **********************************************


    Allocate(A_stiff(size(matrix)))

    DO i=1,SIZE(matrix,1)
	    A_stiff(i)%K = matrix(i)%K
	    A_stiff(i)%JCOL = matrix(i)%JCOL
	    A_stiff(i)%SROW = matrix(i)%SROW
    ENDDO
    DEALLOCATE(matrix)  !> ab.ZWZ 2021/7/9
    !DEALLOCATE(matrix,matrix1,matrix2)
    ! ========== ab.ZWZ 2021/7/9 ============= \\
    Allocate(A_stiff_xt(size(matrix_xt)))
    DO i=1,SIZE(matrix_xt,1)
	    A_stiff_xt(i)%K = matrix_xt(i)%K
        !print*,A_stiff_xt(i)%K 
	    A_stiff_xt(i)%JCOL = matrix_xt(i)%JCOL
	    A_stiff_xt(i)%SROW = matrix_xt(i)%SROW
    ENDDO
    DEALLOCATE(matrix_xt)   
    ! =========== ab.ZWZ 2021/7/9 ============= //

    WRITE(6,*) 'Constructing Sitffness Matrix Done.'
    WRITE(6,*)

    !==============================GLOBAL STIFF DONE==============================

    IF (NSolver<=1) THEN	! Gauss-Seidel or Linear
	    CALL Diagonal_Preconditioner(A_stiff, A_diag)
    ENDIF

    If (.Not.IMPIC_index) StiffUpdateFlag = .False. !ab.ZWZ

Endif !ab.ZWZ

n_nodes_in_elem = 4 ! Retangular element

!ALLOCATE(el_type(size(information_1,2)))
!
!el_type(:)  =information_1(6,:)
!
!ALLOCATE(BETA_SIGN(size(information_1,2),2))

!$ ============== mb.ZWZ 2021/7/9 ============= \\
IF( N_Objects /= 0) THEN
    ALLOCATE(el_type(size(information_1,2)))
    el_type(:)  =information_1(6,:)
    ALLOCATE(BETA_SIGN(size(information_1,2),2))
ELSE
    ALLOCATE(el_type(2))
    el_type(:)  = 0
    ALLOCATE(BETA_SIGN(2,2))
END IF
!$ ============== mb.ZWZ 2021/7/9 ============= //

!------------------------------------------------------------------------------------



!WRITE(6,*) 'Forming Essential BC RHS Vector ....'
!!CALL Global_RHS_EBC_2D (    information_1,information_2, element_index, FUN_2D_ONE, U_fe_full, p_basic, t_c, e_basic,	&
!!					  	    el_type , p_int_x, p_int_y,  Global_Beta,	&
!!						    node_type, bnd_elem_index, G_Stiff_HW, rhs_EBC, delta)
!!DO k = 1,size(RHS_EBC)
!!    IF(rhs_EBC(k)>1D-5)THEN
!!        print*,rhs_EBC(k),'    ',k
!!    ENDIF
!!ENDDO
!WRITE(6,*) 'Forming Essential BC RHS Vector Done.'
!WRITE(6,*)




!!========================NBC====================================================================
WRITE(6,*) 'Forming Neumann BC RHS Vector ....'
CALL Global_RHS_NBC_2D_dj (element_index, el_type, information_1, information_2, HP, HT, e_basic,&
					  	   p_int_x, p_int_y,  Global_Beta,	&
						   node_type, bnd_elem_index, rhs_NBC, delta)	

!========================NBC====================================================================



!ALLOCATE(rhs_FIX(SIZE(rhs_EBC)))
If (Recall_flag == 0) Then          ! mb.ZWZ for avoid allocating rhs_Fix multiple times
    ALLOCATE(rhs_FIX(num_of_unknowns))    ! mb.ZWZ for moving EBC assembling to IFE_Solve
Endif

DO k=1,SIZE(rhs_FIX)
	rhs_FIX(k) = 0.D0
END DO


!LY REVISE, 2021-11-30
!This is FIX=FIX-IDG, then rhs = rhs-FIX = rhs-(FIX-IDG) = rhs + IDG -FIX on later.
DO k=1,SIZE(rhs_FIX)
	rhs_FIX(k) = rhs_FIX(k)+rhs_NBC(k)
  !rhs_FIX(k) = rhs_FIX(k)
END DO


!DO k=1,SIZE(rhs_FIX)
!	!rhs_FIX(k) = rhs_FIX(k)+rhs_EBC(k)+rhs_NBC(k)
!	rhs_FIX(k) = rhs_FIX(k)+rhs_NBC(k)  !$ mb.ZWZ for moving EBC assembling to IFE_Solve
!END DO

DO ni = 1, SIZE(HP,2)
  Phi(ni,1) = U_fe_full(ni)
END DO



!DEALLOCATE(rhs_EBC)
DEALLOCATE(rhs_NBC)
If (IMPIC_index == .False.) Then
    DEALLOCATE(bnd_elem_index)
End If
DEALLOCATE(el_type,BETA_SIGN)
DEALLOCATE(U_fe,U_fe_full)

WRITE(6,*)
WRITE(6,*) 'IFE Start Done!'
WRITE(6,*) '===================== '




END