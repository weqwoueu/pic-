SUBROUTINE IFE_Solve_2D( Rho, Phi, nnx, nny, delta, time)

USE IFE_MAIN_PARAM
USE IFE_RHS_PARAM
USE IFE_Data
USE Gauss_Data
USE IFE_Boundary
USE IFE_INTERFACE, ONLY: Global_RHS_PDE_2D, Global_Mass, Diagonal_Preconditioner, &
                         Generate_Gauss_reference_triangle, generate_Gauss_reference_1D, &
						 generate_load_vector_local_IFE_on_interface, &
						 generate_load_vector_local_IFE_nonhomogeneous, &
						 periodic_boundary_conditions, &
                         MY_JPCG_Solver, & ! ab.ZWZ 2021/7/9
                         Global_RHS_EBC_2D_bjw, generate_load_vector_local_IDG_D ! ab.ZWZ for moving EBC assembling to IFE_Solve

USE IMPIC_Data_2D, only: IMPIC_index ! ab.ZWZ 2021/11/9 if IFE_start added into loop
USE Constant_Variable_2D, ONLY: Phi_ref, time_ref ! wsy 2022 8 18 for Phi-Periodical case
USE TimeControl, Only : Frequency ! wsy 2022 8 18 for Phi-Periodical case
IMPLICIT NONE

INTEGER	nnx, nny, delta
!REAL(8)	Rho(0:nnx+1,0:nny+1), Phi(0:nnx+1,0:nny+1)
REAL(8) :: Rho(Field_Size, 1)     !wsy revise for sidg 2022 7 23
REAL(8) :: Phi(Field_Size, 1)     !wsy revise for sidg 2022 7 23
REAL(8), INTENT(IN) :: time

!====================================ADD===========================
REAL(8),DIMENSION(:,:),	POINTER		  ::	ALQ
INTEGER    num
!===============================================

INTEGER		ni, nj, nk, nindex, num_of_nodes, num_of_unknowns, ITER,&
			i, j, k
REAL(8)		Err_Iter, CPUtime, row_sum, xxx, yyy

REAL(8), DIMENSION(:), POINTER	::		rhs, R_fe_full, U_fe_full, U_fe, U_fe_old, rhs_res
EXTERNAL								FUN_Density_rhs_2D, FUN_Density_rhs_TRIANGLE_2D
!EXTERNAL								FUN_Density_Drhs_2D
INTEGER									MaxIT
REAL(8)									MaxTol, RelResidue, Residue, RHSResidue

INTEGER									max_err_node(1)

!=======================================================================================
EXTERNAL								FUN_sphere_u_ErrAna
REAL(8), DIMENSION(:), POINTER	::		Error_2, Error_0, IError_2
REAL(8), DIMENSION(:,:), POINTER::		Error_1, IError_1
!=========================================================================================

!==============================new add for surface jump CYC===============================
REAL(8), DIMENSION(:), POINTER				   ::		b2, b3, b4
EXTERNAL												GetFlux
REAL(8), DIMENSION(:), POINTER				   ::		Gauss_coefficient_reference_1D, Gauss_point_reference_1D
INTEGER													Gauss_point_number_triangle
REAL(8),DIMENSION(:),POINTER                   ::		Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:),POINTER                 ::		Gauss_point_reference_triangle
REAL(8), DIMENSION(:), POINTER				::	rhs_EBC, b_IDG !$ mb.ZWZ for moving EBC assembling to IFE_Solve
!=========================================================================================

!==============================new add for mass===========================================
!TYPE(SPARSE), DIMENSION(:), POINTER		::	A_mass
EXTERNAL								FUN_Density_Drhs_2D
!=========================================================================================


num_of_nodes	= SIZE(HP,2)
num_of_unknowns = SIZE(rhs_FIX)

ALLOCATE(R_fe_full(num_of_nodes),U_fe_full(num_of_nodes),	&
         U_fe(num_of_unknowns), U_fe_old(num_of_unknowns))

!------wsy revise for sidg 2022 7 23-------
DO ni=1,SIZE(HP,2)
  R_fe_full(ni) = Rho(ni,1)
  U_fe_full(ni) = Phi(ni,1)
END DO
!-------------------------------------------


!! Exract U_fe (unknown nodes) from U_fe_full (all nodes)
DO k=1,num_of_nodes
	IF (node_type(2, k) > 0) THEN	! Unknown node
		U_fe(node_type(2,k)) = U_fe_full(k)
	END IF
END DO


PRINT*, "Max Rho = ",MAXVAL(Rho)
PRINT*, "Min Rho = ",MINVAL(Rho)
PRINT*, "Max U_fe = ",MAXVAL(U_fe)
PRINT*, "Min U_fe = ",MINVAL(U_fe)
PRINT*, "Max U_fe_full = ",MAXVAL(U_fe_full)
PRINT*, "Min U_fe_full = ",MINVAL(U_fe_full)
PRINT*

DO k=1,SIZE(U_fe)
	U_fe_old(k) = U_fe(k)
END DO



!=========LY modification for Phi-Periodical case, 2022-8-11=========
!======Part1: Settle Phi-Periodical Potential======
!In subroutine 'SetObjects_2D', we construct the array 'PhiPeriodicalNode' to finding the node index in Phi-Periodical object.
!Do i = 1, Size(PhiPeriodicalNode)
  !Phi(PhiPeriodicalNode(i),1) = objects(1)%Phi + 200.0720*DSIN(2.0*pi*30.0D6*xt*(5.605381896D-11)) !time_ref = 5.605381896D-11
  !Phi(PhiPeriodicalNode(i),1) = objects(1)%Phi + 200.0720*DSIN(2.0*pi*30.0D6*xt*(5.605381896D-11)) !time_ref = 5.605381896D-11
!End Do

!======Part2: Update Injected Boundary Potential======
!bc_value(1) = 6.3023 + Phi(PhiPeriodicalNode(1),1)
!bc_value(1) = 450/Phi_ref*SIN(2*Pi*Frequency*time*time_ref) !mb.ZWZ for period changing D boundary

!======Part3: Update the E_Dirichlet(7,:) boundary value======
!In subroutine 'Setup_IFE_Mesh_2D', we construct the array 'E_DirichletValue' to store the Dirichlet boundary value.
!So, we can update the using boundary value according to E_DirichletValue and bc_value.
Do i = 1, Size(E_Dirichlet,2)
  If (E_Dirichlet(3,i)==-1 .AND. E_Dirichlet(4,i)==0 .AND. E_Dirichlet(6,i)==1) Then  !left Dirichlet boundary
    E_DirichletValue(i) = bc_value(1)
  End If
End Do
!=========LY modification for Phi-Periodical case, 2022-8-11=========


!WRITE(6,*) 'Forming Essential BC RHS Vector ....'
CALL Global_RHS_EBC_2D_bjw(A_stiff_xt, HP, HT, e_basic, node_type, U_fe_full, bnd_elem_index, rhs_EBC, delta)
!WRITE(6,*) 'Forming Essential BC RHS Vector Done.'


!---------------------------------LY REVISE, 2021-11-29------------------------------
!============There is the vector b_SIDG/DG============
WRITE(6,*) 'Forming SIDG Vector ....'
CALL generate_load_vector_local_IDG_D(Global_Beta, node_type, num_of_unknowns, &
                                      information_1, information_2, information_3_D, HP, element_index, edge_index_D, &
                                      HT, E_Dirichlet, &
                                      Gauss_coefficient_reference_1D_Eight, Gauss_point_reference_1D_Eight, &
                                      1, delta, con_penalty, b_IDG)     !LY modification.
WRITE(6,*) 'Forming SIDG Vector Done.'
WRITE(6,*)
!============There is the vector b_SIDG/DG============



IF (NSolver==0) THEN		! Linear Solver

	    PRINT*, 'Linear Solver :'
        
		!ALLOCATE(rhs_res(SIZE(rhs_FIX)), b2(SIZE(rhs_FIX)))
		ALLOCATE(rhs_res(SIZE(rhs_FIX)))

		!PRINT*, 'Forming PDE RHS Vector ....'

		CALL Global_RHS_PDE_2D(FUN_Density_rhs_2D,FUN_Density_rhs_TRIANGLE_2D, U_fe_full, R_fe_full, HP, HT,	&
							   h_partition, node_type, num_of_unknowns, G_RHS_HW, G_Interpol_HW, rhs, &
							   delta, element_index,information_1,information_2)

		!PRINT*, 'Forming PDE RHS Vector Done.'
  !      PRINT*

!	=================================================CYC ADD for surface jump=====================================================================================
		!Gauss_point_number_triangle = 9		
		!CALL Generate_Gauss_reference_triangle( Gauss_point_number_triangle, Gauss_coefficient_reference_triangle_Nine, &
		!                                        Gauss_point_reference_triangle_Nine)
		!
		!CALL generate_Gauss_reference_1D(4, Gauss_coefficient_reference_1D_Four, Gauss_point_reference_1D_Four)
		!
		!CALL generate_load_vector_local_IFE_on_interface(GetFlux, information_1, information_2,	&
		!												p_basic, element_index, t_c, nnx, nny,		&
		!												Gauss_coefficient_reference_1D_Four, Gauss_point_reference_1D_Four, 1, b2, node_type, node_index)
  !
		!CALL generate_load_vector_local_IFE_nonhomogeneous(GetFlux, Global_Beta, information_1, information_2,	&
		!													p_basic, element_index, t_c, nnx, nny, Gauss_coefficient_reference_1D_Four, Gauss_point_reference_1D_Four,	&
		!													Gauss_coefficient_reference_triangle_Nine, Gauss_point_reference_triangle_Nine, 1, 1, 0, 1, 1, 0, b3, node_type)
  !
  !
		!CALL generate_load_vector_local_IFE_nonhomogeneous(GetFlux, Global_Beta, information_1, information_2,	&
		!													p_basic, element_index, t_c, nnx, nny, Gauss_coefficient_reference_1D_Four, Gauss_point_reference_1D_Four,	&
		!													Gauss_coefficient_reference_triangle_Nine, Gauss_point_reference_triangle_Nine, 1, 0, 1, 1, 0, 1, b4, node_type)
  !
  !
		!DO i=1,SIZE(rhs,1)
		!	rhs(i) = rhs(i)-b2(i)-b3(i)-b4(i)
		!ENDDO
!	=================================================CYC ADD for surface jump=====================================================================================

		PRINT*, "Max rhs = ", MAXVAL(rhs)
		PRINT*, "Min rhs = ", MINVAL(rhs)
		DO j=1,SIZE(rhs)
			! rhs_FIX = rhs_EBC
			rhs(j)		= rhs(j)-rhs_FIX(j)-rhs_EBC(j) ! the ImproveSIDG method do not need to handle Dirichle edge
			rhs_res(j)	= rhs(j)
		END DO
		PRINT*, "Max rhsFIX = ", MAXVAL(rhs)
		PRINT*, "Min rhsFIX = ", MINVAL(rhs)
        PRINT*
        Deallocate(b_IDG,rhs_EBC)
        
!=============================================Solve Periodic Bound ary Conditions============================================
		num=0
		DO i=1,size(bc_index)
			IF(bc_index(i)==-1)THEN
				num=num+1
			ENDIF
		ENDDO
		
		IF(num>0) THEN
			ALLOCATE(ALQ(2,SIZE(bc_index)))
			OPEN(1, ACTION='READ', FILE='ALQ.msh')
			DO i=1,SIZE(bc_index)
				READ(1,*) ALQ(1,i),ALQ(2,i)	
			ENDDO
			CLOSE(1)
			CALL periodic_boundary_conditions(nnx,nny,bc_index,p_basic,node_type,bc_point_1,bc_point_2,ALQ,A_stiff,rhs)
			DEALLOCATE(ALQ)
		ENDIF
			
!!=========================ADD==================================================================================
 
		!PRINT*, 'Solving the Linear System ....'

		IF (Blocking==1) THEN
     !		CALL IMSL_JPCG_Solver_Permuted( A_mass, node_permute, rhs, U_fe, A_diag, ITER )
		ELSE
			IF(SLSolver==1) THEN
!!!				!CALL IMSL_JPCG_Solver( A_stiff, node_permute, rhs, U_fe, A_diag, ITER )
!!!				CALL IMSL_JPCG_Solver(A_stiff, rhs, U_fe, A_diag, SIZE(A_stiff,1), SIZE(rhs,1), ITER)
			ELSEIF(SLSolver==2) THEN
!!!				CALL DLAP_JPCG_Solver(A_stiff, rhs, U_fe, SIZE(A_stiff,1), SIZE(rhs,1), ITER)
            ELSEIF(SLSolver==3) THEN
                !print*,A_stiff
                !print*,''
                !print*,rhs
				CALL MY_JPCG_Solver(A_stiff, rhs, U_fe, A_diag, SIZE(A_stiff,1), SIZE(rhs,1), ITER)				
			END IF
		END IF

		PRINT*, 'Solving the Linear System Done.'
    PRINT*
		
		PRINT*, 'PCCG converged in',ITER

		! Get U_fe_full (all nodes) from U_fe (unknown nodes)
		DO j=1,num_of_nodes
			IF (node_type(2, j) > 0) THEN	! Unknown node
				U_fe_full(j) = U_fe(node_type(2,j))
			END IF
    END DO          
        
		! Calculate the residue
		CALL Sparse_Residue(A_stiff, U_fe, rhs_res, SIZE(A_stiff,1), SIZE(rhs,1), Residue)
		PRINT*, 'Max Absolute Residue (After) = ', Residue
		PRINT*

		DEALLOCATE(rhs, rhs_res)
        !$ ab.ZWZ 2021/11/9 if IFE_start added into loop \\
        IF (IMPIC_Index) THEN
            DEALLOCATE(A_stiff) 
            DEALLOCATE(A_diag) 
        ENDIF
        !$ ab.ZWZ 2021/11/9 if IFE_start added into loop //
ELSEIF (NSolver==1) THEN		! Gauss-Seidel

	PRINT*, 'Gasuss-Seidel Solver :'

	ALLOCATE(rhs_res(SIZE(rhs_FIX)))

	DO k=1,BGS_MaxIT        

		PRINT*, 'Iteration: ',k
			
		PRINT*, 'Forming PDE RHS Vector ....'

		CALL Global_RHS_PDE_2D(FUN_Density_rhs_2D,FUN_Density_rhs_TRIANGLE_2D, U_fe_full, R_fe_full, p_basic, t_c,	&
							   h_partition, node_type, num_of_unknowns, G_RHS_HW, G_Interpol_HW, rhs, delta, element_index,information_1,information_2)

		PRINT*, 'Forming PDE RHS Vector Done.'

		DO j=1,SIZE(rhs)
			! rhs_FIX = rhs_EBC
			rhs(j) = rhs(j)-rhs_FIX(j)
			rhs_res(j) = rhs(j)
		END DO

		RHSResidue = DSQRT(SUM(rhs*rhs))
		
		! Calculate Residue
		CALL Sparse_Residue(A_stiff, U_fe, rhs_res, SIZE(A_stiff,1), SIZE(rhs,1), Residue)		
		RelResidue = Residue/RHSResidue
		
		PRINT*, 'Relative Residue = ', RelResidue

		IF (RelResidue <= BGS_Tol) THEN
			PRINT*, 'Gasuss-Seidel Solver Converged in ', k
			EXIT
		END IF

		PRINT*, 'Solving the Linear System ....'

		IF (Blocking==1) THEN
    !		CALL IMSL_JPCG_Solver_Permuted( A_mass, node_permute, rhs, U_fe, A_diag, ITER )
		ELSE
			!CALL IMSL_JPCG_Solver( A_stiff, node_permute, rhs, U_fe, A_diag, ITER )
			IF(SLSolver==1) THEN
				!!!!CALL IMSL_JPCG_Solver( A_stiff, node_permute, rhs, U_fe, A_diag, ITER )
!!!!!				CALL IMSL_JPCG_Solver(A_stiff, rhs, U_fe, A_diag, SIZE(A_stiff,1), SIZE(rhs,1), ITER)
			ELSEIF(SLSolver==2) THEN
!!!!!				CALL DLAP_JPCG_Solver(A_stiff, rhs, U_fe, SIZE(A_stiff,1), SIZE(rhs,1), ITER)
			ELSEIF(SLSolver==3) THEN
				CALL MY_JPCG_Solver(A_stiff, rhs, U_fe, A_diag, SIZE(A_stiff,1), SIZE(rhs,1), ITER)				
			END IF
		END IF

		PRINT*, 'Solving the Linear System Done.'
		
		PRINT*, 'PCCG converged in',ITER

		! Get U_fe_full (all nodes) from U_fe (unknown nodes)
		DO j=1,num_of_nodes
			IF (node_type(2, j) > 0) THEN	! Unknown node
				U_fe_full(j) = U_fe(node_type(2,j))
			END IF
		END DO

		max_err_node	= MAXLOC(DABS(U_fe_old - U_fe))
		Err_Iter		= MAXVAL(DABS(U_fe_old(max_err_node)-U_fe(max_err_node)))

		PRINT*, 'Max Absolute Change  = ', Err_Iter
		PRINT*, 'Location of node of max error     = ',max_err_node
		PRINT*, 'New solution at node of max error = ',U_fe(max_err_node)
		PRINT*, 'Old solution at node of max error = ',U_fe_old(max_err_node)
		PRINT*, 'MAXVAL(U_fe)=',MAXVAL(U_fe), ' MINVAL(U_fe)=',MINVAL(U_fe)		
        
		IF (Err_Iter <= BGS_Tol) THEN            
			EXIT
		ELSE
			DO j=1,SIZE(U_fe)
				U_fe_old(j) = U_fe(j)
			END DO
		END IF
			
	END DO

	DEALLOCATE(rhs, rhs_res)


ELSEIF (NSolver==2) THEN	! Newton

	PRINT*, 'Newton Solver :'
	
	ALLOCATE(A_mass(SIZE(A_stiff)))
	DO j=1,SIZE(A_stiff)
		A_mass(j)%JCOL = A_stiff(j)%JCOL
		A_mass(j)%SROW = A_stiff(j)%SROW
	END DO

	DO k=1,Newton_MaxIT        

		PRINT*, 'Iteration: ',k
					
		PRINT*, 'Forming PDE RHS Vector ....'

		CALL Global_RHS_PDE_2D(FUN_Density_rhs_2D,FUN_Density_rhs_TRIANGLE_2D, U_fe_full, R_fe_full, p_basic, t_c,	&
							   h_partition, node_type, num_of_unknowns, G_RHS_HW, G_Interpol_HW, rhs, delta, element_index,information_1,information_2)

		PRINT*, 'Forming PDE RHS Vector Done.'
!
		DO j=1,SIZE(rhs)
!			! rhs_FIX = rhs_EBC
			rhs(j) = rhs(j)-rhs_FIX(j)
!			rhs_res(j) = rhs(j)
		END DO

		RHSResidue = DSQRT(SUM(rhs*rhs))

		! Calculate Residue
		DO i=1,SIZE(rhs)
			row_sum = Zero
			DO j=A_stiff(i)%SROW, A_stiff(i+1)%SROW-1
				row_sum = row_sum + A_stiff(j)%K * U_fe(A_stiff(j)%JCOL)
			END DO

			rhs(i) = rhs(i) - row_sum 
		END DO

		Residue = DSQRT(SUM(rhs*rhs))

		RelResidue = Residue/RHSResidue

		PRINT*, 'Relative Residue = ', RelResidue

		IF (RelResidue <= Newton_Tol) THEN
			PRINT*, 'Newton Solver Converged in ', k
			EXIT
		END IF

		PRINT*, 'Forming Mass Matrix ....'
		CALL CPU_TIME(CPUtime); PRINT*, " CPU Time [msec] = ",CPUtime*1000

		CALL Global_Mass( FUN_Density_Drhs_2D, U_fe_full, R_fe_full, p_basic, t_c, &
    					  h_partition, node_type, num_of_unknowns, G_Mass_HW, G_Interpol_HW, A_mass, element_index,information_1,information_2)

		PRINT*, 'Forming Mass Matrix Done.'

		DO j=1,SIZE(A_mass)
!			! dG/dU = K - dF/dU
			A_mass(j)%K = A_stiff(j)%K - A_mass(j)%K
		END DO

		! Prepare preconditioning matrix: daigonal matrix
		CALL Diagonal_Preconditioner(A_mass, A_diag)

		PRINT*, 'Solving the Linear System ....'

		! Initial Guess for Solution-change Vector (dU)
		DO j=1,SIZE(U_fe)
			U_fe(j)		= Zero
!			!U_fe_1(j)	= 0.0
		END DO

		IF (Blocking==1) THEN
!!			CALL IMSL_JPCG_Solver_Permuted( A_mass, node_permute, rhs, U_fe, A_diag, ITER )
		ELSE
!			!CALL PCCG_ICholesky_Solver(A_mass, rhs, U_fe_1, ITER, SIZE(A_mass,1), SIZE(rhs,1))
!			!CALL IMSL_JPCG_Solver( A_mass, node_permute, rhs, U_fe, A_diag, ITER )
!			IF(SLSolver==1) THEN
!				!CALL IMSL_JPCG_Solver( A_stiff, node_permute, rhs, U_fe, A_diag, ITER )
!				CALL IMSL_JPCG_Solver(A_mass, rhs, U_fe, A_diag, SIZE(A_mass,1), SIZE(rhs,1), ITER)
!			ELSEIF(SLSolver==2) THEN
!				CALL DLAP_JPCG_Solver(A_stiff, rhs, U_fe, SIZE(A_stiff,1), SIZE(rhs,1), ITER)
!			ELSEIF(SLSolver==3) THEN
				CALL MY_JPCG_Solver(A_mass, rhs, U_fe, A_diag, SIZE(A_mass,1), SIZE(rhs,1), ITER)
!			END IF
		END IF
		PRINT*, 'PCCG converged in',ITER

		PRINT*, 'Solving the Linear System Done.'

		max_err_node	= MAXLOC(DABS(U_fe))
		Err_Iter		= MAXVAL(DABS(U_fe(max_err_node)))

		DO j=1,SIZE(U_fe)
			U_fe(j) = U_fe(j) + U_fe_old(j)
		END DO
!		! Now U_fe = U_fe_new
!		! Get U_fe_full (all nodes) from U_fe (unknown nodes)
		DO j=1,num_of_nodes
			IF (node_type(2, j) > 0) THEN	! Unknown node
				U_fe_full(j) = U_fe(node_type(2,j))
			END IF
		END DO

		WRITE(6,*) 'Max Absolute Change  = ', Err_Iter
		WRITE(6,*) 'Location of node of max error     = ',max_err_node
		WRITE(6,*) 'New solution at node of max error = ',U_fe(max_err_node)
		WRITE(6,*) 'Old solution at node of max error = ',U_fe_old(max_err_node)
		WRITE(6,*) MAXVAL(U_fe), MINVAL(U_fe)		

		IF (Err_Iter <= Newton_Tol ) THEN
			PRINT*, 'Newton Solver Converged in ', k
			EXIT
		ELSE
			DO j=1,SIZE(U_fe)
				U_fe_old(j) = U_fe(j)
			END DO
		END IF
		
	END DO

	DEALLOCATE(A_mass, rhs)

ENDIF


!=====================neW ADD================================
!CALL Error_Analysis(FUN_sphere_u_ErrAna, U_fe_full, p_basic, t_basic_int,	&
!					t_iel, p_int_x, p_int_y, p_int_z, Global_Beta,		&
!					Error_2, Error_0, Error_1)
!======================dukun=============================================
!CALL Error_Analysis(   Global_Beta, FUN_sphere_u_ErrAna, U_fe_full, p_basic, t_c,	&
!                       h_partition, node_type, num_of_unknowns, delta, &!!!
!					   element_index,information_1,information_2, &
!					   Error_2, Error_0, Error_1)
!=========================================================================

!WRITE(6,*)
!WRITE(6,*) 'Element of Maximum L0-Error    = ',MAXLOC(Error_0)
!WRITE(6,*) '           Maximum L0-Error    = ',MAXVAL(Error_0)
!WRITE(6,*)
!WRITE(6,*) 'Element of Maximum L2-Error    = ',MAXLOC(Error_2)
!WRITE(6,*) '           Maximum L2-Error    = ',MAXVAL(Error_2)
!WRITE(6,*) '                   L2-Error    = ',DSQRT(SUM(Error_2))
!WRITE(6,*)
!WRITE(6,*) 'Element of Maximum H1x-Error   = ',MAXLOC(Error_1(:,1))
!WRITE(6,*) '           Maximum H1x-Error   = ',MAXVAL(Error_1(:,1))
!WRITE(6,*) '                   H1x-Error   = ',DSQRT(SUM(Error_1(:,1)))
!WRITE(6,*) 'Element of Maximum H1y-Error   = ',MAXLOC(Error_1(:,2))
!WRITE(6,*) '           Maximum H1y-Error   = ',MAXVAL(Error_1(:,2))
!WRITE(6,*) '                   H1y-Error   = ',DSQRT(SUM(Error_1(:,2)))
!WRITE(6,*) 'Element of Maximum H1z-Error   = ',MAXLOC(Error_1(:,3))
!WRITE(6,*) '                   H1z-Error   = ',DSQRT(SUM(Error_1(:,3)))
!WRITE(6,*)
!============================================================



!! Copy U_fe_full from IFE onto Phi for PIC
!DO ni=1,nnx
!	DO nj=1,nny
!		nindex = nj+(ni-1)*nny
!		Phi(ni,nj) = U_fe_full(nindex)
!	END DO
!END DO

DO ni=1, SIZE(HP,2)
  Phi(ni,1) = U_fe_full(ni)
END DO


DEALLOCATE(R_fe_full, U_fe_full, U_fe, U_fe_old)

!$ ab.ZWZ 2021/11/9 if IFE_start added into loop
!IF (.Not. IMPIC_INDEX)THEN
!    DEALLOCATE(rhs_FIX) 
!ENDIF
END