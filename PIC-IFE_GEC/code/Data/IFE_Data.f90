!					==================================
!					||      IFE_Data_Interface      ||
!					==================================


MODULE IFE_Data

USE IFE_MAIN_PARAM

!LY Add for Multigrid Store, 2022-1-17
USE Cell_Data_2D
!LY Add for Multigrid Store, 2022-1-17

REAL(8) :: con_penalty

INTEGER, DIMENSION(:,:), POINTER	  ::	t_basic_int, t_iel, node_type, node_permute, t_c, t_basic, &
                                        DGT_IFE_partition, HT, DGT, AHT, DG_flag, HT_flag

!LY ADD FOR SIDG-PIC COUPLE: Phi, Rho, Rho_s. 2022-1-14
INTEGER :: Field_Size
!LY ADD FOR SIDG-PIC COUPLE: Phi, Rho, Rho_s. 2022-1-14

INTEGER, DIMENSION(:,:), POINTER	  ::	e_basic, HE, DGE, AHE, E_Dirichlet

REAL(8), DIMENSION(:,:), POINTER	  ::	p_basic, p_int_x, p_int_y, p_int_z, HP, P_average, DGP, AHP, P_flag

TYPE(SPARSE), DIMENSION(:), POINTER	::	A_stiff, A_mass, A_stiff_xt, A_stiff_Test

REAL(8), DIMENSION(:), POINTER		  ::	A_diag, rhs_FIX, Global_Beta

!==========LY Add for Multigrid Store, 2022-1-17==========
TYPE(CellDataType), DIMENSION(:), POINTER :: CellMesh, CellMesh_Temp
INTEGER :: N_Tier_Old, N_Tier_Pres
!==========LY Add for Multigrid Store, 2022-1-17==========

!=========LY modification for Speeding Particle Positioning, 2022-7-25=========
Integer, Dimension(:,:), Pointer :: ElementParentChild, ElementParentChildTemp
Integer, Dimension(:,:), Pointer :: EdgeElement, EdgeElementTemp
Integer, Dimension(:), Pointer :: EdgeParentNumber
Integer, Dimension(:,:), Pointer :: EdgeParent
!=========LY modification for Speeding Particle Positioning, 2022-7-25=========

!=========LY modification for Periodic Boundary Condition, 2022-7-25=========
Integer, Dimension(:,:), Pointer :: PairNodes
!=========LY modification for Periodic Boundary Condition, 2022-7-25=========

!=========LY modification for dealing with Floating Dirichlet boundary condition, 2022-7-25=========
Real(8), Dimension(:), Allocatable :: E_DirichletValue
!=========LY modification for dealing with Floating Dirichlet boundary condition, 2022-7-25=========

!=========LY modification for Phi-Periodical case, 2022-8-11=========
Integer, Dimension(:), Allocatable :: PhiPeriodicalNode
!=========LY modification for Phi-Periodical case, 2022-8-11=========

!=================MY NEW ADD=======================================================
INTEGER, DIMENSION(:), POINTER				::    element_index, node_index, edge_index, edge_index_D
INTEGER,DIMENSION(:,:),POINTER              ::    information_1
REAL(8),DIMENSION(:,:),POINTER              ::    information_2, information_3, &
                                                    information_3_D
REAL(8)	                                          h_partition(2)
!==================================================================================

!================CYC NEW ADD FOR FILTER============================================
INTEGER, DIMENSION(:), POINTER				::	element_index_for_curve
INTEGER, DIMENSION(:,:), POINTER			::	information_1_for_curve
REAL(8), DIMENSION(:,:), POINTER			::	information_2_for_curve
!==================================================================================

! Hardwire Global Arrays
REAL(8)		G_Stiff_HW(2,4,4), G_RHS_HW(2,4,4), G_Mass_HW(4,4,4), G_Interpol_HW(4,4)

INTEGER, DIMENSION(:), POINTER				::	bnd_elem_index  !$ ab.ZWZ for move EBC assembling to IFE_Solve

!========================== WSY ADD FOR PPR ============================
INTEGER, DIMENSION(:),POINTER           :: node_intrs
REAL(8), DIMENSION(:,:), ALLOCATABLE    :: coef_ppr
REAL(8), DIMENSION(:,:), POINTER        :: efx_1,efy_1
REAL(8), DIMENSION(:,:), POINTER        :: efx_2,efy_2
INTEGER, DIMENSION(:,:), POINTER	    :: HT_PPR
INTEGER(8)                              :: n_node_dg_ppr
INTEGER, DIMENSION(:), ALLOCATABLE      :: node_type_cg
REAL(8), DIMENSION(:,:), ALLOCATABLE    :: P_average_PPR 
INTEGER, DIMENSION(:,:), ALLOCATABLE    :: point_patch
INTEGER, DIMENSION(:), ALLOCATABLE      :: re_times
INTEGER, DIMENSION(:), POINTER          :: PPRPointIndex
!========================================================================

!========================== WSY ADD FOR OUTPUT ============================
Real(8), Allocatable :: node_map(:,:)
!========================================================================

! ===================wsy add for ParticlePositioning ======================
Integer, Allocatable ::  Element_Map(:,:)
Integer              ::  repeat_refinement
Real(8) :: Mapxmin, Mapymin
! ============================================================
END MODULE