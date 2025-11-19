SUBROUTINE FE_Interpol_HW_2D(vert,  hardwire)

USE IFE_Data
USE Gauss_Data
USE IFE_INTERFACE, ONLY: Linear_FE_Basis_Eval_2D, Generate_Gauss_local, Generate_Gauss_reference
IMPLICIT NONE



REAL(8), INTENT(IN)					::	vert(2,4)
REAL(8), INTENT(OUT)				::	hardwire(4,4)

INTEGER	i
REAL(8), DIMENSION(:), ALLOCATABLE		::	basis
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	basis_coef, gnodes
INTEGER, DIMENSION(:), ALLOCATABLE		::	dind(:)

!===============================================================================
INTEGER                                           Gauss_point_number
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_reference
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_reference
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_reference_triangle

REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_local
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_local
!================================================================================

ALLOCATE(basis(4))

ALLOCATE(basis_coef(4,4))
CALL Linear_FE_Basis_Coeff_2D(vert, basis_coef)

ALLOCATE(gnodes(2,4))
    
    ! wsy revise for sidg
    h_partition(1) = vert(1,2) - vert(1,1)
    h_partition(2) = vert(2,4) - vert(2,1)

   CALL Generate_Gauss_local(Gauss_coefficient_reference_Four,Gauss_point_reference_Four, &
                             vert(1:2,1),h_partition, Gauss_coefficient_local_Four, Gauss_point_local_Four)

gnodes=TRANSPOSE(Gauss_point_local_Four)

ALLOCATE(dind(2))
dind = (/0, 0/)

DO i=1,4

	CALL Linear_FE_Basis_Eval_2D( basis_coef, gnodes(1,:), gnodes(2,:), i,		&
							      dind(1), dind(2), basis)

	hardwire(i,:) = basis

END DO

DEALLOCATE(basis, basis_coef,gnodes, dind)

END