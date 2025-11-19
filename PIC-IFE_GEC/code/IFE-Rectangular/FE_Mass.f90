SUBROUTINE	FE_Mass(RHS_FUN, U_val, R_val, el_region, &
					n_nodes_in_elem_1, n_nodes_in_elem_2, mass_hw, interpol_hw, matrix)


USE	IFE_MAIN_PARAM

IMPLICIT NONE

EXTERNAL					RHS_FUN
INTEGER, INTENT(IN)		::	el_region,	&
							n_nodes_in_elem_1, n_nodes_in_elem_2
REAL(8), INTENT(IN)		::	U_val(n_nodes_in_elem_1), R_val(n_nodes_in_elem_1)
REAL(8), INTENT(OUT)	::	matrix(4,4)
REAL(8), INTENT(IN)		::	mass_hw(4,4,4), interpol_hw(4,4)
!REAL(8), INTENT(IN)		::	gnodes(2,3)

REAL(8)									::	V, eint
REAL(8), DIMENSION(:), ALLOCATABLE		::	rhs_val, u_gnodes, r_gnodes

INTEGER									::	i, j

ALLOCATE(u_gnodes(4), r_gnodes(4))
u_gnodes = Zero
r_gnodes = Zero
DO i=1,4
	u_gnodes = u_gnodes + u_val(i)*interpol_hw(i,:)
	r_gnodes = r_gnodes + r_val(i)*interpol_hw(i,:)
END DO
!print *,u_gnodes
!print *,r_gnodes
!stop
ALLOCATE(rhs_val(4))
CALL RHS_FUN(r_gnodes, u_gnodes, el_region, rhs_val, SIZE(u_gnodes, 1))

!print *,rhs_val
!stop


matrix = Zero
DO i=1,n_nodes_in_elem_1    
   		            
   	DO j=1, i
	
		eint = SUM(mass_hw(i,j,:)*rhs_val)
    		        
		matrix(i,j) = matrix(i,j) + eint
				
	END DO

END DO
!print *,matrix
!stop
DO i=1,n_nodes_in_elem_1    
   		            
   	DO j=i+1,n_nodes_in_elem_2
	
		matrix(i,j) = matrix(j,i)
				
	END DO

END DO
!print *,matrix
!stop
DEALLOCATE(rhs_val, u_gnodes, r_gnodes)

END