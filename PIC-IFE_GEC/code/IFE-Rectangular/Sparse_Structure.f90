SUBROUTINE Sparse_Structure(t_basic_int, HE, node_type, K_VROW, SMatrix, NZ)

USE IFE_MAIN_PARAM

IMPLICIT NONE

INTEGER, DIMENSION(:,:), INTENT(IN)		::	t_basic_int, node_type
INTEGER, DIMENSION(:,:), POINTER        ::	HE
INTEGER, DIMENSION(:,:), POINTER		::	K_VROW
TYPE(SPARSE), DIMENSION(:), POINTER		::	SMatrix
INTEGER, INTENT(OUT)					::	NZ

INTEGER		i, j, k, n, n_e, n_n, NC, n_nodes_in_elem, n_edge, brother_element, br, ee
LOGICAL		found

n_e	= SIZE(t_basic_int,2)
n_n = MAXVAL(node_type(2,:))
n_edge = SIZE(HE,2)

n_nodes_in_elem	= 4

ALLOCATE(K_VROW(n_n,Max_Nodal_Connect + 2))

DO i = 1, SIZE(K_VROW,1)
	K_VROW(i,:) = 0
END DO

DO k=1,n_e 

	DO i=1,n_nodes_in_elem

		If ( node_type(2,t_basic_int(i,k)) > 0 ) Then	! This is an UNKNOWN node
    		NC = K_VROW(node_type(2,t_basic_int(i,k)),1)

    		DO j=1,n_nodes_in_elem

    			If ( node_type(2,t_basic_int(j,k)) > 0 ) Then	! This is an UNKNOWN node

					found = .FALSE.

					DO n=2,NC+1

						If ( K_VROW(node_type(2,t_basic_int(i,k)),n)==node_type(2,t_basic_int(j,k)) ) Then
							found = .TRUE.

						END If
					END DO
					If (.NOT. found) Then
						NC = NC+1
						K_VROW(node_type(2,t_basic_int(i,k)),1)		= NC
						K_VROW(node_type(2,t_basic_int(i,k)),NC+1)	= node_type(2,t_basic_int(j,k))

					END If
				END If
			END DO
		END If
	END DO

END DO


DO ee=1,n_edge 

    k = HE(5, ee)
	DO i=1,n_nodes_in_elem

		If ( node_type(2,t_basic_int(i,k)) > 0 ) Then	! This is an UNKNOWN node
    		NC = K_VROW(node_type(2,t_basic_int(i,k)),1)

            If (HE(6, ee) > 0) Then
                    
                brother_element = HE(6, ee)

    		    DO j=1,n_nodes_in_elem

    			    If ( node_type(2,t_basic_int(j,brother_element)) > 0 ) Then	! This is an UNKNOWN node

					    found = .FALSE.

					    DO n=2,NC+1

						    If ( K_VROW(node_type(2,t_basic_int(i,k)),n)==node_type(2,t_basic_int(j,brother_element)) ) Then
							    found = .TRUE.

						    END If
					    END DO
					    If (.NOT. found) Then
						    NC = NC+1
						    K_VROW(node_type(2,t_basic_int(i,k)),1)		= NC
						    K_VROW(node_type(2,t_basic_int(i,k)),NC+1)	= node_type(2,t_basic_int(j,brother_element))

					    END If
				    END If
                END DO
            
            END If
            
		END If
    END DO

    If (HE(6, ee) > 0) Then
    ! this_element is DG element but bro_element is not DG element
        If (t_basic_int(5, HE(5, ee)) >= 1 .AND. t_basic_int(5, HE(6, ee)) == 0) Then
        
            k = HE(6, ee)
            DO i=1,n_nodes_in_elem

		        If ( node_type(2,t_basic_int(i,k)) > 0 ) Then	! This is an UNKNOWN node
    		        NC = K_VROW(node_type(2,t_basic_int(i,k)),1)

                    If (HE(6, ee) > 0) Then
                    
                        brother_element = HE(5, ee)

    		            DO j=1,n_nodes_in_elem

    			            If ( node_type(2,t_basic_int(j,brother_element)) > 0 ) Then	! This is an UNKNOWN node

					            found = .FALSE.

					            DO n=2,NC+1

						            If ( K_VROW(node_type(2,t_basic_int(i,k)),n)==node_type(2,t_basic_int(j,brother_element)) ) Then
							            found = .TRUE.

						            END If
					            END DO
					            If (.NOT. found) Then
						            NC = NC+1
						            K_VROW(node_type(2,t_basic_int(i,k)),1)		= NC
						            K_VROW(node_type(2,t_basic_int(i,k)),NC+1)	= node_type(2,t_basic_int(j,brother_element))

					            END If
				            END If
                        END DO
            
                    END If
            
                END If
            
            END DO
        
        END If
    END If
END DO





NZ = SUM(K_VROW(:,1),1)


ALLOCATE(SMatrix(NZ))

DO i = 1, SIZE(SMatrix)
	SMatrix(i)%K	= Zero
	SMatrix(i)%JCOL	= 0
	SMatrix(i)%SROW	= 0
END DO


SMatrix(1)%SROW = 1
k = 0
DO i=1,n_n
	DO j=2,K_VROW(i,1)+1

		k = k+1
		SMatrix(k)%JCOL = K_VROW(i,j)	

	END DO

	SMatrix(i+1)%SROW = SMatrix(i)%SROW+K_VROW(i,1)		

END DO

! -------------- wsy test SIDG-------------------
!WRITE(6,*) K_VROW(1,:)
!WRITE(6,*) K_VROW(2,:)
!WRITE(6,*) K_VROW(3,:)
!WRITE(6,*) K_VROW(4,:)
!WRITE(6,*) K_VROW(5,:)
!WRITE(6,*) K_VROW(6,:)
!WRITE(6,*) K_VROW(7,:)
!WRITE(6,*) K_VROW(8,:)
!WRITE(6,*) K_VROW(9,:)
!WRITE(6,*) K_VROW(10,:)
!WRITE(6,*) K_VROW(11,:)
!WRITE(6,*) K_VROW(12,:)
! --------------------------------------------------

DEALLOCATE(K_VROW)

END Subroutine Sparse_Structure