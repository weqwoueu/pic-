SUBROUTINE Sparse_Structure_xt(t_basic_int, HE, node_type, K_VROW, SMatrix, NZ)

USE IFE_MAIN_PARAM

IMPLICIT NONE

INTEGER, DIMENSION(:,:), INTENT(IN)		::	t_basic_int, node_type
INTEGER, DIMENSION(:,:), POINTER        ::	HE
INTEGER, DIMENSION(:,:), POINTER		::	K_VROW
TYPE(SPARSE), DIMENSION(:), ALLOCATABLE	::	SMatrix
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

		IF ( node_type(2,t_basic_int(i,k)) > 0 ) THEN	! This is an UNKNOWN node
    		NC = K_VROW(node_type(2,t_basic_int(i,k)),1)

    		DO j=1,n_nodes_in_elem
                     
    			IF ( node_type(2,t_basic_int(j,k)) < 0 ) THEN	! This is an KNOWN node

					found = .FALSE.

					DO n=2,NC+1

						IF ( K_VROW(node_type(2,t_basic_int(i,k)),n)== t_basic_int(j,k) ) THEN
							found = .TRUE.

						END IF
					END DO
					IF (.NOT. found) THEN
						NC = NC+1
						K_VROW(node_type(2,t_basic_int(i,k)),1)		= NC
						K_VROW(node_type(2,t_basic_int(i,k)),NC+1)	= t_basic_int(j,k)

					END IF
				END IF
			END DO
		END IF
	END DO

END DO




DO ee=1,n_edge

    k = HE(5, ee)
	DO i=1,n_nodes_in_elem

		IF ( node_type(2,t_basic_int(i,k)) > 0 ) THEN	! This is an UNKNOWN node
    		NC = K_VROW(node_type(2,t_basic_int(i,k)),1)

            IF (HE(6, ee) > 0) THEN
                    
                brother_element = HE(6, ee)

    		    DO j=1,n_nodes_in_elem

    			    IF ( node_type(2,t_basic_int(j,brother_element)) < 0 ) THEN	! This is an KNOWN node

					    found = .FALSE.

					    DO n=2,NC+1

						    IF ( K_VROW(node_type(2,t_basic_int(i,k)),n)==t_basic_int(j,brother_element) ) THEN
							    found = .TRUE.

						    END IF
					    END DO
					    IF (.NOT. found) THEN
						    NC = NC+1
						    K_VROW(node_type(2,t_basic_int(i,k)),1)		= NC
						    K_VROW(node_type(2,t_basic_int(i,k)),NC+1)	= t_basic_int(j,brother_element)

					    END IF
				    END IF
                END DO
            
            END IF     
            
		END IF
    END DO
    
    ! this_element is DG element but bro_element is not DG element
    IF (HE(6, ee) > 0) THEN
        
        IF (t_basic_int(5, HE(5, ee)) >= 1 .AND. t_basic_int(5, HE(6, ee)) == 0) THEN
              
            k = HE(6, ee)
            DO i=1,n_nodes_in_elem

		        IF ( node_type(2,t_basic_int(i,k)) > 0 ) THEN	! This is an UNKNOWN node
    		        NC = K_VROW(node_type(2,t_basic_int(i,k)),1)

                    IF (HE(6, ee) > 0) THEN
                    
                        brother_element = HE(5, ee)

    		            DO j=1,n_nodes_in_elem

    			            IF ( node_type(2,t_basic_int(j,brother_element)) < 0 ) THEN	! This is an KNOWN node

					            found = .FALSE.

					            DO n=2,NC+1

						            IF ( K_VROW(node_type(2,t_basic_int(i,k)),n)==t_basic_int(j,brother_element) ) THEN
							            found = .TRUE.

						            END IF
					            END DO
					            IF (.NOT. found) THEN
						            NC = NC+1
						            K_VROW(node_type(2,t_basic_int(i,k)),1)		= NC
						            K_VROW(node_type(2,t_basic_int(i,k)),NC+1)	= t_basic_int(j,brother_element)

					            END IF
				            END IF
                        END DO
            
                    END IF
            
                END IF
            
            END DO
            
        END IF
        
    END IF
    
    
END DO


NZ = SUM(K_VROW(:,1),1)
!WRITE(*,*) '========',NZ
IF (.NOT.ALLOCATED(SMatrix)) THEN
    IF(NZ < n_n .OR. NZ == n_n)THEN
        ALLOCATE(SMatrix(n_n+1))
    ELSE
        ALLOCATE(SMatrix(NZ))
    ENDIF    
END IF
!WRITE(*,*) '========',SIZE(SMatrix)

!ALLOCATE(SMatrix(NZ))

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

! -------------- wsy test IDG-------------------
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

END