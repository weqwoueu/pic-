SUBROUTINE periodic_boundary_conditions(nnx,nny,bc_index,p_basic,node_type,bc_point_1,bc_point_2,ALQ,A_stiff,f)

!
!Purpose:
!  removal of constraint periodic boundary condtions
!  the correspond boundary's relation xj=alpha*xi+Q

USE IFE_INTERFACE, ONLY: periodic_boundary_corresponding_nodes,solve_periodic_boundary_conditions
USE IFE_MAIN_PARAM
IMPLICIT NONE 

INTEGER nnx,nny
INTEGER, DIMENSION(:),POINTER	::	bc_index
REAL(8), DIMENSION(:,:),POINTER	      ::	p_basic
INTEGER, DIMENSION(:,:), POINTER  ::node_type
TYPE(SPARSE), DIMENSION(:),POINTER		::	A_stiff
!ALQ: record alpha and Q
REAL(8),DIMENSION(:,:),POINTER  ::	ALQ

REAL(8), DIMENSION(:,:),POINTER	::	bc_point_1, bc_point_2

REAL(8), DIMENSION(:),POINTER	::f

!unknow_nodes(1,2,:):the correspond node and their alpha and Q
REAL(8)	,DIMENSION(:,:),POINTER 	::	unknow_nodes
INTEGER num_unknows

REAL(8) AL,Q
INTEGER m,i,j,num,n

CALL periodic_boundary_corresponding_nodes(nnx,nny,bc_index,p_basic,node_type,ALQ,bc_point_1,bc_point_2,unknow_nodes)


DO m=1,SIZE(unknow_nodes,2)
	IF(unknow_nodes(1,m)<unknow_nodes(2,m))THEN
		j=unknow_nodes(1,m)
		i=unknow_nodes(2,m)
		AL=unknow_nodes(3,m)
		Q=unknow_nodes(4,m)
	ELSEIF(unknow_nodes(1,m)>=unknow_nodes(2,m))THEN
		i=unknow_nodes(1,m)
		j=unknow_nodes(2,m)
		AL=1.0/unknow_nodes(3,m)
		Q=-unknow_nodes(4,m)/unknow_nodes(3,m)
	ENDIF
	num=0
	DO n=1,j
		IF(node_type(2,n)==-1)THEN
			num=num+1
		ENDIF
	ENDDO

	j=j-num

	num=0
	DO n=1,i
		IF(node_type(2,n)==-1)THEN
			num=num+1
		ENDIF
	ENDDO

	i=i-num

	
	num_unknows=SIZE(p_basic,2)

CALL solve_periodic_boundary_conditions(i,j,AL,Q,num_unknows ,A_stiff,f)
ENDDO

DEALLOCATE(unknow_nodes)

END

