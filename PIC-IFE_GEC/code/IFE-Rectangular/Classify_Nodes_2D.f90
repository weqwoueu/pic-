SUBROUTINE Classify_Nodes_2D(	p_basic, node_type,	node_loc, node_permute )

! Purpose:		Classify mesh nodes and generate node_type and node_permute arrays
! Last Update:	10/5/2003 04:15 AM


USE IFE_MAIN_PARAM
USE Object_Data_2D

IMPLICIT NONE

REAL(8), DIMENSION(:,:), POINTER	::	p_basic
INTEGER, DIMENSION(:), POINTER		::	node_loc
INTEGER, DIMENSION(:,:), POINTER	::	node_type, node_permute

INTEGER		n_nodes, n_unknowns, n, count , vnode_count, cnode_count, n_vnodes, n_cnodes 

! SIZE(node_type) = 2 x num_of_nodes

! node_type(1,:) =	1		Drichilet boundary node,	
!					0		Neumann boundary node,	
!					-200	Interior node.
! node_type(2,:) =	>0		Unknown node index for Interior and Neumann boundary nodes,
!				 =	-1		Drichilet boundary node.

n_nodes		= SIZE(p_basic,2);
n_unknowns	= MAXVAL(node_type(2,:))


! * Assign node lcoations to interior nodes
! * Assign Drichilet boundary node index to Neumann nodes which touch conducting object
!!Y.C. Delete temporay: if here undelete, delete Y.C add, keep here, NON-vacuum point is treated as fix point
DO n=1,n_nodes
	IF		(node_type(1,n)<0) THEN		! * Interior node (either known or unknown) 
			node_type(1,n) = node_loc(n)	!    then set node_type(1,n) = node_loc(n) < 0
	ELSEIF  (node_type(1,n)==0) THEN	! * Neumann boundary node
		IF (node_loc(n)<vacuum) THEN		! If it belongs to an object,
			node_type(1,n) = 1				!   then reassign it as a Drichilet node
			node_type(2,n) = node_loc(n)	!   and set node_type(2,n) = node_loc(n) < 0
		ELSE
			node_type(1,n) = node_loc(n)	! else, set node_type(1,n) = vacuum = -1
		END IF
	ELSE								! * Drichilet boundary node
										!     then do nothing
										!     [ We might also set 
										!       node_type(1,n) = node_loc(n) < 0 ]
	ENDIF
END DO
!!Y.C. Delete temporay: if here undelete, delete Y.C add, keep here, NON-vacuum point is treated as fix point


! Now we have
! node_type(1,:) =	1		Drichilet boundary node, or 
!							Neumann boundary node belonging to an object
!					<0		node location for Interior nodes, or
!							node location for Neumann boundary nodes
! node_type(2,:) =	>0		Unknown node index for Interior nodes, or
!							Unknown node index for Neumann boundary nodes
!				 =	-1		Drichilet boundary node.
!				 =  <-1     node location for Neumann boundary nodes belonging
!							to objects


! N.B.: node_type(2,:) may be also used to identify the location of Drichilet nodes 
!		by setting it to a negative integer

!!Y.C. Delete temporay: if here undelete, delete Y.C add, keep here, NON-vacuum point is treated as fix point
count = 0
DO n=1,n_nodes
	IF (node_type(1,n)>=vacuum .AND. node_type(1,n)<0) THEN		! * Neumann or Interior floating node
		count = count + 1		!				increase unknowns count by one
		node_type(2,n) = count	!				set node_type(2,n) = count
	!  RKBegin
	ELSEIF (node_type(1,n)<vacuum) THEN		! * Interior fixed node
			node_type(1,n) = 1				!   then reassign it as a Drichilet node
			node_type(2,n) = node_loc(n)	!   and set node_type(2,n) = node_loc(n) < 0
	!  RKEnd
	ELSE 
		!node_type(2,n) = -1		! * Drichilet node
									!				do nothing
	END IF
END DO
!!Y.C. Delete temporay: if here undelete, delete Y.C add, keep here, NON-vacuum point is treated as fix point

!!Y.C. Add: if here undelete, delete Y.C. Delete temporay, keep here, NON-vacuum point is treated as float point
!count = 0
!DO n=1,n_nodes
!!	IF (node_type(1,n)<=vacuum) THEN		! * Neumann or Interior floating node
!	IF (node_type(1,n)<=0) THEN		! * Neumann or Interior floating node
!		count = count + 1		!				increase unknowns count by one
!		node_type(2,n) = count	!				set node_type(2,n) = count
!!	!  RKBegin
!!	ELSEIF (node_type(1,n)<vacuum) THEN		! * Interior fixed node
!!			node_type(1,n) = 1				!   then reassign it as a Drichilet node
!!			node_type(2,n) = node_loc(n)	!   and set node_type(2,n) = node_loc(n) < 0
!!	!  RKEnd
!!	ELSE 
!!		!node_type(2,n) = -1		! * Drichilet node
!!									!				do nothing
!	END IF
!END DO
!!Y.C. Add: if here undelete, delete Y.C. Delete temporay, keep here, NON-vacuum point is treated as float point


! Now we have
! node_type(1,:)	same as before
! node_type(2,:)	same as before but unknown node indeces have been recounted

n_unknowns	= MAXVAL(node_type(2,:))	! The number of unknowns recalculated

ALLOCATE(node_permute(2,n_unknowns))
node_permute = 0

vnode_count = 0
cnode_count = 0

!!Y.C. Delete
!DO n=1,n_nodes
!!print*,'n=',n, node_type(1,n)
!
!	IF		(node_type(1,n)==vacuum) THEN	! * Neumann boundary node or Interior node
!											!   belonging to vacuum
!		vnode_count  = vnode_count+1		!      increase vnode_count by one
!		node_permute(1,node_type(2,n)) = vnode_count
!											!      set node_permute(1,node_type(2,n)) =
!											!	      vnode_count, where
!											!         node_type(2,n) = unknown node index
!!print*,'v=',vnode_count, node_permute(1,node_type(2,n)) 
!!print*, '  '
!!pause
!	ELSEIF	(node_type(1,n)<vacuum) THEN	! * Interior node belonging to an object 
!		cnode_count = cnode_count+1			!      increase cnode_count by one
!		node_permute(2,node_type(2,n)) = cnode_count
!											!      set node_permute(2,node_type(2,n)) =
!											!	      cnode_count, where
!											!         node_type(2,n) = unknown node index	ENDIF											
!!print*,'c=',cnode_count, node_permute(2,node_type(2,n)) 
!!print*,'  '
!	END IF
!END DO
!!Y.C. Delete

!Y.C. Change
DO n=1,n_nodes

	IF		((node_type(1,n)<=0).AND.(node_loc(n)==vacuum)) THEN	! * Neumann boundary node or Interior node
											!   belonging to vacuum
		vnode_count  = vnode_count+1		!      increase vnode_count by one
		node_permute(1,node_type(2,n)) = vnode_count
											!      set node_permute(1,node_type(2,n)) =
											!	      vnode_count, where
											!         node_type(2,n) = unknown node index

	ELSEIF	((node_type(1,n)<=0).AND.(node_loc(n)<vacuum)) THEN	! * Interior node belonging to an object 
		cnode_count = cnode_count+1			!      increase cnode_count by one
		node_permute(2,node_type(2,n)) = cnode_count
											!      set node_permute(2,node_type(2,n)) =
											!	      cnode_count, where
											!         node_type(2,n) = unknown node index	ENDIF											

	END IF
END DO
!Y.C. Add


n_vnodes = MAXVAL(node_permute(1,:))
n_cnodes = MAXVAL(node_permute(2,:))

!print*, node_permute(2,:)
!print*, n_vnodes, n_cnodes

! On leaving this subroutine, we have

! node_type(1,:) =	1		Drichilet boundary node, or 
!							Neumann boundary node belonging to an object
!					<0		node location for Interior nodes, or
!							node location for Neumann boundary nodes
! node_type(2,:) =	>0		Unknown node index for Interior nodes, or
!							Unknown node index for Neumann boundary nodes
!				 =	-1		Drichilet boundary node.
!				 =  <-1     node location for Neumann boundary nodes belonging
!							to objects


! node_permute(2, num_of_unknown_nodes)
! node_permute(1,:) =	>0	Unknown node index for Interior nodes, or
!							Unknown node index for Neumann boundary nodes in < Vacuum >
!					=	0	Node is NOT in < Vacuum >
! node_permute(2,:)	=	>0	Unknown node index for Interior nodes, or
!							Unknown node index for Neumann boundary nodes in 
!							< Conductor-Object >
!					=	0	Node is NOT in < Conductor-Object >


END